"""
Filter MT-SNVs from per-donor allelic / coverage tables.
"""

import re
import argparse
import numpy as np
import pandas as pd
import pysam


##


def parse_args():

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--allelic_table',  required=True, help='Per-donor allelic table (tsv.gz)')
    p.add_argument('--coverage_table', required=True, help='Per-donor coverage table (tsv.gz)')
    p.add_argument('--ref',            required=True, help='Faidx-indexed reference FASTA')
    p.add_argument('--donor',          required=True, help='Donor identifier')
    p.add_argument('--min_strand_ratio',         type=float, default=0.1)
    p.add_argument('--max_strand_ratio',         type=float, default=0.9)
    p.add_argument('--af_threshold',             type=float, default=0.05)
    p.add_argument('--min_callable_coverage',    type=float, default=10)
    p.add_argument('--max_sb_pval',              type=float, default=0.05)
    p.add_argument('--max_prevalence_low_AF',    type=float, default=0.25)
    p.add_argument('--min_gs_count_low_AF',      type=int,   default=1)
    p.add_argument('--min_overlap_count_low_AF', type=int,   default=1)
    p.add_argument('--min_n_mutations',          type=int,   default=1)
    p.add_argument('--max_prevalence_germline',  type=float, default=0.9)

    return p.parse_args()


##


def infer_mt_contig(path_ref):

    refs = pysam.FastaFile(path_ref).references
    for c in ('chrM', 'MT'):
        if c in refs:
            return c
    raise ValueError("No chrM/MT contig found in reference")


##


def annotate_3nt(df, path_ref, region):
    """
    Annotate 3nt context for each MT-SNV in the allelic table.
    """

    ref = pd.DataFrame(
        [ (i+1, x) for i, x in enumerate(pysam.FastaFile(path_ref).fetch(region)) ],
        columns=['POS', 'REF']
    )
    COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    ref_dict = ref.set_index('POS')['REF'].to_dict()
    chrM_len = max(ref_dict)

    # Strand (H/L) and SBS96 trinucleotide context in pyrimidine convention:
    #   rCRS is the L-strand sequence, so:
    #   L strand: REF is pyrimidine (C/T) -> read directly from reference (L-strand)
    #   H strand: REF is purine (A/G)     -> reverse-complement to pyrimidine (H-strand)
    def _annotate(row):

        pos, ref_base, alt_base = row['POS'], row['REF'], row['ALT']
        prev = ref_dict[pos - 1 if pos > 1 else chrM_len]
        nxt  = ref_dict[pos + 1 if pos < chrM_len else 1]

        if ref_base in ('C', 'T'):
            return 'L', f"{ref_base}>{alt_base}", f"{prev}[{ref_base}>{alt_base}]{nxt}", f"{prev}{ref_base}{nxt}"

        r, a = COMP[ref_base], COMP[alt_base]
        return 'H', f"{r}>{a}", f"{COMP[nxt]}[{r}>{a}]{COMP[prev]}", f"{COMP[nxt]}{r}{COMP[prev]}"

    df[['strand', 'MUT_TYPE', 'SBS96', '3nt_context']] = df.apply(_annotate, axis=1, result_type='expand')

    return df


##


def main():

    args = parse_args()
    donor = re.sub(r'\.v[^.]*\.dupmarked.*$', '', args.donor)

    df  = pd.read_csv(args.allelic_table, sep='\t')
    cov = pd.read_csv(args.coverage_table, sep='\t')

    metrics = {}

    df = df.query('REF!="N"').copy()
    df = df.merge(cov, on=['sample', 'POS', 'REF'], how='left')

    df['donor']  = donor
    df['sample'] = df['sample'].astype(str).str.replace(r'\.v[^.]*\.dupmarked.*$', '', regex=True)
    metrics['n_unfiltered_samples'] = df['sample'].nunique()

    df['MUT'] = df['POS'].astype(str) + '_' + df['REF'] + '_' + df['ALT']
    df['AF']  = df['total_count'] / df['callable_coverage']
    df['sample'] = pd.Categorical(df['sample'])
    metrics['n_unfiltered_calls'] = df.shape[0]
    metrics['n_unfiltered_muts']  = df['MUT'].nunique()

    region = infer_mt_contig(args.ref)
    df = annotate_3nt(df, args.ref, region)

    n_samples = df['sample'].nunique()
    df = df.merge(df.groupby('MUT').size().to_frame('n_samples').reset_index(), on='MUT', how='left')
    df = df.merge((df.groupby('MUT').size() / n_samples).to_frame('prevalence').reset_index(), on='MUT', how='left')

    # Stage I: base filtering on summary stats
    df_high_AF = df.loc[
        (df['strand_ratio']>=args.min_strand_ratio) &
        (df['strand_ratio']<=args.max_strand_ratio) &
        (df['AF']>=args.af_threshold) &
        (df['callable_coverage']>=args.min_callable_coverage) &
        (df['sb_pval']>=args.max_sb_pval)
    ]
    df_low_AF = df.loc[
        (df['strand_ratio']>=args.min_strand_ratio) &
        (df['strand_ratio']<=args.max_strand_ratio) &
        (df['AF']<args.af_threshold) &
        (df['prevalence']<=args.max_prevalence_low_AF) &
        (df['gs>1_count']>=args.min_gs_count_low_AF) &
        (df['overlap_count']>=args.min_overlap_count_low_AF) &
        (df['callable_coverage']>=args.min_callable_coverage) &
        (df['sb_pval']>=args.max_sb_pval)
    ]
    metrics['n_high_AF_muts'] = df_high_AF['MUT'].nunique()
    metrics['n_low_AF_muts']  = df_low_AF['MUT'].nunique()

    df_filtered = pd.concat([df_high_AF, df_low_AF]).drop_duplicates()

    # Stage II: remove outlier samples
    n_muts_per_sample = df_filtered.groupby('sample')['MUT'].nunique()
    median = n_muts_per_sample.median()
    mad = np.median(np.abs(n_muts_per_sample - median))
    chosen_samples = (
        n_muts_per_sample
        .loc[lambda x: (x>=args.min_n_mutations) & (x<=median + 5*mad)]
        .index
    )
    df_filtered = df_filtered.loc[df_filtered['sample'].isin(chosen_samples)].copy()

    n_samples = df_filtered['sample'].nunique()
    df_filtered.drop(columns=['n_samples', 'prevalence'], inplace=True)
    df_filtered = df_filtered.merge(df_filtered.groupby('MUT').size().to_frame('n_samples').reset_index(), on='MUT', how='left')
    df_filtered = df_filtered.merge((df_filtered.groupby('MUT').size() / n_samples).to_frame('prevalence').reset_index(), on='MUT', how='left')

    # Stage II: exclude MT-SNVs with too variable AF in positive samples
    exclude = (
        df_filtered.groupby('MUT')['AF']
        .apply(lambda x: x.max() / x.min()).loc[lambda x: x>5]
        .index
    )
    df_filtered = df_filtered.loc[~df_filtered['MUT'].isin(exclude)].copy()

    df_filtered.to_csv(f'{donor}.confident_allelic_table.tsv.gz', sep='\t', index=False)

    # Stage III: germline
    df_filtered = df_filtered.loc[df_filtered['prevalence']<=args.max_prevalence_germline].copy()
    metrics['n_filtered_samples'] = df_filtered['sample'].nunique()
    metrics['n_filtered_calls']   = df_filtered.shape[0]
    metrics['n_filtered_muts']    = df_filtered['MUT'].nunique()

    df_filtered.to_csv(f'{donor}.final_allelic_table.tsv.gz', sep='\t', index=False)

    (
        pd.DataFrame(
            [(donor, k, v) for k, v in metrics.items()],
            columns=['donor', 'metric', 'value'],
        )
        .to_csv(f'{donor}.metrics.txt', sep='\t', index=False)
    )


##


if __name__ == '__main__':
    main()


##
