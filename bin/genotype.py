"""
Genotype annotated alleles in individual samples.

For each MUT in the confident callset, scan the AF distribution from the
annotated table for a relative gap (>= gap_th). If a gap is found, keep
annotated alleles above it (also filtered by total_count and sb_pval);
otherwise, if the mean AF fold-change between consensus and annotated is
larger than gap_th, keep only the consensus alleles; otherwise discard
the site.
"""

import argparse
import pandas as pd


##


def parse_args():

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--annotated_allelic_table', required=True,
                   help='Annotated per-donor allelic table (tsv.gz)')
    p.add_argument('--confident_allelic_table', required=True,
                   help='Confident per-donor allelic table (tsv.gz)')
    p.add_argument('--donor',           required=True, help='Donor identifier')
    p.add_argument('--gap_th',          type=float, default=0.5,
                   help='Relative AF gap threshold (also used as min fold-change)')
    p.add_argument('--coverage_th',     type=float, default=10,
                   help='Minimum callable_coverage on annotated table')
    p.add_argument('--min_total_counts', type=int,   default=3,
                   help='Minimum total_count when keeping annotated alleles above the gap')
    p.add_argument('--max_sb_pval',     type=float, default=0.05,
                   help='Minimum sb_pval when keeping annotated alleles above the gap')

    return p.parse_args()


##


def main():

    args = parse_args()

    df     = pd.read_csv(args.confident_allelic_table, sep='\t')
    df_all = pd.read_csv(args.annotated_allelic_table, sep='\t')

    CALLSET = df['MUT'].unique()
    df_all  = df_all.query('MUT in @CALLSET and callable_coverage>=@args.coverage_th')

    ALLELES = []
    for MUT in CALLSET:

        df_mut     = df.query('MUT==@MUT')
        df_all_mut = df_all.query('MUT==@MUT')

        if df_all_mut.empty:
            continue

        mean_consensus = df_mut['AF'].mean()
        mean_all       = df_all_mut['AF'].mean()
        mean_fc        = (mean_consensus - mean_all) / mean_all if mean_all > 0 else 0

        # Find first relative gap in AF distribution
        t = 0
        found_gap = False
        x = df_all_mut['AF'].sort_values()
        for i in range(len(x) - 1):
            x_i    = x.iloc[i]
            x_next = x.iloc[i + 1]
            if x_i <= 0:
                continue
            delta = (x_next - x_i) / x_i
            if delta > args.gap_th:
                found_gap = True
                t = x_next - ((x_next - x_i) / 2)
                break

        if found_gap:
            ALLELES.append(
                df_all_mut.loc[
                    (df_all_mut['AF'] >= t) &
                    (df_all_mut['total_count'] >= args.min_total_counts) &
                    (df_all_mut['sb_pval'] >= args.max_sb_pval)
                ].copy()
            )
        else:
            if mean_fc > args.gap_th:
                ALLELES.append(df_mut.copy())
            # else: discard site

    if ALLELES:
        df_genotyped = pd.concat(ALLELES, ignore_index=True)
        df_genotyped['is_consensus'] = (
            (df_genotyped['gs>1_count'] > 0) & (df_genotyped['overlap_count'] > 0)
        )
    else:
        df_genotyped = df.iloc[0:0].copy()
        df_genotyped['is_consensus'] = pd.Series(dtype=bool)

    df_genotyped.to_csv(
        f'{args.donor}.genotyped_allelic_table.tsv.gz',
        sep='\t', index=False
    )


##


if __name__ == '__main__':
    main()
