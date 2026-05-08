#!/usr/bin/env python3

"""
Make per-colony allelic and coverage tables from a mitobam.

Fragments are grouped into eUMI groups defined by
(fwd_start, frag_end, fwd_cigar, rev_cigar). All output counts are
in molecules (eUMI groups), not reads. PCR copies are collapsed by
consensus before any count is incremented.

Consensus rules (per position, per group of N pairs):
  N == 1 : non-overlap always passes; overlap requires 100% R1/R2 agreement
  N == 2 : 100% of all votes must agree
  N  > 2 : strict majority > 50% of votes

Output
------
{sample}.allelic_table.tsv.gz
    One row per (POS, ALT) where ALT != REF and total_count > 0.
    Columns:
      sample, POS, REF, ALT,
      total_count, fw_count, rev_count, overlap_count,
      strand_ratio, sb_pval, pos_end_read,
      mean_group_size, gs>1_count

{sample}.coverage_table.tsv.gz
    One row per position where total_coverage > 0.
    Columns:
      sample, POS, REF,
      total_coverage, overlap_coverage, callable_coverage

Column definitions
------------------
  total_count    = fw_count + rev_count + overlap_count  (molecules)
  fw_count       : molecules where POS is fwd-only AND consensus = ALT
  rev_count      : molecules where POS is rev-only AND consensus = ALT
  overlap_count  : molecules where POS is in R1/R2 overlap AND consensus = ALT
  strand_ratio   : fw_count / (rev_count + 1e-18)
  sb_pval        : binomial test on (fw_count, fw+rev, p=0.5); NaN if fw+rev=0
  pos_end_read   : minimum distance from nearest read end across all reads
                   that voted for the consensus ALT base
  mean_group_size: mean eUMI group size across molecules calling this ALT
  gs>1_count     : molecules with group_size > 1 (PCR-confirmed calls)

  total_coverage   : molecules aligning to this position (no BQ filter)
  overlap_coverage : molecules where this position is in R1/R2 overlap
  callable_coverage: molecules where consensus was successfully called (REF or ALT)
"""

import gzip
import argparse
import numpy as np
import pysam
from collections import defaultdict
from scipy.stats import binomtest


##


my_parser = argparse.ArgumentParser(
    prog='make_strand_tables',
    description='Build per-colony allelic and coverage tables from a mitobam.'
)
my_parser.add_argument('--bam',    type=str, required=True,  help='Path to mitobam.bam')
my_parser.add_argument('--sample', type=str, required=True,  help='Colony / sample name')
my_parser.add_argument('--ref',    type=str, required=True,  help='Path to faidx-indexed reference FASTA')
my_parser.add_argument('--region', type=str, default='chrM', help='MT contig name (default: chrM)')
my_parser.add_argument('--min_mq', type=int, default=30,     help='Min mapping quality (default: 30)')
my_parser.add_argument('--min_bq', type=int, default=30,     help='Min base quality (default: 30)')
my_parser.add_argument('--max_softclip', type=float, default=0.1,
                        help='Max fraction of softclipped bases per read (default: 0.1)')

args          = my_parser.parse_args()
bam_path      = args.bam
sample        = args.sample
ref_path      = args.ref
region        = args.region
min_mq        = args.min_mq
min_bq        = args.min_bq
max_softclip  = args.max_softclip

BASES = ['A', 'C', 'G', 'T']
maxBP = 16569
MIN_DIST_FROM_END = 5


##


def get_ref_dict(ref_path, region):
    fasta = pysam.FastaFile(ref_path)
    seq   = fasta.fetch(region).upper()
    fasta.close()
    return {i + 1: base for i, base in enumerate(seq)}   # 1-based


def dist_from_end(qpos, read_length):
    """Minimum distance of a query position from either end of the read."""
    return min(qpos, read_length - 1 - qpos)


##


def main():

    ref_dict = get_ref_dict(ref_path, region)

    # allelic[pos1][base] = {fw, rev, overlap, pos_end, group_sizes}
    allelic = defaultdict(lambda: {
        b: {'fw': 0, 'rev': 0, 'overlap': 0, 'pos_end': [], 'group_sizes': []}
        for b in BASES
    })
    # coverage[pos1] = {total, overlap, callable}
    coverage = defaultdict(lambda: {'total': 0, 'overlap': 0, 'callable': 0})

    # ------------------------------------------------------------------ #
    # Phase 1: load all reads grouped by read name
    # ------------------------------------------------------------------ #
    bam = pysam.AlignmentFile(bam_path, 'rb', require_index=False)
    read_pairs = defaultdict(list)
    for read in bam:
        if (read.is_unmapped        or
            read.is_supplementary   or
            read.is_secondary       or
            not read.is_proper_pair or
            read.mapping_quality < min_mq):
            continue
        qlen = read.query_length
        if qlen:
            sc = sum(length for op, length in (read.cigartuples or []) if op == 4)
            if sc / qlen > max_softclip:
                continue
        read_pairs[read.query_name].append(read)
    bam.close()

    # ------------------------------------------------------------------ #
    # Phase 2: group pairs into eUMI groups
    # eUMI = (fwd_start, frag_end, fwd_cigar, rev_cigar)
    # ------------------------------------------------------------------ #
    eumi_groups = defaultdict(list)
    for _, reads in read_pairs.items():
        if len(reads) != 2:
            continue
        r0, r1 = reads
        if not r0.is_reverse and r1.is_reverse:
            fwd, rev = r0, r1
        elif r0.is_reverse and not r1.is_reverse:
            fwd, rev = r1, r0
        else:
            continue
        if fwd.cigarstring is None or rev.cigarstring is None:
            continue
        tlen = abs(fwd.template_length)
        if tlen == 0:
            continue
        eumi = (fwd.reference_start, fwd.reference_start + tlen,
                fwd.cigarstring, rev.cigarstring)
        eumi_groups[eumi].append((fwd, rev))

    # ------------------------------------------------------------------ #
    # Phase 3: process each eUMI group
    # ------------------------------------------------------------------ #
    for eumi, pairs in eumi_groups.items():

        group_size = len(pairs)

        # Alignment structure derived from the first pair.
        # All pairs share the same CIGAR so the reference-position →
        # query-position mapping is identical across pairs.
        fwd0, rev0 = pairs[0]
        fwd_map0 = {rp: qp for qp, rp in fwd0.get_aligned_pairs(matches_only=True)
                    if rp is not None and rp < maxBP}
        rev_map0 = {rp: qp for qp, rp in rev0.get_aligned_pairs(matches_only=True)
                    if rp is not None and rp < maxBP}

        overlap_pos  = set(fwd_map0) & set(rev_map0)
        fwd_only_pos = set(fwd_map0) - overlap_pos
        rev_only_pos = set(rev_map0) - overlap_pos
        all_pos      = set(fwd_map0) | set(rev_map0)

        # Coverage: raw molecule depth — no BQ filter
        for rp in all_pos:
            coverage[rp + 1]['total'] += 1
        for rp in overlap_pos:
            coverage[rp + 1]['overlap'] += 1

        # Accumulate BQ-filtered base votes from every read in every pair.
        # group_votes[rp][base] = {'count': N, 'pos_end': [...]}
        # Query positions reused from fwd_map0/rev_map0 (same CIGAR).
        group_votes = defaultdict(
            lambda: defaultdict(lambda: {'count': 0, 'pos_end': []})
        )

        for fwd, rev in pairs:
            fwd_seq = fwd.query_sequence
            fwd_bq  = fwd.query_qualities
            rev_seq = rev.query_sequence
            rev_bq  = rev.query_qualities
            if fwd_seq is None or fwd_bq is None or rev_seq is None or rev_bq is None:
                continue
            fwd_len = fwd.query_length
            rev_len = rev.query_length

            # fwd read: covers fwd_only + overlap positions
            for rp, qp in fwd_map0.items():
                q = fwd_bq[qp]
                if q < min_bq:
                    continue
                b = fwd_seq[qp]
                if b not in BASES:
                    continue
                d = dist_from_end(qp, fwd_len)
                if d < MIN_DIST_FROM_END:
                    continue
                group_votes[rp][b]['count'] += 1
                group_votes[rp][b]['pos_end'].append(d)

            # rev read: covers rev_only + overlap positions
            for rp, qp in rev_map0.items():
                q = rev_bq[qp]
                if q < min_bq:
                    continue
                b = rev_seq[qp]
                if b not in BASES:
                    continue
                d = dist_from_end(qp, rev_len)
                if d < MIN_DIST_FROM_END:
                    continue
                group_votes[rp][b]['count'] += 1
                group_votes[rp][b]['pos_end'].append(d)

        # Consensus pass: one molecule count per position if threshold met
        for rp, base_votes in group_votes.items():
            total_votes = sum(v['count'] for v in base_votes.values())
            if total_votes == 0:
                continue

            consensus_base = max(base_votes, key=lambda b: base_votes[b]['count'])
            css  = base_votes[consensus_base]['count'] / total_votes
            pos1 = rp + 1

            # Apply group-size-dependent consensus threshold
            if group_size == 1:
                # overlap requires both R1 and R2 to agree (css == 1.0)
                if rp in overlap_pos and css < 1.0:
                    continue
            elif group_size == 2:
                if css < 1.0:
                    continue
            else:   # group_size > 2
                if css <= 0.5:
                    continue

            # Callable coverage: consensus met for REF or ALT
            coverage[pos1]['callable'] += 1

            # Strand bucket
            if rp in fwd_only_pos:
                bucket = 'fw'
            elif rp in rev_only_pos:
                bucket = 'rev'
            else:
                bucket = 'overlap'

            # Allelic table: only non-REF calls
            ref_base = ref_dict.get(pos1, 'N')
            if consensus_base != ref_base:
                d = allelic[pos1][consensus_base]
                d[bucket] += 1
                d['pos_end'].extend(base_votes[consensus_base]['pos_end'])
                d['group_sizes'].append(group_size)

    # ------------------------------------------------------------------ #
    # Write allelic table
    # ------------------------------------------------------------------ #
    allelic_header = [
        'sample', 'POS', 'REF', 'ALT',
        'total_count', 'fw_count', 'rev_count', 'overlap_count',
        'strand_ratio', 'sb_pval', 'pos_end_read',
        'mean_group_size', 'gs>1_count'
    ]

    with gzip.open(f'{sample}.allelic_table.tsv.gz', 'wt') as fh:
        fh.write('\t'.join(allelic_header) + '\n')

        for pos1 in sorted(allelic):
            ref_base = ref_dict.get(pos1, 'N')
            for alt in BASES:
                if alt == ref_base:
                    continue
                d     = allelic[pos1][alt]
                fw    = d['fw']
                rv    = d['rev']
                ov    = d['overlap']
                total = fw + rv + ov
                if total == 0:
                    continue

                pos_end = int(np.median(d['pos_end']))
                sr      = fw / (rv + 1e-18)
                ss_n    = fw + rv
                sb_pval = round(
                    binomtest(fw, ss_n, 0.5, alternative='two-sided').pvalue, 6
                ) if ss_n > 0 else 'nan'

                gs         = d['group_sizes']
                mean_gs    = round(float(np.mean(gs)), 2)
                gs1p_count = sum(1 for x in gs if x > 1)

                fh.write('\t'.join(map(str, [
                    sample, pos1, ref_base, alt,
                    total, fw, rv, ov,
                    round(sr, 4), sb_pval, pos_end,
                    mean_gs, gs1p_count
                ])) + '\n')

    # ------------------------------------------------------------------ #
    # Write coverage table
    # ------------------------------------------------------------------ #
    coverage_header = [
        'sample', 'POS', 'REF',
        'total_coverage', 'overlap_coverage', 'callable_coverage'
    ]

    with gzip.open(f'{sample}.coverage_table.tsv.gz', 'wt') as fh:
        fh.write('\t'.join(coverage_header) + '\n')

        for pos1 in sorted(coverage):
            ref_base = ref_dict.get(pos1, 'N')
            cov      = coverage[pos1]
            if cov['total'] == 0:
                continue
            fh.write('\t'.join(map(str, [
                sample, pos1, ref_base,
                cov['total'], cov['overlap'], cov['callable']
            ])) + '\n')


##


if __name__ == '__main__':
    main()
