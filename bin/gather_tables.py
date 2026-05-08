#!/usr/bin/env python3

"""
Gather per-colony allelic and coverage tables into individual-level TSVs.

Expects all {colony}.allelic_table.tsv.gz and {colony}.coverage_table.tsv.gz
files to be present in the working directory (staged by Nextflow).

Output
------
{individual}.allelic_table.tsv.gz
{individual}.coverage_table.tsv.gz
"""

import glob
import argparse
import pandas as pd


##


my_parser = argparse.ArgumentParser(
    prog='gather_strand_concordance',
    description='Concatenate per-colony allelic and coverage tables.'
)
my_parser.add_argument('--individual', type=str, default='unknown',
                        help='Donor / individual identifier')

args       = my_parser.parse_args()
individual = args.individual


##


def concat_tables(pattern, out_path):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise RuntimeError(f'No files found matching: {pattern}')
    df = pd.concat(
        [pd.read_csv(p, sep='\t') for p in paths],
        ignore_index=True
    )
    df.to_csv(out_path, sep='\t', index=False)
    print(f'Written: {out_path}  ({len(df):,} rows, {df["sample"].nunique()} colonies)')


def concat_mt_cn(pattern, out_path):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise RuntimeError(f'No files found matching: {pattern}')
    df = pd.concat(
        [pd.read_csv(p, sep='\t', header=None,
                     names=['colony_name', 'nuclear_coverage', 'mt_coverage', 'mt_cn'])
         for p in paths],
        ignore_index=True
    )
    df.to_csv(out_path, sep='\t', index=False)
    print(f'Written: {out_path}  ({len(df):,} colonies)')


##


if __name__ == '__main__':
    concat_tables('*.allelic_table.tsv.gz',  f'{individual}.allelic_table.tsv.gz')
    concat_tables('*.coverage_table.tsv.gz', f'{individual}.coverage_table.tsv.gz')
    concat_mt_cn('*.mt_cn.txt',              f'{individual}.mt_cn.tsv.gz')
