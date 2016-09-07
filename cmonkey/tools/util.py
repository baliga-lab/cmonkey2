"""util,py - reusable functionality"""
import os
import pandas


def read_ratios(result_dir):
    csvpath = os.path.join(result_dir, 'ratios.tsv.gz')
    df = pandas.read_csv(csvpath, index_col=0, sep='\t')
    df.index = map(str, df.index)
    return df

