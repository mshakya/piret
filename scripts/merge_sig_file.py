#! /usr/bin/env python

"""
Script that joins all csv from EdgeR.

It creates one csv file for easier manual analysis.
"""

from __future__ import print_function

import argparse
import pandas as pd
import os


def cmdline_parser():
    """
    Create an argparse instance.

    for inputting options.
    """
    parser = argparse.ArgumentParser(description="""merges EdgeR comparison files
        to one csv file""")
    parser.add_argument("-i", "--INPUT_FOLDER",
                        help="folder that has the EdgeR results")
    return parser


def main():
    """
    Main function.

    All functions are here.
    """
    parser = cmdline_parser()
    args = parser.parse_args()

    gene_list_filename = []
    for dirname, dirnames, filenames in os.walk(args.INPUT_FOLDER):
        filenames.sort()
        for filename in filenames:
            if "__et.csv" in filename:
                gene_list_filename.append(filename)
    df = pd.read_csv(gene_list_filename[0])
    df.columns = ["Geneid1", "Geneid", "Chr", "Start", "End", "Strand",
                  "Length", "logFC" + gene_list_filename[0].split("__et.csv")[0],
                  "logCPM" + gene_list_filename[0].split("__et.csv")[0],
                  "PValue" + gene_list_filename[0].split("__et.csv")[0],
                  "FDR" + gene_list_filename[0].split("__et.csv")[0]]
    gene_list_filename.pop(0)
    for table in gene_list_filename:
        df1 = pd.read_csv(table)
        df1.columns = ["Geneid1", "Geneid", "Chr", "Start", "End", "Strand",
                       "Length", "logFC" + str(table).split("__et.csv")[0],
                       "logCPM" + str(table).split("__et.csv")[0],
                       "PValue" + str(table).split("__et.csv")[0],
                       "FDR" + str(table).split("__et.csv")[0]]
        df = pd.merge(df, df1, how='left', on=["Geneid1", "Geneid", "Chr",
                                               "Start", "End", "Strand",
                                               "Length"])
    output_file = '/'.join([args.INPUT_FOLDER, "all_fold_changes.csv"])
    df.to_csv(output_file, index=False)


if __name__ == '__main__':
    main()
