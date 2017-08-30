#! /usr/bin/env python

import argparse
import pandas as pd


def cmdline_parser():
    """
    creates an argparse instance
    """
    parser = argparse.ArgumentParser(description="""program that generates histogram
        from contig length distribution""")
    parser.add_argument("-i", "--INPUT_FOLDER",
                        help="folder that has the EdgeR results")
    return parser



def main():
    """
    """
    parser = cmdline_parser()
    args = parser.parse_args()

    in_table = pd.read_csv(args.INPUT_TABLE, sep="\t")
    length_table = in_table[['Length']]
    plt.figure()
    length_table.plot.hist(alpha=0.5)
    plt.yscale('log', nonposy='clip')

    plt.show()
    plt.savefig(args.OUT_FIGURE)

if __name__ == '__main__':
    main()

