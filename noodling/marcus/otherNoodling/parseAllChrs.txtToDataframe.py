"""
parseAllChrs.txtToDataframe.py
Marcus Viscardi     April 13, 2020

Josh already has a fully functional version of this parser within assignReadsToGenes4.py

But I thought it would be good practice to write my own(?)

Either way the real goal here is to parse the *allChrs.txt annotation file into a pandas dataframe, from there we can
    use pandas to split the dataframe up into individual chromosomes. The advantage here will be to do the same with
    the reads, allowing us to parallelize each chromosome. Further splitting could be done within each chr. This avoids
    my major worry of trying to avoid memory conflicts while setting up multiprocessing.
"""

import argparse
import os
import sys
import pandas as pd
# Pandas default would cut off any columns beyond 5 so:
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


def parseArgs():  # Adapted from parseFastqToDataframe.py
    parser = argparse.ArgumentParser(description="Parse a *.allChrs.txt annotations file into a pandas dataframe")
    parser.add_argument('filename', metavar='filename', type=str,
                        default=None, help="Path to *.allChrs.txt file")
    parser.add_argument('-s', '--split_chrs', metavar='split_chrs', type=bool,
                        default=False, help="Option to output several dataframes, one for each chromosome")
    parser.add_argument('-n', '--num_lines', metavar='num_lines', type=int,
                        default=None, help="Option to only read N number of lines of file")

    args = parser.parse_args()

    # Quickly convert Namespace object to dictionary
    arg_dict = {arg: vars(args)[arg] for arg in vars(args)}

    # Print Given arguments, currently for debugging
    print("\nGiven Arguments (ArgParse):")
    for key, arg in arg_dict.items():
        if not arg:
            print(f"\t{key} = {arg} -> (Will not be passed)")
        else:
            print(f"\t{key} = {arg}")
    print("\tDone.\n")

    # Recreate dict without arguments that did not receive any input, print and return final dict
    arg_dict = {k: v for k, v in arg_dict.items() if v}
    return arg_dict


def parseAllChrsToDataframe(filename, num_lines=None, split_chrs=False):
    if num_lines:
        annot_df = pd.read_csv(filename,
                               # index_col=0,
                               delimiter='|',
                               names=['chr_index', 'genes'],
                               nrows=num_lines,
                               )
    else:
        annot_df = pd.read_csv(filename,
                               # index_col=0,
                               delimiter='|',
                               names=['chr_index', 'genes'],
                               )
    # print(annot_df)
    annot_df[['chr', 'chr_pos']] = pd.DataFrame(annot_df['chr_index'].str.split('_').values.tolist(),
                                                index=annot_df.index)
    ################
    print(annot_df)


if __name__ == '__main__':
    arg_dict = parseArgs()
    parseAllChrsToDataframe(**arg_dict)
