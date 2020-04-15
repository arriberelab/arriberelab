"""
parseAllChrstxtToDataframe.py
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
    parser.add_argument('-s', '--split_chrs', action='store_true',
                        help="Boolean flag to split chromosomes and return a dict of dataframes")
    parser.add_argument('-n', '--num_lines', metavar='num_lines', type=int,
                        default=None, help="Option to only read N number of lines of file")
    parser.add_argument('-p', '--print_rows', metavar='print_rows', type=int,
                        default=None, help="Option to print 'n' number of lines of final dataframe")

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


def parseAllChrsToDataframe(filename, num_lines=None, split_chrs=False, print_rows=None, deep_memory=None):
    # Quick check to ensure the passed file path exists
    if not os.path.isfile(filename):
        print(f"\nFile does not exist at: {filename}, Terminating Script\n")
        sys.exit()
    else:
        print(f"\nParsing identified file at: {filename}")
    
    if num_lines:
        annot_df = pd.read_csv(filename,
                               # index_col=0,
                               delimiter='|',
                               names=['chr', 'genes'],
                               nrows=num_lines,
                               )
    else:
        annot_df = pd.read_csv(filename,
                               # index_col=0,
                               delimiter='|',
                               names=['chr', 'genes'],
                               )
    # Split the chr_index column into two
    annot_df[['chr', 'chr_pos']] = pd.DataFrame(annot_df['chr'].str.split('_').values.tolist(),
                                                index=annot_df.index)
    # Sort by Chromosome and index on chromosome
    annot_df = annot_df.sort_values(by=['chr', 'chr_pos'])
    # TODO: Reorganize DF columns
    # TODO: Add a lot of the functionality of the parseSAMToDataframe
    
    print(f'Finished parsing and sorting of file at: {filename}\n')
    
    if print_rows:
        print(f"\nFirst {print_rows} rows of dataframe:\n", annot_df.head(print_rows), '\n')
    if deep_memory:
        print(annot_df.info(memory_usage='deep'))
        
    if split_chrs:
        # Right now this is the only hardcoded piece which will cause major issues, the column containing the chrs
        #  has to be explicitly passed into this. Besides that it is extremely reusable
        annot_df_dict = dict(tuple(annot_df.groupby('chr')))
        if print_rows:
            for chr, df in annot_df_dict.items():
                print(f"\nFirst {print_rows} rows of Chromosome-{chr}:\n", df.head(print_rows))
        print(f"\nSplit dataframe into {len(annot_df_dict.keys())} separate df's:\n{annot_df_dict.keys()}")
        return annot_df_dict
    else:
        return annot_df


if __name__ == '__main__':
    arg_dict = parseArgs()
    parseAllChrsToDataframe(**arg_dict)
