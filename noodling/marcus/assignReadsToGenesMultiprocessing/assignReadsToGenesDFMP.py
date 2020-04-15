"""
assignReadsToGenesDFMP.py
Marcus Viscardi     April 14, 2020

Basically just rewriting Josh's assignReadsToGenes4.py, with the major change of utilizing pandas dataframes and 
    multiprocessing in an effort to speed up read assignment.
"""

import sys, os
import argparse
import multiprocessing
import parseAllChrstxtToDataframe
import parseSAMToDataframe
import pandas as pd
# Pandas default would cut off any columns beyond 5 so:
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


def parseArgs():
    parser = argparse.ArgumentParser(description="Take reads.SAM file + *.allChrs.txt"
                                                 "file and assign genes to the reads")
    parser.add_argument('sam_file', metavar='sam_file',
                        type=str, help="Path to .sam file")
    parser.add_argument('annot_file', metavar='annot_file',
                        type=str, help="Path to .allChrs.txt file")
    parser.add_argument('-n', '--num_lines', metavar='num_lines', type=int,
                        default=None, help="Option to only read 'n' number of lines of each file,"
                                           "mostly if you're doing a quick test")
    parser.add_argument('-p', '--print_rows', metavar='print_rows', type=int,
                        default=None, help="Option to print 'n' number of lines of final dataframe,"
                                           "can get lengthy fast as each split chromosome will print this many lines")
    parser.add_argument('-m', '--deep_memory', action='store_true',
                        help="Boolean flag to print dataframe deep memory info\n"
                             "(this can be very CPU/time intensive, but informative)")
    
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
    arg_dict = {k: v for k, v in arg_dict.items() if v is not None}
    return arg_dict


def parseSamToDF(sam_file, num_lines=None, print_rows=None, deep_memory=False, **kwargs):
    return parseSAMToDataframe.parseSamToDataframe(sam_file, 20, num_lines=num_lines, split_chrs=True,
                                                   print_rows=print_rows, deep_memory=deep_memory)


def parseAllChrsToDF(annot_file, num_lines=None, print_rows=None, deep_memory=False, **kwargs):
    return parseAllChrstxtToDataframe.parseAllChrsToDataframe(annot_file, num_lines=num_lines, split_chrs=True,
                                                              print_rows=print_rows, deep_memory=deep_memory)


if __name__ == '__main__':
    arg_dict = parseArgs()
    SAM_df = parseSamToDF(**arg_dict)
    annot_df = parseAllChrsToDF(**arg_dict)
