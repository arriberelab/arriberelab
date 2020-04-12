"""
parseSAMToDataframe.py
Marcus Viscardi     April 12, 2020

Trying to parse a fastq file into a pandas dataframe, with the goal of being able to utilize read data faster

The goal will be to  parse SAM files into a dataframe, which would (maybe) allow for acceleration of the
    assignReadsToGenes functionality.
Going to initially run testing with data from:
    /data15/joshua/working/200331_yeastTestRun/processedData
"""

import sys, os
import argparse
import numpy as np
import pandas as pd
# Pandas default would cut off basically all the columns so:
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


def parseArgs():
    parser = argparse.ArgumentParser(description="Parse a SAM file into a pandas dataframe")
    parser.add_argument('filename', metavar='filename', type=str, help="Path to .sam file")
    parser.add_argument('headerlines', metavar='headerlines', type=int, help="Number of header lines before data")

    args = parser.parse_args()
    file = args.filename

    if not os.path.isfile(file):
        print("This is not a file path\n\n\nTERMINATING\n\n")
        sys.exit()
    else:
        print(f"Parsing file @: {file}")

    arg_dict = {arg: vars(args)[arg] for arg in vars(args)}
    print(arg_dict)
    return arg_dict


def parseSamToDataframe(arg_dict: dict):

    o = 'object'
    i = 'int64'
    sam_dtypes_list = [o, i, o, i, i, o, o, i, i, o, o, o, o, o, o]
    sam_dtypes_dict = {i: sam_dtypes_list[i] for i in range(15)}

    SAM_df = pd.read_csv(arg_dict['filename'],
                         sep="\t",
                         header=arg_dict['headerlines'],
                         names=range(15),
                         dtype=sam_dtypes_dict)  # The hope is this allows pandas to skip the type calling step
    print(SAM_df.head(10))
    print(SAM_df.info(memory_usage='deep'))
    return SAM_df


if __name__ == '__main__':
    arguments = parseArgs()
    df = parseSamToDataframe(arguments)