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
import pandas as pd
# Pandas default would cut off basically all the columns so:
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


def parseArgs():
    parser = argparse.ArgumentParser(description="Parse a SAM file into a pandas dataframe")
    parser.add_argument('filename', metavar='filename',
                        type=str, help="Path to .sam file")
    parser.add_argument('headerlines', metavar='headerlines',
                        type=int, help="Number of header lines before data")
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
    arg_dict = {k: v for k, v in arg_dict.items() if v is not None}
    return arg_dict


def parseSamToDataframe(filename, headerlines, num_lines=None):
    
    # Quick check to ensure the passed file path exists
    if not os.path.isfile(filename):
        print(f"File does not exist at: {filename}, Terminating Script\n")
        sys.exit()
    else:
        print(f"Parsing identified file at: {filename}\n")
    
    # Assign known order of datatypes for .SAM file, this helps speed up parsing
    o = 'object'
    i = 'int64'
    sam_dtypes_list = [o, i, o, i, i, o, o, i, i, o, o, o, o, o, o]
    sam_dtypes_dict = {i: sam_dtypes_list[i] for i in range(15)}
    
    if num_lines:
        SAM_df = pd.read_csv(filename,
                             sep="\t",
                             header=headerlines,
                             nrows=num_lines,
                             names=range(15),  # Just name columns by an index, maybe better naming in the future?
                             dtype=sam_dtypes_dict,  # The hope is this allows pandas to skip the type calling step
                             )
    else:
        SAM_df = pd.read_csv(filename,
                             sep="\t",
                             header=headerlines,
                             names=range(15),
                             dtype=sam_dtypes_dict,
                             )
    
    print("\nFirst 5 rows of dataframe:\n", SAM_df.head(5), "\n")
    # Print dataframe info. This is currently really intensive with the deep memory usage call.
    # Remove the df.info print for speed as needed
    # print(SAM_df.info(memory_usage='deep'))
    return SAM_df


if __name__ == '__main__':
    arguments = parseArgs()
    df = parseSamToDataframe(**arguments)
