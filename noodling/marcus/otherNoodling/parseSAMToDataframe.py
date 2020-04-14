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
                        default=None, help="Option to only read 'n' number of lines of file")
    parser.add_argument('-p', '--print_rows', metavar='print_rows', type=int,
                        default=None, help="Option to print 'n' number of lines of final dataframe")
    parser.add_argument('-m', '--deep_memory', action='store_true',
                        help="Boolean flag to print dataframe deep memory info\n"
                             "(this can be very CPU/time intensive)")
    
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


def parseSamToDataframe(filename, headerlines, num_lines=None, print_rows=None, deep_memory=False):
    
    # Quick check to ensure the passed file path exists
    if not os.path.isfile(filename):
        print(f"File does not exist at: {filename}, Terminating Script\n")
        sys.exit()
    else:
        print(f"Parsing identified file at: {filename}\n")
    
    # Assign known order of datatypes for .SAM file, this helps speed up parsing
    o = 'object'
    i = 'int64'
    # I can't find a clean way to do this quickly but not explicitly
    yeast_chr_datatype = pd.api.types.CategoricalDtype(categories=['I', 'II', 'III', 'IV', 'V',
                                                                   'VI', 'VII', 'VIII', 'IX', 'X',
                                                                   'XI', 'XII', 'XIII', 'XIV', 'XV',
                                                                   'XVI'], ordered=True)
    chr_cat = yeast_chr_datatype
    sam_dtypes_list = [o, i, chr_cat, i, i, o, o, i, i, o, o, o, o, o, o]
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
    
    SAM_df = SAM_df.sort_values(by=[2, 3], ignore_index=True)
    # This sort_values with the ignore_index option on will currently reset the indexes
    # to the new sorted order, I don't know if this is good or bad...
    
    print(f'Finished parsing and sorting of file at: {filename}\n')
    
    if print_rows:
        print(f"\nFirst {print_rows} rows of dataframe:\n", SAM_df.head(print_rows), "\n")
    if deep_memory:
        # Print dataframe info. This is currently really intensive with the deep memory usage call.
        # Remove the df.info print for speed as needed
        print(SAM_df.info(memory_usage='deep'))
    return SAM_df


if __name__ == '__main__':
    arguments = parseArgs()
    df = parseSamToDataframe(**arguments)
