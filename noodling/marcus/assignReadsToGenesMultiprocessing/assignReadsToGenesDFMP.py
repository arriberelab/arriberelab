"""
assignReadsToGenesDFMP.py
Marcus Viscardi     April 14, 2020

Basically just rewriting Josh's assignReadsToGenes4.py, with the major change of utilizing pandas dataframes and 
    multiprocessing in an effort to speed up read assignment.
"""

import sys, os
import argparse
import multiprocessing
import re
import parseAllChrstxtToDataframe
import parseSAMToDataframe
import pandas as pd
# Pandas default would cut off any columns beyond 5 so:
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)

# TODO: name read.SAM df columns for ease of writing code, as well as reading code. Currently column indexes just make
#       the applied df functions harder to follow. How far down this path do we want to go? Name all columns? Relevant
#       ones only? Just add comments to clarify?


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


def recoverMappedPortion(Cigar, Read):
    # April 15, 2020: Stolen verbatim (+ my comments) from assignReadsToGenes4.py
    
    """Given a Cigar string and a Read, will return the sequence of the read that mapped to the genome."""
    # Edit Oct 10, 2013 to include skipped portions of reference sequence (introns)
    
    # first process the CIGAR string
    cigarSplit = re.findall('(\d+|[a-zA-Z]+)', Cigar)  # RE func: matches sets of digits or letters
    cigarSplit = [[int(cigarSplit[ii]), cigarSplit[ii + 1]] for ii in range(0, len(cigarSplit), 2)]
    
    # Then use that information to parse out nts of the read sequence
    mappedRead = ''
    ii = 0
    N = 0
    for entry in cigarSplit:
        if entry[1] in ['M', 'I']:  # then it's either aligned to the genomic sequence or has an insert relative to it
            mappedRead += Read[ii:ii + entry[0]]
            ii += entry[0]
        elif entry[1] == 'S':
            ii += entry[0]
        elif entry[1] == 'N':
            N += entry[0]
            # N is used for "skipped region from the reference". I keep track of Ns and
            #  return them for calculation of position on the - strand
        # else:
        #     print(entry[1], end='')
    return mappedRead, N


def assignReadsToGenes(sam_file, annot_file, **kwargs):
    sam_df_dict = parseSamToDF(sam_file, **kwargs)
    annot_df_dict = parseAllChrsToDF(annot_file, **kwargs)
    print()
    
    # Going for the df.merge() function for mapping annotations onto reads
    for chr, df in sam_df_dict.items():
        try:
            print(f"\nPreforming alignment for Chr-{chr} containing {len(df.index)} reads")
            sam_df_dict[chr] = df.merge(annot_df_dict[chr], left_on=3, right_on='chr_pos')
            print(sam_df_dict[chr][[0, 2, 'chr', 3, 'chr_pos', 5, 9, 'gene', 'gene_string']].head(5))
            # print(sam_df_dict[chr].head(5))
        except KeyError as key:
            print(f"Chr-{chr:<4} not found in annotations!\n {key}")
    
    # Use df.dropna() to remove rows which did not map in previous step
    print("\nCleaning unmapped rows...\n")
    for chr, df in sam_df_dict.items():
        sam_df_dict[chr] = df.dropna(axis=0, how='any')
        print(f'Chr-{chr:<4}, Rows:{len(sam_df_dict[chr].index):>5}, Columns:{len(sam_df_dict[chr].columns):>3}')
    
    # Lets try to apply the recoverMappedPortion() to dataframe to see how it does
    for chr, df in sam_df_dict.items():
        try:
            sam_df_dict[chr][[15, 16]] = pd.DataFrame(df.apply(lambda x: recoverMappedPortion(x[5], x[9]),
                                                 axis=1).tolist(), index=df.index)
            print(f'Recovery of mapped portion complete for Chr-{chr:->4}, read count={len(sam_df_dict[chr].index)}')
        except AttributeError:
            print(f'No reads for mapped portion recovery in Chr-{chr:->4}, read count={len(sam_df_dict[chr].index)}')
    
    print("\n\nDone?")


if __name__ == '__main__':
    arg_dict = parseArgs()
    assignReadsToGenes(**arg_dict)
