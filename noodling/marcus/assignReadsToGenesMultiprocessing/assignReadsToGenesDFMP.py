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
import numpy as np
import pandas as pd

# Pandas default would cut off any columns beyond 5 so:
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


def parseArgs():
    parser = argparse.ArgumentParser(description="Take reads.SAM file + *.allChrs.txt"
                                                 "file and assign genes to the reads")
    parser.add_argument('sam_file', metavar='sam_file',
                        type=str, help="Path to .sam file")
    parser.add_argument('annot_file', metavar='annot_file',
                        type=str, help="Path to .allChrs.txt file")
    parser.add_argument('output_prefix', metavar='output_prefix',
                        type=str, help="Prefix for output file names")
    parser.add_argument('-n', '--num_lines', metavar='num_lines', type=int,
                        default=None, help="Option to only read 'n' number of lines of each file,"
                                           "mostly if you're doing a quick test")
    parser.add_argument('-p', '--print_rows', metavar='print_rows', type=int,
                        default=None, help="Option to print 'n' number of lines of final dataframe,"
                                           "can get lengthy fast as each split chromosome will print this many lines")
    parser.add_argument('-m', '--deep_memory', action='store_true',
                        help="Boolean flag to print dataframe deep memory info\n"
                             "(this can be very CPU/time intensive, but informative)")
    parser.add_argument('-o', '--concatenate_output', action='store_true',
                        help="Boolean flag to output a single file for all chromosomes")
    parser.add_argument('-u', '--keep_non_unique', action='store_true',
                        help="Boolean flag to allow non-unique reads through the assignment process")
    
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


def recoverMappedPortion_dfWrapper(sam_df_dict, print_rows=None, **kwargs):
    for chr_key, df in sam_df_dict.items():
        try:
            sam_df_dict[chr_key][['map_read_seq', 'N']] = pd.DataFrame(df.apply(lambda x:
                                                                                recoverMappedPortion(x['cigar'],
                                                                                                     x['read_seq']),
                                                                                axis=1).tolist(), index=df.index)
            print(f'Recovery of mapped portion complete for Chr-{chr_key:->4}, '
                  f'read count={len(sam_df_dict[chr_key].index)}')
            if print_rows:
                print(sam_df_dict[chr_key][['read_id',
                                            'chr',
                                            'chr_pos',
                                            'cigar',
                                            'read_seq',
                                            'map_read_seq']].head(print_rows))
        except AttributeError:
            print(f'No reads for mapped portion recovery in Chr-{chr_key:->4}, '
                  f'read count={len(sam_df_dict[chr_key].index)}')
    return sam_df_dict


def assignReadsToGenes(sam_df_dict, annot_df_dict, print_rows=None, keep_non_unique=False, **kwargs):
    # Going for the df.merge() function for mapping annotations onto reads
    print(f"\nAnnotation alignment for {len(sam_df_dict.keys())} chromosomes:")
    for chr_key, df in sam_df_dict.items():
        try:
            # ONLY KEEP UNIQUELY MAPPING READS: >>
            if not keep_non_unique:
                df = sam_df_dict[chr_key][sam_df_dict[chr_key]['NH'].str.endswith(':1')]
            # << ONLY KEEP UNIQUELY MAPPING READS
            
            print(f"\tPreforming alignment for Chr-{chr_key:->4} containing {len(df.index)} reads")
            sam_df_dict[chr_key] = df.merge(annot_df_dict[chr_key], on=['chr_pos', 'chr'])
            print(f"\t\tSuccess, {len(sam_df_dict[chr_key][sam_df_dict[chr_key]['gene'] != np.nan].index):>7} "
                  f"reads assigned to genes in Chr-{chr_key:->4}\n")
            if print_rows:
                print(sam_df_dict[chr_key][['read_id',
                                            'chr',
                                            'chr_pos',
                                            'cigar',
                                            'read_seq',
                                            'gene',
                                            'gene_string']].head(print_rows))
            
            # TODO: some functionality to drop full df if nothing is assigned at all
        except KeyError as key:
            if str(key).strip("'") == chr_key:
                print(f"\t\tChr-{chr_key:->4} not found in annotations! -> Adding empty columns to compensate\n")
                sam_df_dict[chr_key]['gene'], sam_df_dict[chr_key]['gene_string'] = np.nan, np.nan
            else:
                print(f"\tOther KeyError pertaining to:", str(key), chr_key)
    return sam_df_dict


def fixSenseNonsense(sam_df_dict, print_rows=None, **kwargs):
    """Check for strand and sense/antisense, create rev-compliment as needed, and edits gene_string"""
    def revCompl(seq):
        """Will return the reverse complement of a sequence"""
        a = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
             'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
             'N': 'N', 'n': 'n', '.': '.', '-': '-'}
        return ''.join([a[seq[-i]] for i in range(1, len(seq) + 1)])
    
    def S_AS_flag_rework(S_AS, chr_pos, read_seq, map_read_seq, N, gene_string):
        if S_AS & 16 != 0:  # check strand
            chr_pos += len(map_read_seq) + N - 1
            read_strand = '-'
            map_read_seq = revCompl(map_read_seq)
            read_seq = revCompl(read_seq)
        else:
            read_strand = '+'
        # This has some issues if reads without annotations get through >>>
        if str(gene_string)[-1] == read_strand:
            gene_string = str(gene_string)[:-1] + 'S'
        else:
            gene_string = str(gene_string)[:-1] + 'AS'
        # <<< This has some issues if reads without annotations get through
        return read_strand, read_seq, map_read_seq, chr_pos, gene_string
    
    for chr_key, df in sam_df_dict.items():
        sam_df_dict[chr_key][['strand', 'read_seq', 'map_read_seq', 'chr_pos', 'gene_string']] = \
            pd.DataFrame(df.apply(lambda x: S_AS_flag_rework(x['strand'],
                                                             x['chr_pos'],
                                                             x['read_seq'],
                                                             x['map_read_seq'],
                                                             x['N'],
                                                             x['gene_string']),
                                  axis=1).tolist(), index=df.index)
    return sam_df_dict


def outputToCSV(jam_file ):
    pass


def main(sam_file, annot_file, output_prefix, print_rows=None,
         concatenate_output=False, keep_non_unique=False, **kwargs):
    sam_df_dict = parseSamToDF(sam_file, **kwargs)
    annot_df_dict = parseAllChrsToDF(annot_file, **kwargs)
    unassigned_df = pd.DataFrame()
    
    # Going for the df.merge() function for mapping annotations onto reads
    sam_df_dict = assignReadsToGenes(sam_df_dict, annot_df_dict,
                                     print_rows=print_rows, keep_non_unique=keep_non_unique, **kwargs)
    
    # Lets try to apply the recoverMappedPortion() to dataframe to see how it does
    #  Currently doing this after dropping unassigned reads as this is a more time intensive step.
    sam_df_dict = recoverMappedPortion_dfWrapper(sam_df_dict, print_rows=print_rows, **kwargs)
    
    jam_df_dict = fixSenseNonsense(sam_df_dict, print_rows=print_rows, **kwargs)
    
    if not concatenate_output:
        for chr_key, df in jam_df_dict.items():
            jam_df_dict[chr_key].sort_values(by=['read_id', 'chr', 'chr_pos'])
            jam_df_dict[chr_key].to_csv(f"{output_prefix}.chr{chr_key}.jam",
                                        index=False, sep='\t',
                                        columns=['read_id',
                                                 'chr',
                                                 'chr_pos',
                                                 'strand',
                                                 'mapq',
                                                 'cigar',
                                                 'map_read_seq',
                                                 'NH', 'HI',
                                                 'gene',
                                                 'gene_string'])
    else:
        jam_all_chrs = pd.concat(jam_df_dict.values(), ignore_index=True)
        jam_all_chrs.sort_values(by=['read_id', 'chr', 'chr_pos'])
        jam_all_chrs.to_csv(f"{output_prefix}.allChrs.jam",
                            index=False, sep='\t',
                            columns=['read_id',
                                     'chr',
                                     'chr_pos',
                                     'strand',
                                     'mapq',
                                     'cigar',
                                     'map_read_seq',
                                     'NH', 'HI',
                                     'gene',
                                     'gene_string'])
    
    print("\n\nDone?")


if __name__ == '__main__':
    arg_dict = parseArgs()
    main(**arg_dict)
