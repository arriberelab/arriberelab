"""
assignReadsToGenesDFMP.py
Marcus Viscardi     Adapted from assignReadsToGenesDF.py on April 18, 2020

Basically just rewriting Josh's assignReadsToGenes4.py, with the change of
    utilizing pandas dataframes in an effort to speed up read assignment.
Currently full script is self contained. Pieces like the parsers and the
    revCompliment functions were written into here to make it a little more
    portable.
Addition of multiprocessing rather than iterating though dictionary of
    chr:chr_df sets. The hope is that this will be able to speed up the
    processing of larger read sets.
"""

from os import path, getpid
from sys import exit
from argparse import ArgumentParser
from re import findall
from typing import Dict, Any

import multiprocessing
from numpy import nan as numpy_nan
from pandas import set_option, read_csv, DataFrame, concat
# Pandas default would cut off any columns beyond 5 so:
set_option('display.max_rows', 50)
set_option('display.max_columns', 20)
set_option('display.width', 300)

# Some type hints help IDEs like PyCharm read the code a bit better:
CHR_DF_DICT = Dict[str, DataFrame]
ARG_DICT = Dict[str, Any]


def parseArgs() -> ARG_DICT:
    """
    parseArgs
    
    A tool to parse  CLI (Command Line Interface) arguments. Major advantage is that this
        allows users to type -h or --help, which will show all of the required and possible
        arguments along with the help statements written below.
    
    Output -> A dictionary of the argument key (ie. 'sam_file') and the passed value (ie.
        '/my/path/to/a/sam/file.sam'). Currently this system is set up to not pass any
        values which are not assigned by the user in the CLI call. This provides the advantage
        that another dictionary  of key:value defaults can be used to update the argument
        dictionary with any defaults that the script-writers would want.
    """
    
    # Start a argument parser instance which will accept CLI input
    parser = ArgumentParser(description="Take reads.SAM file + *.allChrs.txt"
                                        "file and assign genes to the reads,"
                                        "output *.jam file")
    # Required arguments:
    parser.add_argument('sam_file', metavar='sam_file',
                        type=str, help="Path to .sam file")
    parser.add_argument('annot_file', metavar='annot_file',
                        type=str, help="Path to .allChrs.txt file")
    parser.add_argument('output_prefix', metavar='output_prefix',
                        type=str, help="Prefix for output file names")
    # Optional arguments accepting integer values:
    parser.add_argument('-n', '--num_lines', metavar='num_lines', type=int,
                        default=None, help="Option to only read 'n' number of lines of each file,"
                                           "mostly if you're doing a quick test")
    parser.add_argument('-p', '--print_rows', metavar='print_rows', type=int,
                        default=None, help="Option to print 'n' number of lines of final dataframe,"
                                           "can get lengthy fast as each split chromosome will print this many lines")
    # Flag arguments which turn on functionality if passed:
    parser.add_argument('-m', '--deep_memory', action='store_true',
                        help="Boolean flag to print dataframe deep memory info\n"
                             "(this can be very CPU/time intensive, but informative)")
    parser.add_argument('-o', '--concatenate_output', action='store_true',
                        help="Boolean flag to output a single file for all chromosomes")
    parser.add_argument('-u', '--keep_non_unique', action='store_true',
                        help="Boolean flag to allow non-unique reads through the assignment process")
    
    # This creates a namespace object from the users CLI call
    args = parser.parse_args()
    
    # Quickly convert Namespace object to dictionary
    arg_dict = {arg: vars(args)[arg] for arg in vars(args)}
    
    # Print arguments, specifies arguments which will not be passed in the arg_dict
    print("\nGiven Arguments (ArgParse):")
    for key, arg in arg_dict.items():
        if not arg:
            print(f"\t{key} = {arg} -> (Will not be passed)")
        else:
            print(f"\t{key} = {arg}")
    print("\tDone.\n")
    
    # Recreate dict without arguments that did not receive any input
    arg_dict = {k: v for k, v in arg_dict.items() if v is not None}
    
    return arg_dict


def parseSamToDF(sam_file: str, headerlines: int = 20,
                 split_chrs: bool = True, num_lines: int = None,
                 print_rows: int = None, deep_memory: bool = False,
                 **kwargs) -> DataFrame or CHR_DF_DICT:
    """
    parseSamToDF
    
    This function does a lot of the grunt work in parsing, organizing/sorting and
        splitting up the specified SAM file.
    TODO: Need to pass in the headerlines value somehow. It will differ across sam files!
    """
    # Quick check to ensure the passed file path exists
    if not path.isfile(sam_file):
        print(f"\nFile does not exist at: {sam_file}, Terminating Script\n")
        exit()
    else:
        print(f"\nParsing identified file at: {sam_file}")
    
    # Assign known order of datatypes for .SAM file, this helps to speed up parsing
    o = 'object'
    i = 'int64'
    chr_cat = 'category'
    sam_dtypes_list = [o, i, chr_cat, i, i, o, o, i, i, o, o, o, o, o, o]
    sam_dtypes_dict = {i: sam_dtypes_list[i] for i in range(15)}
    
    # Parse sam file using pandas.read_csv function
    #   There is currently more work that can be front-loaded into this function call:
    #   - 'converters': The use of a the converters parameter could allow us to work with some of the
    #     functionality that is currently in the later steps of this script
    #   - 'low_memory': By default pandas will start breaking the file into chunks if there is not
    #     enough memory for the file to be read in all at once. For the most part this functionality
    #     is advantageous for our use case, but may be worth being aware of. The read_csv processes
    #     is dramatically stalled when this 'chunking' process is implemented.
    #   - 'memory_map': In the opposite direction, memory_map allows the read_csv process to load the
    #     full file to be parsed into memory, then reads it from there. If the used system is able to
    #     handle this much memory usage, this option can improve performance because there is no
    #     longer any I/O overhead.
    #   - 'comment': This parameter would provide the advantage of not specifying how many header
    #     lines are before the data in the file. The issue currently is that the '@' symbol used
    #     as the comment indicator for sam files allows comes up in the Phred alignment score scale.
    if num_lines:
        SAM_df = read_csv(sam_file,
                          sep="\t",
                          header=headerlines,
                          nrows=num_lines,  # If the CLI user specifies a number of lines to use, it comes in here
                          names=range(15),  # Just name columns by an index, to be renamed below
                          dtype=sam_dtypes_dict,  # The hope is this allows pandas to skip the type calling step
                          )
    else:
        SAM_df = read_csv(sam_file,
                          sep="\t",
                          header=headerlines,
                          names=range(15),
                          dtype=sam_dtypes_dict,
                          )
    
    # This sort_values with the ignore_index option on will currently reset the
    # indexes to the new sorted order, I don't know if this is good or bad...
    SAM_df = SAM_df.sort_values(by=[2, 3], ignore_index=True)
    
    # Rename columns, this can be the source of an error if a SAM file is not in
    #   the 15 column format - be wary!
    SAM_df = SAM_df.rename(columns={0: 'read_id',
                                    1: 'strand',
                                    2: 'chr',
                                    3: 'chr_pos',
                                    4: 'mapq',
                                    5: 'cigar',
                                    6: 'rnext',
                                    7: 'pnext',
                                    8: 'tlen',
                                    9: 'read_seq',
                                    10: 'phred_qual',
                                    11: 'NH',
                                    12: 'HI',
                                    13: 'AS',
                                    14: 'nM'})
    
    print(f'Finished parsing and sorting of file at: {sam_file}')
    
    if print_rows:
        print(f"\nFirst {print_rows} rows of dataframe:\n", SAM_df.head(print_rows), "\n")
    if deep_memory:
        print(SAM_df.info(memory_usage='deep'))
    
    # Attempting to split chromosomes
    if split_chrs:
        SAM_df_dict = dict(tuple(SAM_df.groupby('chr')))
        # if print_rows:
        #     for chr, df in SAM_df_dict.items():
        #         print(f"\nFirst {print_rows} rows of Chromosome-{chr}:\n", df.head(print_rows))
        print(f"Split -reads.SAM- dataframe into {len(SAM_df_dict.keys())} separate df's:\n{SAM_df_dict.keys()}")
        return SAM_df_dict
    else:
        return SAM_df


def parseAllChrsToDF(annot_file: str,
                     split_chrs: bool = True, num_lines: int = None,
                     print_rows: int = None, deep_memory: bool = False,
                     **kwargs) -> DataFrame or CHR_DF_DICT:
    """
    parseAllChrsToDF
    
    This function does a lot of the grunt work in parsing, organizing/sorting and
        splitting up the specified annotations file.
    """
    
    # Quick check to ensure the passed file path exists
    if not path.isfile(annot_file):
        print(f"\nFile does not exist at: {annot_file}, Terminating Script\n")
        exit()
    else:
        print(f"\nParsing identified file at: {annot_file}")
    
    # Parse annotations file using pandas.read_csv function
    #   There is currently more work that can be front-loaded into this function call:
    #   - 'converters': The use of a the converters parameter could allow us to work with some of the
    #     functionality that is currently in the later steps of this script
    #   - 'low_memory': By default pandas will start breaking the file into chunks if there is not
    #     enough memory for the file to be read in all at once. For the most part this functionality
    #     is advantageous for our use case, but may be worth being aware of. The read_csv processes
    #     is dramatically stalled when this 'chunking' process is implemented.
    #   - 'memory_map': In the opposite direction, memory_map allows the read_csv process to load the
    #     full file to be parsed into memory, then reads it from there. If the used system is able to
    #     handle this much memory usage, this option can improve performance because there is no
    #     longer any I/O overhead.
    if num_lines:
        annot_df = read_csv(annot_file,
                            delimiter='|',
                            names=['chr', 'gene_string'],
                            nrows=num_lines,
                            )
    else:
        annot_df = read_csv(annot_file,
                            delimiter='|',
                            names=['chr', 'gene_string'],
                            )
    
    # Split the chr_index column into two
    annot_df[['chr', 'chr_pos']] = DataFrame(annot_df['chr'].str.split('_').values.tolist(),
                                             index=annot_df.index)
    
    # Reorganize gene info for final .JAM file:
    #   This first step will currently trip up if an annotation file has multiple genes at a single index...
    #   Which will need to be handled when we get to that.
    annot_df[['gene_string', 'genes']] = DataFrame(
        annot_df['gene_string'].str.strip().str.split('\t').values.tolist(),
        index=annot_df.index)
    annot_df[['gene', 'sense']] = DataFrame(annot_df['genes'].str.strip().str.split(':').values.tolist(),
                                            index=annot_df.index)
    annot_df['gene_string'] = annot_df['gene_string'] + ':' + annot_df['sense']
    
    # Sort by Chromosome and index on chromosome
    annot_df = annot_df.sort_values(by=['chr', 'chr_pos'])
    
    # Change data types:
    annot_df = annot_df.astype({'chr': 'category',
                                'chr_pos': 'int64',
                                'gene': 'object',
                                'gene_string': 'object'})
    
    # Quickly reorder columns... might be completely superficial
    annot_df = annot_df[['chr', 'chr_pos', 'gene', 'gene_string']]
    
    print(f'Finished parsing and sorting of file at: {annot_file}')
    
    if print_rows:
        print(f"\nFirst {print_rows} rows of dataframe:\n",
              annot_df.sort_values(by=['chr', 'chr_pos']).head(print_rows), '\n')
    if deep_memory:
        print(annot_df.info(memory_usage='deep'))
    
    if split_chrs:
        annot_df_dict = dict(tuple(annot_df.groupby('chr')))
        # # Print for each chromosome if the user asked for print_rows
        # if print_rows:
        #     for chr, df in annot_df_dict.items():
        #         print(f"\nFirst {print_rows} rows of Chromosome-{chr}:\n", df.head(print_rows))
        # Short print to show how the dataframe was split up
        print(f"Split annotations dataframe into {len(annot_df_dict.keys())} separate df's:\n{annot_df_dict.keys()}")
        return annot_df_dict
    else:
        return annot_df


def recoverMappedPortion_perRead(Cigar, Read):
    # April 15, 2020: Stolen verbatim (+ my comments) from assignReadsToGenes4.py

    """Given a Cigar string and a Read, will return the sequence of the read that mapped to the genome."""
    # Edit Oct 10, 2013 to include skipped portions of reference sequence (introns)

    # first process the CIGAR string
    cigarSplit = findall('(\d+|[a-zA-Z]+)', Cigar)  # RE func: matches sets of digits or letters
    cigarSplit = [[int(cigarSplit[ii]), cigarSplit[ii + 1]] for ii in range(0, len(cigarSplit), 2)]

    # Then use that information to parse out nts of the read sequence
    mappedRead = ''
    ii = 0
    N = 0
    for entry in cigarSplit:
        if entry[1] in ['M',
                        'I']:  # then it's either aligned to the genomic sequence or has an insert relative to it
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


def recoverMappedPortion_perChr(chr_key, df_dict_in, df_dict_out, print_rows):
    df_in = df_dict_in[chr_key]
    print(f"Process {getpid():>5} working on {chr_key}")
    try:
        df_out = df_in.copy()
        # Apply the recover mapped portion function to each read
        df_out[['map_read_seq', 'N']] = \
            DataFrame(df_in.apply(lambda x: recoverMappedPortion_perRead(x['cigar'], x['read_seq']),
                                  axis=1).tolist(), index=df_in.index)

        print(f'Recovery of mapped portion complete for Chr-{chr_key:->4}, '
              f'read count={len(df_out.index)}/{len(df_in.index)}')
        if print_rows:
            print(df_out[chr_key][['read_id',
                                   'chr',
                                   'chr_pos',
                                   'cigar',
                                   'read_seq',
                                   'map_read_seq']].head(print_rows))
    # Catch error of empty dataframe for a chromosome
    except AttributeError:
        print(f'No reads for mapped portion recovery in Chr-{chr_key:->4}, '
              f'read count={len(df_in.index)}')
        df_out = df_in
    # This is where we will write to the new dictionary
    df_dict_out[chr_key] = df_out


def recoverMappedPortion_dfWrapper_MP(sam_df_dict: CHR_DF_DICT, print_rows: int = None,
                                      **kwargs) -> CHR_DF_DICT:
    """
    recoverMappedPortion_dfWrapper
    
    Function to accept the sam dataframe dictionary, recover the mapped portion of the
        original read for each read in each chromosome's dataframe
    """
    mappedPortionManager = multiprocessing.Manager()
    shared_input_dict = mappedPortionManager.dict(sam_df_dict)
    shared_output_dict = mappedPortionManager.dict()
    process_list = []
    # Loop through each chromosome
    for chr_key, chr_df in shared_input_dict.items():
        process = multiprocessing.Process(target=recoverMappedPortion_perChr, args=[chr_key,
                                                                                    shared_input_dict,
                                                                                    shared_output_dict,
                                                                                    print_rows])
        process_list.append(process)

    for pro in process_list:
        pro.start()
    for pro in process_list:
        pro.join()
    return shared_output_dict


def assignReadsToGenes(sam_df_dict: CHR_DF_DICT, annot_df_dict: CHR_DF_DICT,
                       print_rows: int = None, keep_non_unique: bool = False,
                       **kwargs) -> CHR_DF_DICT:
    """
    assignReadsToGenes
    
    This function utilizes the pandas.DataFrame.merge function which pulls together
        the reads dataframes and the annotations dataframes matching the rows based
        on the chromosome and the chromosome position. This is much faster than
        iterating through the reads dataframe because it is able to utilize some of
        the C backend of pandas (which was built with numpy).
    """
    
    print(f"\nAnnotation alignment for {len(sam_df_dict.keys())} chromosomes:")
    for chr_key, df in sam_df_dict.items():
        try:
            # ONLY KEEP UNIQUELY MAPPING READS: >>
            if not keep_non_unique:
                df = sam_df_dict[chr_key][sam_df_dict[chr_key]['NH'].str.endswith(':1')]
            # << ONLY KEEP UNIQUELY MAPPING READS
            
            print(f"\tPreforming alignment for Chr-{chr_key:->4} containing {len(df.index)} reads")
            # Using the df.merge() function for mapping annotations onto reads
            sam_df_dict[chr_key] = df.merge(annot_df_dict[chr_key], on=['chr_pos', 'chr'])
            print(f"\t\tSuccess, {len(sam_df_dict[chr_key][sam_df_dict[chr_key]['gene'] != numpy_nan].index):>7} "
                  f"reads assigned to genes in Chr-{chr_key:->4}\n")
            if print_rows:
                print(sam_df_dict[chr_key][['read_id',
                                            'chr',
                                            'chr_pos',
                                            'cigar',
                                            'read_seq',
                                            'gene',
                                            'gene_string']].head(print_rows))
            
            # TODO: some functionality to drop full df if nothing is assigned at all, this is only an issue when
            #       trying to work with a small slice of a annotations file for testing
        except KeyError as key:
            if str(key).strip("'") == chr_key:
                print(f"\t\tChr-{chr_key:->4} not found in annotations! -> Adding empty columns to compensate\n")
                sam_df_dict[chr_key]['gene'], sam_df_dict[chr_key]['gene_string'] = numpy_nan, numpy_nan
            else:
                print(f"\tOther KeyError pertaining to:", str(key), chr_key)
    return sam_df_dict


def fixSenseNonsense(sam_df_dict: CHR_DF_DICT,
                     **kwargs) -> CHR_DF_DICT:
    """
    fixSenseNonsense
    
    Check for strand and sense/antisense, create rev-compliment as needed, and edits gene_string
    """
    def revCompl(seq: str):
        """
        Stolen from assignReadsToGenes4.py
        Will return the reverse complement of a sequence
        """
        a = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
             'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
             'N': 'N', 'n': 'n', '.': '.', '-': '-'}
        return ''.join([a[seq[-i]] for i in range(1, len(seq) + 1)])
    
    def S_AS_flag_rework(S_AS: int, chr_pos: int, read_seq: str,
                         map_read_seq: str, N: int, gene_string: str):
        """Stolen from assignReadsToGenes4.py"""
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
            DataFrame(df.apply(lambda x: S_AS_flag_rework(x['strand'],
                                                          x['chr_pos'],
                                                          x['read_seq'],
                                                          x['map_read_seq'],
                                                          x['N'],
                                                          x['gene_string']),
                               axis=1).tolist(), index=df.index)
    return sam_df_dict


def outputToCSV(jam_file):
    pass


def main(sam_file: str, annot_file: str, output_prefix: str,
         print_rows: int = None, concatenate_output: bool = False,
         keep_non_unique: bool = False, **kwargs) -> None:
    sam_df_dict = parseSamToDF(sam_file, **kwargs)
    annot_df_dict = parseAllChrsToDF(annot_file, **kwargs)
    unassigned_df = DataFrame()
    
    # Going for the df.merge() function for mapping annotations onto reads
    annotated_sam_df_dict = assignReadsToGenes(sam_df_dict, annot_df_dict,
                                     print_rows=print_rows, keep_non_unique=keep_non_unique, **kwargs)
    
    # Apply the recoverMappedPortion() to dataframe to see how it does
    #  Currently doing this after dropping unassigned reads as this is a more time intensive step.
    fixed_annotated_sam_df_dict = recoverMappedPortion_dfWrapper(annotated_sam_df_dict, print_rows=print_rows, **kwargs)
    
    jam_df_dict = fixSenseNonsense(fixed_annotated_sam_df_dict, print_rows=print_rows, **kwargs)
    
    # Output to file
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
        jam_all_chrs = concat(jam_df_dict.values(), ignore_index=True)
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
    
    print("\n\nDone?!")


if __name__ == '__main__':
    arg_dict = parseArgs()
    main(**arg_dict)
