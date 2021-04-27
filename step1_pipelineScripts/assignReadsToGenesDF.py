"""
**newest**
assignReadsToGenesDF.py
Marcus Viscardi     April 14, 2020

Basically just rewriting Josh's assignReadsToGenes4.py, with the change of
    utilizing pandas dataframes in an effort to speed up read assignment.
Currently full script is self contained. Pieces like the parsers and the
    revCompliment functions were written into here to make it a little more
    portable.
"""

from os import path
from sys import exit
from argparse import ArgumentParser
from re import findall
from typing import Dict, Any
from csv import QUOTE_NONE

from timeit import default_timer
from numpy import nan as numpy_nan
from pandas import set_option, read_csv, DataFrame, concat, errors

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
    parser.add_argument('-m', '--minLength', metavar='minLength', type=int,
                        default=None, help="Minimum post read recovery length")
    parser.add_argument('-M', '--maxLength', metavar='minLength', type=int,
                        default=None, help="Maximum post read recovery length")
    # Flag arguments which turn on functionality if passed:
    parser.add_argument('--deep_memory', action='store_true',
                        help="Boolean flag to print dataframe deep memory info\n"
                             "(this can be very CPU/time intensive, but informative)")
    parser.add_argument('-u', '--keep_non_unique', action='store_true',
                        help="Boolean flag to allow non-unique reads through the assignment process")
    parser.add_argument('-j', '--output_joshSAM', action='store_true',
                        help="Boolean flag to output a .joshSAM file instead")
    
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
    
    # Recreate dict without arguments that did not receive any input
    arg_dict = {k: v for k, v in arg_dict.items() if v is not None}
    
    return arg_dict


def parseSamToDF(sam_file: str, keep_non_unique: bool = False,
                 split_chrs: bool = True, deep_memory: bool = False,
                 **kwargs) -> DataFrame or CHR_DF_DICT:
    """
    parseSamToDF
    
    This function does a lot of the grunt work in parsing, organizing/sorting and
        splitting up the specified SAM file.
    """
    # Quick check to ensure the passed file path exists
    if not path.isfile(sam_file):
        print(f"\033[31;1m\nFile does not exist at: {sam_file}, Terminating Script\n\033[0m\n")
        exit()
    else:
        print(f"\nParsing read file at: {sam_file}")
    
    # Assign known order of datatypes for .SAM file, this helps to speed up parsing
    o = 'object'
    i = 'int64'
    chr_cat = 'category'
    sam_dtypes_list = [o, i, chr_cat, i, i, o, o, i, i, o, o, o, o, o, o]
    sam_dtypes_dict = {i: sam_dtypes_list[i] for i in range(15)}
    
    # First figure out how many header lines are in the file,
    #   Pandas has a hard time parsing @comments without messing up the Phred scores in the file
    headerlines = 0
    with open(sam_file, 'r')as sam_file_quick:
        for line in sam_file_quick:
            if line.startswith('@'):
                headerlines += 1
            else:
                break
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
    try:
        SAM_df = read_csv(sam_file,
                          sep="\t",
                          header=headerlines,
                          names=range(15),
                          dtype=sam_dtypes_dict,
                          )
    except errors.ParserError as error_message:
        print(f'\033[31;1m\nPandas Error: {error_message}\nThe SAM files passed to assignReadsToGenes has no reads, '
              f'please check previous steps.\n\033[0m')
        exit()
    # Sort my chr (2) and chr_pos (3)
    SAM_df = SAM_df.sort_values(by=[2, 3])
    
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
    
    print(f'Finished parsing and sorting of read file')
    
    if deep_memory:
        print(SAM_df.info(memory_usage='deep'))
    # ONLY KEEP UNIQUELY MAPPING READS:
    if not keep_non_unique:
        print(f"Filtering out non-uniquely mapped reads...")
        SAM_df = SAM_df[SAM_df['NH'].str.endswith(':1')]
    else:
        print(f"Keeping non-uniquely mapped reads...")
    
    if split_chrs:
        SAM_df_dict = dict(tuple(SAM_df.groupby('chr')))
        print(f"Split reads file dataframe into {len(SAM_df_dict.keys())} separate df's:\n{SAM_df_dict.keys()}")
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
        print(f"\033[31;1m\nFile does not exist at: {annot_file}, Terminating Script\n\033[0m")
        exit()
    else:
        print(f"\nParsing annotation file at: {annot_file}")
    
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
    annot_df = read_csv(annot_file,
                        delimiter='\t',
                        names=['chr', 'gene', 'gene_string'])
    
    # Split the chr_index column into two
    try:
        annot_df[['chr', 'chr_pos']] = DataFrame(annot_df['chr'].str.split('_').values.tolist(),
                                                 index=annot_df.index)
    except ValueError as error:
        with open(annot_file, 'r') as file:
            line_one = file.readline()
        print(f"\n\033[31;1m\nError: \"{error}\"\n"
              f"\tLikely that genome annotation file is in incorrect format (file: {annot_file})\n"
              f"\tPlease ensure that format is:\tchr_chr-pos\tgene\ttranscript(s)(separated by '|')\n"
              f"\tFirst line of passed file:   \t{line_one}\n"
              f"Terminating...\n\033[0m")
        exit()
    
    # Sort by Chromosome and index on chromosome
    annot_df = annot_df.sort_values(by=['chr', 'chr_pos'])
    
    # Change data types:
    annot_df = annot_df.astype({'chr': 'category',
                                'chr_pos': 'int64',
                                'gene': 'object',
                                'gene_string': 'object'})
    
    # Quickly reorder columns (completely superficial)
    annot_df = annot_df[['chr', 'chr_pos', 'gene', 'gene_string']]
    
    print(f'Finished parsing and sorting of annotations file')
    
    if print_rows and print_rows == int:
        print(f"\nFirst {print_rows} rows of dataframe:\n",
              annot_df.sort_values(by=['chr', 'chr_pos']), '\n')
    elif print_rows:
        print(f"PrintError: print_rows parameter must be an integer, you passed: {print_rows},"
              f"which is a {type(print_rows)}.")
    if deep_memory:
        print(annot_df.info(memory_usage='deep'))
    
    if split_chrs:
        annot_df_dict = dict(tuple(annot_df.groupby('chr')))
        # Short print to show how the dataframe was split up
        print(f"Split annotations dataframe into {len(annot_df_dict.keys())} separate df's"
              f":\n{annot_df_dict.keys()}")
        return annot_df_dict
    else:
        return annot_df


def recoverMappedPortion_dfWrapper(sam_df_dict: CHR_DF_DICT, print_rows: int = None,
                                   **kwargs) -> CHR_DF_DICT:
    """
    recoverMappedPortion_dfWrapper
    
    Function to accept the sam dataframe dictionary, recover the mapped portion of the
        original read for each read in each chromosome's dataframe
    """
    
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
    
    print(f"\nRecovering mapped portion of reads for {len(sam_df_dict.keys())} chromosomes:")
    # Loop through each chromosome
    for chr_key, df in sam_df_dict.items():
        try:
            # Apply the recover mapped portion function to each read
            sam_df_dict[chr_key][['map_read_seq', 'N']] = \
                DataFrame(df.apply(lambda x: recoverMappedPortion_perRead(x['cigar'], x['read_seq']),
                                   axis=1).tolist(), index=df.index)
            
            print(f'Chr-{chr_key:->4} mapped portion recovery complete, '
                  f'read count={len(sam_df_dict[chr_key].index)}')
            if print_rows:
                print(sam_df_dict[chr_key][['read_id',
                                            'chr',
                                            'chr_pos',
                                            'cigar',
                                            'read_seq',
                                            'map_read_seq']].head(print_rows))
        # Catch error of empty dataframe for a chromosome
        except AttributeError:
            print(f'No reads for mapped portion recovery in Chr-{chr_key:->4}, '
                  f'read count={len(sam_df_dict[chr_key].index)}')
    return sam_df_dict


def assignReadsToGenes(sam_df_dict: CHR_DF_DICT, annot_df_dict: CHR_DF_DICT,
                       print_rows: int = None, **kwargs) -> (CHR_DF_DICT, CHR_DF_DICT):
    """
    assignReadsToGenes
    
    This function utilizes the pandas.DataFrame.merge function which pulls together
        the reads dataframes and the annotations dataframes matching the rows based
        on the chromosome and the chromosome position. This is much faster than
        iterating through the reads dataframe because it is able to utilize some of
        the C backend of pandas (which was built with numpy).
    """
    unassigned_df_dict = {}
    print(f"\nAnnotation alignment for {len(sam_df_dict.keys())} chromosomes:")
    for chr_key, df in sam_df_dict.items():
        try:
            # Using the df.merge() function for mapping annotations onto reads
            # TODO: In the future we could utilize the parameter: "on='left'" in the merge call. This
            #       will provide the advantage(?) of allowing unmapped reads through, meaning we can
            #       pass these into analysis or QC scripts as needed.
            sam_df_dict[chr_key] = df.merge(annot_df_dict[chr_key], how='left', on=['chr', 'chr_pos'])
            unassigned_df_dict[chr_key] = sam_df_dict[chr_key][sam_df_dict[chr_key]['gene'].isnull()]
            sam_df_dict[chr_key] = sam_df_dict[chr_key][sam_df_dict[chr_key]['gene'].notnull()]
            print(f"Chr-{chr_key:->4} genes assigned, read count="
                  f"{len(sam_df_dict[chr_key].index):>8}", end="")
            # If there are any unassigned reads, print this in the output
            if len(unassigned_df_dict[chr_key].index) > 0:
                print(f", unassigned read count="
                      f"{len(unassigned_df_dict[chr_key].index):>5}")
            else:
                print()
            if print_rows:
                print(sam_df_dict[chr_key][['read_id',
                                            'chr',
                                            'chr_pos',
                                            'cigar',
                                            'read_seq',
                                            'gene',
                                            'gene_string']].head(print_rows))
        except KeyError as key:
            if str(key).strip("'") == chr_key:
                print(f"\t\tChr-{chr_key:->4} not found in annotations! -> Adding empty columns to compensate\n")
                sam_df_dict[chr_key]['gene'], sam_df_dict[chr_key]['gene_string'] = numpy_nan, numpy_nan
            else:
                print(f"\tOther KeyError pertaining to:", str(key), chr_key)
    return sam_df_dict, unassigned_df_dict


def fixSenseNonsense(sam_df_dict: CHR_DF_DICT,
                     **kwargs) -> CHR_DF_DICT:
    """
    fixSenseNonsense
    
    Check for strand and sense/antisense, create rev-compliment as needed, and edits gene_string column
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
                         map_read_seq: str, N: int):
        """Stolen from assignReadsToGenes4.py"""
        if S_AS & 16 != 0:  # check strand
            chr_pos += len(map_read_seq) + N - 1
            read_strand = '-'
            map_read_seq = revCompl(map_read_seq)
            read_seq = revCompl(read_seq)
        else:
            read_strand = '+'
        return read_strand, read_seq, map_read_seq, chr_pos
    
    print(f"\nApplying sense and rev-comp changes to {len(sam_df_dict.keys())} chromosomes:")
    for chr_key, df in sam_df_dict.items():
        sam_df_dict[chr_key][['strand', 'read_seq', 'map_read_seq', 'chr_pos']] = \
            DataFrame(df.apply(lambda x: S_AS_flag_rework(x['strand'],
                                                          x['chr_pos'],
                                                          x['read_seq'],
                                                          x['map_read_seq'],
                                                          x['N'],
                                                          ),
                               axis=1).tolist(), index=df.index)
        print(f'Chr-{chr_key:->4} sense and rev-comp changes applied, '
              f'read count={len(sam_df_dict[chr_key].index)}')
    return sam_df_dict


def finalFixers(sam_df_dict: CHR_DF_DICT, **kwargs) -> CHR_DF_DICT:
    """
    finalFixers
    
    Function to apply last touches on reads:
        -Add sense or antisense call to gene
        -Put HI (hit index) and NH (number of hits) together
    """
    def perReadFinalFix(gene: str, strand: str,
                        HI: str, NH: str):
        # This has some issues if reads without annotations get through >>>
        if str(gene)[-1] == strand:
            new_gene = str(gene)[:-1] + 'S'
        else:
            new_gene = str(gene)[:-1] + 'AS'
        # Pull out everything after the last ':'
        hit_index = HI.split(':')[-1]
        number_of_hits = NH.split(':')[-1]
        # This is marginally faster than f-string concatenation
        HINH = hit_index + ':' + number_of_hits
        return HINH, new_gene
    
    print(f"\nFinalizing formatting for {len(sam_df_dict.keys())} chromosomes:")
    # Apply it to each chr DF:
    for chr_key, df in sam_df_dict.items():
        if not df.empty:
            sam_df_dict[chr_key][['HI:NH', 'gene']] = DataFrame(df.apply(lambda x: perReadFinalFix(x['gene'],
                                                                                                   x['strand'],
                                                                                                   x['HI'],
                                                                                                   x['NH']),
                                                                         axis=1).tolist(), index=df.index)
            # Add read_length column for filtering and joshSAM output (if flagged for)
            sam_df_dict[chr_key]['read_length'] = sam_df_dict[chr_key]['map_read_seq'].str.len()
            print(f'Chr-{chr_key:->4} finalized, '
                  f'read count={len(sam_df_dict[chr_key].index)}')
        else:
            print(f'Chr-{chr_key:->4} is empty!')
    return sam_df_dict


def filterByReadLength(jam_all_chrs: DataFrame, jam_columns: list, maxLength: int, minLength: int, output_prefix: str):
    if minLength and maxLength:
        # Create a dataframe copy with all of the too short/long reads
        filter_out = jam_all_chrs[(jam_all_chrs['read_length'] > maxLength) | (jam_all_chrs['read_length'] < minLength)]
        if len(filter_out.index) > 0:
            print(f"\nFiltering out {len(filter_out.index)} reads shorter than {minLength}nts or longer "
                  f"than {maxLength}nts to file: "
                  f"{output_prefix}.tooShortOrLong.jelly (this will be in a jam format)")
            # Output the jelly format for reads that are too short/long:
            filter_out.to_csv(f"{output_prefix}.tooShortOrLong.jelly",
                              index=False, sep='\t',
                              columns=jam_columns)
            # Overwrite the main dataframe to remove all the reads that were just written out to the jelly file:
            jam_all_chrs = jam_all_chrs[(jam_all_chrs['read_length'] >= minLength) & \
                                        (jam_all_chrs['read_length'] <= maxLength)]
        else:
            print(f"\nNo reads < {minLength}nts or > {maxLength}nts identified. "
                  f"{output_prefix}.tooShortOrLong.jelly will not be generated")
    else:
        print(f"\nSkipping filtering of final read lengths...")
    return jam_all_chrs


def adjustForJoshSAM(jam_all_chrs: DataFrame) -> DataFrame:
    """
    4/26/2021     Parissa noticed that downstream files didn't work with the new joshSAM outputs
    To fix this I am adjusting the data in the main dataframe as it's being outputted so that it
    will look EXACTLY like the old joshSAM files
    
    The main issue seems to be with the gene and gene_string columns
        old style:
            gene column:        gene_id
            gene_string col:    hit1:S/AS\thit2:S/AS
        new (broken style):
            gene column:        gene_id:S/AS
            gene_string col:    hit1|hit2|hit3
            
    Plan:
        Replace Pipe characters with SorAS+\t
        
    *Note: I think Matt also noticed this issue at one point... I said I would fix it and forgot to...
    """
    jam_all_chrs["gene_string"] = jam_all_chrs.apply(lambda row: row["gene_string"].replace("|", f":{row['gene'].split(':')[-1]}\t") + f":{row['gene'].split(':')[-1]}", axis=1)
    jam_all_chrs["gene"] = jam_all_chrs.apply(lambda row: row["gene"].split(":")[0], axis=1)
    return jam_all_chrs


def main(sam_file: str, annot_file: str, output_prefix: str,
         print_rows: int = None, keep_non_unique: bool = False,
         output_joshSAM: bool = False, minLength: int = None,
         maxLength: int = None, **kwargs) -> None:
    
    start_time = default_timer()  # Timer
    
    # Load the read and annotation files into pandas DFs
    annot_df_dict = parseAllChrsToDF(annot_file, **kwargs)
    sam_df_dict = parseSamToDF(sam_file, keep_non_unique=keep_non_unique, **kwargs)
    
    ####################################################################################################################
    """Treatment + Assignment of reads from SAM file"""
    ####################################################################################################################
    # Apply the recoverMappedPortion() to dataframe
    post_map_df_dict = recoverMappedPortion_dfWrapper(sam_df_dict, print_rows=print_rows, **kwargs)
    
    # Handle +/- and Sense/Antisense issues from SAM format
    post_sense_antisense_df_dict = fixSenseNonsense(post_map_df_dict, print_rows=print_rows, **kwargs)
    
    # Use df.merge() function for mapping annotations onto reads, also spits out unassigned reads
    assigned_df_dict, unassigned_df_dict = assignReadsToGenes(post_sense_antisense_df_dict, annot_df_dict,
                                                              print_rows=print_rows, **kwargs)
    
    # Add the HitIndex:NumberOfHits column from the SAM HI and NH columns
    jam_df_dict = finalFixers(assigned_df_dict, **kwargs)
    
    # Columns and order for final output(s):
    jam_columns = ['read_id',
                   'chr',
                   'chr_pos',
                   'strand',
                   'mapq',
                   'cigar',
                   'map_read_seq',
                   'HI:NH',
                   'gene',
                   'gene_string']
    joshSAM_columns = ['chr',
                       'chr_pos',
                       'strand',
                       'map_read_seq',
                       'read_length',
                       'gene',
                       'gene_string']
    
    # Bring all chromosomes together
    jam_all_chrs = concat(jam_df_dict.values(), ignore_index=True)
    unassigned_all_chrs = concat(unassigned_df_dict.values(), ignore_index=True)
    # Sort by chromosome and chr_pos
    jam_all_chrs.sort_values(by=['chr', 'chr_pos'], inplace=True)
    ####################################################################################################################
    
    ####################################################################################################################
    """Output any reads that were not mapped to the annotation file"""
    ####################################################################################################################
    if not unassigned_all_chrs.empty:
        unassigned_all_chrs.sort_values(by=['chr', 'chr_pos'], inplace=True)
        unassigned_all_chrs['HI:NH'] = unassigned_all_chrs.apply(lambda row:
                                                                 row['HI'].split(':')[-1] + ':' +
                                                                 row['NH'].split(':')[-1],
                                                                 axis=1)
        unassigned_all_chrs.to_csv(f"{output_prefix}.unassignedReads.jelly",
                                   index=False, sep='\t',
                                   columns=jam_columns)
    else:
        print(f"\nNo unassigned reads. {output_prefix}.unassignedReads.jelly will not be generated.")
    ####################################################################################################################
    """Filter out and output reads that are too long or too short:"""
    ####################################################################################################################
    jam_all_chrs = filterByReadLength(jam_all_chrs, jam_columns, maxLength, minLength, output_prefix)
    ####################################################################################################################
    
    ####################################################################################################################
    """Output of jam file (and other files by request/argument-flag)"""
    ####################################################################################################################
    # Write unique reads to .jam file:
    #   Everything before the '.to_csv' is the filter action to only output unique reads
    jam_all_chrs[jam_all_chrs['NH'].str.endswith(':1')].to_csv(f"{output_prefix}.allChrs.jam",
                                                               index=False, sep='\t',
                                                               columns=jam_columns)
    # Below is to handle the other cases of wanting to keep multiply-mapping reads or a joshSAM output
    if keep_non_unique:
        # If flag is called, write non-unique reads and unique reads to another .jam file
        jam_all_chrs.to_csv(f"{output_prefix}.redundantAndUnique.allChrs.jam",
                            index=False, sep='\t',
                            columns=jam_columns)
        if output_joshSAM:
            # joshSAM files are organized based on the original SAM file (this info should be retained by the indexes)
            jam_all_chrs.sort_index(inplace=True)
            # Write joshSAM file with non-unique reads and unique reads >>>
            adjustForJoshSAM(jam_all_chrs).to_csv(f"{output_prefix}.redundantAndUnique.allChrs.joshSAM",
                                                  index=False, sep='\t',
                                                  columns=joshSAM_columns, header=False,
                                                  quoting=QUOTE_NONE, quotechar="", escapechar=" ")
    if output_joshSAM:
        # joshSAM files are organized based on the original SAM file (this info should be retained by the indexes)
        jam_all_chrs.sort_index(inplace=True)
        # Write joshSAM file with only unique reads:
        #   Everything before the '.to_csv' is the filter action to only output unique reads
        adjustForJoshSAM(jam_all_chrs[jam_all_chrs['NH'].str.endswith(':1')]).to_csv(f"{output_prefix}.allChrs.joshSAM",
                                                                                     index=False, sep='\t',
                                                                                     columns=joshSAM_columns,
                                                                                     header=False,
                                                                                     quoting=QUOTE_NONE, quotechar="",
                                                                                     escapechar=" ")
    ####################################################################################################################
    
    end_time = default_timer()  # Timer
    total_seconds = end_time - start_time
    print(f"\nTotal Time for assignReadsToGenesDF:{total_seconds // 60:3.0f}min{total_seconds % 60:3.0f}sec\n")


if __name__ == '__main__':
    arg_dict = parseArgs()
    main(**arg_dict)
