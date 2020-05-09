#!/usr/bin/env python3
"""
Joshua Arribere, July 25, 2014
Mar 23, 2020: Converted to python 3
April 2, 2020: Now accepts settings via a line-delimited txt file

Input: settings.txt - a line-delimited settings file in the format:
    adaptorSeq|(raw sequence of adaptor)
    minimumReadLength|(min length after adaptor and UMI trimming)
    maximumReadLength|(max length after adaptor and UMI trimming)
    UMI5|(5' UMI length in nts)
    UMI3|(3' UMI length in nts)
    genomeDir|(full path to genome directory)
    genomeAnnots|(full path to genome annotation file in gtf format)
    cores|(number of cores to use--ground control has 16 cores total)
    misMatchMax|(number of allowed mismatches)
    optString|(parameters for STAR run)
    optional:
        genomeDir2|(full path to genome directory for filter round of mapping)
        genomeAnnots2|(full path to genome annotation file for filter mapping)
        optString2|(parameters for filter mapping STAR run)

    inputReads.fastq - a fastq file of reads

run as python3 pipelineWrapper9.py inputReads.fastq settings.txt outPrefix
"""

# import sys, common
import sys
import os, readCollapser4, filterJamByReadLength, thecountReads
import argparse
import assignReadsToGenesDF
import infoGraphQC
# import metaStartStop
from logJosh import Tee

# Absolute defaults are overwritten by the given settings file and any command line arguments given
ABSOLUTE_DEFAULT_DICT = {'cores': 7, 'misMatchMax': 0,
                         'optString': '--outFilterScoreMinOverLread 1 '
                                      '--outFilterMatchNminOverLread 1 '
                                      '--outReadsUnmapped Fastx '
                                      '--outSJfilterOverhangMin 6 6 6 6',
                         'regenerate': False,
                         'misMatchMax2': 3,
                         'optString2': f'--outFilterScoreMin 14 '
                                       f'--outFilterScoreMinOverLread 0.3 '
                                       f'--outFilterMatchNmin 14 '
                                       f'--outFilterMatchNminOverLread 0.3 '
                                       f'--outReadsUnmapped Fastx '
                                       f'--outSJfilterOverhangMin 1000 1000 1000 1000 ',
                         'genomeDir2': False, 'genomeAnnots2': False}


def main(fastqFile, settings, outPrefix, adaptorSeq, minimumReadLength,
         maximumReadLength, genomeDir, genomeAnnots, cores, misMatchMax,
         umi5, umi3, optString, filterMap, optString2, genomeDir2, genomeAnnots2, misMatchMax2,
         keepNonUnique, outputJoshSAM, regenerate, **otherkwargs):
    """
    Main: Does the work of the pipeline. Large number of parameters accepts input from combineSettingsAndArguments as a
    keyword dictionary
    
    WARNING: if changing parameter names here, they will need to be changed:
        - in ABSOLUTE_DEFAULT_DICT,
        - in the settings files
        - in the argParse options
    """
    
    ############################################################################################################
    """First check to ensure fastq file exists, cutadapt WILL NOT throw an error if a non-existent fastq is passed"""
    ############################################################################################################
    if not os.path.isfile(fastqFile):
        print(f"\033[31;1m\nThe fastq file does not exist at: {fastqFile}, Terminating Script\n\033[0m\n")
        exit()
    
    ############################################################################################################
    """Trim adaptor from reads and sort by desired length"""
    ############################################################################################################
    # print('skipping trimming of any kind, so min/maximumReadLength restrictions ignored.')
    # os.system(f'cp {fastqFile} {outPrefix}.trimmed')
    cutadaptOutput = outPrefix + ".trimmed.selfDestruct.fastq"
    if not os.path.isfile(cutadaptOutput) or regenerate:
        print(f'Read length restriction does not include {umi5}Ns and {umi3}Ns.')
        print(f'To accommodate UMI length, this program will add {umi5 + umi3}nts to acceptable length.')
        os.system(f'cutadapt -a {adaptorSeq} '
                  f'-j {cores} '
                  f'-m {minimumReadLength+umi5+umi3} '
                  f'-M {maximumReadLength+umi5+umi3} '
                  f'--too-short-output {outPrefix + ".trimmed.selfDestruct.tooShort.fastq"} '
                  f'--too-long-output {outPrefix + ".trimmed.selfDestruct.tooLong.fastq"} '
                  f'{fastqFile} > {cutadaptOutput} '
                  f'2>/dev/null'
                  )
        regenerate = True
    else:
        print(f"Reusing cutadapt output from {os.stat(cutadaptOutput).st_mtime}\n\t"
              f"(file: {cutadaptOutput})\n\t"
              f"If this is not intended: use -r or --regenerate flag to regenerate all files.\n\t"
              f"\033[1mThis functionality does not take into account changes in run parameters!!\033[0m")
    ############################################################################################################
    """Collapse reads and trim off UMIs"""
    ############################################################################################################
    # It doesn't make sense to run readCollapser if there's no UMI
    if 0<umi5+umi3<=6:
        print(f"\nYour combined UMI length is {umi5 + umi3}, which is pretty short.\n"
              f"I'm going to try and collapse based on it, assuming you know what\n"
              f"you are doing. But if you do not understand this message, please\n"
              f"go find someone who can help you.\n")
    # Need UMIs in order to run readCollapser
    if umi5+umi3!=0:
        readCollapsedOutput = outPrefix+".trimmed.collapsed.selfDestruct.fastq"
        if not os.path.isfile(readCollapsedOutput) or regenerate:
            readCollapser4.main([outPrefix+'.trimmed.selfDestruct.fastq',
                                 umi5, umi3, readCollapsedOutput])
            regenerate = True
        else:
            print(f"Reusing readCollapser output from {os.stat(readCollapsedOutput).st_mtime}\n\t"
                  f"(file: {readCollapsedOutput})\n\t"
                  f"If this is not intended: use -r or --regenerate flag to regenerate all files.\n\t"
                  f"\033[1mThis functionality does not take into account changes in run parameters!!\033[0m")
    else:
        print('No UMI length given: Skipping collapsing...')
        # The next line creates a symbolic link for the .collapsed file location
        # instead of the prior way, which copied them.
        os.system(f'ln -s {outPrefix}.trimmed.selfDestruct.fastq '
                  f'{outPrefix}.trimmed.collapsed.selfDestruct.fastq')
    
    ############################################################################################################
    """Introduce a variable to make reading code easier"""
    ############################################################################################################
    readFile = f'{outPrefix}.trimmed.collapsed.selfDestruct.fastq'
    
    ############################################################################################################
    """Perform a filter round of mapping. e.g. to rDNA or RNAi trigger"""
    ############################################################################################################
    if filterMap:
        if genomeDir2 and genomeAnnots2:
            filterMapOutput = outPrefix+'.trimmed.collapsed.mapped.filterUnmapped.out.mate1'
            if not os.path.isfile(filterMapOutput) or regenerate:
                print(f'Performing filter round of mapping to {genomeDir2}')
                print(f'Only running on {cores} cores.')
                print(f'{misMatchMax2} mismatch max!')
                print(f'Length/Score parameters: --outFilterMismatchNmax {misMatchMax2} {optString2}')
                os.system(f'STAR {optString2} '
                          f'--outFilterMismatchNmax {misMatchMax2} '
                          f'--alignIntronMax 1 '
                          f'--sjdbGTFfile {genomeAnnots2} '
                          f'--genomeDir {genomeDir2} '
                          f'--readFilesIn {readFile} '
                          f'--runThreadN {cores} '
                          f'--outFileNamePrefix {outPrefix}.trimmed.collapsed.mapped.filter'
                          )
                regenerate = True
            else:
                print(f"Reusing filterMap output from {os.stat(filterMapOutput).st_mtime}\n\t"
                      f"(file: {filterMapOutput})\n\t"
                      f"If this is not intended: use -r or --regenerate flag to regenerate all files.\n\t"
                      f"\033[1mThis functionality does not take into account changes in run parameters!!\033[0m")
            # Now rewrite the read file to map from the unmapped reads (regen or not!)
            readFile = outPrefix + '.trimmed.collapsed.mapped.filterUnmapped.out.mate1'
        else:
            print('No file given for genomeDir2 or genomeAnnots2, skipping filter round of mapping...')
    else:
        print('skipping filter round of mapping...')
    
    ############################################################################################################
    """Commence read mapping"""
    ############################################################################################################
    starCheckFile = outPrefix + ".finalMapped.Aligned.out.sam"
    if not os.path.isfile(starCheckFile) or regenerate:
        print(f'Only running on {cores} cores.')
        print(f'{misMatchMax} mismatch max!')
        print(f'Length/Score parameters: {optString} --outFilterMismatchNmax {misMatchMax}')
        os.system(f'STAR {optString} '
                  f'--outFilterMismatchNmax {misMatchMax} '
                  f'--alignIntronMax 1 '
                  f'--sjdbGTFfile {genomeAnnots} '
                  f'--genomeDir {genomeDir} '
                  f'--readFilesIn {readFile} '
                  f'--runThreadN {cores} '
                  f'--outFileNamePrefix {outPrefix}.finalMapped.')
        regenerate = True
    else:
        print(f"Reusing STAR run output from {os.stat(starCheckFile).st_mtime}\n\t"
              f"(file checked: {starCheckFile})\n\t"
              f"If this is not intended: use -r or --regenerate flag to regenerate all files.\n\t"
              f"\033[1mThis functionality does not take into account changes in run parameters!!\033[0m")
    
    ############################################################################################################
    """Assign reads to genes"""
    ############################################################################################################
    assignReadsOutput = outPrefix + ".allChrs.jam"
    if not os.path.isfile(assignReadsOutput) or regenerate:
        print('\nAssigning reads to genes', end=' ')
        if keepNonUnique:
            print('allowing for multiply-mapping reads...')
        else:
            print('only allowing uniquely-mapping reads...')
        
        genomeAnnotProcessed=genomeAnnots.strip('gtf')+'allChrs.txt'
        assignReadsToGenesDF.main(outPrefix +'.finalMapped.Aligned.out.sam',
                                  genomeAnnotProcessed,
                                  outPrefix,
                                  keep_non_unique=keepNonUnique,  # If this flag is passed from argParse
                                  #                                   (or settings.txt) assignReadsToGenesDF will *also*
                                  #                                   output a file titled:
                                  #                                   outputPrefix.redundantAndUnique.allChrs.jam
                                  output_joshSAM=outputJoshSAM,  # Output old format joshSAM file in addition to .jam if
                                  #                                 this flag is passed. Also stacks with keepNonUnique
                                  )
        print('Done with read assignment!')
        regenerate = True
    else:
        print(f"Reusing assignReadsToGenes output from {os.stat(assignReadsOutput).st_mtime}\n\t"
              f"(file checked: {assignReadsOutput})\n\t"
              f"If this is not intended: use -r or --regenerate flag to regenerate all files.\n\t"
              f"\033[1mThis functionality does not take into account changes in run parameters!!\033[0m")
    
    ############################################################################################################
    """Additional filtering of reads by length"""
    ############################################################################################################
    print('filtering read lengths again...')
    filterJamByReadLength.main([outPrefix+'.allChrs.jam',
                                minimumReadLength,
                                maximumReadLength,
                                outPrefix+'.filtered_%s-%snt.jam'%(minimumReadLength,maximumReadLength)])
    
    ############################################################################################################
    """Creating infographic"""
    ############################################################################################################
    print('Making infographic')
    infoGraphQC.main([outPrefix+'.allChrs.jam',minimumReadLength,maximumReadLength,-21,21,-30,12,outPrefix+'.qc'])
    thecountReads.main([fastqFile, outPrefix])

def parseArguments():
    """
    Outputs a dictionary of argument keys (taken either from the metavar parameter below or the --name for flags)
        and the passed argument's item. This can be used with
    """
    parser=argparse.ArgumentParser(description='Wrapper for the handling of libraries starting from fastq files.')
    # Required Arguments:
    parser.add_argument('fastqFile', metavar='fastqFile',
                        type=str, help='The fastq file input')
    parser.add_argument('settings', metavar='settings', type=str,
                        help='A text file containing a line-delimited input of adaptorSeq, genomeDir, and genomeAnnots.'
                             ' Plus any of the following optional arguments. All in arg|value format')
    parser.add_argument('outPrefix', metavar='outPrefix', type=str,
                        help='The outPrefix that all files produced will be named as.')
    # Optional Arguments: (None default here allows us to not pass anything that isn't given by user.
    #                      This helps to simplify settings.txt/default overwrites further down the line.)
    parser.add_argument('--umi5', metavar='umi5', type=int, default=None,
                        help='The number of nucleotides to be trimmed from the 5prime end of reads.')
    parser.add_argument('--umi3', metavar='umi3', type=int, default=None,
                        help='The number of nucleotides to be trimmed from the 3prime end of reads.')
    parser.add_argument('--minimumReadLength', '--min', type=int, default=None,
                        help='The shortest reads to be mapped to the genome.')
    parser.add_argument('--maximumReadLength', '--max', type=int, default=None,
                        help='The longest reads to be mapped to the genome.')
    parser.add_argument('--adaptorSeq', '--adaptor', metavar='adaptor', type=str, default=None,
                        help='3\' adaptor to be trimmed off.')
    parser.add_argument('--genomeDir', metavar='genomeDir', type=str, default=None,
                        help='Genome directory where STAR index can be found.')
    parser.add_argument('--genomeAnnots', metavar='genomeAnnots', type=str, default=None,
                        help='Genome annotations (gtf format).')
    parser.add_argument('--cores', metavar='cores', type=int, default=None,
                        help='Number of cores to use.')
    parser.add_argument('--misMatchMax', metavar='misMatchMax', type=int, default=None,
                        help='Number of mismatches to tolerate during mapping.')
    # Flag Arguments: (just add these as tags to change pipeline functionality)
    parser.add_argument('-u', '--keepNonUnique', action='store_true',
                        help="Boolean flag to allow non-unique reads in the final .jam file(s).")
    parser.add_argument('-j', '--outputJoshSAM', action='store_true',
                        help="Boolean flag to also output joshSAM format files in addition to the jam")
    parser.add_argument('-p', '--print_arguments', action='store_true',
                        help="Boolean flag to show how arguments are overwritten/accepted")
    parser.add_argument('-f', '--filterMap', action='store_true',
                        help="Boolean flag to perform filter mapping")
    parser.add_argument('-r', '--regenerate', action='store_true',
                        help="Boolean flag to ignore previously produced files and generate all files anew")
    
    # Spit out namespace object from argParse
    args = parser.parse_args()
    # Quickly convert Namespace object to dictionary
    arg_dict = {arg: vars(args)[arg] for arg in vars(args)}
    
    if arg_dict['print_arguments']:
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

def parseSettings(settings,print_arguments=False,**other_kwargs):
    """
    Will loop through and replace variables that are None
    """
    ##first, parse the settings file to a dictionary called settingsDict
    settingsDict={}
    with open(settings,'r') as f:
        for line in f:
            if not line.startswith('#'):
                line=line.strip()
                if line!='':
                    line=line.split('|')
                    if len(line) == 2:
                        settingsDict[line[0]] = line[1]
                    else:
                        print("\033[31;1m\nRemove pipes ('|') from settings file arguments (or rewrite parser)\n\033[0m")
                        raise ImportError
    if print_arguments:
        print(f"\nSettings Arguments (file: '{settings}')")
        for key, arg in settingsDict.items():
            if not arg:
                print(f"\t{key} = {arg} -> (Will not be passed)")
            else:
                print(f"\t{key} = {arg}")
    settingsDict = {k: v for k, v in settingsDict.items() if v is not None and v is not ''}
    return settingsDict

def combineSettingsAndArguments():
    absoluteDefDict = ABSOLUTE_DEFAULT_DICT
    argDict = parseArguments()
    settingsDict = parseSettings(**argDict)
    finalArgDict = {}
    
    finalArgDict.update(absoluteDefDict)
    finalArgDict.update(settingsDict)
    finalArgDict.update(argDict)
    print("\033[1m\nPipeline Arguments:")
    # Absolute defaults overwritten by settings.txt then overwritten by CLI args
    for key, arg in finalArgDict.items():
        print(f"\t{key} = {arg}")
        try:
            finalArgDict[key] = int(arg)
        except ValueError:
            finalArgDict[key] = str(arg)
        except TypeError:
            finalArgDict[key] = str(arg)
    print(f"cutadapt version:")
    os.system('cutadapt --version')
    print(f"STAR version:")
    os.system('STAR --version')
    print('\n\033[0m\n')
    return finalArgDict

if __name__ == '__main__':
    Tee()
    argument_dict = combineSettingsAndArguments()
    main(**argument_dict)
