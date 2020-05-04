#!/usr/bin/env python3
"""
Joshua Arribere, July 25, 2014
Mar 23, 2020: Converted to python 3
April 2, 2020: Now accepts settings via a line-delimited txt file

Input: settings.txt - a line-delimited settings file in the format:
    adaptorSeq (raw sequence of adaptor)
    minReadLength (min length after adaptor and UMI trimming)
    maxReadLength (max length after adaptor and UMI trimming)
    UMI5 (5' UMI length in nts)
    UMI3 (3' UMI length in nts)
    genomeDir (full path to genome directory)
    genomeAnnots (full path to genome annotation file in gtf format)
    cores (number of cores to use--ground control has 16 cores total)
    misMatchMax (number of allowed mismatches)
    optString (parameters for STAR run)
    optional:
        genomeDir2 (full path to genome directory for filter round of mapping)
        genomeAnnots2 (full path to genome annotation file for filter mapping)
        optString2 (parameters for filter mapping STAR run)

    inputReads.fastq - a fastq file of reads

run as python3 pipelineWrapper9.py settings.txt inputReads.fastq outPrefix
"""

# import sys, common
import os, readCollapser4, filterJoshSAMByReadLength, thecountReads
import argparse
import assignReadsToGenesDF
import infoGraphQC
# import metaStartStop
from logJosh import Tee

# Absolute defaults are overwritten by the given settings file and any command line arguments given
ABSOLUTE_DEFAULT_DICT = {'cores': 7, 'misMatchMax': 0,
                         'optString': '--outFilterScoreMin 14'
                                      '--outFilterScoreMinOverLread 0.3'
                                      '--outFilterMatchNmin 14'
                                      '--outFilterMatchNminOverLread 0.3'
                                      '--outReadsUnmapped Fastx'
                                      '--outFilterMismatchNmax 0'  # Need an edit here to allow us to drop in the arg if passed
                                      '--outSJfilterOverhangMin 6 6 6 6'}


def main(fastqFile,settings,outPrefix,adaptorSeq,minimumReadLength,
    maximumReadLength,genomeDir,genomeAnnots,cores,misMatchMax,
    umi5,umi3,optString,**otherkwargs):
    
    #Uncomment the following for filter round of mapping:
    """
    genomeDir2=setList[10]
    genomeAnnots2=setList[11]
    optString2=setList[12]
    """
    ############################################################################################################
    """First set some parameters-- all of this can be deleted if we're good""" 
    ############################################################################################################
    #adaptorSeq='CTGTAGGCACCATCAAT'#adaptor for Komili loc1 dataset (looks same as rachel's adaptor)
    #adaptorSeq='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'#oJA126
    #adaptorSeq='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'#revCompl of forward adaptor
    #adaptorSeq='CACTCGGGCACCAAGGAC'#from BZ boris oligos
    #adaptorSeq='CTGTAGGCACCATCAATC'#Jonathan Gent
    #adaptorSeq='TGGAATTCTCGGGTGCCAAGG'#hendriks 3'adaptor for Ribo-seq time course
    #adaptorSeq='AGATGACGGGAAATAAAAGACACGTGCTGAAGTCAA'#another possibility for nextSeq adaptor from NEB
    #adaptorSeq='AGATCGGAAGAGCACACGTCTGAACTCC'#possibly the nextSeq adaptor from NEB
    # Atail
    #adaptorSeq='AAAAAAAAAAAAAAAAAA'
    #adaptorSeq='CTGTAGGCACCATCAA'#this is rachel's adaptor
    
    ############################################################################################################
    
    print(f'adaptorseq: {adaptorSeq}\n'
          f'minimumReadLength (not including UMI): {minimumReadLength}\n'
          f'maximumReadLength (not including UMI): {maximumReadLength}\n'
          f'5\' UMI length: {umi5}\n'
          f'3\' UMI length: {umi3}\n'
          f'genomeDir: {genomeDir}\n'
          f'genomeAnnots: {genomeAnnots}\n'
          f'cores To Use: {cores}\n'
          f'misMatchMax: {misMatchMax}\n'
          )
    
    ############################################################################################################
    """Trim adaptor from reads and sort by desired length"""
    ############################################################################################################
    # print('skipping trimming of any kind, so min/maximumReadLength restrictions ignored.')
    # os.system(f'cp {fastqFile} {outPrefix}.trimmed')
    print('read length restriction does not include %sNs and %sNs.'%(umi5,umi3))
    print('to accomodate UMI length, this program will add %snts to acceptable length.'%(umi5+umi3))
    
    os.system(f'cutadapt -a {adaptorSeq} '
              f'-m {minimumReadLength+umi5+umi3} '
              f'-M {maximumReadLength+umi5+umi3} '
              f'--too-short-output {outPrefix + ".trimmed.selfDestruct.tooShort.fastq"} '
              f'--too-long-output {outPrefix + ".trimmed.selfDestruct.tooLong.fastq"} '
              f'{fastqFile} > {outPrefix + ".trimmed.selfDestruct.fastq"} '
              f'2>/dev/null'
              )
    
    ############################################################################################################
    """Collapse reads and trim off UMIs"""
    ############################################################################################################
    ##it doesn't make sense to run readCollapser if there's no UMI
    if 0<umi5+umi3<=6:
        print(f'Your combined UMI length is {umi5+umi3}, which is pretty short.\
            I\'m going to try and collapse based on it, assuming you know what \
            you are doing. But if you do not understand this message, please \
            go find someone who can help you.')
    ##
    if umi5+umi3!=0:
        readCollapser4.main([outPrefix+'.trimmed.selfDestruct.fastq', 
                         umi5, umi3, 
                         outPrefix+'.trimmed.collapsed.selfDestruct.fastq'])
    else:
        print('Skipping collapsing...')
        ##the next line creates a symbolic link for the .collapsed file location
        ##instead of the prior way, which copied them.
        os.system(f'ln -s {outPrefix}.trimmed.selfDestruct.fastq '
                    f'{outPrefix}.trimmed.collapsed.selfDestruct.fastq')
    
    ############################################################################################################
    """Introduce a variable to make reading code easier"""
    ############################################################################################################
    readFile = f'{outPrefix}.trimmed.collapsed.selfDestruct.fastq'
    
    ############################################################################################################
    """Perform a filter round of mapping. e.g. to rDNA or RNAi trigger"""
    ############################################################################################################
    print('skipping filter round of mapping...')
    
    ##uncomment to do filter round of mapping
    """
    misMatchMax2=0
    print(f'performing filter round of mapping to {genomeDir2}')
    print(f'Only running on {cores} cores.')
    print(f'{misMatchMax2} mismatch max!')
    optString2= f'--outFilterScoreMin 14 ' \
        f'--outFilterScoreMinOverLread 0.3 ' \
        f'--outFilterMatchNmin 14 ' \
        f'--outFilterMatchNminOverLread 0.3 ' \
        f'--outReadsUnmapped Fastx ' \
        f'--outFilterMismatchNmax {misMatchMax2} ' \
        f'--outSJfilterOverhangMin 1000 1000 1000 1000 '
    print(f'Length/Score parameters: {optString2}')
    os.system(f'STAR {optString2} '
              f'--alignIntronMax 1 '
              f'--sjdbGTFfile {genomeAnnots2} '
              f'--genomeDir {genomeDir2} '
              f'--readFilesIn {readFile} '
              f'--runThreadN {cores} '
              f'--outFileNamePrefix {outPrefix}.trimmed.collapsed.mapped.filter'
              )
    #Now rewrite the read file to map from the unmapped reads
    #readFile=outPrefix+'.trimmed.collapsed.mapped.filterUnmapped.out.mate1'
    """
    
    ############################################################################################################
    """Commence read mapping"""
    ############################################################################################################
    print(f'Only running on {cores} cores.')
    print(f'{misMatchMax} mismatch max!')
    """
    optString= f'--outFilterScoreMin 14 ' \
        f'--outFilterScoreMinOverLread 0.3 ' \
        f'--outFilterMatchNmin 14 ' \
        f'--outFilterMatchNminOverLread 0.3 ' \
        f'--outReadsUnmapped Fastx ' \
        f'--outFilterMismatchNmax {misMatchMax} ' \
        f'--outSJfilterOverhangMin 6 6 6 6'
    """
    print(f'Length/Score parameters: {optString}')
    os.system(f'STAR {optString} '
              f'--alignIntronMax 1 '
              f'--sjdbGTFfile {genomeAnnots} '
              f'--genomeDir {genomeDir} '
              f'--readFilesIn {readFile} '
              f'--runThreadN {cores} '
              f'--outFileNamePrefix {outPrefix}.finalMapped.')
    
    # print(f'Printing file {outPrefix}Log.final.out')
    # os.system(f'lpr -p {outPrefix}Log.final.out')
    ############################################################################################################
    """Assign reads to genes"""
    ############################################################################################################
    print('Assigning reads to genes...')
    genomeAnnotProcessed=genomeAnnots.strip('gtf')+'allChrs.txt'
    # assignReadsToGenes4.main([genomeAnnotProcessed,
    #                          outPrefix+'.finalMapped.Aligned.out.sam',
    #                          outPrefix])
    assignReadsToGenesDF.main(outPrefix+'finalMapped.Aligned.out.sam',
                              genomeAnnotProcessed,
                              outPrefix,
                              keep_non_unique=False,
                              concatenate_output=True)
    print('Assigning reads to genes allowing for multiply-mapping reads...')
    # assignReadsToGenes5.main([genomeAnnotProcessed,
    #                          outPrefix+'.finalMapped.Aligned.out.sam',
    #                          outPrefix+'.redundantAndUnique'])
    assignReadsToGenesDF.main(outPrefix + 'finalMapped.Aligned.out.sam',
                              genomeAnnotProcessed,
                              outPrefix+'.redundantAndUnique',
                              keep_non_unique=True,
                              concatenate_output=True)
    print('Done with read assignment!')
    ############################################################################################################
    """Additional filtering of reads by length"""
    ############################################################################################################
    # print('Quitting early!!!'), sys.exit()
    
    # TODO: Definitely broken:
    print('filtering read lengths again...')
    filterJoshSAMByReadLength.main([outPrefix+'.jam',
                                minimumReadLength,
                                maximumReadLength,
                                outPrefix+'.filtered_%s-%snt.jam'%(minimumReadLength,maximumReadLength)])
    #print('Quitting early!!!'), sys.exit()
    
    ############################################################################################################
    """Creating infographic"""
    ############################################################################################################
    print('Making infographic')
    infoGraphQC.main([outPrefix+'.jam',minimumReadLength,maximumReadLength,-21,21,-30,12,outPrefix+'.qc'])
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
                        help='A text file containing a line-delimited input of adaptorSeq, genomeDir, and genomeAnnots.')
    parser.add_argument('outPrefix', metavar='outPrefix', type=str,
                        help='The outPrefix that all files produced will be named as.')
    # Optional Arguments: (None default here allows us to not pass anything that isn't given by user.
    #                      This helps to simplify settings.txt/default overwrites further down the line.)
    parser.add_argument('--umi5', metavar='umi5', type=int, default=None,
                        help='The number of nucleotides to be trimmed from the 5prime end of reads.')
    parser.add_argument('--umi3', metavar='umi3', type=int, default=None,
                        help='The number of nucleotides to be trimmed from the 3prime end of reads.')
    parser.add_argument('--minimumReadLength', '--min', metavar='min', type=int, default=None,
                        help='The shortest reads to be mapped to the genome.')
    parser.add_argument('--maximumReadLength', '--max', metavar='max', type=int, default=None,
                        help='The longest reads to be mapped to the genome.')
    parser.add_argument('--adaptorSeq', '--adaptor', metavar='adaptor', type=str, default=None,
                        help='3\' adaptor to be trimmed off.')
    parser.add_argument('--genomeDir', metavar='genomeDir', type=str, default=None,
                        help='Genome directory where STAR index can be found.')
    parser.add_argument('--genomeAnnots', metavar='genomeAnnots', type=str, default=None,
                        help='Genome annotations (gtf format).')
    parser.add_argument('--cores', metavar='cores', type=int, default=None,  # 7 cores as default!
                        help='Number of cores to use.')
    parser.add_argument('--misMatchMax', metavar='misMatchMax', type=int, default=None,  # 0 mismatchmax as default
                        help='Number of mismatches to tolerate during mapping.')
    # Flag Arguments: (just add these as tags to change pipeline functionality)
    parser.add_argument('-u', '--keep_non_unique', action='store_true',  # This is currently unused!
                        help="Boolean flag to allow non-unique reads in the final .jam file(s).")
    parser.add_argument('-p', '--print_arguments', action='store_true',
                        help="Boolean flag to show how arguments are overwritten/accepted")
    
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
                        print("Remove pipes ('|') from settings file arguments (or rewrite parser)")
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
    if finalArgDict['print_arguments']:
        print("\nFinal Arguments: (absolute defaults overwritten by settings.txt overwritten by CLI arguments)")
    for key, arg in finalArgDict.items():
        if finalArgDict['print_arguments']:
            print(f"\t{key} = {arg}")
        try:
            finalArgDict[key] = int(arg)
        except ValueError:
            finalArgDict[key] = str(arg)
    return finalArgDict

if __name__ == '__main__':
    Tee()
    arguments = combineSettingsAndArguments()
    main(**arguments)