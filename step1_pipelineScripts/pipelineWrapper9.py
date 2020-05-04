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

import sys, common, os, assignReadsToGenes4, readCollapser4, filterJoshSAMByReadLength, thecountReads
import argparse
import assignReadsToGenes5
import assignReadsToGenesDF
import infoGraphQC
#import metaStartStop
from logJosh import Tee

def parseSettings(settingsFile,adaptorSeq,minimumReadLength,\
        maximumReadLength,genomeDir,genomeAnnots,cores,\
        misMatchMax,umi5,umi3):
    """
    Will loop through and replace variables that are None
    """
    ##first, parse the settings file to a dictionary called settingsDict
    settingsDict={}
    with open(settingsFile,'r') as f:
        for line in f:
            if not line.startswith('#'):
                line=line.strip()
                if line!='':
                    line=line.split()
                    ##the next line is to allow for the optString
                    ##parameter. You could also do this by
                    #changing the delimiter. May revisit this.
                    settingsDict[line[0]]=' '.join(line[1:])
    ##now loop through the variables and update any that are not None
    if adaptorSeq==None:
        adaptorSeq=settingsDict['adaptorSeq']
    if minimumReadLength==None:
        minimumReadLength=int(settingsDict['minimumReadLength'])
    if maximumReadLength==None:
        maximumReadLength=int(settingsDict['maximumReadLength'])
    if genomeDir==None:
        genomeDir=settingsDict['genomeDir']
    if genomeAnnots==None:
        genomeAnnots=settingsDict['genomeAnnots']
    if umi5==None:
        umi5=int(settingsDict['umi5'])
    if umi3==None:
        umi3=int(settingsDict['umi3'])
    ##
    return adaptorSeq,minimumReadLength,maximumReadLength,\
        genomeDir,genomeAnnots,cores,misMatchMax,umi5,umi3,\
        settingsDict['optString']

def main(fastqFile,settings,outPrefix,adaptorSeq,minimumReadLength,
    maximumReadLength,genomeDir,genomeAnnots,cores,misMatchMax,
    umi5,umi3):
    ##parse the settings file w/ override from the command line
    adaptorSeq,minimumReadLength,maximumReadLength,genomeDir,\
        genomeAnnots,cores,misMatchMax,umi5,umi3,optString=\
        parseSettings(settings,adaptorSeq,minimumReadLength,\
        maximumReadLength,genomeDir,genomeAnnots,cores,\
        misMatchMax,umi5,umi3)
    ##
    
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
    Still learning how to do argparse.
    """
    parser=argparse.ArgumentParser(description='Wrapper for the handling of libraries starting from fastq files.')
    parser.add_argument('fastqFile', help='The fastq file input')
    parser.add_argument('settings', help='A text file containing a line-delimited input of adaptorSeq, genomeDir, and genomeAnnots.')
    parser.add_argument('outPrefix', type=str, help='The outPrefix that all files produced will be named as.')
    parser.add_argument('--umi5', type=int, default=None, help='The number of nucleotides to be trimmed from the 5prime end of reads.')
    parser.add_argument('--umi3', type=int, default=None, help='The number of nucleotides to be trimmed from the 3prime end of reads.')
    parser.add_argument('--minimumReadLength','--min', type=int, default=None, help='The shortest reads to be mapped to the genome.')
    parser.add_argument('--maximumReadLength','--max', type=int, default=None, help='The longest reads to be mapped to the genome.')
    parser.add_argument('--adaptorSeq','--adaptor', type=str, default=None, help='3\' adaptor to be trimmed off.')
    parser.add_argument('--genomeDir', type=str, default=None, help='Genome directory where STAR index can be found.')
    parser.add_argument('--genomeAnnots', type=str, default=None, help='Genome annotations (gtf format).')
    parser.add_argument('--cores', type=int, default=7, help='Number of cores to use.')
    parser.add_argument('--misMatchMax', type=int, default=0, help='Number of mismatches to tolerate during mapping.')
    return parser.parse_args()

if __name__ == '__main__':
    Tee()
    arguments=parseArguments()
    main(**vars(arguments))
