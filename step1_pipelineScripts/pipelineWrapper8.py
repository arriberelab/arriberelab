#!/usr/bin/env python3
"""
Joshua Arribere, July 25, 2014
Mar 23, 2020: Converted to python 3
April 2, 2020: Now accepts settings via a txt file

Input:
    inputReads.fastq - fastq file of reads
    settings.txt - a line-delimited file with labels of the format:
        adaptorSeq sequence
        minimumReadLength length
        maximumReadLength length
        genomeDir /path/to/genome/directory/
        gneomeAnnots /path/to/genome/directory/annots.gtf
        umi5 length
        umi3 length
        optString STARparameters
    outPrefix - prefix for all output files

run as python3 pipelineWrapper8.py inputReads.fastq settings.txt [options] outPrefix
"""

import sys, common, os, argparse, assignReadsToGenes4, assignReadsToGenes5
import readCollapser4, filterJoshSAMByReadLength, thecountReads
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
    
    ############################################################################################################
    """Some common adaptor sequences to use in settings file"""
    ############################################################################################################
    #adaptorSeq='CTGTAGGCACCATCAAT'#adaptor for Komili loc1 dataset (looks same as rachel's adaptor)
    #adaptorSeq='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'#oJA126
    #adaptorSeq='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'#revCompl of forward adaptor
    #adaptorSeq='CACTCGGGCACCAAGGAC'#from BZ boris oligos
    #adaptorSeq='CTGTAGGCACCATCAATC'#Jonathan Gent
    #adaptorSeq='TGGAATTCTCGGGTGCCAAGG'#hendriks 3'adaptor for Ribo-seq time course
    #adaptorSeq='AGATGACGGGAAATAAAAGACACGTGCTGAAGTCAA'#another possibility for nextSeq adaptor from NEB
    #adaptorSeq='AGATCGGAAGAGCACACGTCTGAACTCC'#possibly the nextSeq adaptor from NEB
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
    genomeDir2 = #path to filter genome dir
    genomeAnnots2 = #path to filter genome annots
    
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
    readFile=outPrefix+'.trimmed.collapsed.mapped.filterUnmapped.out.mate1'
    """
    
    ############################################################################################################
    """Commence read mapping"""
    ############################################################################################################
    print(f'Only running on {cores} cores.')
    print(f'{misMatchMax} mismatch max!')
    """
    #optString for normal pipeline:
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
    assignReadsToGenes4.main([genomeAnnotProcessed,
                             outPrefix+'.finalMapped.Aligned.out.sam',
                             outPrefix])
    print('Assigning reads to genes allowing for multiply-mapping reads...')
    assignReadsToGenes5.main([genomeAnnotProcessed,
                             outPrefix+'.finalMapped.Aligned.out.sam',
                             outPrefix+'.redundantAndUnique'])
    print('Done with read assignment!')
    ############################################################################################################
    """Additional filtering of reads by length"""
    ############################################################################################################
    # print('Quitting early!!!'), sys.exit()
    print('filtering read lengths again...')
    filterJoshSAMByReadLength.main([outPrefix+'.jam',
                                minimumReadLength,
                                maximumReadLength,
                                outPrefix+'.filtered_%s-%snt.jam'%(minimumReadLength,maximumReadLength)])
    #print('Quitting early!!!'), sys.exit()
    
    ############################################################################################################
    """Create .bam and .bai files"""
    ############################################################################################################
    picardPath = '/home/annie/programs/picard.jar'
    print('Making BAM file')
    #converting to .bam file
    os.system('samtools view -S -b %s.finalMapped.Aligned.out.sam > %s.finalMapped.Aligned.out.bam' % (outPrefix, outPrefix))
    #sort
    os.system('samtools sort %s.finalMapped.Aligned.out.bam -o %s.finalMapped.Aligned.out.sorted.bam' % (outPrefix, outPrefix))
    #index .bam file
    print('Indexing BAM file')
    os.system('samtools index %s.finalMapped.Aligned.out.sorted.bam' % (outPrefix))
    
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
    parser=argparse.ArgumentParser(description='Wrapper for handling libraries starting from fastq files.')
    parser.add_argument('fastqFile', help='Input reads in fastq format.')
    parser.add_argument('settings', help='Line-delimited text file with adaptorSeq, genomeDir, genomeAnnots, 5\' UMI, and 3\' UMI.')
    parser.add_argument('outPrefix', type=str, help='Prefix to name all output files.')
    parser.add_argument('--umi5', type=int, default=None, help='Number of nts to trim from the 5\' end of reads.')
    parser.add_argument('--umi3', type=int, default=None, help='Number of nts to trim from 3\' end of reads.')
    parser.add_argument('--minimumReadLength','--min', type=int, default=None, help='Shortest reads to map to the genome.')
    parser.add_argument('--maximumReadLength','--max', type=int, default=None, help='Longest reads to map to the genome.')
    parser.add_argument('--adaptorSeq','--adaptor', type=str, default=None, help='3\' adaptor to trim off.')
    parser.add_argument('--genomeDir', type=str, default=None, help='Genome directory where STAR index can be found.')
    parser.add_argument('--genomeAnnots', type=str, default=None, help='Genome annotations in gtf format.')
    parser.add_argument('--cores', type=int, default=7, help='Number of cores to use.')
    parser.add_argument('--misMatchMax', type=int, default=0, help='Number of mismatches to tolerate during mapping.')
    return parser.parse_args()

if __name__ == '__main__':
    Tee()
    arguments=parseArguments()
    main(**vars(arguments))
