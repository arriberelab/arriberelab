#!/usr/bin/env python3

"""
Joshua Arribere, July 25, 2014
Converted to python 3: Mar 23, 2020

Parissa edited to make use of UMI lengths given in command line and revised readCollapser.

run as python3 pipelineWrapper8.py inputReads.fastq outPrefix --options
"""

import sys, common, os, assignReadsToGenes4, metaStartStop, readCollapser4, filterJoshSAMByReadLength, argparse
from logJosh import Tee
parser = argparse.ArgumentParser(description='Wrapper for the handling of libraries starting from fastq files.')
parser.add_argument('fastqFile', help='The fastq file input')
parser.add_argument('settings', help='A text file containing a line-delimited input of adaptorSeq, genomeDir, and genomeAnnots.')
parser.add_argument('outPrefix', type=str, help='The outPrefix that all files produced will be named as.')
parser.add_argument('--umi5', type=int, default=0, help='The number of nucleotides to be trimmed from the 5prime end of reads.')
parser.add_argument('--umi3', type=int, default=6, help='The number of nucleotides to be trimmed from the 3prime end of reads.')
parser.add_argument('--minimumReadLength','--min', type=int, default=15, help='The shortest reads to be mapped to the genome.')
parser.add_argument('--maximumReadLength','--max', type=int, default=30, help='The longest reads to be mapped to the genome.')
args = parser.parse_args()
def parser(fastqFile,settings,outPrefix,umi5,umi3,minimumReadLength,maximumReadLength):
    with open(settings, 'r') as s:
        p=[]
        for line in s:
            line=line.strip()
            if not line.startswith("#"):
                p.append(line)
        adaptorSeq=p[0]
        genomeDir=p[1]
        genomeAnnots=p[2]
    main(fastqFile,outPrefix,umi5,umi3,minimumReadLength,maximumReadLength,adaptorSeq,genomeDir,genomeAnnots)
def main(fastqFile,outPrefix,umi5,umi3,minimumReadLength,maximumReadLength,adaptorSeq,genomeDir,genomeAnnots):
    
    #fastqFile, minimumReadLength, maximumReadLength, umi5, umi3, outPrefix = args[0:]
    
    ############################################################################################################
    """First set some parameters"""
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
    
    #minimumReadLength = int(minimumReadLength)
    #maximumReadLength = int(maximumReadLength)
    
    #genomeDir='/data3/genomes/170622_yeastWithUTRs/'
    #genomeDir='/data1/genomes/161002_yeast/'
    #genomeDir='/data1/genomes/160212_Celegans_Ce10with_unc-54Degenerate/'
    #genomeDir='/home/josh/genomes/150226_Scer/'
    #genomeDir='/homAGATCGGAAGAGCACACGTCTGAACTCCAGTCACe/josh/genomes/150321_GFP_with_Ce_Chr/150325_PD8362/'
    #genomeDir='/home/josh/genomes/131217_Ty1/140127_Ty1/'
    #genomeDir='/home/josh/genomes/150518_Celegans/'
    #genomeDir='/data3/genomes/170620_unc-54GFPwithT2A/'
    #genomeDir='/data3/genomes/170719_unc-54e1092/'
    #genomeDir='/data4/genomes/171207_Celegans_release90/'
    #genomeDir='/data5/marissa/180411_genomes/180417_ensembl/'
    #genomeDir='/data4/genomes/171218_historicalGenomeT2A/'
    #genomeDir='/data1/genomes/170322_genomeWithUnc-54GFPNonStopDegenerate2/'
    #genomeDir='/data1/genomes/160110_Celegans_rel83/'
    #genomeDir='/data7/180330_backups/160520_joshGallifreyBackup/data1/genomes/161002_yeast/'
    #genomeDir='/data12/joshua/genomes/180402_clone_160110_Celegans_rel83/160110_Celegans_rel83/'
    #genomeDir='/data12/joshua/genomes/171218_historicalGenomeT2A/'
    #genomeDir='/data12/joshua/genomes/191125_srf0788FromParissa/'
    #genomeDir='/data15/joshua/genomes/200329_cerevisiae/'
    
    #genomeAnnots='/home/josh/genomes/131217_Ty1/140127_Ty1/131217_M18706_revised.gtf'
    #genomeAnnots='/home/josh/working/141117_working_newGTF/Caenorhabditis_elegans.WBcel215.70.sansBJA7_40_77.gtf'
    #genomeAnnots='/data3/genomes/170622_yeastWithUTRs/170622_Saccharomyces_cerevisiae.R64-1-1.85.gtf'
    #genomeAnnots='/data1/genomes/161002_yeast/Saccharomyces_cerevisiae.R64-1-1.85.gtf'
    #genomeAnnots='/data1/genomes/160212_Celegans_Ce10with_unc-54Degenerate/Caenorhabditis_elegans.WBcel235.83.gtf'
    #genomeAnnots='/home/josh/genomes/150226_Scer/Saccharomyces_cerevisiae.R64-1-1.78.gtf'
    #genomeAnnots='/home/josh/genomes/150321_GFP_with_Ce_Chr/150325_WBcel215.70_chrPD2874-2876and8362.gtf'
    #genomeAnnots='/data1/genomes/160110_Celegans_rel83/Caenorhabditis_elegans.WBcel235.83.gtf'
    #genomeAnnots='/data1/genomes/170322_genomeWithUnc-54GFPNonStopDegenerate2/170322_genomeWithExtraUnc-54Chr.gtf'
    #genomeAnnots='/data3/genomes/170620_unc-54GFPwithT2A/170612_genomeWithUnc-54ChrAsItAppearsInPD4092.gtf'
    #genomeAnnots='/data3/genomes/170719_unc-54e1092/Caenorhabditis_elegans.WBcel235.83.gtf'
    #genomeAnnots='/data4/genomes/171207_Celegans_release90/Caenorhabditis_elegans.WBcel235.90.gtf'
    #genomeAnnots='/data5/marissa/180411_genomes/180417_ensembl/Schizosaccharomyces_pombe.ASM294v2.22.gtf'
    #genomeAnnots='/data4/genomes/171218_historicalGenomeT2A/170612_genomeWithUnc-54ChrAsItAppearsInPD4092.gtf'
    #genomeAnnots='/data7/180330_backups/160520_joshGallifreyBackup/data1/genomes/161002_yeast/Saccharomyces_cerevisiae.R64-1-1.85.gtf'
    #genomeAnnots='/data12/joshua/genomes/180402_clone_160110_Celegans_rel83/160110_Celegans_rel83/Caenorhabditis_elegans.WBcel235.83.gtf'
    #genomeAnnots='/data12/joshua/genomes/171218_historicalGenomeT2A/170612_genomeWithUnc-54ChrAsItAppearsInPD4092.gtf'
    #genomeAnnots='/data12/joshua/genomes/191125_srf0788FromParissa/191122_genomeWithUnc-54AsItAppearsInWJA0788.gtf'
    #genomeAnnots='/data15/joshua/genomes/200329_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.99.gtf'
    
    cores = 10  #groundcontrol has 16 cores total: cat /proc/cpuinfo | grep processor | wc -l
    misMatchMax = 0
    
    print(f'adaptorSeq: {adaptorSeq}\n'
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
    umi5=int(umi5)
    umi3=int(umi3)
    print('read length restriction does not include %sNs and %sNs.'%(umi5,umi3))
    print('to accomodate UMI length, this program will add %snts to acceptable length.'%(umi5+umi3))
    #print('read length restriction does not include 6Ns. This program will NOT add 6')
    #print('read length restriction does not include 6Ns and 4Ns. This program WILL add 10')
    #print('read length restriction does not include 6Ns. This program WILL add 6')
    
    os.system(f'cutadapt -a {adaptorSeq} '
              f'-m {minimumReadLength+umi5+umi3} '
              f'-M {maximumReadLength+umi5+umi3} '
              f'--too-short-output {outPrefix + ".trimmed.tooShort"} '
              f'--too-long-output {outPrefix + ".trimmed.tooLong"} '
              f'{fastqFile} > {outPrefix + ".trimmed"} '
              f'2>/dev/null'
              )
    
    ############################################################################################################
    """Collapse reads and trim off UMIs"""
    ############################################################################################################
    readCollapser4.main([outPrefix+'.trimmed', 
                         umi5, umi3, 
                         outPrefix+'.trimmed.collapsed.fastq'])
    ##print('skipping collapsing...')
    #os.system(f'cp {outPrefix}.trimmed '
    #          f'{outPrefix}.trimmed.collapsed.fastq')
    
    ############################################################################################################
    """Introduce a variable to make reading code easier"""
    ############################################################################################################
    readFile = f'{outPrefix}.trimmed.collapsed.fastq'
    
    ############################################################################################################
    """Perform a filter round of mapping. e.g. to rDNA or RNAi trigger"""
    ############################################################################################################
    
    #genomeDir2='/data1/genomes/170320_unc-54GFPNonStop/'
    #genomeAnnots2='/data1/genomes/170320_unc-54GFPNonStop/170320_unc-54GFPNonStop.gtf'
    #genomeDir2='/data1/genomes/170320_unc-54BJA40TriggerChr/'
    #genomeAnnots2='/data1/genomes/170320_unc-54BJA40TriggerChr/170320_dummyGTF.gtf'
    #genomeDir2='/data1/genomes/161031_Ty1/'
    #genomeAnnots2='/data1/genomes/161031_Ty1/131217_M18706_revised.gtf'
    #genomeDir2='/home/josh/genomes/150519_triggerChr4/'
    #genomeDir2='/home/josh/genomes/150608_triggerChr5_onlyunc-22andunc-54/'
    #genomeDir2='/data3/genomes/170626_rDNAcerevisiae/'
    #genomeAnnots2='/data3/genomes/170626_rDNAcerevisiae/blah.gtf'
    #genomeDir2='/data4/genomes/171207_BJA40BJA7chr/'
    #genomeDir2='/data4/genomes/180103_BJA7_40_77_chr/'
    #genomeAnnots2=genomeDir2+'blah.gtf'
    #genomeDir2='/data4/genomes/180331_triggerChr_pJA7_pJA40_pJA77/'
    #genomeDir2='/data4/genomes/180514_triggerChr_pJA7_pJA77_pJA151_pJA153/'
    #genomeAnnots2=genomeDir2+'blah.gtf'
    """
    genomeDir2='/data8/genomes/181106_pMPmismatchFeeding/'
    genomeAnnots2=genomeDir2+'blah.gtf'
    misMatchMax2=0
    print(f'performing filter round of mapping to {genomeDir2}')
    print(f'Only running on {cores} cores.')
    print(f'{misMatchMax2} mismatch max!')
    optString= f'--outFilterScoreMin 14 ' \
        f'--outFilterScoreMinOverLread 0.3 ' \
        f'--outFilterMatchNmin 14 ' \
        f'--outFilterMatchNminOverLADAPTORSEQread 0.3 ' \
        f'--outReadsUnmapped Fastx ' \
        f'--outFilterMismatchNmax {misMatchMax2} ' \
        f'--outSJfilterOverhangMin 1000 1000 1000 1000 '
    print(f'Length/Score parameters: {optString}')
    os.system(f'STAR {optString} '
              f'--alignIntronMax 1 '
              f'--sjdbGTFfile {genomeAnnots2} '
              f'--genomeDir {genomeDir2} '
              f'--readFilesIn {readFile} '
              f'--runThreadN {cores} '
              f'--outFileNamePrefix {outPrefix}.trimmed.collapsed.mapped.filter'
              )
    """
    #"""Now rewrite the read file to map from the unmapped reads"""
    print('skipping filter round of mapping...')
    #readFile=outPrefix+'.trimmed.collapsed.mapped.filterUnmapped.out.mate1'
    
    ############################################################################################################
    """Commence read mapping"""
    ############################################################################################################
    print(f'Only running on {cores} cores.')
    print(f'{misMatchMax} mismatch max!')
    # optString= f'--outFilterScoreMin 14 ' \
    #     f'--outFilterScoreMinOverLread 0.3 ' \
    #     f'--outFilterMatchNmin 14 ' \
    #     f'--outFilterMatchNminOverLread 0.3 ' \
    #     f'--outFilterMismatchNmax {misMatchMax} ' \
    #     f'--outSJfilterOverhangMin 20 10 10 10'
    # os.system(f'STAR {optString} '
    #           f'--alignIntronMax 1 '
    #           f'--sjdbGTFfile {genomeAnnots} '
    #           f'--genomeDir {genomeDir} '
    #     AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC      f'--readFilesIn {readFile} '
    #           f'--runThreadN {cores} '
    #           f'--outFileNamePrefix {outPrefix}.finalMapped.')
    #optString = f'--outFilterMatchNmin 70 ' \
        #f'--outReadsUnmapped Fastx ' \
        #f'--outFilterMismatchNmax {misMatchMax} ' \
        #f'--outSJfilterOverhangMin 6 6 6 6'
    
    #Use the next optString for the normal pipeline
    optString= f'--outFilterScoreMin 14 ' \
        f'--outFilterScoreMinOverLread 0.3 ' \
        f'--outFilterMatchNmin 14 ' \
        f'--outFilterMatchNminOverLread 0.3 ' \
        f'--outReadsUnmapped Fastx ' \
        f'--outFilterMismatchNmax {misMatchMax} ' \
        f'--outSJfilterOverhangMin 6 6 6 6'
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
    
    ############################################################################################################
    """Additional filtering of reads by length"""
    ############################################################################################################
    # print('Quitting early!!!', sys.exit())
    print('Filtering read lengths again...')
    filterJoshSAMByReadLength.main([outPrefix+'.joshSAM',
                                minimumReadLength,
                                maximumReadLength,
                                outPrefix+'.joshSAM.filtered_%s-%snt'%(minimumReadLength,maximumReadLength)])
    print('Quitting early!!!', sys.exit())
    
    ############################################################################################################
    """Make a metagene plot of start/stop codon"""
    ############################################################################################################
    # print('skipping the output metagene plot...')
    print('Plotting metagene...')
    metaStartStop.main([genomeAnnots,
                        f'{outPrefix}.plot',
                        f'{outPrefix}.joshSAM.filtered_{minimumReadLength}-{maximumReadLength}nt',
                        'Library'])
    print(f'Done! {outPrefix}')

if __name__ == '__main__':
    Tee()
    parser(**vars(args))
