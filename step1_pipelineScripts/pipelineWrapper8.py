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

run as python3 pipelineWrapper8.py settings.txt inputReads.fastq outPrefix
"""

import sys, common, os, assignReadsToGenes4, readCollapser4, filterJoshSAMByReadLength, thecountReads
#import metaStartStop
from logJosh import Tee

def main(args):
    settings, fastqFile, outPrefix = args[0:]
    
    with open (settings,'r') as s:
        setList=[]
        for line in s:
            line=line.strip()
            if not line.startswith("#"):
                setList.append(line)
    adaptorSeq=setList[0]
    minimumReadLength=setList[1]
    maximumReadLength=setList[2]
    umi5=setList[3]
    umi3=setList[4]
    genomeDir=setList[5]
    genomeAnnots=setList[6]
    cores=setList[7]
    misMatchMax=setList[8]
    optString=setList[9]
    
    #Uncomment the following for filter round of mapping:
    """
    genomeDir2=setList[10]
    genomeAnnots2=setList[11]
    optString2=setList[12]
    """
    
    minimumReadLength = int(minimumReadLength)
    maximumReadLength = int(maximumReadLength)
    
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
    
    #genomeDir='/data3/genomes/170622_yeastWithUTRs/'
    #genomeDir='/data1/genomes/161002_yeast/'
    #genomeDir='/data1/genomes/160212_Celegans_Ce10with_unc-54Degenerate/'
    #genomeDir='/home/josh/genomes/150226_Scer/'
    #genomeDir='/home/josh/genomes/150321_GFP_with_Ce_Chr/150325_PD8362/'
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
    
    #cores = 10  #groundcontrol has 16 cores total: cat /proc/cpuinfo | grep processor | wc -l
    #misMatchMax = 0
    
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
    umi5=int(umi5)
    umi3=int(umi3)
    print('read length restriction does not include %sNs and %sNs.'%(umi5,umi3))
    print('to accomodate UMI length, this program will add %snts to acceptable length.'%(umi5+umi3))
    
    os.system(f'cutadapt -a {adaptorSeq} '
              f'-m {minimumReadLength+umi5+umi3} '
              f'-M {maximumReadLength+umi5+umi3} '
              f'--too-short-output {outPrefix + ".trimmed.tooShort.fastq"} '
              f'--too-long-output {outPrefix + ".trimmed.tooLong.fastq"} '
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
    print('skipping filter round of mapping...')
    
    #Can delete this if we're good:
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
    #genomeDir2='/data8/genomes/181106_pMPmismatchFeeding/'
    #genomeAnnots2=genomeDir2+'blah.gtf'
    
    #Uncomment the following to commence filter round of mapping:
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
    #Can delete these/put in other settings files:
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
    #           f'--readFilesIn {readFile} '
    #           f'--runThreadN {cores} '
    #           f'--outFileNamePrefix {outPrefix}.finalMapped.')
    #optString = f'--outFilterMatchNmin 70 ' \
        #f'--outReadsUnmapped Fastx ' \
        #f'--outFilterMismatchNmax {misMatchMax} ' \
        #f'--outSJfilterOverhangMin 6 6 6 6'
    
    #Use the next optString for the normal pipeline
    #Can delete this bc it is now in the settings file:
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
    assignReadsToGenes4.main([genomeAnnotProcessed,
                             outPrefix+'.finalMapped.Aligned.out.sam',
                             outPrefix])
    
    ############################################################################################################
    """Additional filtering of reads by length"""
    ############################################################################################################
    # print('Quitting early!!!'), sys.exit()
    print('filtering read lengths again...')
    filterJoshSAMByReadLength.main([outPrefix+'.joshSAM',
                                minimumReadLength,
                                maximumReadLength,
                                outPrefix+'.joshSAM.filtered_%s-%snt'%(minimumReadLength,maximumReadLength)])
    #print('Quitting early!!!'), sys.exit()
    
    ############################################################################################################
    """Make a metagene plot of start/stop codon"""
    ############################################################################################################
    print('skipping the output metagene plot...')
    """
    print('Plotting metagene...')
    metaStartStop.main([genomeAnnots,
                        f'{outPrefix}.plot',
                        f'{outPrefix}.joshSAM.filtered_{minimumReadLength}-{maximumReadLength}nt',
                        'Library'])
	print(f'Done! {outPrefix}')
	
	############################################################################################################
	"""Run number and percentages of reads/riboinforgraphic"""
	############################################################################################################
	print('Calculating read percentages')
	thecountReads.main(fastq, outPrefix)
	"""
if __name__ == '__main__':
    Tee()
    main(sys.argv[1:])
