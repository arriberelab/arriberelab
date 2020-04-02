"""
Script to parse a line-delimited settings file for pipelineWrapper

Input: settings.txt - a line-delimited settings file in the format:
    adaptorSeq (raw sequence of adaptor)
    minReadLength (min length after adaptor and UMI trimming)
    maxReadLength (max length after adaptor and UMI trimming)
    UMI5 (5' UMI length in nts)
    UMI3 (3' UMI length in nts)
    genomeDir (full path to genome directory)
    genomeAnnots (full path to genome annotation file in gtf format)

Output: strings to feed into pipelineWrapper
"""

import sys

def parseSettingsFile(settings):
    with open (settings,'r') as s:
        setList=[]
        for line in s:
            line=line.strip()
            setList.append(line)
    adaptorSeq=setList[0]
    minimumReadLength=setList[1]
    maximumReadLength=setList[2]
    umi5=setList[3]
    umi3=setList[4]
    genomeDir=setList[5]
    genomeAnnots=setList[6]

def main(args):
    parsedSettings=parseSettingsFile(args)

if __name__=='__main__':
    main(sys.argv[1])
