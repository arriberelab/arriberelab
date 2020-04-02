"""
Script to parse a line-delimited settings file for pipelineWrapper

Input: settings.txt - a line-delimited settings file in the format:
    header WITHOUT info
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
    adaptorSeq=setList[1]
    minimumReadLength=setList[2]
    maximumReadLength=setList[3]
    umi5=setList[4]
    umi3=setList[5]
    genomeDir=setList[6]
    genomeAnnots=setList[7]

def main(args):
    parsedSettings=parseSettingsFile(args)

if __name__=='__main__':
    main(sys.argv[1])
