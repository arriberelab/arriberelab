"""
Script to parse a line-delimited settings file for pipelineWrapper

Input: settings.txt - a line-delimited settings file in the format:
    adaptorSeq=raw sequence of adaptor
    minReadLength=min length after adaptor and UMI trimming
    maxReadLength=max length after adaptor and UMI trimming
    UMI5=5' UMI length in nts
    UMI3=3' UMI length in nts
    genomeDir=full path to genome directory
    genomeAnnots=full path to genome annotation file in gtf format

Output: a dict of strings to feed into pipelineWrapper
"""

import sys

def parseSettingsFile(settings):
    with open (settings,'r') as s:
        setDict={}
        for line in s:
            if not line.startswith("#"):
                line=line.strip().split("=")
                for item in line:
                    setDict[line[0]]=line[1]
        print(setDict)

def main(args):
    parsedSettings=parseSettingsFile(args)

if __name__=='__main__':
    main(sys.argv[1])

