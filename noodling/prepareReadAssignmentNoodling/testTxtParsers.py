"""
Joshua Arribere, March 26, 2020

Script to test how we're parsing information out of the .txt files of
    prepareReadAssignmentFiles and see what's most RAM-efficient.

Input: chr.txt - some chr.txt file as output from prepReadAssFile.

Output: print-to-screen of memory usage.

run as python testTxtParsers.py inFile.txt
"""
import sys
from logJosh import Tee
import os, psutil
import numpy

def parseInFile1(inFile):
    aa=[]
    with open(inFile,'r') as f:
        for line in f:
            line=line.strip()
            aa.append('\t'.join(line.split('\t')[1:]))
    return aa

def parseInFile2(inFile):
    #this approach considerably more expensive than method 1
    aa={}
    with open(inFile,'r') as f:
        for line in f:
            line=line.strip()
            key=int(line.split('\t')[0])
            aa[key]='\t'.join(line.split('\t')[1:])
    return aa

def parseInFile3(inFile):
    #prohibitively slow
    aa=numpy.array([])
    ct=0
    with open(inFile,'r') as f:
        for line in f:
            line=line.strip()
            key=int(line.split('\t')[0])
            aa=numpy.append(aa,'\t'.join(line.split('\t')[1:]))
            ct+=1
            if ct%10000==0:
                print 'Working on %s...'%(ct)
    return aa

def parseInFile4(inFile):
    #scrapped this approach when I realized that I have to know the size
    #of each element in the array
    sys.exit()
    bb=int(os.popen('wc -l %s'%(inFile)).read().split(  )[0])
    aa=numpy.array(['']*bb)
    with open(inFile,'r') as f:
        for line in f:
            line=line.strip()
            key=int(line.split('\t')[0])
            aa=numpy.append(aa,'\t'.join(line.split('\t')[1:]))
    return aa

def main(args):
    #first pass the input file name
    inFile=args[0]
    #
    #aa=parseInFile1(inFile)
    #del(aa)
    #aa=numpy.array(aa)#throws a memory error, even for small chrs
    #
    #aa=parseInFile2(inFile)
    #
    aa=parseInFile3(inFile)
    #
    #aa=parseInFile4(inFile)
    #
    process = psutil.Process(os.getpid())
    print(process.memory_info().rss)
    #this output matches up with what I see by top (gc has 32GB RAM)
    #1GB = 1e9 bytes


if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
