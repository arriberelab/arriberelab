"""
Joshua Arribere July 23, 2017

Script to filter .sam file by read lengths

Input - inFile.sam
    X - min length
    Y - max length

Output: a .sam file only with reads in length [X,Y] inclusive

run as python filterSAMByReadLength.py inFile.sam 15 18 outPrefix
"""
import sys, common
from logJosh import Tee
from assignReadsToGenes3 import recoverMappedPortion

def main(args):
    inFile,X,Y,outPrefix=args[0:]
    #
    X=int(X)
    Y=int(Y)
    #
    with open(inFile,'r') as f:
        with open(outPrefix+'.sam','w') as g:
            for line in f:
                if line.startswith('@'):
                    g.write(line)
                else:
                    line2=line.strip().split('\t')
                    Cigar=line2[5]
                    Read=line2[9]
                    #necessary to remove soft clipping for calculated
                    #read lengths
                    mappedRead,N=recoverMappedPortion(Cigar,Read)
                    if X<=len(mappedRead)<=Y:
                        g.write(line)
    

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
