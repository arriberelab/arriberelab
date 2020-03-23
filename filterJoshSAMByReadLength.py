"""
Joshua Arribere June 26, 2017

Script to filter a joshSAM file by readlength

run as python filterJoshSAMByReadLength inFile.joshSAM min max outFile.joshSAM
"""
import sys, common
from logJosh import Tee

def main(args):
    inFile,minL,maxL,outPrefix=args[0:]
    minL,maxL=int(minL),int(maxL)
    with open(inFile,'r') as f:
        with open(outPrefix,'w') as g:
            for line in f:
                if line.startswith('@'):
                    g.write(line)
                else:
                    line2=line.strip().split('\t')
                    """
                    if len(line2[3])!=int(line2[4]):
                        print line
                    """
                    if minL<=len(line2[3])<=maxL:
                        g.write(line)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
