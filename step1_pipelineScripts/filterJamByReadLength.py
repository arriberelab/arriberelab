"""
Joshua Arribere May 4, 2020

Script to filter a jam file by readlength

run as python filterJamByReadLength inFile.jam min max outFile.jam
"""
import sys, common
from logJosh import Tee

def main(args):
    inFile,minL,maxL,outPrefix=args[0:]
    minL,maxL=int(minL),int(maxL)
    ct=0
    with open(inFile,'r') as f:
        with open(outPrefix,'w') as g:
            for line in f:
                if ct==0:
                    g.write(line)
                    ct+=1
                else:
                    line2=line.strip().split('\t')
                    if minL<=len(line2[6])<=maxL:
                        g.write(line)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
