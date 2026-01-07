"""
Joshua Arribere, Sept 2, 2024

Script to plot output of DESeq2, simply highlighting genes that change
    (and that are up, probably).

Input: inFile.DESeqct - output from medianNormalizer.py

Output: MAplot

run as python3 plotDESeq2HighlightChanges.py inFile outPrefix
"""
import sys, common
from logJosh import Tee
from pyx import *

def parse(inFile):
    """
    Will parse output of DESeq2 to list of format:
    [(x,y,1/0),...] where 1/0 is:
        1 if padj<=0.01
        0 else
    """
    aa=[]
    pvalCut=0.01
    print('Flipping L2FC...')
    with open(inFile,'r') as f:
        f.readline()#skip header
        for line in f:
            line=line.strip()
            line=line.split(',')
            pval=line[-1]
            if pval!='NA':
                x=float(line[1])
                y=-float(line[2])
                temp=[x,y]
                if float(pval)<=pvalCut and y>=0:
                    temp.append(1)
                else:
                    temp.append(0)
                aa.append(temp)
    ##
    return aa

def mkPlot(inList,outPrefix):
    """
    inList=[(x,y,1/0),...]
    Will plot points as scatter, and then re-plot highlighting points w/
    1 in the 1/0 position.
    """
    g=graph.graphxy(width=4,height=2,
        x=graph.axis.log(min=10,max=10**5,title='Avg Expression'),
        y=graph.axis.linear(min=-2.5,max=5,title='Log2FC'))
    ##plot everything
    g.plot(graph.data.points(inList,x=1,y=2,title=None),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk.black,deco.filled],
            size=0.01)])
    ##highlight points
    temp=[entry for entry in inList if entry[-1]==1]
    g.plot(graph.data.points(temp,x=1,y=2,title=None),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk(0.85, 0.35, 0, 0),deco.filled],
            size=0.04)])
    ##
    g.writePDFfile(outPrefix)



def main(args):
    inFile,outPrefix=args[0:]
    ##
    inList=parse(inFile)
    ##
    mkPlot(inList,outPrefix)


if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
