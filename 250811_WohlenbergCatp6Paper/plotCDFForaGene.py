"""
Joshua Arribere, Jan 18, 2024

Script to plot a eCDF of reads across a txtID's CDS.

Input: inFiles.txt - tab-delimited of the format
    name1\tinFile1.jam
    name2\tinFile2.jam
    ...
    features.txt - text file containing the beginning
        and end of interesting features, along with
        names. For example
        polyPro\t500\t509
        ^this will draw lines at x=500,x=509.
    txtID - the txtID you wish to plot reads of.

Output: pdf

run as python3 plotCDFForAGeneAndMultipleJams.py 
    inFiles.txt features.txt txtID outPrefix
"""
import sys, common, collections, pickle
from logJosh import Tee
from pyx import *

def getPositions(jamFile,txtID):
    """
    will look for lines in jamFile with txtID, and 
    then extract their position, create a dict of
    {position:ct}. Will then convert that dict to
    a CDF
    """
    print('Offsetting for P-site...')
    aa=collections.defaultdict(int)
    with open(jamFile,'r') as f:
        for line in f:
            if txtID in line:
                line=line.strip().split('\t')
                ##
                if ':S' in line[8]:##check for sense
                    for entry in line[9:]:
                        if txtID in entry:
                            posit=int(entry.split(':')[1])
                            aa[posit+12]+=1
                            break#one count per read
                            #even though each txt should
                            #only be present once per
                            #line
    ##
    total=sum(aa.values())
    bb=[k for k in aa.keys()]
    bb.sort()
    cc=[]
    tally=0
    for k in bb:
        tally+=aa[k]
        cc.append((k,tally/total))
    ##
    return cc

def parseReads(inFile,txtID):
    """
    Will input a text file w/ names and jamFile names
    and count reads across txtID
    """
    aa={}
    ##
    with open(inFile,'r') as f:
        for line in f:
            name,theFile=line.strip().split('\t')
            ##
            aa[name]=getPositions(theFile,txtID)
    ##
    return aa

def mkPlot(reads,featureFile,outPrefix):
    ##
    g=graph.graphxy(width=8,height=3,
        key=graph.key.key(pos='br',hinside=0),
        x=graph.axis.linear(title='Position (nt)'),
        y=graph.axis.linear(min=0,max=1,title='CDF'))
    ##
    jj=0
    for name,data in reads.items():
        print(name)
        g.plot(graph.data.points(data,x=1,y=2,
            title=name),
            [graph.style.line([style.linestyle.solid,
                common.colors(jj)])])
        jj+=1
    ##
    with open(featureFile,'r') as f:
        for line in f:
            name,left,right=line.strip().split('\t')
            g.plot(graph.data.points(
                [(left,0),(left,1)],x=1,y=2,
                title=name),
                [graph.style.line([style.linestyle.solid,
                    color.cmyk.black])])
            g.plot(graph.data.points(
                [(right,0),(right,1)],x=1,y=2,
                title=None),
                [graph.style.line([style.linestyle.solid,
                    color.cmyk.black])])
    ##
    ##
    g.writePDFfile(outPrefix)

def main(args):
    inFile,featureFile,txtID,outPrefix=args[0:]
    ##
    reads=parseReads(inFile,txtID)
    ##
    common.rePickle(reads,outPrefix+'%s.p'%(txtID))
    with open(outPrefix+'%s.p'%(txtID),'rb') as f:
        reads=pickle.load(f)
    ##
    mkPlot(reads,featureFile,outPrefix+'%s'%(txtID))

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
