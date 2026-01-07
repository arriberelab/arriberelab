"""
Joshua Arribere, Jan 22, 2025

Script to calculate and plot the polarity scores for jam files in elegans.

Input: inFiles.jam - on the command line
    cds.fa - for CDS lengths
    outPrefix

Output: Will plot the polarity scores of all genes

run as python3 polarityScore.py inFile1.jam ... inFileN.jam 
    cds.fa outPrefix
"""
import sys, common, collections
from logJosh import Tee
from pyx import *

def parse(inFiles):
    """
    Will parse a list of jam files to a list of dictionaries, each of
    the format:
    {wbgene:{position:ct}}
    Will only consider positions that are uniquely assignable to a single
    position within all transcripts that that read could correspond to.
    """
    theMin=28
    theMax=31
    print('Restricting to read lengths in [%s,%s].'%(theMin,theMax))
    aa=[]
    for inFile in inFiles:
        temp=collections.defaultdict(lambda:collections.defaultdict(int))
        with open(inFile,'r') as f:
            f.readline()#skip a line, the header
            for line in f:
                ##
                line=line.strip().split('\t')
                if line[-2].endswith(':S') and theMin<=len(line[6])<=theMax:
                    wbgene=line[-2].split(':')[0]
                    ##
                    txts=line[-1].split('|')
                    startPositions=list(set([int(txt.split(':')[1])+12 for txt in txts]))##P-site offset of 12
                    stopPositions=list(set([int(txt.split(':')[2])+12 for txt in txts]))##P-site offset of 12
                    if len(startPositions)==1 and len(stopPositions)==1:
                        if startPositions[0]>=15 and stopPositions[0]<=-15:
                            position=startPositions[0]
                            temp[wbgene][position]+=1
        aa.append(temp)
    return aa

def restrict(inList,N):
    """
    inList is a list of dicts, each dict of the format:
    {wbgene:{position:ct}}.
    Will ID wbgenes with a total of at least N cts in all list entries.
    Will return a list of dicts with those genes.
    """
    aa=collections.defaultdict(list)
    for inDict in inList:
        for wbgene,temp in inDict.items():
            aa[wbgene].append(sum(temp.values()))
    bb=[]
    for k,v in aa.items():
        if len(v)==len(inList) and min(v)>=N:
            bb.append(k)
    cc=[]
    for inDict in inList:
        tempDict={}
        for k,v in inDict.items():
            if k in bb:
                tempDict[k]=inDict[k]
        cc.append(tempDict)
    return cc

def parseCDS(cdsDict):
    """
    Will return a dict of {wbgene:length} for all wbgenes
    that have a single length
    """
    aa=collections.defaultdict(list)
    for k,v in cdsDict.items():
        wbgene=k.split('gene:')[1].split('.')[0]
        aa[wbgene].append(len(v))
    bb={}
    for k,v in aa.items():
        if len(list(set(v)))==1:
            bb[k]=v[0]
    return bb

def getPolScore(i,l,d,totalCt):
    wi=(2*i-l-1)/(l-1)
    pi=d*wi/totalCt
    return pi

def getPolScores(inList2,cdsDict):
    """
    inList2 is a list of dicts, each dict of the format:
    {wbgene:{position:ct}}
    cdsDict is of the format: {wbgene:length}
    Will compute the polarity score for each wbgene, and
    return a list where each entry is a dict of format
    {wbgene:pol_score}
    """
    aa=[]
    for inDict in inList2:
        temp=collections.defaultdict(list)
        for wbgene,positionDict in inDict.items():
            if wbgene in cdsDict:
                l=cdsDict[wbgene]
                totalCt=sum(positionDict.values())
                for i,ct in positionDict.items():
                    temp[wbgene].append(getPolScore(i,l,ct,totalCt))
                temp[wbgene]=sum(temp[wbgene])
        aa.append(temp)
    return aa

def mkCDF(list1):
    list1.sort()
    k=len(list1)
    aa=[]
    for ii in range(k):
        aa.append((list1[ii],(ii+1)/k))
    return aa

def mkPlot(someTuples,outPrefix):
    g=graph.graphxy(width=4,height=4,key=graph.key.key(pos='tr',hinside=0),
        x=graph.axis.linear(min=-0.25,max=0.25,
            title='Polarity Score'),
        y=graph.axis.linear(title='CDF'))
    ##
    ii=0
    for someTuple in someTuples:
        theKey=someTuple[0].replace('_','-')
        theData=someTuple[1]
        for k,v in theData.items():
            if v==0:
                print(k,v)
        theVals=list(theData.values())
        theData=mkCDF(theVals)
        ##
        g.plot(graph.data.points(theData,x=1,y=2,title=theKey),
            [graph.style.line([style.linestyle.solid,
                common.colors(ii)])])
        ii+=1
    ##
    g.writePDFfile(outPrefix)

def main(args):
    inFiles=args[:-2]
    cdsFile=args[-2]
    outPrefix=args[-1]
    ##
    inList=parse(inFiles)
    ##
    N=50#min read ct
    print('Restricting to genes with at least %s cts in all libs...'%(N))
    inList2=restrict(inList,N)
    ##
    cdsDict=common.parseFasta(cdsFile)
    cdsDict=parseCDS(cdsDict)
    ##
    polarityScoresPerLib=getPolScores(inList2,cdsDict)
    ##
    mkPlot(zip(inFiles,polarityScoresPerLib),outPrefix)

if __name__ == '__main__':
    Tee()
    main(sys.argv[1:])
