"""
Joshua Arribere, Mar 21, 2024

Script to plot the distribution of a bunch of Ribo-seq reads about AA
    motifs.

Input: seqs.fa - fasta-formatted transcript sequences of CDSes only
    inFiles.txt - tab-delimited file containing name and location of
        jamFiles to plot, like this:
            riboSeq1\twhereTheFileIs.jam
    N - how long of a window to examine, in AAs

Output: For ii in [0-N], will plot the meta Ribo-seq read distribution
    about positions with at least ii of a given AA in the next N AAs.
    Will be a ton of plots.

run as python3 polyAARuns.py seqs.fa inFiles.txt N outPrefix
"""
import sys, common, collections, numpy
from logJosh import Tee
from pyx import *
print('This is a RAM-hungry script. It will use tens of Gb for C.elegans.')

def convertToProteins(inDict):
    aa={}
    for k,v in inDict.items():
        k2=k.strip().split()[0]
        v2=common.translate(v)[1]
        aa[k2]=v2
    return aa

def getRuns(protSeqs,N):
    aa={}
    AAs=['A','C','D','E','F','G','H','I','K','L',
        'M','N','P','Q','R','S','T','V','W','Y']
    for k,v in protSeqs.items():
        aa[k]={}
        for ii in range(len(v)-N):
            aa[k][ii]={}
            temp=v[ii:ii+10]
            for AA in AAs:
                aa[k][ii][AA]=temp.count(AA)
    return aa

def getPosition(txtList):
    starts=list(set([int(entry.split(':')[1]) for entry in txtList]))
    stops=list(set([int(entry.split(':')[2]) for entry in txtList]))
    if len(starts)==1 and len(stops)==1:
        return starts[0]
    else:
        return 'na'


def parseReads(jamFile):
    aa=collections.defaultdict(lambda:collections.defaultdict(int))
    with open(jamFile,'r') as f:
        f.readline()#skips header
        for line in f:
            line=line.strip().split('\t')
            if line[8].endswith(':S'):
                txtList=line[-1].split('|')
                position=getPosition(txtList)
                if position!='na':
                    for entry in txtList:
                        aa[entry.split(':')[0]][position]+=1
    return aa

def getReads(jamFiles):
    """
    jamFiles is a tab-delimited list of name\tfile. Will parse read
    locations as a dict of {txt:{position:ct}} according to 5'end.
    """
    aa=[]
    with open(jamFiles,'r') as f:
        for line in f:
            line=line.strip().split('\t')
            aa.append([line[0],parseReads(line[1])])
    return aa

def getMetaData(ii,AA,runDict,readDicts,X,N):
    ##
    bb=[]
    count=0
    for name,readDict in readDicts:
        aa=collections.defaultdict(float)
        for txtID,ctDictPerPosition in runDict.items():
            if txtID in readDict:
                for position,ctDictPerAA in ctDictPerPosition.items():
                    position2=position*3
                    if ctDictPerAA[AA]>=ii:##then there are at least ii AAs
                        tot=0.
                        for jj in range(position2-X,position2+X+N):
                            tot+=readDict[txtID][jj]
                        if tot>0:
                            count+=1
                            for jj in range(position2-X-12,position2+X+N):
                                aa[jj-position2]+=readDict[txtID][jj]/tot
        cc={}
        for k,v in aa.items():
            cc[k]=v/count
        ##do another round s.t. everything sums to one.
        theMetaTotal=sum(cc.values())
        dd={}
        for k,v in cc.items():
            dd[k+12]=v/theMetaTotal##add 12 nt offset for P site
        bb.append((name,dd))
    return bb

def getMetaDataPrepared(metaData,X,N):
    aa=[]
    counter=0
    for name,metaDict in metaData:
        temp=[]
        for ii in range(-X,X+N):
            if ii in metaDict:
                temp.append((ii,metaDict[ii]))
                counter+=1
            else:
                temp.append((ii,0))
        aa.append((name,temp))
    return aa, counter

def mkPlot(ii,AA,runDict,readDicts,N):
    ##
    X=100
    metaData=getMetaData(ii,AA,runDict,readDicts,X,N)
    metaDataPrepped,counter=getMetaDataPrepared(metaData,X,N)
    ##
    ##
    if counter!=0:
        g=graph.graphxy(width=16,height=2,ypos=ii*4,
            key=graph.key.key(pos='tr',hinside=0),
            x=graph.axis.linear(min=-X,max=X+N,
                title='Pos Rel Start of >=%s of %s are %s'%(ii,N,AA)),
            y=graph.axis.linear(min=0,
                title='Avg Read Density'))
    else:##to deal w/ y-axis range zero error.
        g=graph.graphxy(width=16,height=2,ypos=ii*4,
            key=graph.key.key(pos='tr',hinside=0),
            x=graph.axis.linear(min=-X,max=X+N,
                title='Pos Rel Start of >=%s of %s are %s'%(ii,N,AA)),
            y=graph.axis.linear(min=0,max=0.01,
                title='Avg Read Density'))
    ##
    ct=0
    for name,theData in metaDataPrepped:
        ct+=1
        ##
        g.plot(graph.data.points(theData,x=1,y=2,
            title=name),
            [graph.style.line([style.linestyle.solid,
                common.colors(ct)])])
        ##
        ct1,ct2=[],[]
        for entry in theData:
            if 0<=entry[0]<=30:
                ct1.append(entry[1])
            elif entry[0]<=-10 or entry[0]>=30+10:
                ct2.append(entry[1])
        ct1=numpy.average(ct1)
        ct2=numpy.average(ct2)
        if ct2>0:
            print(name,ct1/ct2)
    ##
    return g

def main(args):
    seqFile,jamFiles,N,outPrefix=args[0:]
    ##
    N=int(N)
    ##parse input fasta file
    cdsSeqs=common.parseFasta(seqFile)
    ##translate to proteins
    protSeqs=convertToProteins(cdsSeqs)
    ##scan through windows
    runDict=getRuns(protSeqs,N)
    ##parse the reads
    readDicts=getReads(jamFiles)
    ##summary
    """
    At this point,
    runDict={txtID:{position:{AA:ct_in_next_N_AAs}}}
    readDicts=[[name,{txtID:{position:ct}}],...]
    I want to loop through each AA, and then through each number from
    1 through N, and make a meta plot about all
    the positions. It would probably be good to do this for a minimum
    number of read counts.
    """
    AAs=['A','C','D','E','F','G','H','I','K','L',
        'M','N','P','Q','R','S','T','V','W','Y']
    for AA in AAs:
        print('Working on %s...'%(AA))
        c=canvas.canvas()
        for ii in range(3,7):
            print('Working on runs of length %s.'%(ii))
            g=mkPlot(ii,AA,runDict,readDicts,N)
            c.insert(g)
        c.writePDFfile(outPrefix+'.%s'%(AA))



if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
