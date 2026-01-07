"""
Joshua Arribere, June 4, 2021

Script to make a heat map of stop codons.

Input: file.jam - will be used to identify reads.
    wbgenes.txt - line-delimited list of genes. Will highlight these.

Output: heatmap of CDF of reads about some window of the stop codon. Will
    have a bunch of filters in this code.

run as python3 stopCodonHeatMap.py inFile.jam wbgenes.txt outPrefix
"""
import sys, common, collections
from logJosh import Tee
from pyx import *

def getGoodAndPosition(txts):
    """
    For a list of txts, will return two things:
    (1) 1 or 0 depending on whether all txts have the same positRelStop
    (2) the positRelStop
    """
    temp=[]
    for txt in txts:
        txtName,startPosit,stopPosit=txt.split(':')
        temp.append(int(stopPosit))
    ##
    temp=list(set(temp))
    if len(temp)==1:
        return 1, temp[0]
    else:
        return 0, 'na'

def getReads(jamFile):
    """
    Will input a jam file, and do the following:
    1. Identify reads w/in X nts of a stop codon (based on 3'end).
    2. Throw out stop codons with fewer than some number of total reads.
    3. Throw out genes w/ multiple possible stop codons.
    After all that, will return a dict of format:
    {wbgene:{position:ct}}
    Note that all of this will be for the 3'end of the read.
    """
    ##set some parameters
    window=50
    readCtCutoff=50
    print('Using a window of %s nts about stop, w/ at least %s reads.'%(
        window,readCtCutoff))
    ##
    aa=collections.defaultdict(lambda:collections.defaultdict(int))
    badGenes={}
    ##
    with open(jamFile,'r') as f:
        f.readline()
        for line in f:
            line=line.strip().split('\t')
            ##
            wbgene,strand=line[8].split(':')
            ##
            if wbgene not in badGenes and strand=='S' and line[7]=='1:1':
                ##
                readSeq=line[6]
                txts=line[9].split('|')
                ##
                good,position=getGoodAndPosition(txts)
                if good:
                    position+=len(readSeq)-2
                    ##position=0 > T of stop codon
                    ##position=1 > A/G of stop codon
                    ##position=2 > A/G of stop codon
                    ##position=3 > nt after stop codon
                    if position in range(-window,window+1):
                        aa[wbgene][position]+=1
                else:
                    badGenes[wbgene]=1
    ##
    cc={}
    for wbgene in aa:
        if wbgene not in badGenes:
            total=sum(aa[wbgene].values())
            if total>=readCtCutoff:
                cc[wbgene]=aa[wbgene]
    ##
    print('%s genes survived the cutoffs.'%(len(cc)))
    return cc, window

def process(inDict,window):
    """
    inDict of the format
    {wbgene:{position:ct}}
    and window an integer.
    Will convert to a list where each entry is
    [wbgene,[(position,cdf)]]
    where all positions in range(-window,window+1) are covered and y-value
    is the cdf
    """
    ##
    sortSite=0.5
    print('Will sort stops by when they pass %s%s of read cts.'%(
        int(sortSite*100),"%"))
    ##
    aa={}
    bb={}
    for wbgene,v in inDict.items():
        temp=[]
        ##
        total=sum(v.values())
        runningTally=0
        ##
        for ii in range(-window,window+1):
            runningTally+=v[ii]
            temp.append((ii,runningTally/total))
            ##
            if runningTally/total>=sortSite:
                if wbgene not in bb:
                    bb[wbgene]=ii
        ##
        aa[wbgene]=temp
    ##
    cc=collections.defaultdict(list)
    for k,v in bb.items():
        cc[v].append(k)
    ##
    dd=[]
    for ii in range(-window,window+1):
        for wbgene in cc[ii]:
            dd.append([wbgene,aa[wbgene]])
    ##
    return dd

def getColor(y):
    """
    Given 0<=y<=1, will return:
    black>blue if y<0.25
    blue>green if 0.25<=y<0.5
    green>yellow if 0.5<=y<0.75
    yellow>red if 0.75<=y
    """
    if y<0.25:
        x=y/0.25
        return color.rgb(0,114*x/255,178*x/255)
    elif 0.25<=y<0.5:
        x=(y-.25)/.25
        return color.rgb(0,(114+44*x)/255,(178-63*x)/255)
    elif 0.5<=y<0.75:
        x=(y-.5)/.25
        return color.rgb(240*x/255,(158+70*x)/255,(115-49*x)/255)
    else:
        x=(y-.75)/.25
        return color.rgb((240-27*x)/255,(228-134*x)/255,(66-66*x)/255)

def mkHeatMap(stopPositions,window,wbgenes,outPrefix):
    """
    stopPositions is of the format:
    {wbgene:{position:ct}} where pos:ct is a defaultdict
    window is an integer
    wbgenes is a dict of genes
    Will draw a heatmap cdf for each wbgene in stopPositions. Will
    draw a box if that wbgene is also present in wbgenes. Will probably
    order the genes, but not yet sure how.
    """
    stopPositions=process(stopPositions,window)
    ##
    c=canvas.canvas()
    unit.set(defaultunit="pt")
    pixelSize=1
    ##
    ##the next line restricts to user-provided list of genes
    #stopPositions=[entry for entry in stopPositions if entry[0] in wbgenes]
    ##
    for ii in range(len(stopPositions)):
        entry=stopPositions[ii]
        wbgene,data=entry
        ##
        for position in data:
            x,y=position
            ##
            theColor=getColor(y)
            c.fill(path.rect(x,-ii,pixelSize,pixelSize),
                [theColor])
        ##
        if wbgene in wbgenes:
            c.fill(path.rect(x+pixelSize*2,-ii,4*pixelSize,pixelSize),
                [color.rgb(0,0,0)])
            print('%s special'%(wbgene))
        else:
            print(wbgene)
        ##
    ##add the axis
    for ii in range(-window,window+1):
        if ii%3==0:
            c.fill(path.rect(ii,pixelSize,pixelSize,pixelSize),
                [color.rgb(0,0,0)])
    ##label the axis
    c.text(0+pixelSize/2.,pixelSize*3,"0",
        [text.halign.boxcenter,text.valign.bottom])
    ##add a key
    scale=10
    pixelSize=pixelSize*3
    for ii in range(0,scale+1):
        temp=ii/scale
        theColor=getColor(temp)
        x=-window-pixelSize*10
        y=-pixelSize*ii
        c.fill(path.rect(x,y,pixelSize,pixelSize),[theColor])
        ##
        if ii%5==0:
            c.text(x-pixelSize*2,y+pixelSize/2.,'%0.2f'%(temp),
                [text.halign.boxright,text.valign.middle])
    ##
    c.writePDFfile(outPrefix)

def main(args):
    inFile,wbgeneFile,outPrefix=args[0:]
    ##
    wbgenes=common.parseGeneList(wbgeneFile)
    ##
    stopPositions,window=getReads(inFile)
    ##
    mkHeatMap(stopPositions,window,wbgenes,outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
