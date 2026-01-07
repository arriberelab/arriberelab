"""
Joshua Arribere, May 10, 2022

Script to make the infograph in a new way, using my home-rolled
    nested dictionaries of doom. Now with statistics!

Input: inFile.jam - jam-formatted file
    startOrStop - whether to focus about the start or stop codon
    leftBound rightBound - how many nts to look up/downstream
    myFavGenes.txt - line-delimited list of genes to restrict to
    cutoff - min number of reads

Output: heatmap of read densities about start, stop codons showing:
    rpm and permutation p-values

run as python3 infoGraphPermutation.py inFile.jam myFavGenes.txt
    cutoff startOrStop leftB rightB outPrefix
"""
import sys, common, collections, csv, numpy, random, scipy.stats
import math
from logJosh import Tee
from pyx import *

def getPosition(txtList,startOrStop):
    """
    txtList is a list of [txt:posRelATG:posRelStop, ...]
    startOrStop is start/stop. Dep on whether start or stop, will
    return the positionRelStart/Stop if it's the same across all
    txts. Otherwise return 'na'
    """
    if startOrStop=='start':
        theVals=[int(entry.split(':')[1]) for entry in txtList]
    elif startOrStop=='stop':
        theVals=[int(entry.split(':')[2]) for entry in txtList]
    else:
        print('startOrStop must be one of start or stop. Duh. Exiting...')
        sys.exit()
    ##
    theVals=list(set(theVals))
    ##
    if len(theVals)==1:
        return theVals[0]
    else:
        return 'na'

def parseJam(inFile,startOrStop):
    """
    inFile is a jam file. Will parse to two dicts:
    readCts:{gene:ct}
    and
    dataDict:{readLen:{gene:{position:ct}}}
    Will only consider reads where all isoforms have the same position
    relative to startOrStop
    """
    ##initialize the dicts
    readCts=collections.defaultdict(int)
    dataDict=collections.defaultdict(lambda:
        collections.defaultdict(lambda:
        collections.defaultdict(float)))
    ##
    with open(inFile,'r') as f:
        f.readline()#skip the first line (header)
        for line in f:
            line=line.strip().split('\t')
            if len(line)==10:
                ##
                readLength=len(line[6])
                gene,strand=line[8].split(':')
                txtList=line[9].split('|')
                ##now figure out if the read maps to the same site in all
                ##the txts
                position=getPosition(txtList,startOrStop)
                ##
                if position!='na':#then it's the same site across all txts
                    readCts[gene]+=1
                    ##
                    dataDict[readLength][gene][position]+=1
    ##
    return readCts, dataDict

def filterByReadCtsAndGeneList(readCts,dataDict,readCutoff,myFavGenes):
    """
    readCts={gene:ct}. Will only include genes that have at least
    readCutoff number of counts and are in myFavGenes.
    dataDict={readLen:{gene:{position:ct}}}
    Will return a dict just like dataDict, but restricted to genes
    passing the above criteria.
    """
    ##first figure out who passed
    #print('Not restricting to your favorite genes...')
    passed={}
    for gene,ct in readCts.items():
        if ct>=readCutoff and gene in myFavGenes:
            passed[gene]=1
    ##
    print('%s genes passed cutoffs.'%(len(passed)))
    ##
    aa=collections.defaultdict(lambda:
        collections.defaultdict(dict))
    ##
    for readLength,geneDict in dataDict.items():
        for gene,positionDict in geneDict.items():
            if gene in passed:
                aa[readLength][gene]=positionDict
    ##
    return aa

def mkMeta(readDict,leftB,rightB):
    """
    readDict={readLen:{gene:{position:ct}}}
    Will tally up reads within [leftB,rightB], inclusive
    across all genes, and then take the average and return it
    in a dict as
    {readLen:[(position,ct)]}
    where positions are in order
    """
    ##aa will be of format
    ##{readLen:{position:[ctList]}}
    aa=collections.defaultdict(lambda:
        collections.defaultdict(list))
    ##
    for readLen,geneDict in readDict.items():
        for gene,positionDict in geneDict.items():
            for position in range(leftB,rightB+1):
                if position in positionDict:
                    aa[readLen][position].append(positionDict[position])
                else:
                    aa[readLen][position].append(0)
    ##
    ##now get the average
    bb=collections.defaultdict(list)
    for readLen,positionDict in aa.items():
        for position in range(leftB,rightB+1):
            bb[readLen].append((position,
                numpy.average(positionDict[position])))
    ##
    return bb

def scramble(dataDict):
    """
    dataDict={readLen:{gene:{position:ct}}}
    Will scramble {position:ct} and return
    """
    aa=collections.defaultdict(lambda:
        collections.defaultdict(lambda:
            collections.defaultdict(float)))
    ##
    for readLen,geneDict in dataDict.items():
        for gene,positionDict in geneDict.items():
            ##first get non-zero read counts
            tempDict=dict((k,v) for k,v in positionDict.items() if v!=0)
            ##
            theVals=list(tempDict.values())
            random.shuffle(theVals)
            shuffleDict=dict(zip(tempDict,theVals))
            ##
            aa[readLen][gene]=shuffleDict
    ##
    return aa

def tallyUp(dataDictRandMeta,dataDictMetaMeta):
    """
    dataDictRandMeta={readLen:[(position,ct)]}
    dataDictMetaMeta={readLen:{position:[ctList]}}
    will append ct to ctList
    """
    for readLen,entryList in dataDictRandMeta.items():
        for entry in entryList:
            position,ct=entry[0],entry[1]
            dataDictMetaMeta[readLen][position].append(ct)
    return dataDictMetaMeta

def mkMetaPermutation(dataDict,dataDictMeta,leftB,rightB):
    """
    This is going to be a massive function-of-functions.
    dataDict={readLen:{gene:{position:ct}}}
    dataDictMeta={readLen:[(position,ct)]} where position is in
    [leftB,rightB], inclusive.
    For each of N permutations, will:
    (1) Scramble the {position:ct} in dataDict.
    (2) Compute metaDict of format
        {readLen:[(position,ct)]}
    (3) Tally up:
        {readLen:{position:[randCt1,randCt2,...]}}
    (4) Compute a z-score for how many stdevs the real metaCt is
        away from the distribution of randCts at that position.
    (5) Output a dict of
        {readLen:[(position,pval)]}
        where pval is the one-tailed p-value from the z-score.
        If pval is less than some cutoff, will make it 'na'
    """
    ##hard code the number of iterations
    N=30#blah blah XXXX
    ##
    dataDictMetaMeta=collections.defaultdict(lambda:
                        collections.defaultdict(list))
    ##
    for ii in range(N):
        ##print info
        if ii%10==0:
            print('Working on %s of %s...'%(ii,N))
        ##(1) scramble the dataDict positions
        dataDictRand=scramble(dataDict)
        ##(2) compute metaDict
        dataDictRandMeta=mkMeta(dataDictRand,leftB,rightB)
        ##(3) Tally up
        dataDictMetaMeta=tallyUp(dataDictRandMeta,dataDictMetaMeta)
    ##randomizations completed.
    ##now loop through the readLengths and positions and
    ##compute z-scores
    outDict=collections.defaultdict(list)
    for readLen,entryList in dataDictMeta.items():
        for entry in entryList:
            position,ct=entry[0],entry[1]
            ##(4) Compute a z-score
            temp=dataDictMetaMeta[readLen][position]
            tempAvg=numpy.average(temp)
            tempStDev=numpy.std(temp)
            ##
            if tempStDev==0:
                pval='na'
            else:
                zScore=(ct-tempAvg)/tempStDev
                pval=scipy.stats.norm.sf(zScore)
                if pval>0:
                    pval=-math.log(pval,10)
                else:#rounds
                    pval=10
                pval=min([pval,10])
                ##pval cutoff
                if pval<=-math.log(0.01,10):#pval cutoff
                    pval='na'
            ##(5) Output a dict
            outDict[readLen].append((position,pval))
    ##
    return outDict

def mkHeatMap(metaDict,leftB,rightB,ii):
    """
    metaDict={readLen:[(position,ct)]}
    will make a heatmap with min=0, max=max(cts)
    """
    ##get an ordered list of read lengths
    readLengths=list(metaDict.keys())
    readLengths.sort()
    ##get the max value to set the scale
    theVals=[]
    for readLen,entryList in metaDict.items():
        theVals+=[entry[1] for entry in entryList]
    theMax=max([entry for entry in theVals if entry!='na'])
    ##initialize the canvas
    g=canvas.canvas()
    ##
    for jj in range(len(readLengths)):
        readLen=readLengths[jj]
        g.text(0,jj+0.5-ii*20,readLen,
            [text.size(4),text.halign.right,text.valign.middle])
        for kk in range(len(metaDict[readLen])):
            entry=metaDict[readLen][kk]
            ##retrieve the value
            position,theVal=entry[0],entry[1]
            ##add text
            if jj==0 and position%3==0:
                g.text(kk-0.5,jj-0.5-ii*20,position,
                [text.size(4),text.halign.center,text.valign.top])
            ##
            if theVal!='na':
                theVal/=theMax
                ##convert that to a color
                theColor=color.cmyk(1*theVal,0.5*theVal,0,0)
            else:
                theColor=color.cmyk(0,0,0,0.1)#light grey
            ##
            g.stroke(path.rect(kk,jj-ii*20,1,1),[style.linewidth.Thin,
                                 color.cmyk(0,0,0,0),
                                 deco.filled([theColor])])
            ##
    ##add key
    N=4#number of steps in the key
    for ll in range(N):
        ##get the value
        theFrac=ll/N
        theVal=theFrac*theMax
        theColor=color.cmyk(1*theFrac,0.5*theFrac,0,0)
        ##draw the rect
        g.stroke(path.rect(kk+2,jj-ii*20-ll,1,1),
            [style.linewidth.Thin,color.cmyk(0,0,0,0),
                deco.filled([theColor])])
        ##write the text
        g.text(kk+2+1,jj+0.5-ii*20-ll,round(theVal,2),
            [text.size(4),text.halign.left,text.valign.middle])
    ##
    return g

def mkPlot(metaDictList,outPrefix,leftB,rightB):
    """
    metaDictList is a list of dicts of the format
    {readLen:[(position,ct)]}
    where positions are ordered between [leftB,rightB],
    inclusive. readLens are not ordered.
    Will plot heatmaps arrayed vertically.
    """
    ##
    c=canvas.canvas()
    ##
    for ii in range(len(metaDictList)):
        ##
        metaDict=metaDictList[ii]
        ##
        g=mkHeatMap(metaDict,leftB,rightB,ii)
        ##
        c.insert(g)
    ##
    c.writePDFfile(outPrefix)

def main(args):
    ##input the files and parameters
    inFile,annotFile,myFavGeneFile,readCutoff,startOrStop,leftB,rightB,\
        outPrefix=args[0:]
    #
    leftB,rightB=int(leftB),int(rightB)
    #
    ##input the jam file
    #inFile, startOrStop
    #out: readCts, dataDict
    readCts,dataDict=parseJam(inFile,startOrStop)
    #
    ##filter by read counts
    #readCutoff, myFavGeneFile, readCts, dataDict
    #out: readDict
    dataDict=filterByReadCtsAndGeneList(readCts,dataDict,int(readCutoff),\
        common.parseGeneList(myFavGeneFile))
    #
    ##make ready for plotting
    #readDict
    #out: readDictMeta
    dataDictMeta=mkMeta(dataDict,leftB,rightB)
    #
    ##Do permutations to get p-values
    #readDict & readDictmeta
    #out: readDictPvalForPlotting
    dataDictMetaPval=mkMetaPermutation(dataDict,dataDictMeta,leftB,rightB)
    #
    ##plot everything
    #[startDictMeta,startDictPvalForPlotting]
    #out: graph, saved to file outPrefix.pdf
    mkPlot([dataDictMeta,dataDictMetaPval],
        outPrefix,leftB,rightB)
    #

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
