"""
Joshua Arribere, May 5, 2022

Script to analyze smg-1 and smg-6 targets w/ DESeq2.

Input: cts_S - sense-stranded cts
    cts_AS - antisense-stranded cts
    conditions - file of format:
        libName\tuntreated/treated1/treated6\tsingle-read

Output: Will make a series of files:
 - combine S and AS cts into a single file
 - run diffExpr on untreated v treated1, output results
 - run diffExpr on untreated v treated6, output results
 - plots of untreated v treated1 and untreated v treated6 with:
     - treated1 targets highlighted
     - treated6 targets highlighted
     - treated1 but not 6 targets highlighted
     - treated6 but not 1 targets highlighted

run as python3 plotWtSmg1Smg6.py cts_S cts_AS conditions outPrefix
"""
import sys, common, collections, os, math
from logJosh import Tee
from pyx import *

def parse(senseCtsFile):
    """
    Header is of format "\tlibName1\tlibName2\t...libNameN"
    Will return a dict of format:
    {gene:{libNameI:ct}
    """
    aa=collections.defaultdict(lambda:collections.defaultdict(int))
    start=True
    with open(senseCtsFile,'r') as f:
        for line in f:
            line=line.strip().split('\t')
            if start==True:
                header=line
                start=False
            else:
                gene=line[0]
                for ii in range(len(line[1:])):
                    aa[gene][header[ii]]=float(line[ii+1])
    return aa

def combine(senseCts,antisenseCts):
    """
    Will combine, and discard anything w/ less than some number of cts
    """
    aa=list(senseCts.keys())
    bb=list(antisenseCts.keys())
    cc=list(set(aa+bb))
    ##
    dd=collections.defaultdict(lambda:collections.defaultdict(int))
    for gene in cc:
        aa=list(senseCts[gene].keys())
        bb=list(antisenseCts[gene].keys())
        temp=list(set(aa+bb))
        for key in temp:
            dd[gene][key]=senseCts[gene][key]+antisenseCts[gene][key]
    ##now filter by minimum counts
    ee={}
    minCt=5.
    for gene,libDict in dd.items():
        if min(libDict.values())>=minCt:
            ee[gene]=libDict
    ##
    return ee

def writeCtDict(ctDict,outPrefix):
    header=list(ctDict['WBGene00022277'].keys())
    header.sort()
    ##
    with open(outPrefix,'w') as f:
        f.write('\t'+'\t'.join(header))
        for k,v in ctDict.items():
            temp='\n%s'%(k)
            for libName in header:
                temp+='\t%s'%(ctDict[k][libName])
            f.write(temp)

def getTargets(deSeq2File,pvalCut):
    """
    Will return a list of wbgenes for deSeq2 output file and a pval
    cutoff (e.g., 0.01)
    """
    aa=[]
    with open(deSeq2File,'r') as f:
        f.readline()#skips the header
        for line in f:
            line=line.strip().split(',')
            ##
            wbgene=line[0].strip('"')
            log2FC=float(line[2])
            try:
                pval=float(line[-1])
                ##
                if pval<=pvalCut and log2FC>0:
                    aa.append(wbgene)
            except ValueError:
                pass
    ##
    """
    print('%s targets extracted at cutoff %s.'%(len(aa),pvalCut))
    with open(deSeq2File+'.targets','w') as f:
        for gene in aa:
            f.write(gene+'\n')
    """
    with open(deSeq2File+'.targets','w') as f:
        for gene in aa:
            f.write(gene+'\n')
    ##
    return aa

def mkPlot(deSeq2File,targetName,targetList,theColor):
    """
    deSeq2File will be plotting and saved in deSeq2File+targetName.
    Will highlight genes in targetList.
    """
    aa={}
    with open(deSeq2File,'r') as f:
        f.readline()
        for line in f:
            line=line.strip().split(',')
            ##
            wbgene=line[0].strip('"')
            x=float(line[1])
            y=float(line[2])
            ##
            aa[wbgene]=(x,y)
    ##
    g=graph.graphxy(width=8,height=4,
        key=graph.key.key(pos='tr'),
        x=graph.axis.log(min=5,title='Avg Expression'),
        y=graph.axis.linear(min=-2.5,max=5,title='Log2FC'))
    ##
    g.plot(graph.data.points(list(aa.values()),x=1,y=2,title=None),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk.black,deco.filled],
            size=0.01)])
    ##now plot the targets
    bb=[aa[wbgene] for wbgene in aa if wbgene in targetList]
    g.plot(graph.data.points(bb,x=1,y=2,title=targetName),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[theColor,deco.filled],
            size=0.04)])
    ##the following lines can be commented in/out to print
    ##targets of interest
    for wbgene,val in aa.items():
        if wbgene in targetList:
            if val[1]<0:
                print(wbgene,val)
    print('Done printing genes of interest...')
    ##
    g.writePDFfile(deSeq2File+'_%s'%(targetName))

def getPvals(deSeq2File):
    """
    Will return a dict of wbgene:log_pval for all genes w/ a positive L2FC
    """
    aa={}
    with open(deSeq2File,'r') as f:
        f.readline()#skip header
        for line in f:
            line=line.strip().split(',')
            ##
            wbgene=line[0].strip('"')
            ##
            LFC=float(line[2])
            if LFC>0:
                try:
                    pval=-math.log(float(line[-1]),10)
                    aa[wbgene]=pval
                except ValueError:
                    pass
    return aa

def mkPvalPlot(outPrefix1,outPrefix6,bothTargets,smg1Not6Targets,
        smg6Not1Targets,outPrefix):
    """
    outPrefix1 and 6 are deSeq2 output files. targets are lists of wbgenes
    Will plot the log of the adjusted pvalue against each other only for
    genes w/ a LFC>0
    """
    pvals1=getPvals(outPrefix1)
    pvals6=getPvals(outPrefix6)
    ##
    pvals={}
    for gene,v in pvals1.items():
        if gene in pvals6:
            pvals[gene]=(v,pvals6[gene])
    ##use to print some things
    ct=0
    for k,v in pvals.items():
        if v[0]>5 and v[1]<2:
            print(k,v)
            ct+=1
    print('%s genes satisfied the requirements.'%(ct))
    ##
    g=graph.graphxy(width=4,height=4,key=graph.key.key(pos='br',
        hinside=0),
        x=graph.axis.linear(max=20,title='Smg1 v WT'),
        y=graph.axis.linear(max=20,title='Smg6 v WT'))
    ##
    g.plot(graph.data.points(list(pvals.values()),x=1,y=2,title=None),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk.black,deco.filled],
            size=0.01)])
    ##now plot the targets
    ##bothTargets,smg1Not6Targets,smg6Not1Targets
    g.plot(graph.data.points(
        [pvals[wbgene] for wbgene in pvals if wbgene in 
        bothTargets],
        x=1,y=2,title='Smg16 Targets'),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk(0.97,0,0.75,0),deco.filled],
            size=0.02)])
    g.plot(graph.data.points(
        [pvals[wbgene] for wbgene in pvals if wbgene in 
        smg1Not6Targets],
        x=1,y=2,title='Smg1Not6 Targets'),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk(0.8,0,0,0),deco.filled],
            size=0.02)])
    g.plot(graph.data.points(
        [pvals[wbgene] for wbgene in pvals if wbgene in 
        smg6Not1Targets],
        x=1,y=2,title='Smg6Not1 Targets'),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk(0,0.8,1,0),deco.filled],
            size=0.02)])
    ##
    g.writePDFfile(outPrefix)

def main(args):
    senseCtsFile,antisenseCtsFile,conditionsFile,outPrefix=args[0:]
    ##
    senseCts=parse(senseCtsFile)
    antisenseCts=parse(antisenseCtsFile)
    ##
    ctDict=combine(senseCts,antisenseCts)
    writeCtDict(ctDict,outPrefix+'.cts')
    ##Call DESeq2
    outPrefix1=outPrefix+'treated1_v_untreated'
    os.system('Rscript %s %s %s %s %s %s'%('/data16/joshua/scripts/median_normalize_DESeq2.r',
        outPrefix+'.cts',conditionsFile,'treated1','untreated',
        outPrefix1))
    outPrefix6=outPrefix+'treated6_v_untreated'
    os.system('Rscript %s %s %s %s %s %s'%('/data16/joshua/scripts/median_normalize_DESeq2.r',
        outPrefix+'.cts',conditionsFile,'treated6','untreated',
        outPrefix6))
    ##
    ##ID targets
    print('#################################################')
    smg1Targets=getTargets(outPrefix1,0.01)
    smg6Targets=getTargets(outPrefix6,0.01)
    ##
    smg1Not6Targets=[wbgene for wbgene in smg1Targets if 
        wbgene not in smg6Targets]
    print('%s of %s targets were unique to smg-1.'%(len(smg1Not6Targets),
        len(smg1Targets)))
    smg6Not1Targets=[wbgene for wbgene in smg6Targets if 
        wbgene not in smg1Targets]
    print('%s of %s targets were unique to smg-6.'%(len(smg6Not1Targets),
        len(smg6Targets)))
    print('#################################################')
    ##
    ##Now plot everything
    mkPlot(outPrefix1,'smg1Targets',smg1Targets,color.cmyk(1,0.5,0,0))
    mkPlot(outPrefix1,'smg6Targets',smg6Targets,color.cmyk(0,0.5,1,0))
    mkPlot(outPrefix6,'smg6Targets',smg6Targets,color.cmyk(0,0.5,1,0))
    mkPlot(outPrefix6,'smg1Targets',smg1Targets,color.cmyk(1,0.5,0,0))
    #Now plot the ones that are unique to one or the other
    mkPlot(outPrefix1,'smg1Not6Targets',smg1Not6Targets,
        color.cmyk(0.8,0,0,0))
    mkPlot(outPrefix1,'smg6Not1Targets',smg6Not1Targets,
        color.cmyk(0,0.8,1,0))
    mkPlot(outPrefix6,'smg1Not6Targets',smg1Not6Targets,
        color.cmyk(0.8,0,0,0))
    mkPlot(outPrefix6,'smg6Not1Targets',smg6Not1Targets,
        color.cmyk(0,0.8,1,0))
    ##
    ##Now plot the p-values against each other
    bothTargets=[entry for entry in smg1Targets if entry in smg6Targets]
    mkPvalPlot(outPrefix1,outPrefix6,bothTargets,smg1Not6Targets,
        smg6Not1Targets,outPrefix+'_pval')

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])

