"""
Joshua Arribere, Sept 2, 2024

Script to plot L2FC for genes that significantly change in subsets of
    smg-5d(smg2 in script), smg-5m, and smg-6d.

Input: DESeqFileSmg5d
    DESeqFileSmg6d
    DESeqFileSmg5m

Output: L2FC plots with
    smg-5d only, smg-6d only, smg-5d+6d only, smg-5d+6d+m only

run as python3 plotDESeq2_smg5dSmg6dSmg5m.py Smg5d Smg6d Smg5m outPrefix
"""
import sys, common, scipy.stats
from logJosh import Tee
from pyx import *

def parse(smgFile):
    print('Flipping L2FC.')
    aa={}
    with open(smgFile,'r') as f:
        f.readline()#skip header
        for line in f:
            line=line.strip().split(',')
            gene=line[0].strip('"')
            l2fc=line[2]
            padj=line[-1]
            if padj!='NA':
                aa[gene]=(-float(l2fc),float(padj))
    ##
    return aa

def getGeneLists(smg2Dict,smg5Dict,smg6Dict):
    """
    Each smgDict is of the format {wbgene:(l2fc,pval)}
    Will sort genes into whether they change in
    smg2 only
    smg5 only
    smg2+5 but not smg6
    smg2+5+6
    """
    smg2,smg5,smg25,smg256=[],[],[],[]
    pvalCut=0.05
    for k,v in smg2Dict.items():
        if k in smg5Dict and k in smg6Dict:
            smg2l2fc=v[0]
            smg2p=v[1]
            smg5l2fc=smg5Dict[k][0]
            smg5p=smg5Dict[k][1]
            smg6l2fc=smg6Dict[k][0]
            smg6p=smg6Dict[k][1]
            ##
            if smg2p<=pvalCut:
                if smg5p>=pvalCut and smg6p>=pvalCut:
                    if smg2l2fc>=0:
                        smg2.append(k)
                elif smg5p<=pvalCut and smg6p>=pvalCut:
                    if smg2l2fc>=0:
                        smg25.append(k)
                elif smg5p<=pvalCut and smg6p<=pvalCut:
                    if smg2l2fc>=0:
                        smg256.append(k)
            elif smg5p<=pvalCut and smg6p>=pvalCut:
                if smg5l2fc>=0:
                    smg5.append(k)
    ##
    print('Counts for diff expr\'ed genes:',
        '\nsmg-2: %s'%(len(smg2)),
        '\nsmg-5: %s'%(len(smg5)),
        '\nsmg-25: %s'%(len(smg25)),
        '\nsmg-256: %s'%(len(smg256)))
    ##
    return smg2, smg5, smg25, smg256

def mkCDF(inList):
    inList.sort()
    return [(entry,inList.index(entry)/len(inList)) for entry in inList]

def mkPlot(inDict,smg2,smg5,smg25,smg256,ct):
    """
    Will plot a cdf, highlighting the various gene lists. Will
    place the plot at ct position.
    inDict={wbgene:[l2fc,pvalAdj])
    """
    ##
    gridpainter=graph.axis.painter.regular(gridattrs=[
        attr.changelist([style.linestyle.dashed,None])])
    ##initialize the graph
    g=graph.graphxy(width=2.5,height=2,xpos=ct*3,
        x=graph.axis.linear(min=-1,max=4,title='L2FC',
            painter=gridpainter),
        y=graph.axis.linear(title='CDF',painter=gridpainter))
    ##add the data
    temp1=[v[0] for k,v in inDict.items() if v[1]!='NA']
    g.plot(graph.data.points(mkCDF(temp1),x=1,y=2),
        [graph.style.line()])
    ##now highlight the gene lists
    ii=0
    colors=[color.cmyk(1,0.5,0,0),#blue
        color.cmyk(0.97,0,0.75,0),#green
        color.cmyk(0,0.8,1,0),#red
        color.cmyk(0.1,0.7,0,0)]#redPurple
    for geneList in [smg2,smg5,smg25,smg256]:
        theColor=colors[ii]
        temp=[v[0] for k,v in inDict.items() if v[1]!='NA' and k in geneList]
        g.plot(graph.data.points(mkCDF(temp),x=1,y=2),
            [graph.style.line([theColor])])
        ii+=1
        ##run some stats
        print('Working on dataset %s, analyzing gene list %s'%(ct,ii))
        print(scipy.stats.ks_2samp(temp1,temp))
    ##
    return g

def main(args):
    smg2File,smg5File,smg6File,outPrefix=args[0:]
    ##
    smg2Dict=parse(smg2File)
    smg5Dict=parse(smg5File)
    smg6Dict=parse(smg6File)
    #smg7Dict=parse(smg7File)
    ##
    ##get the gene lists
    smg2,smg5,smg25,smg256=getGeneLists(smg2Dict,smg5Dict,smg6Dict)
    ##initialize the canvas, then make the plots
    c=canvas.canvas()
    ct=0
    for entry in [smg2Dict,smg5Dict,smg6Dict]:
        g=mkPlot(entry,smg2,smg5,smg25,smg256,ct)
        c.insert(g)
        ct+=1
    ##
    c.writePDFfile(outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
