"""
Joshua Arribere, Jan 4, 2018
Converted to python3: April 15, 2020

Script to plot select columns of a geneCt file, as output
    from median_normalizer

Input: inFile.geneCt - tab-delimited file of format:
    Gene\tCol1\tCol2...
    gene\tct\tct\t...
    myFavoriteGenes.txt - line-delimited list of genes\tgeneLabel
    colNames - space-delimited list of column names you
        want to plot. By default, the first column name
        will be the y-axis for all plots

Output: scatter plot with myFavoriteGenes highlighted

run as python3 plotGeneCts.py inFile.geneCt myFavoriteGenes.txt
    outPrefix colName colName 
"""
import sys, common, csv
from logJosh import Tee
from pyx import *

def parseGenes(myFavoriteGenes):
    aa={}
    with open(myFavoriteGenes,'r') as f:
        for line in f:
            line=line.strip().split('\t')
            aa[line[0]]=common.parseGeneList(line[1])
    return aa

def parseGeneCtFile(inFile,columns):
    """Will parse a geneCt file to a dict of format:
    {gene:{col:ct}} where col is one of columns'
    entries"""
    aa={}
    with open(inFile,'r') as f:
        freader=csv.DictReader(f,delimiter='\t')
        for row in freader:
            aa[row['Gene']]={}
            for col in columns:
                aa[row['Gene']][col]=float(row[col])
    return aa

def getDataPoints(geneCts,xAxis,yAxis):
    aa=[]
    for gene in geneCts:
        point=(geneCts[gene][xAxis],geneCts[gene][yAxis])
        if point[0]*point[1]>0:
            if point not in aa:
                aa.append(point)
    return aa

def mkScatterPlot(geneCts,counter,xAxis,yAxis,genesToHighlight):
    #initialize parameters of graph dimensions
    Width=4
    Height=4
    #set the axis parters
    logParter=graph.axis.parter.log([graph.axis.parter.log.pre1exp,
                                    graph.axis.parter.log.pre1to9exp])
    #initialize the graph
    g=graph.graphxy(xpos=counter*1.5*Width,
                    ypos=0,
                    width=Width,height=Height,
                    key=graph.key.key(pos='tr',
                                        vinside=0),
                    x=graph.axis.log(min=1,max=10**5,title=xAxis,
                                    parter=logParter),
                    y=graph.axis.log(min=1,max=10**5,title=yAxis))
    #plot bulk data
    tempData=getDataPoints(geneCts,xAxis,yAxis)
    g.plot(graph.data.points(tempData,x=1,y=2,
        title='all genes'),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk.black,deco.filled],
            size=0.01)])
    #highlight myFavoriteGenes
    ct=0
    for geneTuple in genesToHighlight:
        #gene=geneTuple[0]
        #point=[(geneCts[gene][xAxis],geneCts[gene][yAxis])]
        point=[(geneCts[gene][xAxis],geneCts[gene][yAxis]) for gene in geneTuple[1] if gene in geneCts]
        g.plot(graph.data.points(point,x=1,y=2,
            title=geneTuple[0]),
            [graph.style.symbol(graph.style.symbol.circle,
                symbolattrs=[common.colors(ct),deco.filled],
                size=0.05)])
        ct+=1
    #return the graph
    return g

def mkScatterPlots(geneCts,columns,genesToHighlight,outPrefix):
    """will make a row of scatter plots as ordered by columns.
    Will highlight those genes in genesToHighlight with labels
    being values. NOTE: To make managing this file easier, this
    script will remove points that are redundant. For this
    reason, please please please be careful with doing stats
    within this function because you may not be working with
    the full dataset"""
    
    #initialize the canvas
    c=canvas.canvas()
    #make the plots recursively
    for ii in range(1,len(columns)):
        yAxis=columns[0]
        xAxis=columns[ii]
        g=mkScatterPlot(geneCts,ii,xAxis,yAxis,genesToHighlight)
        c.insert(g)
    #write output file
    c.writePDFfile(outPrefix)
    
    #if you want to print genes in a specific read count range:
    """
    yLib=columns[0]
    xLib=columns[1]
    for gene in geneCts:
        yval=geneCts[gene][yLib]
        xval=geneCts[gene][xLib]
        if yval>200:
            if xval<40:
                #print(gene,xval,yval)
                if gene in genesToHighlight[0][1]:
                    print (gene,yval,xval)
    """

def main(args):
    inFile,myFavoriteGenes,outPrefix=args[0:3]
    columns=args[3:]
    #parse the gene names
    genesToHighlight=parseGenes(myFavoriteGenes)
    genesToHighlight=[(k,v) for (k,v) in genesToHighlight.items()]
    #parse the input data
    geneCts=parseGeneCtFile(inFile,columns)
    #plot the data
    mkScatterPlots(geneCts,columns,genesToHighlight,outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
