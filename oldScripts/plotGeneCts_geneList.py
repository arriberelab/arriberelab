"""
Joshua Arribere, Jan 4, 2018

Script to plot select columns of a geneCt file, as output
    from median_normalizer

Input: inFile.geneCt - tab-delimited file of format:
    Gene\tCol1\tCol2...
    gene\tct\tct\t...
    myFavoriteGenes.txt - line-delimited list of genes
    colNames - space-delimited list of column names you
        want to plot. By default, the first column name
        will be the y-axis for all plots

Output: scatter plot with myFavoriteGenes highlighted

run as python plotGeneCts.py inFile.geneCt myFavoriteGenes.txt
    outPrefix SJA57RZ SJZ60RZ SJA61RZ 
"""
import sys, common, csv
from logJosh import Tee
from pyx import *

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
    highlightPoints=[]
    for geneTuple in genesToHighlight:
        gene=geneTuple[0]
        if gene in geneCts:
            point=(geneCts[gene][xAxis],geneCts[gene][yAxis])
            print(gene, point)
            highlightPoints.append(point)
    g.plot(graph.data.points(highlightPoints,x=1,y=2,
        title='selectGenes'),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[common.colors(1),deco.filled],
            size=0.05)])
    #print some stuff
    """
    for gene in geneCts:
        xVal=geneCts[gene][xAxis]
        yVal=geneCts[gene][yAxis]
        if xVal>40:
            if yVal/float(xVal)>3.5:
                print gene#, xVal,yVal
    """
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

def main(args):
    inFile,myFavoriteGenes,outPrefix=args[0:3]
    columns=args[3:]
    #parse the gene names
    genesToHighlight=common.parseGeneList(myFavoriteGenes)
    genesToHighlight=[(k,v) for (k,v) in genesToHighlight.items()]
    #parse the input data
    geneCts=parseGeneCtFile(inFile,columns)
    #plot the data
    mkScatterPlots(geneCts,columns,genesToHighlight,outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
