"""
Joshua Arribere, Jan 4, 2018
Converted to python3: April 15, 2020
Marcus Viscardi, July 16, 2021
    Added plotly support, so that you can just highlight genes to see differences

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

run as python3 plotGeneCts_plotly.py inFile.geneCt outPrefix colName1 colName2
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
    yLib=columns[0]
    xLib=columns[1]
    for gene in geneCts:
        yval=geneCts[gene][yLib]
        xval=geneCts[gene][xLib]
        if yval>300:
            if xval<300:
                print (gene,yval,xval)

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


def pdParseGeneCtFile(inFile):
    import pandas as pd
    dataframe = pd.read_csv(inFile, sep="\t")
    dataframe.rename(columns={"Unnamed: 0": "gene_id"}, inplace=True)
    gene_name_df = pd.read_csv("../geneNames_and_WBGenes.tsv", sep="\t")
    dataframe = dataframe.merge(gene_name_df[["gene_id", "gene_name"]], how="left")
    dataframe["identity"] = dataframe["gene_id"] + " (" + dataframe["gene_name"] + ")"
    return dataframe


def plotlyMkScatterplots(geneCtDF, columns_to_plot):
    import plotly.express as px
    fig = px.scatter(geneCtDF,
                     x=columns_to_plot[0], y=columns_to_plot[1],
                     hover_name="identity")
    fig.update_xaxes(type="log")
    fig.update_yaxes(type="log")
    fig.show()


def main_plotly_and_pandas(args):
    inFile,outPrefix=args[0:2]
    columns_to_plot = args[2:]
    df = pdParseGeneCtFile(inFile)
    print(df.info())
    plotlyMkScatterplots(df, columns_to_plot)


if __name__=='__main__':
    # Tee()
    # main(sys.argv[1:])
    fake_args = ["/data16/anniec/working/200419_SJA246-248/200427_SJA246-248_S.geneCt",
                 "test_out_prefix", "200210_SJA247Nugen", "200210_SJA246Nugen"]
    main_plotly_and_pandas(fake_args)
    
