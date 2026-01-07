"""
Joshua Arribere, July 24, 2025

Script to quantify RFP signal normalized to GFP.

Input: image_c1ORG.tif - red signal (c1)
    image_c2ORG.tif - green signal (c2)
    in two images

Output: will make a series of plots useful to diagnosing
    the behavior of this script. Will also print-to-screen
    the final RFP/GFP value.

run as python3 redGreenWormsC12.py i1_c1ORG.tif i2_c2ORG.tif
    outPrefix
"""
import sys, common, random
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from logJosh import Tee
from math import log, dist
from scipy.stats import linregress
import tifffile
from pyx import *

#def dist(p, q):
    # drop-in replacement for math.dist (Python < 3.8)
    #return sqrt(sum((pi - qi) ** 2 for pi, qi in zip(p, q)))

def getTheVal(xx,yy,inData,r):
    """
    Will get a list of values of all positions within a circle
    drawn by r, centered about point xx,yy. Will return a single
    summary metric of this list (median or avg)
    """
    temp=[]
    pt1=(xx,yy)
    for ii in range(xx-r,xx+r+1):
        for jj in range(yy-r,yy+r+1):
            pt2=(ii,jj)
            euclDist=dist(pt1,pt2)
            if euclDist<=r:
                temp.append(inData[ii][jj])
    ##
    return np.average(temp)

def processFile(inFile1,inFile2,r,outPrefix):
    """
    Will input a tif file. Will examine pixels using a radius of
    r
    """
    with open(inFile1,'rb') as f:
        inData1=tifffile.imread(f)
    with open(inFile2,'rb') as f:
        inData2=tifffile.imread(f)
    ##
    xDim=inData1.shape[0]#number of x-axis pixels
    yDim=inData1.shape[1]#number of y-axis pixels
    ##
    processedPixels=np.zeros((xDim,yDim))##initialize
    processedPixelsRed=np.zeros((xDim,yDim))##initialize
    ##object for the processed pixels vals
    for xx in range(r,xDim-r):
        for yy in range(r,yDim-r):
            theVal=getTheVal(xx,yy,inData2,r)
            ##
            processedPixels[xx][yy]=theVal
            theRedVal=getTheVal(xx,yy,inData1,r)
            ##
            processedPixelsRed[xx][yy]=theRedVal
    #print(processedPixels.tolist())#for a list of lists
    #print(processedPixels.flatten())# for just a single list
    allProcessedPixelVals=processedPixels.flatten()
    print('Fuzzied Pixel Intensity Vals:')
    print('Median:%s'%(np.median(allProcessedPixelVals)))
    print('Avg:%s'%(np.average(allProcessedPixelVals)))
    print('Currently the script doesn\'t use these vals...')
    ##
    plt.imshow(processedPixels)
    plt.colorbar(label='Fuzzied Green Pixel Value')
    plt.savefig(outPrefix+'.png')
    plt.close()
    ##
    ##
    theFilteredPixels=np.zeros((xDim,yDim))##initialize
    theCut=5000
    for xx in range(r,xDim-r):
        for yy in range(r,yDim-r):
            theVal=processedPixels[xx][yy]
            ##
            if theVal>=theCut:
                theFilteredPixels[xx][yy]=1
    plt.imshow(theFilteredPixels)
    plt.colorbar(label='Passed Fuzzied Green Pixel After Cutoff')
    plt.savefig(outPrefix+'.passed.png')
    plt.close()
    ##
    ##now extract the (r,g) vals for pixels that have the minimum
    ##g intensities
    thePixels=[]
    justTheReds=np.zeros((xDim,yDim))##initialize
    for xx in range(r,xDim-r):
        for yy in range(r,yDim-r):
            if theFilteredPixels[xx][yy]==1:
                #thePixels.append((processedPixelsRed[xx][yy],
                #                  processedPixels[xx][yy]))
                thePixels.append((inData1[xx][yy],
                                  inData2[xx][yy]))
                justTheReds[xx][yy]=inData1[xx][yy]
    ##
    plt.imshow(justTheReds)
    plt.colorbar(label='Passed Fuzzied Green Pixel but with Red Pixel Val')
    plt.savefig(outPrefix+'.passed.redPixelVal.png')
    plt.close()
    ##
    return thePixels

def mkScatterPlot(xyVals,outPrefix):
    """
    xyVals=[(x,y),...]
    Will perform linear regression
    """
    yVals=[entry[0] for entry in xyVals]#red
    xVals=[entry[1] for entry in xyVals]#green
    ##
    slope, intercept, r, p, se = linregress(xVals,yVals)
    ##
    #X = np.vstack([np.ones(len(xVals)), xVals]).T
    #y = np.array(yVals)
    #beta, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
    #intercept, slope = beta
    ##
    print(f"Slope: {slope}")
    print(f"Intercept: {intercept}")
    print(f"R-squared: {r**2}")
    print('Avg Red Pixel:%s'%(np.average(yVals)))
    ##
    #ratio=[]
    #for x,y in zip(xVals,yVals):
    #    if x>0:
    #        ratio.append(y/x)
    #print('######\n%s\t%s\n######'%(np.mean(ratio),np.median(ratio)))
    ##
    g=graph.graphxy(width=4,height=4,
        x=graph.axis.linear(title='GFP'),
        y=graph.axis.linear(title='RFP'),
        key=graph.key.key(pos='tr',hinside=0))
    ##
    tempData=list(zip(xVals,yVals))
    random.shuffle(tempData)
    g.plot(graph.data.points(tempData[:1000],x=1,y=2),
        [graph.style.symbol(graph.style.symbol.circle,
            symbolattrs=[color.cmyk.black,deco.filled],
            size=0.01)])
    ##
    g.plot(graph.data.function("y(x)=%s*x+%s"%(slope,intercept),
        title="y(x)=%s*x+%s"%(slope,intercept)))
    ##
    g.writePDFfile(outPrefix)

def main(args):
    inFile1,inFile2,outPrefix=args[0:]
    ##
    inData=processFile(inFile1,inFile2,5,outPrefix+'.pixels')
    ##
    mkScatterPlot(inData,outPrefix+'.scatter')



if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
