"""
Joshua Arribere Dec 2, 2013

Script to make a pdf plot of an arbitrary number of libraries, focusing on one
    continuous region of the genome.

Input: files.txt - tab-delimited file of the format
        name1\tfile1.sam
        ....
        namen\tfilen.sam
            where filei.sam is a sam-formatted file of mapped reads
    chrFile - chr file containing nothing (no header) except the sequence of the chromosomes. No line breaks.
    chr - chr, must be the exact string as appears in the chr column of the sam
           file
    left - left boundary of region to plot
    right - right boundary of region to plot

Output: plot.pdf - plot of libraries over the region of interest

In addition, there are several options. For more info, call this script without
    any arguments.

run as python genePlotter.py [options] inFiles.txt chr left right outPrefix
EDIT: Added option to plot as bar graph
EDIT: 10/25/2016 - JOSH added option to highlight some features
EDIT: 11/8/2016 - JOSH added option to plot heat map of read frames
EDIT: Go to line 518 to change Y-axis
"""
import sys, common, csv, re, collections, os, pyx
from random import shuffle
from logJosh import Tee
from optparse import OptionParser
from pyx import *

def initOptions(parser):
    """Given OptionParser(), will initialize several options and return the parser"""
    #read display types
    parser.add_option("-d", "--display", type="str", nargs=1, dest="displayType",
        default="pileUp",
        help="will display reads as pileUp, barGraph, or frameHeatMap")
   
    #color options (multi v single) 
    parser.add_option("-c", "--color", type="str", nargs=1, dest="color",
        default="multiColor",
        help="will display reads as multiColor or singleColor")
    
    #gtf-formatted gene annotation file
    parser.add_option("-g", "--genes", type="str", nargs=1, dest="geneFile",
        help="will draw genes overlapping with boundaries on bottom of plots")
    
    #gtf-formatted features annotation file
    parser.add_option("-f", "--features", type="str", nargs=1, dest="featuresFile",
        help="will draw features overlapping boundaries on bottom of plots")
    
    #unique reads or not
    parser.add_option("-u", "--unique", type="str", nargs=1, dest="unique",
                      default=False,
        help="will only show unique reads if set to True")
    
    #plotwidth
    parser.add_option("-p", "--plotwidth", type="int", nargs=1, dest="plotWidth",
        default=800,
        help="will draw the plot this many pixels wide")
    
    #verticalSpacing
    parser.add_option("-a", "--axisheight", type="int", nargs=1, dest="verticalSpacing",
        default=10,
        help="will draw the plots this many pixels apart.")
    
    #sequence
    parser.add_option("-s", "--sequence", type="str", nargs=1, dest="sequence",
        default=False,
        help="will draw sequence on bottom if set to True.")
    
    return parser

def getOverlap(exon1,exon2):
    """Given two tuples exoni=[lefti,righti], will return the intersection"""
    if exon1[0]>exon2[1]:
        return None
    elif exon1[0]<=exon2[1]:
        if exon1[1]<exon2[0]:
            return None
        elif exon1[1]>=exon2[0]:
            return ([max([exon1[0],exon2[0]]),min([exon1[1],exon2[1]])])

def getGenes(geneFile,chr,bounds):
    """Will return all genes that overlap in any way on chr within bounds"""
    aa={}
    with open(geneFile,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row)>1:
                currChr=row[0]
                currLeft=int(row[3])
                currRight=int(row[4])
                if currChr==chr:
                    overlap=getOverlap([currLeft,currRight],bounds)
                    if overlap:
                        gene_id=row[8].split('"')[1]
                        if gene_id not in aa:
                            #get common name
                            if 'gene_name' in row[8]:
                                name=row[8].split('gene_name')[1].split('"')[1]
                            else:
                                name=gene_id
                            aa[gene_id]={'name':name,'strand':row[6],'CDS':[],'exon':[]}
                        
                        if row[2] in aa[gene_id]:
                            aa[gene_id][row[2]].append(overlap)
    return aa

def getFirstNonConsecutive(list1):
    """given a list of numbers, will return the number after the lowest consecutive stretch of numbers"""
    if list1==[]:
        return 1
    list1.sort()
    for ii in range(len(list1)-1):
        if list1[ii+1]>list1[ii]+1:
            return list1[ii]+1
    return list1[-1]+1

def getGeneRows(genes,scaleFactor,leftBound):
    """given genes={gene_id:{...,'exon':[exon_list]}}, will add a key 'row' to each gene_id
    that should put each gene far enough away from other genes on the same row"""
    for gene in genes:
        left,right=min([entry[0] for entry in genes[gene]['exon']]),\
                    max([entry[1] for entry in genes[gene]['exon']])
        genes[gene]['space']=[(left-leftBound)*scaleFactor-100,(right-leftBound)*scaleFactor]#subtract 100 for names
        placed=dict((gene_id,genes[gene_id]) for gene_id in genes if 'row' in genes[gene_id])
        #now compare to all other placed
        othersRows=[]
        for gene_id in placed:
            if getOverlap(placed[gene_id]['space'],genes[gene]['space']):
                othersRows.append(placed[gene_id]['row'])
        genes[gene]['row']=getFirstNonConsecutive(othersRows)

def drawGenes(c,geneFile,plotWidth,chr,bounds,rowFactor=0):
    """Given canvas c, gtf-formatted geneFile, width of plot, chr, and bounds, will draw 
    gene models across the bottom, with arrows indicating strand"""
    #will isolate and return all genes overlapping this region
    genes=getGenes(geneFile,chr,bounds)
    print '%s genes are being drawn...'%len(genes)
    
    #figure out nt multiplier
    scaleFactor=plotWidth/(float(bounds[1]-bounds[0]))
    rowScaleFactor=20
    rowSpaceFactor=25
    
    #now figure out where to put everyone
    getGeneRows(genes,scaleFactor,bounds[0])
    
    #now do the drawing
    for gene in genes:
        space=genes[gene]['space']
        row=(1+genes[gene]['row'])*-rowSpaceFactor+rowFactor#rowFactor is used if I want to plot something later, and shift it down/up relative to other stuff
        c.stroke(path.line(space[0]+100,row,space[1],row))
        c.text(space[0]+75,row,genes[gene]['name'],[text.halign.boxright,text.valign.middle])
        #the -0.5 and +1 on the next line are meant to make the exon/CDS boundaries line up exactly with the sequence
        for exon in genes[gene]['exon']:
            c.fill(path.rect((exon[0]-0.5-bounds[0])*scaleFactor,row-rowScaleFactor/4.,(exon[1]+1-exon[0])*scaleFactor,rowScaleFactor/2.),[color.rgb.black])
        for CDS in genes[gene]['CDS']:
            c.fill(path.rect((CDS[0]-0.5-bounds[0])*scaleFactor,row-rowScaleFactor/2.,(CDS[1]+1-CDS[0])*scaleFactor,rowScaleFactor/1.),[color.rgb.black])

def addBlock(chrFile,currPosit,blockLength,readBlock):
    """Will go through readSeq from readPosit to blockLength+readPosit and check string ID compared to
    currPosit in chrFile"""
    aa=[]
    with open(chrFile,'r') as f:
        f.seek(currPosit-1)
        genomeBlock=f.read(blockLength)
    
    for ii in range(blockLength):
        if readBlock[ii].lower()==genomeBlock[ii].lower():
            match=1
        else:
            match=0
        aa.append([currPosit+ii,match])
    return aa

def recoverBlocksFromCIGAR(readStart,Cigar,bounds,chrFile,readSeq):
    """Given the location a read starts, its Cigar string, and some bounds, will:
    (1) check if the read or any portion of it falls within bounds
    and if so,
    (2) return a list of [[position,match]] where position is the absolute chromosomal
    position and match indicates a match to that position in chrFile"""
    #first process the CIGAR string
    cigarSplit=re.findall('(\d+|[a-zA-Z]+)', Cigar)
    cigarSplit=[[int(cigarSplit[ii]),cigarSplit[ii+1]] for ii in range(0,len(cigarSplit),2)]
    
    #quick and dirty check
    total=sum([entry[0] for entry in cigarSplit])
    if (readStart>=bounds[1]) or (readStart+total<=bounds[0]):
        return 'na','na'
    
    #Joshua added this line circa 6/1/2018 to filter out reads that are more unmapped than mapped
    tempDict=collections.defaultdict(int)
    for entry in cigarSplit:
        tempDict[entry[1]]+=entry[0]
    sVal=tempDict['S']
    mVal=tempDict['M']
    totVal=float(sVal+mVal)
    if totVal>0:
        if sVal/totVal>0.1:
            return 'na','na'
    
    #if you're still here, that means your read overlaps with the region of interest. Time to do more work...
    #blockStarts=[0]
    #blockSizes=[0]
    currPosit=readStart
    readPosit=0
    blocks=[]
    #print Cigar
    #print cigarSplit
    for entry in cigarSplit:
        if entry[1]=='M':#then it's aligned to the genomic sequence
            blocks+=addBlock(chrFile,currPosit,entry[0],readSeq[readPosit:readPosit+entry[0]])
            currPosit+=entry[0]
            readPosit+=entry[0]
        elif entry[1]=='I':#I will ignore this for the block portion--the read contains a sequence that's not present in the reference
            readPosit+=entry[0]#1
        elif entry[1]=='S':#seq not present in reference
            readPosit+=entry[0]
        elif entry[1] in ['N','D']:#there is a gap in the read relative to the reference sequence
            currPosit+=entry[0]
    
    return 1,blocks

def checkUnique(NMstring):
    NMlist=NMstring.split(':')
    if NMlist[-1]=='1':
        return True
    else:
        return False

def getReads(fileName,chrFile,chr,bounds,uniQue=False):
    """Will return all reads on chr within range bounds. If uniQue is set to true, will only
    include uniquely mapping reads"""
    aa={'+':[],'-':[]}
    ##the following variable will keep track of the number of unique reads
    #readCt=0.
    #print 'Filtering read lengths to be [15,19).'
    with open(fileName,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row)>2:
                if (row[2]==chr):# and (len(row[9]) in range(15,19)):
                    isItInsideBounds,blocks=recoverBlocksFromCIGAR(int(row[3]),row[5],bounds,chrFile,row[9])
                    isItInsideBounds2,blocks2=0,False
                    if row[6]=='=':#then it's a paired end read file
                        row2=next(freader)#this should read in the next line of the file
                        if row2[0]==row[0]:
                            isItInsideBounds2,blocks2=recoverBlocksFromCIGAR(int(row2[3]),row2[5],bounds,chrFile,row2[9])
                            #now combine the blocks
                            if blocks=='na' or blocks2=='na':
                                if blocks!='na':
                                    combined=blocks
                                elif blocks2!='na':
                                    combined=blocks2
                                else:
                                    combined='na'
                            else:
                                combined=blocks+blocks2
                            if combined!='na':
                                blockPositions=dict((entry[0],1) for entry in combined)
                                blockPositions=blockPositions.keys()
                                blockPositions.sort()
                                
                                if blocks!='na':
                                    blocks=dict((entry[0],entry[1]) for entry in blocks)
                                if blocks2!='na':
                                    blocks2=dict((entry[0],entry[1]) for entry in blocks2)
                                tempBlocks=[]
                                for position in blockPositions:
                                    if blocks!='na':
                                        if position in blocks:
                                            tempBlocks.append([position,blocks[position]])
                                        elif position in blocks2:
                                            tempBlocks.append([position,blocks2[position]])
                                        else:
                                            print 'error: position not found'+row, sys.exit()
                                blocks=tempBlocks
                    
                    #if row[12].split(':')[-1]=='1':
                    #    readCt+=1
                    if isItInsideBounds==1 or isItInsideBounds2==1:
                        yesItsUnique=True
                        if uniQue:
                            yesItsUnique=checkUnique(row[11])
                        #print yesItsUnique
                        if yesItsUnique:
                            if int(row[1])&16!=0:
                                strand='-'
                            else:
                                strand='+'
                                #print row[9]
                                if row[6]=='=':
                                    strand='-'
                            if blocks!=[]:#one reason it might be empty is due to lack of a cigar string
                                aa[strand].append(blocks)
                                
                                #the following commented out lines can be useful to finding errors in the blocks--where I've incorrectly parsed the cigar string
                                """
                                if float(sum([entry[1] for entry in blocks]))/len(blocks)<0.5:
                                    print row
                                    print blocks
                                    print float(sum([entry[1] for entry in blocks]))/len(blocks)
                                    sys.exit()
                                """
    #return aa, readCt
    """
    N=100
    print 'Subsetting to look at %s random reads.'%(N)
    shuffle(aa['+'])
    shuffle(aa['-'])
    aa['+']=aa['+'][:N]
    aa['-']=aa['-'][:N]
    """
    return aa

def blockOverlap(left,right,placed):
    """given left and right, bounds of a new block, will return the lowest nonnegative number
    different than any block from placed=[[left1,right1,num1]] that overlaps with left,right"""
    overlap=[0]#this way it will always default back to row=1 if there are reads place higher up.
    for block in placed:#I add an extra nt of leeway on either side--that way reads don't abut next to each other
        if block[0]>=right+1:
            pass
        elif block[1]<=left-1:
            pass
        else:
            overlap.append(block[2])
    return getFirstNonConsecutive(overlap)

def assignReadsToRows(readList):
    """given readList=[[[position,1/0],...],...], will give each read block a row s.t. two read
    blocks will not overlap"""
    aa=[]
    placedBlocks=[]
    for block in readList:
        left,right=block[0][0],block[-1][0]
        row=blockOverlap(left,right,placedBlocks)
        aa.append([row,block])
        placedBlocks.append([left,right,row])
    return aa

def drawMinus(c,blocks,leftBound,rightBound,verticalPosition,pixelPerNt):
    """given blocks=[[row,[position,num],...,[position,num]],[row,...]], will first find the largest row,
    then it will assume row=1 is verticalPosition+largest_row. Will draw each block in the vertical spot
    given by row, and its horizontal spot given by positions. If num=0, will draw a black box. Otherwise
    red."""
    try:
        numRows=max([entry[0] for entry in blocks])
    except ValueError:
        return verticalPosition
    readHeight=3.#pixels
    verticalCeiling=verticalPosition+numRows*readHeight
    
    for block in blocks:
        row=verticalCeiling-block[0]*readHeight
        prevNt=block[1][0][0]-1
        currCount=0
        for nt in block[1]:
            if nt[0]>=leftBound and nt[0]<=rightBound:
                graphPosition=(nt[0]-leftBound)*pixelPerNt
                if prevNt!=nt[0]-1:
                    if prevNt>=leftBound:
                        c.stroke(path.line((prevNt-leftBound)*pixelPerNt-pixelPerNt/2.,row+readHeight/2.,
                                    graphPosition+pixelPerNt/2.,row+readHeight/2.),
                                    [style.linewidth.THIN,color.cmyk(0,0.8,1.,0)])
                        prevNt=nt[0]
                else:
                    if nt[1]==1:
                        theColor=color.cmyk(0,0.8,1,0)
                    elif nt[1]==0:
                        theColor=color.cmyk(0,0,0,1)
                    c.fill(path.rect(graphPosition-pixelPerNt/2.,row,
                        pixelPerNt,readHeight-1),[theColor])
                    currCount+=1
                    prevNt=nt[0]
    return verticalCeiling

def drawPlus(c,blocks,leftBound,rightBound,verticalPosition,pixelPerNt):
    """given blocks=[[row,[position,num],...,[position,num]],[row,...]], will first find the largest row,
    then it will assume row=1 is verticalPosition+largest_row. Will draw each block in the vertical spot
    given by row, and its horizontal spot given by positions. If num=0, will draw a black box. Otherwise
    blue."""
    try:
        numRows=max([entry[0] for entry in blocks])
    except ValueError:
        return verticalPosition
    readHeight=3.#pixels
    verticalCeiling=verticalPosition+numRows*readHeight
   
    for block in blocks:
        row=verticalPosition+block[0]*readHeight
        prevNt=block[1][0][0]-1
        for nt in block[1]:
            if nt[0]>=leftBound and nt[0]<=rightBound:
                graphPosition=(nt[0]-leftBound)*pixelPerNt
                if prevNt!=nt[0]-1:
                    if prevNt>=leftBound:
                        c.stroke(path.line((prevNt-leftBound)*pixelPerNt-pixelPerNt/2.,row+readHeight/2.,
                                    graphPosition+pixelPerNt/2.,row+readHeight/2.),
                                    [style.linewidth.THIN,color.cmyk(1.,0.5,0,0)])
                        prevNt=nt[0]
                else:
                    if nt[1]==1:
                        theColor=color.cmyk(1.,0.5,0,0)
                    elif nt[1]==0:
                        theColor=color.cmyk(0,0,0,1)
                    c.fill(path.rect(graphPosition-pixelPerNt/2.,row,
                        pixelPerNt,readHeight-1),[theColor])
                    prevNt=nt[0]
    return verticalCeiling

def getNtColot(nt):
    if nt.lower()=='a':
        return color.cmyk(0,0.8,1,0)
    elif nt.lower()=='c':
        return color.cmyk(0,0.5,1,0)
    elif nt.lower()=='g':
        return color.cmyk(0.97,0,0.75,0)
    elif nt.lower()=='t':
        return color.cmyk(1,0.5,0,0)
    else:
        return color.cmyk(0,0,0,1)

def drawSequence(c,chrFile,leftBound,rightBound,plotWidth):
    """will draw the nt sequence, one nt at a time"""
    ntPerPixel=float(plotWidth)/(rightBound-leftBound)
    with open(chrFile,'r') as f:
        for ii in range(leftBound,rightBound+1):
            f.seek(ii-1)
            currNt=f.read(1)
            
            c.text(ntPerPixel*(ii-leftBound),-10,currNt,[text.halign.boxcenter,text.valign.middle,getNtColot(currNt)])

def prepareReadsForBarGraph(blockList,leftBound,rightBound,readCt,strand='+'):
    """Given a blockList=[[[position1,match],...,[positionn,match]],[...],...[]],
    will return a list of [(position,ct)], where ct yields the number of reads
    whose 5'end lies at that position. If strand is set to reverse, will look at the
    3'end of a read.
    will return a list of points, where every points from leftBound to rightBound,
    inclusive, has some value"""
    aa=collections.defaultdict(int)
    for block in blockList:
        ct=1.
        if strand=='+':
            position=block[0][0]
        elif strand=='-':
            position=block[-1][0]
            ct=-1
        aa[position]+=ct
        #the below two lines will count once for every position in the block
        #for subBlock in block:
        #    aa[subBlock[0]]+=ct/len(subBlock)
    
    bb=[]
    normFactor=readCt/1000000.
    #written={}
    for ii in range(leftBound,rightBound+1):
        if True:#aa[ii]!=0:
            bb.append([ii,aa[ii]/normFactor])
            #written[ii]=1
    
    #the following lines are to reduce the number of points plotted for large views
    """
    for ii in range(leftBound+1,rightBound):
        if ii in written:
            if ii-1 not in written:
                bb.append([ii,0])
            if ii+1 not in written:
                bb.append([ii,0])
    """
    
    return bb

def plotBarGraph(barGraphData,verticalPosition,Height,leftBound,rightBound,plotWidth):
    """Given barGraphData, output of prepareForBarGraph,will make a bargraph starting at
    verticalPosition, of height Height, with minimum leftBound and rightBound on the
    abscissa"""
    g=graph.graphxy(width=plotWidth,ypos=verticalPosition,height=Height,
        x=graph.axis.bar(painter=None),
        y=graph.axis.linear(title='rpm',fallbackrange=5))
    
    unit.set(defaultunit="pt")
    
    """
    #the following lines were for looking at only the plus strand, quick and dirty
    plusData=[]
    for pt in barGraphData['+']:
        if pt[1]>0:
            plusData.append(pt)
        else:
            plusData.append([pt[0],0.00001])
    """
    """
    print barGraphData['-'][:10]
    g.plot(graph.data.points(barGraphData['-'],x=1,y=2),
        [graph.style.histogram()])
    """
    g.plot(graph.data.points(barGraphData['+'],xname=0,y=2),
            [graph.style.barpos(fromvalue=0),
            graph.style.bar(barattrs=[color.cmyk(1,0.5,0,0),
                deco.stroked([color.cmyk(1,0.5,0,0)])])])
    
    g.plot(graph.data.points(barGraphData['-'],xname=0,y=2),
            [graph.style.barpos(fromvalue=0),
            graph.style.bar(barattrs=[color.cmyk(0,0.8,1,0),
                deco.stroked([color.cmyk(0,0.8,1,0)])])])
    
    return g

def plotBarGraph2(barGraphData,verticalPosition,Height,pixelPerNt):
    """second version plots bar graph by drawing rectangles"""
    theMin=min([entry[1] for entry in barGraphData['-']])
    theMax=max([entry[1] for entry in barGraphData['+']])
    #comment out the below three lines to remove this
    theMin=0
    theMax=35
    print 'Min and Max ordinate scale hard-coded as %s and %s.'%(theMin,theMax)
    
    c=canvas.canvas()
    unit.set(defaultunit="pt")
    
    #setting reverse=1 and switching start/end causes the scale to be written to the left of the
    #axis instead of the right (potentially overlapping with the reads)
    #g=graph.axis.pathaxis(path.line(0,verticalPosition,0,verticalPosition+Height),
    #                    graph.axis.linear(min=theMin,max=theMax))
    g=graph.axis.pathaxis(path.line(0,verticalPosition+Height,0,verticalPosition),
                        graph.axis.linear(reverse=1,min=theMin,max=theMax))
    
    c.insert(g)
    
    #determine where zero is
    if theMin>0 or theMax<0:
        print 'Uh-oh. Your rpm doesn\'t seem to be centered around 0...'
    axisSize=float(theMax-theMin)
    middleLine=-(theMin/axisSize)*Height
    verticalScaler=Height/axisSize
    
    #draw the + strand
    #print 'Suppressing the + strand'
    
    for ii in range(len(barGraphData['+'])):
        entry=barGraphData['+'][ii]
        if entry[1]!=0.:#only draw nonzero points
            c.stroke(path.rect(pixelPerNt*ii,verticalPosition+middleLine,pixelPerNt,entry[1]*verticalScaler),
                [color.rgb(0.5,0.5,0.5),
                deco.filled([color.rgb(0.5,0.5,0.5)])])
    #            [color.cmyk(1,0.5,0,0),
    #            deco.filled([color.cmyk(1,0.5,0,0)])])
    
    #draw the - strand
    #print 'Suppressing the - strand...'
    
    for ii in range(len(barGraphData['-'])):
        entry=barGraphData['-'][ii]
        if entry[1]!=0.:#only draw nonzero points
            c.stroke(path.rect(pixelPerNt*ii,verticalPosition+middleLine,pixelPerNt,entry[1]*verticalScaler),
                [color.rgb(0/255.,114/255.,178/255.)])
    #            [color.cmyk(0,0.8,1,0)])
    
    #draw the zeroLine
    c.stroke(path.line(1,verticalPosition+middleLine,ii*pixelPerNt,verticalPosition+middleLine))
    
    return c

def getReadsInRegion(fileName,chr,leftBound,rightBound):
    """Given a fileName (string), a chromosome, and a left and right bound,
    will extract all reads mapping in that region, stored in fileName.temp.region.chr.leftBound.rightBound"""
    
    bamFile='.'.join(fileName.split('.')[:-1]+['bam'])
    fileCore='.'.join(fileName.split('.')[:-1])
    
    #first, look for the .bam file and make it if it isn't there
    if not os.path.isfile(bamFile):
        print 'Could not file .bam file for %s, now making it in %s.'%(fileName,bamFile)
        #convert to sorted bam file
        os.system('samtools view -bS %s | samtools sort - %s'%(fileName,fileCore))#no need to add .bam extension b/c samtools will anyways
        #index the bam file
        os.system('samtools index %s'%(bamFile))
    
    #now extract reads in each region
    os.system('samtools view %s %s:%s-%s > %s.temp.region.%s.%s.%s'%(bamFile,chr,str(leftBound),str(rightBound),fileCore,chr,str(leftBound),str(rightBound)))

def getReadCtFromSam(fileName):
    """Will count the number of reads in a .sam file"""
    lineCt=os.popen('wc -l %s'%fileName).read()#total number of lines
    lineCt=int(lineCt.split()[0])
    
    #all lines starting with '@', followed by any combination of these two lists of letters counts as a header
    headerCt=int(os.popen("grep '^\@[H,S,R,P,C][D,Q,D,O]' %s|wc -l"%(fileName)).read().strip())
    
    lineCt2=os.popen("grep '\t255\t' %s | wc -l "%fileName).read()
    lineCt2=int(lineCt2.split()[0])
    print lineCt2
    print 'rpm calc now on mapq 255 only'
    return lineCt2
    #return lineCt-headerCt#readCt is the number of lines minus the header

def prepareReadsForFrameHeatMap(reads,leftBound,rightBound,strand):
    """reads is a list of lists, where each list is a tuple of [position,match]
    where position denotes a bp that the read overlapped with, and match denotes
    a match (or not, 1 or 0) between the read and the genome. This function
    will take the left-most position in each block as the 5'end. It will return
    a dictionary with keys as positions, and values a list of [0,1,2] where
    the value at each position is the number of reads in each frame."""
    aa={}
    for block in reads:
        if len(block)==28:
            if strand=='+':
                position=min([entry[0] for entry in block])
            elif strand=='-':
                position=max([entry[0] for entry in block])
            position2=3*(position/3)
            if leftBound<=position2<=rightBound:
                if position2 not in aa:
                    aa[position2]=[0,0,0]
                aa[position2][position%3]+=1
    return aa

def plotHeatMap(heatMapData,leftBound,verticalPosition,Height,pixelPerNt):
    """will plot heat map as a series of rectangles, given heatMapData=
    {position:[0,1,2]} where 0,1, and 2 positions are read counts."""
    c=canvas.canvas()
    unit.set(defaultunit="pt")
    #
    #write out the frames
    c.text(-10*pixelPerNt,verticalPosition+Height*(-1/6.+1/3.),'Frame 0',[text.halign.boxright,text.valign.middle])
    c.text(-10*pixelPerNt,verticalPosition+Height*(-1/6.+2/3.),'Frame +1',[text.halign.boxright,text.valign.middle])
    c.text(-10*pixelPerNt,verticalPosition+Height*(-1/6.+3/3.),'Frame -1',[text.halign.boxright,text.valign.middle])
    N=5
    print 'Using a cutoff of at least %s reads per codon.'%(N)
    for position in heatMapData:
        tempData=heatMapData[position]
        total=float(sum(tempData))
        if total>N:
            for ii in range(3):
                theFrac=tempData[ii]/total
                c.fill(path.rect(pixelPerNt*(position-leftBound),verticalPosition+ii*Height/3.,
                        pixelPerNt*3,Height/3.),
                        [color.cmyk(theFrac,0.5*theFrac,0,0)])
    
    #draw a line across the bottom
    c.stroke(path.line(0,verticalPosition-Height/6.,800,verticalPosition-Height/6.))
    return c

def main(args):
    #first initialize options
    parser=initOptions(OptionParser())
    
    (options,args)=parser.parse_args()
    #Note: the above function will parse directly from sys.argv! Which means if you're calling this from
    #within another python script, it will show you the args from the first script's call, which is unlikely
    #to be viable input for genePlotter!
    
    #parse the options
    displayType,color,geneFile=options.displayType,options.color,options.geneFile
    featuresFile=options.featuresFile
    plotWidth,verticalSpacing=options.plotWidth,options.verticalSpacing
    unique=options.unique
    sequence=options.sequence
    
    readFiles,chrFile,chr,leftBound,rightBound,outPrefix=args[0:]
    
    #pass bounds
    leftBound,rightBound=int(leftBound),int(rightBound)
    
    #initialize the canvas
    c=canvas.canvas()
    unit.set(defaultunit="pt")
    
    #draw gene models if present
    if geneFile:
        drawGenes(c,geneFile,plotWidth,chr,(leftBound,rightBound))
    
    #draw features (mutations?) if present
    if featuresFile:
        drawGenes(c,featuresFile,plotWidth,chr,(leftBound,rightBound),rowFactor=-100)
    
    #write sequence if desired
    if sequence:
        drawSequence(c,chrFile,leftBound,rightBound,plotWidth)
    
    #plot key for heat map
    if displayType=='frameHeatMap':
        for ii in range(0,5):
            c.fill(path.rect(-50,-50-10*ii,10,10),[pyx.color.cmyk(1.*ii/4.,0.5*ii/4.,0,0)])
            c.text(-55,-50-10*(ii-0.5),str(100*ii/4),[text.halign.boxright,text.valign.middle])
    
    #now commence the plotting
    verticalPosition=0.
    c.stroke(path.line(0,verticalPosition-0.5,plotWidth,verticalPosition-0.5),
                        [style.linestyle.dashed])
    pixelPerNt=float(plotWidth)/(rightBound-leftBound)
    #draw a scale bar
    N=10#bps
    c.fill(path.rect(plotWidth-N*pixelPerNt,0,pixelPerNt*N,5),[pyx.color.cmyk.black])
    c.text(plotWidth-(N/2.)*pixelPerNt,-5,'%s bp'%(N),[text.halign.boxcenter,text.valign.top])
    if os.path.isfile(readFiles):
        with open(readFiles,'r') as f:
            freader=csv.reader(f,delimiter='\t')
            for row in freader:
                name,fileName=row[0],row[1]
                
                #extract reads in region
                #getReadsInRegion(fileName,chr,leftBound,rightBound)
                #tmpFileName='.'.join(fileName.split('.')[:-1])+'.temp.region.%s.%s.%s'%(chr,str(leftBound),str(rightBound))
                tmpFileName=fileName
                
                #let the user know what's going on...
                print 'Working on %s'%name
                
                #first pull out all the reads in the region
                reads=getReads(tmpFileName,chrFile,chr,(leftBound,rightBound),uniQue=unique)
                
                #now get readCt
                print 'It is assumed only uniquely-mapping reads are included in the input file. If not, the rpm calculation is inaccurate.'
                readCt=getReadCtFromSam(fileName)
                
                #reads are now in the format {strand:[[position,match]]}
                
                #let's put the reads in rows now.
                #first prepareReadsFor the bar graph
                if displayType=='barGraph':
                    barGraphData={'+':prepareReadsForBarGraph(reads['+'],leftBound,rightBound,readCt,strand='+'),
                              '-':prepareReadsForBarGraph(reads['-'],leftBound,rightBound,readCt,strand='-')}
                
                #now it's time to do some drawing!
                if displayType=='pileUp':
                    for strand in reads:
                        reads[strand]=assignReadsToRows(reads[strand])
                    
                    #first draw the minus strand
                    c.stroke(path.line(0,verticalPosition-0.5,plotWidth,verticalPosition-0.5),
                                [style.linestyle.dashed])
                    verticalPosition=drawMinus(c,reads['-'],leftBound,rightBound,verticalPosition,pixelPerNt)
                    #as long as we're here, let's draw the title
                    c.text(-10,verticalPosition,name,[text.halign.boxright,text.valign.middle,text.size(2)])
                    
                    verticalPosition=drawPlus(c,reads['+'],leftBound,rightBound,verticalPosition,pixelPerNt)
                    verticalPosition+=verticalSpacing
                    
                    #now draw a dotted line to separate the libraries
                    #c.stroke(path.line(0,verticalPosition,plotWidth,verticalPosition),
                    #            [style.linestyle.dashed])
                elif displayType=='barGraph':
                    #barGraph=plotBarGraph(barGraphData,verticalPosition,verticalSpacing*20,leftBound,rightBound,plotWidth)
                    barGraph=plotBarGraph2(barGraphData,verticalPosition,verticalSpacing*20,pixelPerNt)
                    c.text(-70,verticalPosition+20*verticalSpacing/2.,name,[text.halign.boxright,text.valign.middle,text.size(2)])
                    verticalPosition+=verticalSpacing*20+verticalSpacing
                    c.insert(barGraph)
                
                elif displayType=='frameHeatMap':
                    #This will display the data as the fraction of reads at each position in each frame.
                    #will arbitrarily start the frame relative to the left end of the chromosome.
                    #more complicated things (e.g. changing frames with exons) will be left for future
                    #work
                    pickAStrand='+'
                    print 'Going with the %s strand for frame heat map'%(pickAStrand)
                    heatMapData=prepareReadsForFrameHeatMap(reads[pickAStrand],leftBound,rightBound,pickAStrand)
                    #heatMapData is a dict of {position:[0,1,2]} where cts in each position indicate the number of
                    #reads mapping in that frame. position is rounded down every third base.
                    heatMap=plotHeatMap(heatMapData,leftBound,verticalPosition,verticalSpacing*3,pixelPerNt)
                    c.text(-70,verticalPosition+3*verticalSpacing/2.,name,[text.halign.boxright,text.valign.middle,text.size(2)])
                    verticalPosition+=verticalSpacing*4
                    c.insert(heatMap)
            
            #get rid of the intermediate file
            #os.remove(tmpFileName)
    else:
        pass
    c.writePDFfile(outPrefix)
    #c.writeGSfile(filename=outPrefix+'.jpg')

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
