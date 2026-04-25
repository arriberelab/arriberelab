"""
Joshua Arribere Sept 27, 2013

Script to make a metaGene Plot of an arbitrary number of libraries

Input: annots.gtf - gtf formatted annotation file (used to get txt lengths)
    outPrefix - will be used to name the outfile
    infile.joshSAM - as output from assignReadsToGenes.py
    geneList.txt - tab-delimited file of format
        namei   filei.txt
        where filei contains a line-delimited list of genes. namei will appear as the label

Output: Plot of metaStart, metaStop

run as python metaStartStop.py annots infile.joshSAM geneList.txt outPrefix
EDIT: Oct 16, 2013 - JOSH edited to input a list of your favorite genes
EDIT: Oct 20, 2019 - JOSH edited to do for an arbitrary number of joshSAM files
EDIT: Jun 1, 2021 - JOSH edited to bring into python3, one jam file, one gene file
"""
import sys, common, csv, collections, numpy, senseAboutAS
from logJosh import Tee
from pyx import *

def getLengths(annots):
    """Will return two dicts of total txt length and cds length"""
    a={}
    b=collections.defaultdict(list)
    with open(annots,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row)>1:
                featureType=row[2]
                if featureType in ['exon','CDS']:# and row[1]=='protein_coding':
                    transcript_id = row[8].split('transcript_id')[1].split('"')[1]
                    if transcript_id not in a:
                        a[transcript_id]={'strand':row[6],'exon':[],'CDS':[]}
                    a[transcript_id][featureType].append([int(row[3]),int(row[4])])
                    
                    gene_id=row[8].split('"')[1]
                    b[gene_id].append(transcript_id)
    
    exon,cds={},{}
    for txt in a:
        exon[txt]=sum([entry[1]-entry[0] for entry in a[txt]['exon']])
        cds[txt]=sum([entry[1]-entry[0] for entry in a[txt]['CDS']])
    return exon, cds, b

def getGeneReadCts(inFile):
    """Will sum up the total number of reads mapping to each gene_id (row[5])"""
    a=collections.defaultdict(float)
    with open(inFile,'r') as f:
        f.readline()
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row)>5:
                if row[8].endswith(':S'):
                    a[row[8].split(':')[0]]+=1.
    b={}
    N=10
    print('Filtering for at least %s read counts...'%(N))
    for key in a:
        if a[key]>=N:#edit this line to include a cutoff for number of reads
            b[key]=a[key]
    return b

def getNormFactor(geneReadCts,txtLengths,txtGroups):
    """For each gene in geneReadCts, will look up its associated txts in txtGroups, and then their lengths in txtLengths.
    Will return the the read count divided by the average txt length"""
    a=dict((gene_id,geneReadCts[gene_id]/numpy.average([txtLengths[txt] for txt in txtGroups[gene_id]])) for gene_id in geneReadCts)
    return a

def metaAverage(metaDict,N):
    """metaDict={key:{position:[list]}} will loop through metaDict and change to
    {key:{position:sum(list)/N}}"""
    #print 'averaging over the number of genes/txts overlapping a given position, not the total number of genes/txts'
    for key in metaDict:
        for position in metaDict[key]:
            metaDict[key][position]=sum(metaDict[key][position])/float(N)
            #metaDict[key][position]=sum(metaDict[key][position])/float(len(metaDict[key][position]))

def getMetaStartStop(inFile,normFactors,selectGenes=None):
    """Will get distribution of reads about start/stop codon. To do this it will first:
    (1) Add up the total number of reads mapping to each gene
    (2) Will find the normalization factor s.t. the mean read count across a gene per kb is 1
        b/c txt size is potentially variable across a group, I take the average of all the txts within a gene. Can't think of a better thing to do right now...
    (3) Will apply that normalization factor to each read
    (4) Concomitant with #3, will add read count to a position-specific read count dict
    I distinguish between genes (keys to txtGroups) and txts (vales in txtGroups)
    """
    if not selectGenes:
        selectGenes=normFactors
    
    ##print "Counting reads by their 3'ends"
    ##print "Adding a few nts at 3'end s.t. 0 is now the 1st nt of stop"
    #print 'Restricting to 26G'
    #print 'No read length/seq restrictions'
    print('Plotting by P-site')##P-site
    
    totalCt={}
    #Apply the normalization factor and record read positions
    metaStart,metaStop={'S':collections.defaultdict(list),'AS':collections.defaultdict(list)},{'S':collections.defaultdict(list),'AS':collections.defaultdict(list)}
    with open(inFile,'r') as f:
        f.readline()
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row)>5:
                assocTxts=row[9].split('|')
                numMapped=len(assocTxts)
                gene_id,SorAS=row[8].split(':')
                readSeq=row[6]
                readLength=len(readSeq)-2#-2 is done s.t. 0 is the 1st nt of the stop codon
                
                if gene_id in selectGenes and gene_id in normFactors:# and len(assocTxts)==1:# and readSeq[0]=='G' and readLength==26:
                    totalCt[gene_id]=1
                    for txt in assocTxts:
                        curr=txt.split(':')
                        #the following adds the readLength
                        #this counts reads by their 3'end
                        relStart=int(curr[1])+12#+readLength#P-site
                        relStop=int(curr[2])+12#+readLength#P-site
                        
                        #reads are indexed by their 5'end, whether they're antisense or sense. Thus if I count the entire length of
                        #the read, I need to move forward from the 5'end on the sense strand. But on the antisense strand, need to
                        #go upstream. Hence I introduce the variable flip.
                        if SorAS=='S':
                            flip=1
                        elif SorAS=='AS':
                            flip=-1
                        
                        #these next few lines are to round by codon
                        #relStart=(int(relStart)/3)*3
                        #relStop=(int(relStop)/3)*3
                        
                        #switch to these two lines if you don't want to count the entire length
                        metaStart[SorAS][relStart].append(1./normFactors[gene_id]/numMapped)
                        metaStop[SorAS][relStop].append(1./normFactors[gene_id]/numMapped)
                        #adding a few lines to print out reads around the stop codon to be certain of where they end
                        #if relStop==0 and readLength>20:
                        #    print readSeq
                        
                        #otherwise use these lines
                        #for ii in range(readLength):#since reads are indexed by the 5'end, this will 
                        #    metaStart[SorAS][relStart+ii*flip].append(1./normFactors[gene_id]/numMapped/readLength)
                        #    metaStop[SorAS][relStop+ii*flip].append(1./normFactors[gene_id]/numMapped/readLength)
    print(len(totalCt))
    metaAverage(metaStart,len(totalCt))
    metaAverage(metaStop,len(totalCt))
    
    return metaStart,metaStop

def processMeta(dict1,flip=0):
    """Will return dict1 as an ordered list of tuples. If flip will flip the sign"""
    a=list(dict1.keys())
    a.sort()
    if not flip:
        return [(key,dict1[key]) for key in a]
    else:
        return [(key,-dict1[key]) for key in a]


def mkStartStopPlot(metaData,metaDataMFG,outPrefix):
    """metaData={name:({'S':metaStartSense,'AS':metaStartAS},{'S':metaStopSense,'AS':metaStopAS})}.
    Will plot sense above x=0, antisense below x=0, Start plot on left and Stop plot on right"""
    xgridpainter = graph.axis.painter.regular(gridattrs=[attr.changelist([style.linestyle.dashed,
                                                                    style.linestyle.dotted])])
    if False:
        start=graph.graphxy(width=8,height=8,
                        x=graph.axis.linear(painter=xgridpainter,min=-30,max=30,
                                            title='Position Relative to Start Codon'),
                        y=graph.axis.linear(min=0,
                                            title='Normalized Ribo-seq Read Density'))
        stop=graph.graphxy(width=8,height=8,xpos=start.width*1.1,
                        x=graph.axis.linear(#painter=xgridpainter,
                                            min=-30,max=30,
                                            title='Position Relative to Stop Codon'),
                        y=graph.axis.linkedaxis(start.axes["y"]))
    if True:
        start=graph.graphxy(width=8,height=8,
                        x=graph.axis.linear(painter=xgridpainter,min=-30,max=30,
                                            title='Position Relative to Start Codon'),
                        y=graph.axis.linear(min=0,max=4,
                                            title='Normalized Ribo-seq Read Density'))
        stop=graph.graphxy(width=8,height=8,xpos=start.width*1.1,
                        x=graph.axis.linear(#painter=xgridpainter,
                                            min=-30,max=30,
                                            title='Position Relative to Stop Codon'),
                        y=graph.axis.linkedaxis(start.axes["y"]))
        
    #manually put lines every 3nts
    for ii in range(-60,40,3):
        stop.plot(graph.data.function("x(y)=%s"%(ii),min=-10,
                                    title=None),
                                    [graph.style.line([style.linestyle.dashed])])
    #
    ii=0
    metaStartSense=processMeta(metaData[0]['S'])
    start.plot(graph.data.points(metaStartSense,x=1,y=2),[graph.style.line([common.colors(ii)])])
    metaStartAS=processMeta(metaData[0]['AS'],flip=1)
    start.plot(graph.data.points(metaStartAS,x=1,y=2),[graph.style.line([common.colors(ii)])])
    ##Note: the naming of variables here is screwy...a relic of many different versions. But the output is correct.
    ##I hesitate to try and 'fix' the variable names at this point for fear of screwing up the output.
    metaStartSense=processMeta(metaData[1]['S'])
    #print(metaStartSense)
    if True:
        norm=[entry[1] for entry in metaStartSense if entry[0] in range(-30,-3)]
        stopSite=[entry[1] for entry in metaStartSense if entry[0] in range(-3,0)]
        print('Average read density in the 27 nts upstream of the stop codon: %s'%(numpy.average(norm)))
        print('Average read density in the 3 nts upstream of the stop codon: %s'%(numpy.average(stopSite)))
        print('Ratio of read density in the 3 nts upstream of the stop codon to the 27 nts upstream: %s'%(numpy.average(stopSite)/numpy.average(norm)))
            
    stop.plot(graph.data.points(metaStartSense,x=1,y=2),[graph.style.line([common.colors(ii)])])
    metaStartAS=processMeta(metaData[1]['AS'],flip=1)
    stop.plot(graph.data.points(metaStartAS,x=1,y=2),[graph.style.line([common.colors(ii)])])
    ##
    #print('Plotting all genes in color %s...'%(common.colors(ii)))
    ii+=1
    #print('Highlighting genes in color %s...'%(common.colors(ii)))
    metaStopSense=processMeta(metaDataMFG[0]['S'])
    start.plot(graph.data.points(metaStopSense,x=1,y=2),[graph.style.line([common.colors(ii)])])
    metaStopAS=processMeta(metaDataMFG[0]['AS'],flip=1)
    start.plot(graph.data.points(metaStopAS,x=1,y=2),[graph.style.line([common.colors(ii)])])
    ##
    metaStopSense=processMeta(metaDataMFG[1]['S'])
    stop.plot(graph.data.points(metaStopSense,x=1,y=2),[graph.style.line([common.colors(ii)])])
    metaStopAS=processMeta(metaDataMFG[1]['AS'],flip=1)
    stop.plot(graph.data.points(metaStopAS,x=1,y=2),[graph.style.line([common.colors(ii)])])
    ##
    c=canvas.canvas()
    c.insert(start)
    c.insert(stop)
    c.writePDFfile(outPrefix)

def main(args):
    annots,inFile,geneList,outPrefix=args[0:]
    
    #first get txt lengths
    txtLengths,cdsLengths,txtGroups=getLengths(annots)
    print('Parsed annotation file '+annots)
    #keys are txts like D46A8.9, C17B7.8a
    
    ##Parse gene list file
    #myFavoriteGenes=None
    myFavoriteGenes=common.parseGeneList(geneList)
    print('Parsed gene list '+geneList)
    ##keys are txts like B0272.2, C56C10.7b.3
    
    #Get read counts/gene
    geneReadCts=getGeneReadCts(inFile)
    #keys are gene names, not txts, like WBgene00001115, WBGene00000618
    
    #Get normalization factor
    normFactors=getNormFactor(geneReadCts,txtLengths,txtGroups)
    ##
    metaData=getMetaStartStop(inFile,normFactors)
    metaDataMFG=getMetaStartStop(inFile,normFactors,selectGenes=myFavoriteGenes)
    ##
    print('Plotting...')
    if myFavoriteGenes:
        outPrefix=outPrefix+'_%s'%(len(myFavoriteGenes))
    mkStartStopPlot(metaData,metaDataMFG,outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
