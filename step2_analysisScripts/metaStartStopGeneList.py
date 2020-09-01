#!/usr/bin/env python3
"""
Joshua Arribere Sept 27, 2013
Converted to python 3: Mar 24, 2020

Script to make a metaGene Plot of an arbitrary number of libraries

Input: annots.gtf - gtf formatted annotation file (used to get txt lengths)
    outPrefix - will be used to name the outfile
    infilei.jam - as output from assignReadsToGenes.py
    namei - will be used as a label for data in infilei.jam

Output: Plot of metaStart, metaStop

run as python3 metaStartStop.py annots.gtf outPrefix infile1.jam name1 GeneList.txt
Matt edited 8/31/20
"""
import sys, common, csv, collections, numpy
from logJosh import Tee
import pyx


def getLengths(annots):
    """Will return two dicts of total txt length and cds length"""
    a = {}
    b = collections.defaultdict(list)
    with open(annots, 'r') as f:
        freader = csv.reader(f, delimiter='\t')
        for row in freader:
            if len(row) > 1:
                featureType = row[2]
                if featureType in ['exon', 'CDS']:  # and row[1]=='protein_coding':
                    transcript_id = row[8].split('transcript_id')[1].split('"')[1]
                    if transcript_id not in a:
                        a[transcript_id] = {'strand': row[6], 'exon': [], 'CDS': []}
                    a[transcript_id][featureType].append([int(row[3]), int(row[4])])
                    
                    gene_id = row[8].split('"')[1]
                    b[gene_id].append(transcript_id)
                    
    exon, cds = {}, {}
    for txt in a:
        exon[txt] = sum([entry[1] - entry[0] for entry in a[txt]['exon']])
        cds[txt] = sum([entry[1] - entry[0] for entry in a[txt]['CDS']])
    return exon, cds, b


def getGeneReadCts(inFile):
    """Will sum up the total number of reads mapping to each gene_id (row[5])"""
    a = collections.defaultdict(float)
    with open(inFile, 'r') as f:
        freader = csv.reader(f, delimiter='\t')
        for row in freader:
            if len(row) > 8: #updated for jam files
                gene = row[8]
                genesplit = gene.split(':') #Matt updated for new format, had to get rid of S/AS information here too.
                a[genesplit[0]] += 1.
    return a


def getNormFactor(geneReadCts, txtLengths, txtGroups):
    """For each gene in geneReadCts, will look up its associated txts in txtGroups, and then their lengths in txtLengths.
    Will return the the read count divided by the average txt length"""
    a = dict((gene_id,
              geneReadCts[gene_id] / numpy.average([txtLengths[txt] for txt in txtGroups[gene_id]])) for gene_id in
             geneReadCts)
    return a


def metaAverage(metaDict, N):
    """metaDict={key:{position:[list]}} will loop through metaDict and change to
    {key:{position:sum(list)/N}}"""
    for key in metaDict:
        for position in metaDict[key]:
            metaDict[key][position] = sum(metaDict[key][position]) / float(N)


def getMetaStartStop(inFile, txtLengths, txtGroups, nmdReadLength=None, readLengthRestriction=None):
    """Will get distribution of reads about start/stop codon. To do this it will first:
    (1) Add up the total number of reads mapping to each gene
    (2) Will find the normalization factor s.t. the mean read count across a gene per kb is 1
        b/c txt size is potentially variable across a group, I take the average of all the txts within a gene.
            Can't think of a better thing to do right now...
    (3) Will apply that normalization factor to each read
    (4) Concomitant with #3, will add read count to a position-specific read count dict
    I distinguish between genes (keys to txtGroups) and txts (vales in txtGroups)
    """
    # Get read counts/gene
    geneReadCts = getGeneReadCts(inFile)
    
    # Get normalization factor
    normFactors = getNormFactor(geneReadCts, txtLengths, txtGroups)
    #print(normFactors)
    totalCt = {}
    # Apply the normalization factor and record read positions
    metaStart, metaStop = {'S': collections.defaultdict(list), 'AS': collections.defaultdict(list)}, \
                          {'S': collections.defaultdict(list), 'AS': collections.defaultdict(list)}
    inFrame = 0.
    outOfFrame = 0.
    with open(inFile, 'r') as f:
        freader = csv.reader(f, delimiter='\t')
        for row in freader:
            if len(row) > 8 and row[1] != 'chr': #Matt did this, the new jam files have a header. 08/12/20
                assocTxts = row[9:]
                numMapped = len(assocTxts)
                gene = row[8]                    #Matt added to fix compatability with new jsams. 8/12/20
                genesplit = gene.split(':')
                SorAS = genesplit[1]
                gene_id = genesplit[0]
                readLength = len(row[6])
                totalCt[gene_id] = 1
                if readLengthRestriction:
                    if readLength in readLengthRestriction:
                        passed = True
                    else:
                        passed = False
                else:
                    passed = True
                if nmdReadLength:
                    if gene_id in nmdReadLength:
                        passtwo = True
                    else:
                        passtwo = False
                else:
                    passtwo = True
                #print(passtwo)
                for txt in assocTxts:
                    pipe = '|' #New joshSAMs have all isoforms listed. Chopping off everything except first listed isoform.
                    newcurr = txt.split(pipe,1)[0]
                    curr = newcurr.split(':')
                    #print(assocTxts)
                    relStart = int(curr[1])
                    relStop = int(curr[2])
                    #SorAS = curr[3]
                    if SorAS == 'S' and passed and passtwo:
                        if relStart % 3. == 0:
                            inFrame += 1
                        else:
                            outOfFrame += 1
                            
                    if passed and passtwo:
                        # EDIT: Mar 30, 2015 Josh changed from plotting over entire read length to just 5' most end
                        """
                        for ii in range(readLength):
                            metaStart[SorAS][relStart + ii].append(1. / normFactors[gene_id] / numMapped / readLength)
                            metaStop[SorAS][relStop + ii].append(1. / normFactors[gene_id] / numMapped / readLength)
                        """
                        metaStart[SorAS][relStart].append(1./normFactors[gene_id]/numMapped)
                        metaStop[SorAS][relStop].append(1./normFactors[gene_id]/numMapped)
                        #print(metaStop)
                        
    print('For metaORF, averaged over whole read')
    print(inFrame, outOfFrame)
    print(inFrame / (outOfFrame + inFrame))
    metaAverage(metaStart, len(totalCt))
    metaAverage(metaStop, len(totalCt))
    
    return metaStart, metaStop


def processMeta(dict1, flip=0):
    """Will return dict1 as an ordered list of tuples. If flip will flip the sign"""
    a = dict1.keys()
    sorted(a) #Matt, was a.sort() which is python 2.7 syntax
    if not flip:
        return [(key, dict1[key]) for key in a]
    else:
        return [(key, -dict1[key]) for key in a]


def mkStartStopPlot(names, metaData, nmdMetaData, outPrefix):
    """metaData={name:({'S':metaStartSense,'AS':metaStartAS},{'S':metaStopSense,'AS':metaStopAS})}.
    Will plot sense above x=0, antisense below x=0, Start plot on left and Stop plot on right"""
    #print(nmdMetaData)
    #print(names)
    start = pyx.graph.graphxy(width=8, height=8,
                              x=pyx.graph.axis.linear(min=-100, max=100,
                                                      title='Position Relative to Start Codon'),
                              y=pyx.graph.axis.linear(title='Normalized Read Density'))
    stop = pyx.graph.graphxy(width=8, height=8, xpos=start.width * 1.1,
                             key=pyx.graph.key.key(pos='tr'),
                             x=pyx.graph.axis.linear(min=-100, max=100,
                                                     title='Position Relative to Stop Codon'),
                             y=pyx.graph.axis.linkedaxis(start.axes["y"]))
    
    for ii in range(len(names)):
        name = names[ii]
        metaStartSense = processMeta(metaData[name][0]['S'])
        metaStartSense.sort() #plotting error if not sorted
        start.plot(pyx.graph.data.points(metaStartSense, x=1, y=2),
                   [pyx.graph.style.line([common.colors(ii)])])
        metaStartAS = processMeta(metaData[name][0]['AS'], flip=1)
        metaStartAS.sort()
        start.plot(pyx.graph.data.points(metaStartAS, x=1, y=2),
                   [pyx.graph.style.line([common.colors(ii)])])
        
        metaStopSense = processMeta(metaData[name][1]['S'])
        metaStopSense.sort()
        #print(metaStopSense)
        stop.plot(pyx.graph.data.points(metaStopSense, x=1, y=2, title=name),
                  [pyx.graph.style.line([common.colors(ii)])])
        metaStopAS = processMeta(metaData[name][1]['AS'], flip=1)
        metaStopAS.sort()
        stop.plot(pyx.graph.data.points(metaStopAS, x=1, y=2, title=name),
                  [pyx.graph.style.line([common.colors(ii)])])
        
        nmdMetaStopSense = processMeta(nmdMetaData[name][1]['S'])
        nmdMetaStopSense.sort()
        print(nmdMetaStopSense)
        stop.plot(pyx.graph.data.points(nmdMetaStopSense, x=1, y=2, title=name),
                  [pyx.graph.style.line([common.colors(ii)])])
        
    c = pyx.canvas.canvas()
    c.insert(start)
    c.insert(stop)
    c.writePDFfile(outPrefix)

def parseNMDgeneList(nmdList):
    GeneList = []
    with open(nmdList, 'r') as f:
        freader = csv.reader(f)
        for row in freader:
            if not row in GeneList:
                GeneList.append(row)
    flattened = [val for sublist in GeneList for val in sublist]
    WBGene = dict.fromkeys(flattened, 1)
    return WBGene
def main(args):
    print('WARNING: Will break with multiple inFiles')
    annots, outPrefix = args[:2]
    inFiles = args[2::2]
    names = args[3::2]
    nmdList = args[4]
    nmdReadLength = parseNMDgeneList(nmdList)  
    #print(nmdReadLength)
    # add a variable for read length restriction
    exactLength = None
    # exactLength=[28,29,30]
    # exactLength=[15,16,17,18]
    print(f'Length restriction on reads: {exactLength}')
    
    txtLengths, cdsLengths, txtGroups = getLengths(annots)
    
    metaData = {}
    nmdMetaData = {}
    for i in range(len(names)):
        name = names[i]
        print(name)
        inFile = inFiles[i]
        metaData[name] = getMetaStartStop(inFile, txtLengths, txtGroups, readLengthRestriction=exactLength)
        nmdMetaData[name] = getMetaStartStop(inFile, txtLengths, txtGroups, nmdReadLength, readLengthRestriction=exactLength)
        #print(nmdMetaData)
    '''for i in range(len(names)):
        name = names[i]
        inFile = inFiles[i]
        nmdMetaData[name] = getMetaStartStop(inFile, txtLengths, txtGroups, nmdRestriction, readLengthRestriction=exactLength) 
    '''
    mkStartStopPlot(names, metaData, nmdMetaData, outPrefix)
    

if __name__ == '__main__':
    Tee()
    main(sys.argv[1:])
