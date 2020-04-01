"""
Joshua Arribere Sept 10, 2013
Updated to python 3: Mar 23, 2020

Script to assign reads to genes. Will first assign sense, then antisense

Input: annots.allChrs.txt - output of prepareReadAssignmentFile2.py
    reads.bowtie - as output from bowtie
    BOTH OF THESE FILES ARE ASSUMED TO BE 1-INDEXED

Output: outPrefix.joshSAM - Output will be just like input, but will include a
    last column with the gene at that position, if applicable.

run as python assignReadsToGenes.py annots.gtf reads.bowtie outPrefix
"""
import sys, os, collections, csv, common, re, time, pickle, copy, linecache
import pandas
from logJosh import Tee
csv.field_size_limit(sys.maxsize)

def parseAnnotationDataFrame(annotFile):
    """Will take as input a file of the format:
    chrName_Number|txtInfo\ttxtInfo...\tGene\n
    And return a pandas DataFrame where the indexes are the first
    column (chrName_Number) and the only column is everything after
    the pipe (|)"""
    df=pandas.read_csv(annotFile,index_col=0,delimiter='|',
        names=['genes'])
    return df

def recoverMappedPortion(Cigar,Read):
    """Given a Cigar string and a Read, will return the sequence of the read that mapped to the genome."""
    #Edit Oct 10, 2013 to include skipped portions of reference sequence (introns)
    
    #first process the CIGAR string
    cigarSplit=re.findall('(\d+|[a-zA-Z]+)', Cigar)
    cigarSplit=[[int(cigarSplit[ii]),cigarSplit[ii+1]] for ii in range(0,len(cigarSplit),2)]
    
    #Then use that information to parse out nts of the read sequence
    mappedRead=''
    ii=0
    N=0
    for entry in cigarSplit:
        if entry[1] in ['M','I']:#then it's either aligned to the genomic sequence or has an insert relative to it
            mappedRead+=Read[ii:ii+entry[0]]
            ii+=entry[0]
        elif entry[1]=='S':
            ii+=entry[0]
        elif entry[1]=='N':
            N+=entry[0]
            #N is used for "skipped region from the reference". I keep track of Ns and return them for calculation of position on the - strand
    
    return mappedRead,N

def loadReadPositions(reads):
    """Will load up all read positions as a dictionary of {chr:{position:{}}}"""
    t=time.time()
    a=collections.defaultdict(dict)
    with open(reads,'rU') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            #print row
            if row[0][0]!='@':
                position=int(row[3])
                Chr=row[2]
                if int(row[1])&16!=0:#check strand
                    readSeq,N=recoverMappedPortion(row[5],row[9])
                    position+=len(readSeq)+N-1
                if row[11].split(':')[-1]=='1':#restriction for uniquely mapping
                    a[Chr][position]={}
                    if 'NH' not in row[11]:
                        print('ERROR: Restriction for uniquely mapping failed due to unexpected formatting.', sys.exit())
    print('Restricted to single mapping reads.')
    print(time.time()-t)
    return a

def parseAnnots(annots,readPositions):
    """Will ID which read positions overlap with a protein_coding gene's exon, regardless of strand"""
    t=time.time()
    a={}
    with open(annots,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            featureType=row[2]
            if featureType=='exon':# and row[1]=='protein_coding':
                Chr=row[0]
                strand=row[6]
                gene_id=row[8].split('"')[1]
                transcript_id = row[8].split('transcript_id')[1].split('"')[1]
                for ii in range(int(row[3]),int(row[4])):
                    if ii in readPositions[Chr]:
                        continue  # This is skipping the following code, intentional?
                        if gene_id not in readPositions[Chr][ii]:
                            readPositions[Chr][ii][gene_id]={'strand':strand,'transcript_id':[]}
                        readPositions[Chr][ii][gene_id]['transcript_id'].append(transcript_id)
            
            if featureType in ['exon','CDS']:# and row[1]=='protein_coding':
                transcript_id = row[8].split('transcript_id')[1].split('"')[1]
                if transcript_id not in a:
                    a[transcript_id]={'strand':row[6],'exon':[],'CDS':[]}
                a[transcript_id][featureType].append([int(row[3]),int(row[4])])
    
    b=dict((txt,a[txt]) for txt in a if a[txt]['CDS']!=[])
    #print len(b), ' number of annotated CDSs'
    print('hello')
    print(time.time() -t)
    print('goodbye')
    return b,readPositions

def assignReads(reads,annotDF,outPrefix):
    """Given annotDF, a pandas DataFrame with indexes as chr_position and
    first (and only column) as the txt information,
    and a reads file (from STAR, sam format), will make a joshSAM file,
    which contains all the info in the samFile along with the txt information
    """
    
    #print len(annots), ' number of annotated txts.'
    #print sum([len(readPositions[Chr]) for Chr in readPositions]), ' number of reads'
    readCt=0.
    unassignedCt=0
    with open(reads,'r') as f:
        with open(outPrefix+'.joshSAM','w') as g:
            freader=csv.reader(f,delimiter='\t')
            gwriter=csv.writer(g,delimiter='\t')
            for row in freader:
                if row[0][0]!='@':#then it's presumably a read line
                    ########parse the read information
                    position=int(row[3])
                    Chr=row[2]
                    
                    alignScore=int(row[13].split(':')[-1])
                    misMatches=int(row[14].split(':')[-1])
                    
                    readSeq,N=recoverMappedPortion(row[5],row[9])
                    alignLength=len(readSeq)
                    
                    if int(row[1])&16!=0:#check strand
                        position+=alignLength+N-1
                        readStrand='-'
                        readSeq=common.revCompl(readSeq)
                    else:
                        readStrand='+'
                    
                    readInfo=[Chr,position,readStrand,readSeq[:alignLength],alignLength]
                    ######done parsing readInformation
                    
                    #figure out if the position is present in the annotations
                    chrPosition=f'{Chr}_{position}'
                    try:
                        txtInfo=annotDF.loc[chrPosition]['genes']
                        annotFileLine=txtInfo.strip().split('\t')
                        
                        if len(annotFileLine)>1:
                            if row[11].split(':')[-1]=='1':#restriction for uniquely mapping#added July 27, 2014
                                geneInfo=copy.copy(annotFileLine)
                                ######using pop below makes the above copy necessary
                                gene,geneStrand=geneInfo.pop().split(':')
                                if geneStrand==readStrand:
                                    SorAS='S'
                                else:
                                    SorAS='AS'
                                
                                readInfo.append(gene)
                                for txt in geneInfo:
                                    readInfo.append(':'.join([txt,SorAS]))
                                gwriter.writerow(readInfo)  # I really think that this may be the time limiting step,
                                # having the script write to the csv during each loop is going to slow the process down.
                                # It may be more functional to have the function hold the reads in memory,
                                # and than batch write all at once?
                                readCt+=1
                        else:
                            #hwriter.writerow(readInfo+readPositions[Chr][position].keys())
                            if row[12].split(':')[-1]=='1':#only count a multiply-mapping read the first time it appears
                                unassignedCt+=1#0.066million in UCSC v 0.334million in ENS
                    except KeyError:
                        #then the chr_position is not in the annotations,
                        #meaning no txt is annotated as overlapping with it
                        pass
    print('We need to fix the next two lines to make more \
        meaningful.')
    print('%s reads (%s of reads) were uniquely assigned to a gene with file %s.'%(readCt,readCt/(unassignedCt+readCt),reads))
    print('%s reads (%s of reads) were unassigned.'%(unassignedCt,unassignedCt/(unassignedCt+readCt)))

def main(args):
    t=time.time()
    annotFile,reads,outPrefix=args[0:]
    
    #first parse in the output of prepareReadAssignmentFile
    annotDF=parseAnnotationDataFrame(annotFile)
    
    #with open(annots,'r') as f:
    #    annots=cPickle.load(f)
    
    assignReads(reads,annotDF,outPrefix)
    print(time.time()-t)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
