"""
Joshua Arribere Sept 10, 2013

Script to assign reads to genes. Will first assign sense, then antisense

Input: annotGTFFile - GTF file used for input for prepareReadAssignmentFile2.py
    reads.bowtie - as output from bowtie
    BOTH OF THESE FILES ARE ASSUMED TO BE 1-INDEXED

Output: outPrefix.joshSAM - Output will be just like input, but will include a last column with the gene at that position, if applicable.
    #outPrefix.geneCount - These output files will have "gene   ct"

run as python assignReadsToGenes.py annots.gtf reads.bowtie outPrefix
EDIT: Sept 27, 2013 JOSH revised to put sense/antisense in one file, also omit .geneCount file
EDIT: March 6, 2014 JOSH revised to input output of prepareReadAssignmentFile.py for read
    assignment instead of recalculating read positions from annotations
EDIT: Sept 19, 2014 JOSH revised to be take advantage of faster prepareReadAssignmentFile2.py output
"""
import sys, os, collections, csv, common, re, time, cPickle, copy, linecache
from logJosh import Tee
csv.field_size_limit(sys.maxsize)

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
                        print 'ERROR: Restriction for uniquely mapping failed due to unexpected formatting.', sys.exit()
    print 'Restricted to single mapping reads.'
    print time.time()-t
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
                        continue
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
    print 'hello'
    print time.time() -t
    print 'goodbye'
    return b,readPositions

def getDist(exon_list,position):
    """Will return distance from the first position of the first exon to position, in mRNA space"""
    dist=0
    for exon in exon_list:
        if position in range(exon[0],exon[1]+1):
            return dist+position-exon[0]
        else:
            dist+=exon[1]-exon[0]+1#Edit July 28, 2014 added +1.
    print exon_list, position, dist
    print 'Error: position not found in txt', sys.exit()

def getTxtRelPositions(transcript_id,txt_annot,readStrand,readPosition):
    """Given txt_annot={strand,exon,CDS},readStrand,readPosiiton, will return a colon-separated string of
    transcript_id, readPosition_rel_to_CDSstart, readPosition_rel_to_CDSEnd, S/AS where S/AS indicates whether
    the read is on the same strand or not"""
    #figure out sense/antisense. Easy
    txtStrand=txt_annot['strand']
    if txtStrand==readStrand:
        SorAS='S'
    else:
        SorAS='AS'
    
    #Great, now figure out position relative to CDSStart
    #To do this, I will calculate the distance from the txtStart to readPosition, and txtStart to CDSStart. I'll subtract the two.
    cdsStart=min([entry[0] for entry in txt_annot['CDS']])
    exons=txt_annot['exon']
    exonStarts=[exon[0] for exon in exons]
    exonEnds=[exon[1] for exon in exons]
    exonStarts.sort(),exonEnds.sort()
    exons=zip(exonStarts,exonEnds)
    #Edit: Apparently the exons are not necessarily orderd in the gtf file.
    
    txtStart_cdsStart_dist=getDist(exons,cdsStart)
    txtStart_readPosition_dist=getDist(exons,readPosition)
    readPosition_rel_to_CDSstart=txtStart_readPosition_dist-txtStart_cdsStart_dist
    
    #now do the same thing with the cdsEnd
    cdsEnd=max([entry[1] for entry in txt_annot['CDS']])
    txtStart_cdsEnd_dist=getDist(exons,cdsEnd)
    #already determined txtStart_readPosition_dist
    readPosition_rel_to_CDSend=txtStart_readPosition_dist-txtStart_cdsEnd_dist
    
    #stranded issues. Although ensembl defines start_codon and stop_codon as those exact locations on the - and + strand, here I find the start/stop by taking min/max of CDS exon boundaries, so I need to flip it
    if txtStrand=='+':
        return ':'.join([transcript_id,str(readPosition_rel_to_CDSstart),str(readPosition_rel_to_CDSend),SorAS])
    else:
        return ':'.join([transcript_id,str(-readPosition_rel_to_CDSend),str(-readPosition_rel_to_CDSstart),SorAS])

def assignReads(reads,annots,outPrefix):
    """Given annots={txt_ID:{strand,CDS,exon}}, readPositions={chr:position:gene_id:{strand,txt_ids}},
    and a reads file (from STAR, sam format), will make a file with
    chr position strand read_seq Aligned_length gene_id txt_id:relToCDSStart:relToCDSEnd:S/AS
    where the last column has potentially multiple comma-separated entries"""
    
    #print len(annots), ' number of annotated txts.'
    #print sum([len(readPositions[Chr]) for Chr in readPositions]), ' number of reads'
    readCt=0.
    unassignedCt=0
    with open(reads,'r') as f:
        with open(outPrefix+'.joshSAM','w') as g:
            #with open(outPrefix+'.unassigned','w') as h:
            freader=csv.reader(f,delimiter='\t')
            gwriter=csv.writer(g,delimiter='\t')
            #hwriter=csv.writer(h,delimiter='\t')
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
                    
                    if os.path.isfile('.'.join(annots.split('.')[:-1])+'.%s.txt'%Chr):
                        annotFileLine=linecache.getline('.'.join(annots.split('.')[:-1])+'.%s.txt'%Chr,position)
                        annotFileLine=annotFileLine.strip().split('\t')
                        
                        if len(annotFileLine)>1:
                            #print row
                            #print row[11]
                            #print row[11].split(':')[-1]
                            if row[11].split(':')[-1]=='1':#restriction for uniquely mapping#added July 27, 2014
                                geneInfo=copy.copy(annotFileLine[1:])
                                ######using pop below makes the above copy necessary
                                gene,geneStrand=geneInfo.pop().split(':')
                                if geneStrand==readStrand:
                                    SorAS='S'
                                else:
                                    SorAS='AS'
                                
                                readInfo.append(gene)
                                for txt in geneInfo:
                                    readInfo.append(':'.join([txt,SorAS]))
                                #if 'T25F10.5' in readInfo:#left over from fixing a bug after doing an analysis for Drew Nager
                                #    print readInfo
                                gwriter.writerow(readInfo)
                                readCt+=1
                        else:
                            #hwriter.writerow(readInfo+readPositions[Chr][position].keys())
                            if row[12].split(':')[-1]=='1':#only count a multiply-mapping read the first time it appears
                                unassignedCt+=1#0.066million in UCSC v 0.334million in ENS
    print '%s reads (%s of reads) were uniquely assigned to a gene with file %s.'%(readCt,readCt/(unassignedCt+readCt),reads)
    print '%s reads (%s of reads) were unassigned.'%(unassignedCt,unassignedCt/(unassignedCt+readCt))

def main(args):
    t=time.time()
    annotGTFFile,reads,outPrefix=args[0:]
    
    #with open(annots,'r') as f:
    #    annots=cPickle.load(f)
    
    assignReads(reads,annotGTFFile,outPrefix)
    print time.time()-t


if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
