"""
Joshua Arribere, March 6, 2014

EDIT Nov 26 2019 - Made it so prepareReadAssignmentFile creates a .txt output for each chr one at a time. Doing all at once was maxxing out memory.

Read assignment with assignReadsToGenes.py is taking an interminable amount of time. Part
    of this is because for each read, I calculate whether a gene overlaps with it, whether
    one and only one gene overlaps with it, and then the position of that read in all
    transcripts for that gene. This is a lot of calculations, for each read.
    
    An alternative that I'm going to try out here, is to calculate this for every possible
    position in the worm genome. Then when I have my reads, loop through them and simply
    perform a look up in that data structure.
    
    This script will create that data structure. It will, for every position in a provided
    genome annotation, return all positions which overlap with exactly one gene. For each
    such position, it will provide the position relative to annotated start/stop codons
    for each transcript overlapping that position.

Input: annots.gtf - gtf-formatted annotations. Will also get output name from this

Output: format TBD

run as python prepareReadAssignmentFile.py annots.gtf
EDIT Sept 19 2014 - JOSH made it output a .txt file for each chr, which should speed up
    assignment b/c it doesn't have to load the entire thing into memory at once.
"""
import sys, common, time, collections, csv, subprocess, cPickle
from logJosh import Tee
import assignReadsToGenes3, multiprocessing

def parseAnnots(annots):#borrowed and adapted from parseAnnots
    """Will output readPositions, a dictionary of {chr:position:{}} where the innermost dictionary
    contains information about all genes/transcripts overlapping with that position."""
    
    ########the below lines are just to estimate time the script will take...
    t=time.time()
    ct=0.
    p=subprocess.Popen(['wc','-l',annots],shell=False,stdout=subprocess.PIPE)
    aa=p.stdout.read()
    aa=aa.split()
    total=int(aa[0])
    ########the above lines are just to estimate time the script will take...
    
    a={}
    readPositions=collections.defaultdict(lambda:collections.defaultdict(dict))
    with open(annots,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row)>1:
                ct+=1
                featureType=row[2]
                if featureType=='exon':# and row[1]=='protein_coding':
                    Chr=row[0]
                    strand=row[6]
                    gene_id=row[8].split('"')[1]
                    transcript_id = row[8].split('transcript_id')[1].split('"')[1]
                    for ii in range(int(row[3]),int(row[4])+1):#edit July 29, 2014 josh added +1 to include right exon boundary
                        if gene_id not in readPositions[Chr][ii]:
                            readPositions[Chr][ii][gene_id]={'strand':strand,'transcript_id':[]}
                        readPositions[Chr][ii][gene_id]['transcript_id'].append(transcript_id)
                
                if featureType in ['exon','CDS']:# and row[1]=='protein_coding':
                    transcript_id = row[8].split('transcript_id')[1].split('"')[1]
                    if transcript_id not in a:
                        a[transcript_id]={'strand':row[6],'exon':[],'CDS':[]}
                    a[transcript_id][featureType].append([int(row[3]),int(row[4])])
                if ct%25000==0:
                    print '%s fraction of the way through the input file...'%str(ct/total)
    
    b=dict((txt,a[txt]) for txt in a if a[txt]['CDS']!=[])
    #print len(b), ' number of annotated CDSs'
    print 'hello'
    print time.time() -t
    print 'goodbye'
    return b,readPositions

def getTxtRelPositions(transcript_id,txt_annot,readPosition):#modified so that it returns strand +/- instead of S/AS
    """Given txt_annot={strand,exon,CDS},readStrand,readPosiiton, will return a colon-separated string of
    transcript_id, readPosition_rel_to_CDSstart, readPosition_rel_to_CDSEnd, S/AS where S/AS indicates whether
    the read is on the same strand or not"""
    #Great, now figure out position relative to CDSStart
    #To do this, I will calculate the distance from the txtStart to readPosition, and txtStart to CDSStart. I'll subtract the two.
    cdsStart=min([entry[0] for entry in txt_annot['CDS']])
    exons=txt_annot['exon']
    exonStarts=[exon[0] for exon in exons]
    exonEnds=[exon[1] for exon in exons]
    exonStarts.sort(),exonEnds.sort()
    exons=zip(exonStarts,exonEnds)
    #Edit: Apparently the exons are not necessarily orderd in the gtf file.
    
    txtStart_cdsStart_dist=assignReadsToGenes3.getDist(exons,cdsStart)
    txtStart_readPosition_dist=assignReadsToGenes3.getDist(exons,readPosition)
    readPosition_rel_to_CDSstart=txtStart_readPosition_dist-txtStart_cdsStart_dist
    
    #now do the same thing with the cdsEnd
    cdsEnd=max([entry[1] for entry in txt_annot['CDS']])
    txtStart_cdsEnd_dist=assignReadsToGenes3.getDist(exons,cdsEnd)
    #already determined txtStart_readPosition_dist
    readPosition_rel_to_CDSend=txtStart_readPosition_dist-txtStart_cdsEnd_dist
    
    #stranded issues. Although ensembl defines start_codon and stop_codon as those exact locations on the - and + strand, here I find the start/stop by taking min/max of CDS exon boundaries, so I need to flip it
    if txt_annot['strand']=='+':
        return ':'.join([transcript_id,str(readPosition_rel_to_CDSstart),str(readPosition_rel_to_CDSend)])
    else:
        return ':'.join([transcript_id,str(-readPosition_rel_to_CDSend),str(-readPosition_rel_to_CDSstart)])

def writeOutput(someTuple):
    """
    Will write the output file. For every position from 1 to the largest key in dict1,
    will write a line. If the key is in dict1, will tab, then write the value.
    """
    Chr,dict1,annotFileName=someTuple[0],someTuple[1],someTuple[2]
    with open('.'.join(annotFileName.split('.')[:-1])+'.%s.txt'%Chr,'w') as f:
        for ii in range(1,max(dict1)+1):
            f.write('%s'%ii)
            if ii not in dict1:
                f.write('\n')
            else:
                for entry in dict1[ii]:
                    f.write('\t%s'%entry)
                f.write('\n')

def positionReadsOnTxts(annots,readPositions,annotFileName):
    """Given annots={txt_ID:{strand,CDS,exon}}, readPositions={chr:position:gene_id:{strand,txt_ids}},
    will make a file with
    chr position strand read_seq Aligned_length gene_id txt_id:relToCDSStart:relToCDSEnd:S/AS
    where the last column has potentially multiple comma-separated entries"""
    #print len(annots), ' number of annotated txts.'
    #print sum([len(readPositions[Chr]) for Chr in readPositions]), ' number of reads'
    
    outputPositions={}
    
    for Chr in readPositions:
        outputPositions[Chr]=collections.defaultdict(list)
        
        for position in readPositions[Chr]:
            if readPositions[Chr][position]=={}:#then the read never mapped within to an annotated feature, Should not exist in this version.
                pass#for now I'll ignore these reads
            elif len(readPositions[Chr][position])==1 and readPositions[Chr][position]!={}:#omiting reads that are ambiguously assignable to genes b/c genes overlap
                gene_id=readPositions[Chr][position].keys()[0]
                positionInfo=[]
                for transcript_id in readPositions[Chr][position][gene_id]['transcript_id']:
                    if transcript_id in annots:
                        if annots[transcript_id]['CDS']!=[]:
                            positionInfo+=[getTxtRelPositions(transcript_id,annots[transcript_id],position)]
                            txtStrand=annots[transcript_id]['strand']
                if len(positionInfo)==0:#then it's not assignable to a single txt--maybe b/c no txts are annotated for that gene, or b/c they all lack a CDS annotation
                    pass#don't care about it if there's no CDS annotated for the txts
                else:#then we've got a gene and at least one transcript
                    positionInfo.append(gene_id+':'+txtStrand)
                    outputPositions[Chr][position]=positionInfo
            elif len(readPositions[Chr][position])>1:#then it's multiply mapping
                pass
    
    with open('.'.join(annotFileName.split('.')[:-1])+'.processed.p','w') as f:
        cPickle.dump(outputPositions,f,protocol=2)
    
    #If you already have a processed.p file. comment out 137 to 157, and uncomment 160 to 162.

    #print 'Unpickling...'
    #outputPositions=common.unPickle('.'.join(annotFileName.split('.')[:-1])+'.processed.p')
    #print 'Done unpickling!'
    
    ##now output .txt file
    #p=multiprocessing.Pool(len(outputPositions)) #do not do multiprocessing because it will max out memory
    #p.map(writeOutput,[(Chr,outputPositions[Chr],annotFileName) for Chr in outputPositions])
    
    for Chr in outputPositions:
        print 'Working on %s...'%(Chr)
        writeOutput((Chr,outputPositions[Chr],annotFileName))
        print 'Finished %s.'%(Chr)

def main(args):
    annotFileName=args[0]
    
    annots,readPositions=parseAnnots(annotFileName)
    #annots='blag'
    #readPositions='abas'
    
    positionReadsOnTxts(annots,readPositions,annotFileName)



if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
