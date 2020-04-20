"""
Converted to python 3 on 04/17/2020

Joshua Arribere Sept 30, 2013

Script to parse a .joshSAM-formatted file and return a python dict of {txt:{position:ct}}

Input: inFile.joshSAM (as output from assignReadsToGenes.py)
    annots.gtf - Used to get txtLength

Output: python dict

run as python parse_joshSAM.py
"""
import sys, common, csv, collections, numpy
from logJosh import Tee

def getLength(annots):
    """Will return dict of total txt length as well as {gene:[txt_list]}"""
    a={}
    b=collections.defaultdict(dict)
    with open(annots,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            featureType=row[2]
            if featureType in ['exon','CDS']:# and row[1]=='protein_coding':
                transcript_id = row[8].split('transcript_id')[1].split('"')[1]
                if transcript_id not in a:
                    a[transcript_id]={'strand':row[6],'exon':[],'CDS':[]}
                a[transcript_id][featureType].append([int(row[3]),int(row[4])])
                
                gene_id=row[8].split('"')[1]
                b[gene_id][transcript_id]=1
    
    exon={}
    for txt in a:
        exon[txt]=sum([entry[1]-entry[0] for entry in a[txt]['exon']])
    return exon, dict(b)

def getLigationCorrection(readSeq,side,ligCorrect):
    """Will look on the side(false=5',true=3') of readSeq and return the value of ligCorrect corresponding to
    that sequence"""
    if side:#then it's looking at the 3' end
        seq=readSeq[-2:]
    else:#then it's looking at the 5'end
        seq=readSeq[:2]
    
    if seq in ligCorrect:
        return ligCorrect[seq]
    else:
        return 1

def determineIfInCDS(cdsBounds,txts):
    """Given cdsBounds=[x,y], will loop through txts=[txt:relStart:relStop:SorAS] and return True if all the
    relStarts are >x and all the relStops are <y"""
    success=True
    for txt in txts:
        relStart,relStop=list(map(int,txt.split(':')[1:3]))
        if relStart>=cdsBounds[0] and relStop<=cdsBounds[1]:
            continue
        else:
            return False
    return True

def getGeneReadCts(inFile,side=False,ligCorrect=False,antisenseOnly=False,senseOnly=False,cdsBounds=False):
    """Will sum up the total number of reads mapping to each gene_id (row[5]). side is assumed 5', and will
    only be used if ligCorrect passed"""
    a=collections.defaultdict(float)
    with open(inFile,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row)>5:# and row[3].startswith('G') and row[4]=='26':
                ct=1.
                #This part will correct for ligation bias based on 3'/5'end dint content
                if ligCorrect:
                    ct*=getLigationCorrection(row[3],side,ligCorrect)
                
                if cdsBounds:
                    withinCDS=determineIfInCDS(cdsBounds,row[6:])
                else:
                    withinCDS=True#assumed true
                
                if withinCDS:
                    if not antisenseOnly and not senseOnly:
                        a[row[5]]+=ct
                    else:
                        if antisenseOnly:
                            if row[-1].endswith('AS'):
                                a[row[5]]+=ct
                        elif senseOnly:
                            if row[-1].endswith('S') and (not row[-1].endswith('AS')):
                                a[row[5]]+=ct
    return dict(a)

def getNormFactor(geneReadCts,txtLengths,txtGroups):
    """For each gene in geneReadCts, will look up its associated txts in txtGroups,
    and then their lengths in txtLengths. Will return the the read count divided by
    the average txt length"""
    a=dict((gene_id,geneReadCts[gene_id]/numpy.average([txtLengths[txt] for txt in txtGroups[gene_id]])) for gene_id in geneReadCts)
    return a

def getFreq(dict1):
    """Given a dict of {key:value}, will switch to {key:freq}, where freq=value/sum(dict.values())"""
    total=float(sum(dict1.values()))
    for key in dict1:
        dict1[key]/=total

def determineLigBias(inFile,threePrime):
    """Given a .joshSAM file, will determine the dint frequency within reads, and the dint frequency at the
    threePrime end (if threePrime is set to True. Else, fivePrime end"""
    nucs=['A','T','C','G']
    endDict=dict((nuc1+nuc2,0.) for nuc1 in nucs for nuc2 in nucs)
    seqCtDict=dict((nuc1+nuc2,0.) for nuc1 in nucs for nuc2 in nucs)
    with open(inFile,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            readSeq=row[3]
            for seq in seqCtDict:
                seqCtDict[seq]+=readSeq.count(seq)
            
            if threePrime:
                seq=readSeq[-2:]
                if seq in endDict:
                    endDict[seq]+=1.
            else:
                seq=readSeq[:2]
                if seq in endDict:
                    endDict[readSeq[:2]]+=1.
    
    getFreq(endDict)
    getFreq(seqCtDict)
    
    return dict((seq,seqCtDict[seq]/endDict[seq]) for seq in endDict)

def norMalize(inFile,txtLengths,txtGroups,threePrime,lengthRestriction,firstNt,withLigCorrection=False):
    """Will get distribution of reads on a txt-by-txt basis relative to start codon. To do this it will first:
    (1) Add up the total number of reads mapping to each gene
    (2) Will find the normalization factor s.t. the mean read count across a gene per kb is 1
        b/c txt size is potentially variable across a group, I take the average of all the txts within a gene. Can't think of a better thing to do right now...
    (3) Will apply that normalization factor to each read
    (4) Concomitant with #3, will add read count to a position-specific read count dict
    I distinguish between genes (keys to txtGroups) and txts (values in txtGroups)
    """
    #First, determine correction for ligation, if required.
    if withLigCorrection:
        ligCorrection=determineLigBias(inFile,threePrime)
        print(ligCorrection)
    else:
        ligCorrection=None
    
    #Get read counts/gene
    geneReadCts=getGeneReadCts(inFile,side=threePrime,ligCorrect=ligCorrection)
    
    #Get normalization factor
    normFactors=getNormFactor(geneReadCts,txtLengths,txtGroups)
    
    #Apply the normalization factor and record read positions
    normStart={'S':collections.defaultdict(lambda:collections.defaultdict(float)),
               'AS':collections.defaultdict(lambda:collections.defaultdict(float))}
    with open(inFile,'r') as f:
        freader=csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row)>5:
                assocTxts=row[6:]
                numMapped=len(assocTxts)
                gene_id=row[5]
                
                readSeq=row[3]
                readLength=len(row[3])
                
                #Now I'll filter by length and firstNt
                if lengthRestriction:
                    if readLength==lengthRestriction:
                        passed=1
                    else:
                        passed=0
                else:
                    passed=1
                
                if firstNt:
                    if readSeq.startswith(firstNt):
                        ntMatch=1
                    else:
                        ntMatch=0
                else:
                    ntMatch=1
                
                if passed and ntMatch:
                    for txt in assocTxts:
                        curr=txt.split(':')
                        relStart=int(curr[1])
                        SorAS=curr[3]
                        if threePrime:
                            if SorAS=='S':
                                relStart+=readLength
                            elif SorAS=='AS':
                                relStart-=readLength
                            else:
                                print('Strand not recognized '+SorAS, sys.exit())
                        
                        ct=1.
                        if ligCorrection:
                            ct*=getLigationCorrection(readSeq,threePrime,ligCorrection)
                        
                        normStart[SorAS][curr[0]][relStart]+=ct/normFactors[gene_id]/numMapped
    
    return normStart

def main(args,threePrime=None,lengthRestriction=None,firstNt=None):
    """Will return a dict of {txt:{position:ct}} and {gene:ct}. The second dict has the total counts for a group of
    txts from a gene. The first dict is normalized s.t. the average read distribution across a gene will be 1/nt"""
    inFile,annots=args[0:]
    
    txtLengths,txtGroups=getLength(annots)
    
    metaStart=norMalize(inFile,txtLengths,txtGroups,threePrime,lengthRestriction,firstNt)
    
    return metaStart

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
