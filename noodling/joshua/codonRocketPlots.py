"""
Joshua Arribere, April 2, 2020

Script to look at the occurrence of a codon at different in-frame positions
    of a length of read.

Input: reads.joshSAM - joshSAM-formatted file
    readLengthLower - shortest length of read to consider
    readLengthUpper - longest length of read to consider

Output: tiled plots where each row will be a read length. Each column
    within a row will be a different offset from the 5'end of the read.
    Each plot will have the observed in-frame frequency of the codon across
    the entire transcriptome on the y-axis, and the x-axis will be the
    counts of that codon in that position.

run as python codonRocketPlots.py inFile.joshSAM 15 34 outPrefix
"""
import sys, common, pickle
import pandas as pd
from logJosh import Tee

def passed(txtList):
    SorAS=list(set([entry.split(':')[-1] for entry in txtList]))
    relStart=list(set([entry.split(':')[1] for entry in txtList]))
    relStop=list(set([entry.split(':')[2] for entry in txtList]))
    #still not getting how python3 does iterators v not...
    relStart=[entry for entry in map(int,relStart)]
    relStop=[entry for entry in map(int,relStop)]
    if min(relStart)>0 and max(relStop)<=-12:
        #then it's in the ORF for all txts
        frames=list(set([entry%3 for entry in relStart]))
        if len(frames)==1:
            return frames[0]
    return 'na'

def getCodonCounts(inFile,X,Y):
    """
    Will go through the reads in the joshSAM file. Will:
    (1) Identify reads that are in ORFs for all txts that overlap that
        position.
    (2) For reads passing #1, will identify the in-frame codons at each
        position from the 5'end. Also determine readLength.
    (3) Will store the codon counts for each read length and each position
        in a DataFrame.
    (4) Will return that DataFrame.
    (5) Will also return a DataFrame of codon:count for all codons,
        irrespective of read length or position (but will be in-frame)
    Will do this for read length in [X,Y] inclusive.
    """
    #make the codons
    nucs=['A','T','G','C']
    codons=[nuc1+nuc2+nuc3 for nuc1 in nucs \
                        for nuc2 in nucs \
                        for nuc3 in nucs]
    #initialize the DataFrame
    #first make the indexes
    indexes=[]
    for ii in range(X,Y+1):
        for jj in range(ii-3):
            for codon in codons:
                indexes.append((ii,jj,codon))
    index=pd.MultiIndex.from_tuples(indexes,
        names=['readLength','offset','codon'])
    #now make the DataFrame using the multiindexes
    df=pd.DataFrame(data={'counts':0},index=index)
    #make a txtome codon count DataFrame
    codonTotalsDF=pd.DataFrame(data={'counts':0},
                                    index=codons)
    #############################################################
    #ok should now have a blank dataframe with multiindex column.
    #the indexes are in order readLength, offset, codon.
    #it's now safe to loop through the reads and get counts.
    cntr=0
    with open(inFile,'r') as f:
        for line in f:
            cntr+=1
            line=line.strip().split('\t')
            readSeq=line[3]
            txts=line[6:]
            #
            frame=passed(txts)
            #frame is 0, 1, 2, or 'na'. If it's 'na', then that
            #means that the frame was different for txts, it
            #was out of at least one ORF, or that it was AS
            if frame!='na':
                readLength=len(readSeq)
                for ii in range(3-frame,readLength-3,3):
                    codon=readSeq[ii:ii+3]
                    if readLength in range(X,Y+1) and codon in codons:
                        df.loc[readLength,ii,codon]['counts']+=1
                        codonTotalsDF.loc[codon]['counts']+=1
            if cntr%10000==0:
                print(f'Just finished line {cntr}.')
    
    return df, codonTotalsDF

def main(args):
    #parse the input
    inFile,readLengthLower,readLengthUpper,outPrefix=args[0:]
    #convert strings to ints
    readLengthLower=int(readLengthLower)
    readLengthUpper=int(readLengthUpper)
    #now parse the joshSAM file
    df,codonTotalsDF=getCodonCounts(inFile,readLengthLower,readLengthUpper)
    #pickle the output
    with open(outPrefix+'.p','wb') as f:
        pickle.dump(df,f)
    with open(outPrefix+'.control.p','wb') as f:
        pickle.dump(codonTotalsDF,f)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
