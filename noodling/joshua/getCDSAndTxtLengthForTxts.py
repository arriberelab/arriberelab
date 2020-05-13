"""
Joshua Arribere, April 8, 2020

Script to pull transcript names and annotations from a gtf file. Will
    also determine CDS and txt length

Input: annots.gtf - gtf-formatted annotations.

Output: File of the format:
    txtName\tgeneName\tcdsLength\ttxtLength

run as python getCDSAndTxtLengthForTxts.py annots.gtf outPrefix
"""
import sys, common
from logJosh import Tee
import pandas

def parseAnnotationFile(annotFile):
    """
    Will go through a gtf-formatted file and identify all genes and
    associated transcripts. Will return a pandas DataFrame of
    txtName, geneName,cdsLength,txtLength
    """
    aa={}
    with open(annotFile,'r') as f:
        for line in f:
            if not line.startswith('#'):
                if 'protein_coding' in line:
                    line=line.strip().split('\t')
                    featureType=line[2]
                    if featureType in ['exon','CDS']:
                        ##I'm going to add the version number to the
                        ##geneID and txtID as an added safety check
                        ##for later
                        geneID=line[8].split('gene_id "')[1].\
                                    split('"')[0]
                        #geneVersion=line[8].\
                        #            split('gene_version "')[1].\
                        #            split('"')[0]
                        #geneID=geneID+'.'+geneVersion
                        ##
                        txtID=line[8].split('transcript_id "')[1].\
                                    split('"')[0]
                        #txtVersion=line[8].\
                        #            split('transcript_version "')[1].\
                        #            split('"')[0]
                        #txtID=txtID+'.'+txtVersion
                        ##
                        if geneID not in aa:
                            aa[geneID]={}
                        if txtID not in aa[geneID]:
                            aa[geneID][txtID]={'CDS':0,'exon':0}
                        ##
                        length=int(line[4])-int(line[3])
                        ##this next line is needed b/c the bounds of
                        ##the features are done by the coding nts,
                        ##inclusive
                        length+=1
                        if length<=0:
                            print('Error, you\'re in non-positive \
                                genomic space! Here be dragons!')
                            print(line)
                            print(length,line[3],line[4])
                            sys.exit()
                        aa[geneID][txtID][featureType]+=length
    ##identify the monoTxt genes
    genes,txts,CDSlengths,txtLengths=[],[],[],[]
    for gene in aa:
        for txt in aa[gene]:
            genes.append(gene)
            txts.append(txt)
            CDSlengths.append(aa[gene][txt]['CDS'])
            txtLengths.append(aa[gene][txt]['exon'])
    ##make a DataFrame from that
    data={'gene':genes,
        'txt':txts,
        'CDSlength':CDSlengths,
        'txtLength':txtLengths}
    df=pandas.DataFrame.from_dict(data)
    ##pass the DataFrame back
    return df

def writeOutput(df,outPrefix):
    """
    Will write a dataframe to a file.
    """
    df.to_csv(outPrefix,sep='\t',index=False)

def main(args):
    ##pass the command line inputs
    annotFile,outPrefix=args[0:]
    ##parse the annotation file
    annotInfo=parseAnnotationFile(annotFile)
    ##write the output
    writeOutput(annotInfo,outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
