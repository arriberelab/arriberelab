"""
Joshua Arribere March 4, 2014

Script to extract lists of genes deemed to be differentially expressed from
   output of medianNormalizerWithDESeq.py

Input: inFile_S/AS.DESeqgeneCts_diffExpression - output from R of nbinomTest, written as a .csv file
   p_val - p-value cutoff. Will pull out everything less than or equal to this p-val (padj)
   FC - fold change cutoff. Will only pull out genes that change at least this much.
    If you don't want to deal with this option, just use 1 (one fold change).

run as python extractDifferentiallExpressedGenesFromDESeqTable.py inFile_S.FESeqgeneCts_diffExpression p_val FC outPrefix
"""
import sys
from logJosh import Tee

def parseInFile(inFile,FC):
   """Will parse a .csv file, skipping the first line, outputting second col value and last col value a key,value
   pairs in a dict
   EDIT: Will split into genes that are up/down"""
   aa,bb={},{}
   with open(inFile,'r') as f:
      f.readline()
      for line in f:
         line=line.strip().split(',')
         foldChange=float(line[5])
         if foldChange<1./FC:
            aa[line[1].strip('"')]=float(line[-2])
         elif foldChange>FC:
            bb[line[1].strip('"')]=float(line[-2])
         else:
            pass
   return aa,bb

def writeTable(pvalDict,pval,outFile):
   passedDict=dict((key,pvalDict[key]) for key in pvalDict if pvalDict[key]<=pval)
   
   with open(outFile,'w') as f:
      for gene in passedDict:
         f.write(gene+'\n')

def main(args):
   inFile,pval,FC,outPrefix=args[0:]
   
   FC=float(FC)
   pvalDictUp,pvalDictDown=parseInFile(inFile,FC)
   
   pval=float(pval)
   writeTable(pvalDictUp,pval,outPrefix+'Up.txt')
   writeTable(pvalDictDown,pval,outPrefix+'Down.txt')

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
