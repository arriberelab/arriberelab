"""$
  2 Joshua Arribere Feb 24, 2014$
  3 $
  4 Script to collapse reads containing the exact same sequence. Will pick the read information$
  5     from the first occurence of that read for the outFile.$
  6 $
  7 Input: reads.fastq - fastq-formatted reads$
  8 $
  9 Output: reads.collapsed.fastq - same as input, with redundant reads removed$
 10 $
 11 run as python readCollapser.py inFile.fastq outPrefix$
 12 """$
 13 import sys, common, collections, numpy$
 14 from logJosh import Tee$
 15 $
 16 def main(args):$
 17     inFile,outPrefix=args[0:]$
 18     print 'Trimming of 6 nucleotides from 3\'end and 4 nts from 5\'end'$
 19     #print 'NOT TRIMMING 6 Nts FROM THE 3\'end!!!'$
 20     #print 'trimming 6 Ns from the 3\'end'$
 21     #print 'NOT COLLAPSING READS!'$
 22     reads=collections.defaultdict(lambda:collections.defaultdict(int))$
 23     with open(inFile,'r') as f:$
 24         with open(outPrefix,'w') as g:$
 25             currRead=[]$
 26             for line in f:$
 27                 line=line.strip()$
 28                 if len(currRead)==4:$
 29                     if reads[currRead[1][:-6]][currRead[1][-6:]]==0:$
 30                         #the next line has been moved ~6 lines down$
 31                         #reads[currRead[1]]+=1#also need to change this line for 3 or 6EDIT: DO NOT NEED TO CHANGE: This keeps track of read+barcode, and only writes if read appears with new barcode$
 32                         #g.write('%s\n%s\n%s\n%s\n'%(currRead[0],currRead[1],currRead[2],currRead[3]))$
 33                         g.write('%s\n%s\n%s\n%s\n'%(currRead[0],currRead[1][4:][:-6],currRead[2],currRead[3][4:][:-6]))$
 34                         #g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1][:-6],currRead[2],currRead[3][:-6]))$
 35                         if len(currRead[1])!=len(currRead[3]):$
 36                             print currRead, sys.exit()$
 37                         #must trim quality scores and reads simultaneously. Hence the double -6$
 38                     reads[currRead[1][:-6]][currRead[1][-6:]]+=1$
 39                     currRead=[]$
 40                     currRead.append(line)$
 41                 else:$
 42                     currRead.append(line)$
 43     #now look at complexity$
 44     cntr=[]$
 45     for read in reads:$
 46         if 2<=sum(reads[read].values())<=5:$
 47             #cntr.append(float(len(reads[read]))/sum(reads[read].values()))$
 48             cntr.append(sum([1 for barcode in reads[read] if reads[read][barcode]==1])/float(len(reads[read])))$
 49     print '%s average +/-SD fraction of reads with 2-5 counts with unique barcodes: % s+/- %s'%(inFile,numpy.average(cntr), numpy.std(cntr))$
 50 $
 51 if __name__=='__main__':$
 52     """run as python readCollapser.py inFile.fastq outPrefix"""$
 53     Tee()$
 54     main(sys.argv[1:])
