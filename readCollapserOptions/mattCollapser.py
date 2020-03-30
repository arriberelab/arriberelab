 """$
  2 Matt Modena Mar30, 2020$
  3 $
  4 Script to collapse reads containing the exact same sequence. Will pick the read information$
  5     from the first occurence of that read for the outFile.$
  6 $
  7 Input: reads.fastq - fastq-formatted reads$
  8 $
  9 Output: reads.collapsed.fastq - same as input, with redundant reads removed$
 10 $
 11 run as python readCollapser.py inFile.fastq 5'UMI 3'UMI outPrefix$
 12 """$
 13 import sys, common, collections, numpy$
 14 from logJosh import Tee$
 15 $
 16 def main(args):$
 17     inFile,x,y,outPrefix=args[0:]$
 18     print 'Trimming of y nucleotides from 3\'end and x nts from 5\'end'$
 19     #print 'NOT TRIMMING 6 Nts FROM THE 3\'end!!!'$
 20     #print 'trimming 6 Ns from the 3\'end'$
 21     #print 'NOT COLLAPSING READS!'$
 22     x=int(x)$
 23     y=int(y)$
 24     reads=collections.defaultdict(lambda:collections.defaultdict(int))$
 25     with open(inFile,'r') as f:$
 26         with open(outPrefix,'w') as g:$
 27             currRead=[]$
 28             for line in f:$
 29                 line=line.strip()$
 30                 if len(currRead)==4:$
 31                     if reads[currRead[1][:x]][currRead[1][-y:]]==0:$
 32                         #the next line has been moved ~6 lines down$
 33                         #reads[currRead[1]]+=1#also need to change this line for 3 or 6EDIT: DO NOT NEED TO CHANGE: This keeps track of read+barcode, and only writes if read appears with new barcode$
 34                         #g.write('%s\n%s\n%s\n%s\n'%(currRead[0],currRead[1],currRead[2],currRead[3]))$
 35                         g.write('%s\n%s\n%s\n%s\n'%(currRead[0],currRead[1][:x][-y:],currRead[2],currRead[3][:x][-y:]))$
 36                         #g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1][:-6],currRead[2],currRead[3][:-6]))$
 37                         if len(currRead[1])!=len(currRead[3]):$
 38                             print currRead, sys.exit()$
 39                         #must trim quality scores and reads simultaneously. Hence the double -6$
 40                     reads[currRead[1][:x]][currRead[1][-y:]]+=1$
 41                     currRead=[]$
 42                     currRead.append(line)$
 43                 else:$
 44                     currRead.append(line)$
 45     #now look at complexity$
 46     cntr=[]$
 47     for read in reads:$
 48         if 2<=sum(reads[read].values())<=5:$
 49             #cntr.append(float(len(reads[read]))/sum(reads[read].values()))$
 50             cntr.append(sum([1 for barcode in reads[read] if reads[read][barcode]==1])/float(len(reads[read])))$
 51     print '%s average +/-SD fraction of reads with 2-5 counts with unique barcodes: % s+/- %s'%(inFile,numpy.average(cntr), numpy.std(cntr))$
 52 $
 53 if __name__=='__main__':$
 54     """run as python readCollapser.py inFile.fastq outPrefix"""$
 55     Tee()$
 56     main(sys.argv[1:])$
