"""
Joshua Arribere Feb 24, 2014
Still python2 2

Script to collapse reads containing the exact same sequence. Will pick the read information
    from the first occurence of that read for the outFile.

Input: reads.fastq - fastq-formatted reads

Output: reads.collapsed.fastq - same as input, with redundant reads removed

run as python readCollapser.py inFile.fastq outPrefix
"""
import sys, common, collections, numpy
from logJosh import Tee

def main(args):
    inFile,outPrefix=args[0:]
    #print 'Trimming of 6 nucleotides from 3\'end and 3 nts from 5\'end'
    #print 'NOT TRIMMING 6 Nts FROM THE 3\'end!!!'
    print 'trimming 6 Ns from the 3\'end and 4nts from 5\'end'
    #print 'trimming 6 Ns from the 3\'end'
    #print 'NOT COLLAPSING READS!'
    reads=collections.defaultdict(lambda:collections.defaultdict(int))
    with open(inFile,'r') as f:
        with open(outPrefix,'w') as g:
            currRead=[]
            for line in f:
                line=line.strip()
                if len(currRead)==4:
                    #UMI=currRead[1][:4]+currRead[1][-6:]
                    UMI=currRead[1][-6:]
                    #if reads[currRead[1][4:-6]][UMI]==0:
                    if reads[currRead[1][:-6]][UMI]==0:
                        #the next line has been moved ~6 lines down
                        #reads[currRead[1]]+=1#also need to change this line for 3 or 6EDIT: DO NOT NEED TO CHANGE: This keeps track of read+barcode, and only writes if read appears with new barcode
                        #g.write('%s\n%s\n%s\n%s\n'%(currRead[0],currRead[1],currRead[2],currRead[3]))
                        #g.write('%s\n%s\n%s\n%s\n'%(currRead[0],currRead[1][:-6],currRead[2],currRead[3][:-6]))
                        g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1][4:-6],currRead[2],currRead[3][4:-6]))
                        if len(currRead[1])!=len(currRead[3]):
                            print currRead, sys.exit()
                        #must trim quality scores and reads simultaneously. Hence the double -6
                    #reads[currRead[1][4:-6]][UMI]+=1
                    reads[currRead[1][:-6]][UMI]+=1
                    currRead=[]
                    currRead.append(line)
                else:
                    currRead.append(line)
    #now look at complexity
    cntr=[]
    for read in reads:
        if 2<=sum(reads[read].values())<=5:
            #cntr.append(float(len(reads[read]))/sum(reads[read].values()))
            cntr.append(sum([1 for barcode in reads[read] if reads[read][barcode]==1])/float(len(reads[read])))
    print '%s average +/-SD fraction of reads with 2-5 counts with unique barcodes: % s+/- %s'%(inFile,numpy.average(cntr), numpy.std(cntr))

if __name__=='__main__':
    """run as python readCollapser.py inFile.fastq outPrefix"""
    Tee()
    main(sys.argv[1:])
