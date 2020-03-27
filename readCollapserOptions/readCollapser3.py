"""
Joshua Arribere Feb 24, 2014
Converted to python 3: Mar 23, 2020

Script to collapse reads containing the exact same sequence. 
    Will pick the read information from the first occurence of that read for the outFile.
    
Mar 24, 2020: 
    Contains main functions to toggle between for various UMI situations.
    Make sure to comment out the nonapplicable sections!!

Input: reads.fastq - fastq-formatted reads

Output: reads.collapsed.fastq - same as input, with redundant reads removed

run as python readCollapser3.py inFile.fastq outPrefix
"""

import sys, common, collections, numpy
from logJosh import Tee

"""

############################################################################################################
#Use the following to trim 6N from the 3'end
############################################################################################################
def main(args):
    inFile,outPrefix=args[0:]
    print('trimming 6 Ns from the 3\'end')
    reads=collections.defaultdict(lambda:collections.defaultdict(int))
    with open(inFile,'r') as f:
        with open(outPrefix,'w') as g:
            currRead=[]
            for line in f:
                line=line.strip()
                if len(currRead)==4:
                    #change the bracketed numbers for currReads 1 through 3 below based on the UMI used
                    UMI=currRead[1][-6:]
                    if reads[currRead[1][:-6]][UMI]==0:
                        g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1][:-6],currRead[2],currRead[3][:-6]))
                        if len(currRead[1])!=len(currRead[3]):
                            print(currRead, sys.exit())
                        #must trim quality scores and reads simultaneously. Hence the double -6
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
    print('%s average +/-SD fraction of reads with 2-5 counts with unique barcodes: % s+/- %s'%(inFile,numpy.average(cntr), numpy.std(cntr)))

"""

############################################################################################################
#Use the following to trim 6N from the 3'end and 4N from the 5'end
############################################################################################################
def main(args):
    inFile,outPrefix=args[0:]                                           #order of stuff given in command line. outPrefix is second thing given
    print('trimming 6 Ns from the 3\'end and 4nts from 5\'end')
    reads=collections.defaultdict(lambda:collections.defaultdict(int))  
    #create dictionary of reads? will automatically add keys as they're called, with value as an integer
    with open(inFile,'r') as f:                                         #read inFile as f
        with open(outPrefix,'w') as g:                                  #write outFile as g
            currRead=[]                                                 #create list
            for line in f:                                              #lines in inFile
                line=line.strip()                                       #remove first character
                if len(currRead)==4:                                    #has 4 lines?
                    #UMI=currRead[1][:4]+currRead[1][-6:]               #if you use this, you must use other if reads... below
                    UMI=currRead[1][-6:]                                #UMI = last 6 characters of line 1. line 1 is 1-index in currRead
                    #line 1 contains the read. line 3 contains the quality score.
                    #if reads[currRead[1][4:-6]][UMI]==0:
                    if reads[currRead[1][:-6]][UMI]==0:                 #check if we've seen this read+UMI before
                        g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1][4:-6],currRead[2],currRead[3][4:-6]))
                        #write 0th w/o spaces, 1st from index 4 to 6th from last (don't include -6), 2nd, and 3rd same as 1st
                        if len(currRead[1])!=len(currRead[3]):          #only want reads the same length as quality score
                            print(currRead, sys.exit())
                    #reads[currRead[1][4:-6]][UMI]+=1
                    reads[currRead[1][:-6]][UMI]+=1                     #mark this one as done?
                    currRead=[]                                         #clear list
                    currRead.append(line)                               #add new line to list and start over
                else:
                    currRead.append(line)                               #handles first case

    #now look at complexity
    cntr=[]
    for read in reads:
        if 2<=sum(reads[read].values())<=5:
            #cntr.append(float(len(reads[read]))/sum(reads[read].values()))
            cntr.append(sum([1 for barcode in reads[read] if reads[read][barcode]==1])/float(len(reads[read])))
    print('%s average +/-SD fraction of reads with 2-5 counts with unique barcodes: % s+/- %s'%(inFile,numpy.average(cntr), numpy.std(cntr)))

"""

############################################################################################################
#Use the following to trim NOTHING and NOT COLLAPSE READS
############################################################################################################
def main(args):
    inFile,outPrefix=args[0:]
    print('NOT TRIMMING 6 Nts FROM THE 3\'end!!!')
    print ('NOT COLLAPSING READS!')
    reads=collections.defaultdict(lambda:collections.defaultdict(int))
    with open(inFile,'r') as f:
        with open(outPrefix,'w') as g:
            currRead=[]
            for line in f:
                line=line.strip()
                if len(currRead)==4:
                    UMI=currRead[1][-6:]
                    if True:
                        g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1],currRead[2],currRead[3]))
                        if len(currRead[1])!=len(currRead[3]):
                            print(currRead, sys.exit())
                        #must trim quality scores and reads simultaneously. Hence the double -6
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
    print('%s average +/-SD fraction of reads with 2-5 counts with unique barcodes: % s+/- %s'%(inFile,numpy.average(cntr), numpy.std(cntr)))

"""

############################################################################################################

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
