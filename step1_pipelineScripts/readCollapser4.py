"""
Joshua Arribere Feb 24, 2014
Converted to python 3: Mar 23, 2020

Script to collapse reads containing the exact same sequence. 
    Will pick the read information from the first occurence of that read for the outFile.
    Will trim off UMI nts as specified by user.
    
Mar 30, 2020: Parissa and Matt updated to make use of UMI lengths given in command line.

Input: reads.fastq - fastq-formatted reads

Output: reads.collapsed.fastq - same as input, with redundant reads removed

run as python3 readCollapser3.py inFile.fastq number5'N number3'N outPrefix
"""

import sys, common, collections, numpy
from logJosh import Tee

def main(args):
    inFile,x,y,outPrefix=args[0:]
    x=int(x)        #5'N
    y=int(y)*-1     #3'N, made negative
    print('trimming %sNs from the 3\'end and %sNs from 5\'end'%(y*-1,x))

    reads=collections.defaultdict(lambda:collections.defaultdict(int))  
    #create nested dict, automatically adding keys as they're called and only if they haven't been seen before.
    #key = read. value = dict with key = UMI, value = count
    
    with open(inFile,'r') as f:                                            #read inFile as f
        with open(outPrefix,'w') as g:                                     #write outFile as g
            currRead=[]                                                    #create list
            for line in f:                                                 #lines in inFile
                line=line.strip()                                          #make everything one line
                if len(currRead)==4:                                       #stop when currRead has 4 items in list
                    if y==0:
                        UMI=currRead[1][:x]
                        if reads[currRead[1][x:]][UMI]==0:
                            g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1][x:],currRead[2],currRead[3][x:]))
                            if len(currRead[1])!=len(currRead[3]):
                                print(currRead), sys.exit()
                        reads[currRead[1][x:]][UMI]+=1
                        currRead=[]
                        currRead.append(line)
                    else:
                        UMI=currRead[1][:x]+currRead[1][y:]
                        if reads[currRead[1][x:y]][UMI]==0:                #check if we've seen this read+UMI before
                            g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1][x:y],currRead[2],currRead[3][x:y]))
                            if len(currRead[1])!=len(currRead[3]):         #only want reads the same length as quality score
                                print(currRead), sys.exit()
                        reads[currRead[1][x:y]][UMI]+=1                    #mark this one as done
                        currRead=[]                                        #clear list
                        currRead.append(line)                              #add line to list
                else:
                    currRead.append(line)                                  #add line by line to currRead list
    
    #now look at complexity
    cntr=[]
    for read in reads:
        if 2<=sum(reads[read].values())<=5: #use range 2-5 because they're likely PCR duplicates. beyond 5 are probably rRNA
            #cntr.append(float(len(reads[read]))/sum(reads[read].values()))
            cntr.append(sum([1 for barcode in reads[read] if reads[read][barcode]==1])/float(len(reads[read])))
    print('%s average +/-SD fraction of reads with 2-5 counts with unique barcodes: % s+/- %s'%(inFile,numpy.average(cntr), numpy.std(cntr)))

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
