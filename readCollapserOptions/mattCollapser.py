 """$
 2 Joshua Arribere Feb 24, 2014$
  3 Converted to python 3: Mar 23, 2020$
  4 Script to collapse reads containing the exact same sequence. $
  5     Will pick the read information from the first occurence of that read for the outFile.$
  6     Will trim off UMI nts as specified by user.$
  7     $
  8 Mar 30, 2020: Matt made some edits to Parissa's edits.$
  9 Input: reads.fastq - fastq-formatted reads$
 10 Output: reads.collapsed.fastq - same as input, with redundant reads removed$
 11 run as python readCollapser3.py inFile.fastq number5'N number3'N outPrefix$
 12 """$
 13 $
 14 import sys, common, collections, numpy$
 15 from logJosh import Tee$
 16 $
 17 def main(args):$
 18     inFile,x,y,outPrefix=args[0:]$
 19     x=int(x) #5'N$
 20     y=int(y)*-1 #3'N makes the y negative for the 3' end$
 21     print(y)$
 22     if y == 0:  #replacing a 0 with blank space for slicing later$
 23         y = ''$
 24     print(y)$
 25     print('trimming %sNs from the 3\'end and %sNs from 5\'end'%(y,x))$
 26     reads=collections.defaultdict(lambda:collections.defaultdict(int))  $
 27     #create nested dict, automatically adding keys as they're called and only if they haven't been seen before.$
 28     #key = read. value = dict with key = UMI, value = count$
 29     with open(inFile,'r') as f:                                         #read inFile as f$
 30         with open(outPrefix,'w') as g:                                  #write outFile as g$
 31             currRead=[]                                                 #create list$
 32             for line in f:                                              #lines in inFile$
 33                 line=line.strip()                                       #make everything one line$
 34                 if len(currRead)==4:                                    #stop when currRead has 4 items in list$
 35                     if y == '':                                         #loop to deal with trimming 0 from the 3' end$
 36                         UMI=currRead[1][:x]$
 37                         #print(UMI)$
 38                     else:$
 39                         UMI=currRead[1][:x]+currRead[1][y:]            #UMI = first 4 and last 6 chars of line 1 aka 1-$
 40                         if reads[currRead[1][x:y]][UMI]==0:            #check if we've seen this read+UMI before$
 41                             g.write('%s\n%s\n%s\n%s\n'%(currRead[0].replace(' ','-'),currRead[1][x:y],currRead[2],currRead[3][x:y]))$
 42                             #write 0th w/o spaces, 1st from index 4 to 6th from last, 2nd, and 3rd same as 1st$
 43                             if len(currRead[1])!=len(currRead[3]):          #only want reads the same length as quality score$
 44                                 print(currRead, sys.exit())$
 45                         reads[currRead[1][x:y]][UMI]+=1          #mark this one as done$
 46                         currRead=[]                                         #clear list$
 47                         currRead.append(line)                               #add line to list$
 48                 else:$
 49                     currRead.append(line)                               #add line by line to currRead list$
 50 $
 51     #now look at complexity$
 52     cntr=[]$
 53     for read in reads:$
 54         if 2<=sum(reads[read].values())<=5: #use range 2-5 because they're likely PCR duplicates. beyond 5 are probably rRNA$
 55             #cntr.append(float(len(reads[read]))/sum(reads[read].values()))$
 56             cntr.append(sum([1 for barcode in reads[read] if reads[read][barcode]==1])/float(len(reads[read])))$
 57     print('%s average +/-SD fraction of reads with 2-5 counts with unique barcodes: % s+/- %s'%(inFile,numpy.average(cntr), numpy.std(cntr)))$
 
