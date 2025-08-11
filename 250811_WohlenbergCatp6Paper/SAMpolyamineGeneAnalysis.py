"""
Chloe Wohlenberg 4/10/24
Python3

Script to calculate a score representing the prediction of low vs high polyamine concentrations based on ribosome occupancy in different regions of polyamine synthesis genes
THIS IS AN ALTERNATIVE VERSION OF polyamineGeneAnalysis.py, where input is .SAM files rather than .JAM files

Input:
readsFile - path to library sam-formatted file
readSizes - as a list of comma-delimited sizes. Will only consider reads of these sizes. By default, align by the 3'end.
offset - ribosome offset size needed to adjust location of to desired ribosome site (in nts) translating.
chromsome - what chromosome of the genomic region you're looking at 
genomic bounds of (two) comparable regions for gene of choice (ex. ORF1 vs ORF2)
        start_region1 - start position of region 1 of gene
        end_region1 - end position of region 1 of gene
        start_region2 - start position of region 2 of gene
        end_region2 - end position of region 2 of gene
NOTE: genomic bounds must be listed from min to max for each region. So if your region is on opposite strand, swap your start/end region positions order

Output: numerical ratio score for gene of interest in library of interest

Input: path to library sam file, desired FP size (Psite = 16), desired offset, chromosome of region of interest, genomic bounds of (two!) comparable regions for gene of choice (ex. ORF1 vs ORF2)
"""
import sys
import collections

#first we take the sam file and make a dictionary based on desired read size and offset. Note that genomic position is used, rather than transcript position. This is modified from metaGeneForaBunchofPositions3.py       
def getReads(readsFile, readSizes, chromosome_name, offset):
        """
        Will convert a sam file, restricting to a list of readSizes to a dict of the format:
        {txtName:{position:ct}} where position is relATG and of the 3'end of the read
        """
        print('Using offset such that we record the modified positions of the desired offset.')
        reads_dict= {}
        with open(readsFile,'r') as f:
                       # Skips the first 11 lines to account for the header of file.
                        for i in range(11):
                        	f.readline()
                        for line in f:
                                line=line.strip().split('\t')
                                ##
                                readLen=len(line[9])
                                if readLen == readSizes:
                                        ##
                                        chromosome = line[2]
                                        ##
                                        if chromosome == chromosome_name:
                                                readStartPos = int(line[3]) + offset
                                                if readStartPos not in reads_dict:
                                                        reads_dict[readStartPos] = 1
                                                else:
                                                        reads_dict[readStartPos] += 1

        return reads_dict

#now we want to sum the total number of reads in region1 vs region2. This is also modified from metaGeneForaBunchofPositions3.py
def getRegionReads(reads_dict, start_region, end_region):
        """
        Will get the read counts for all positions w/in a window of position. Will then calculate
        the total number of reads. If this is lower than some number (perhaps
        zero), will return a dict of zeros.
        """
        print('Summing total reads per region')

        region_reads = 0
        for readStartPos in reads_dict:
                if (readStartPos >= (start_region) and readStartPos <= (end_region)):
                        region_reads += reads_dict[readStartPos]

        return region_reads

if __name__ == "__main__":
        readsFile = sys.argv[1]
        readSize = int(sys.argv[2])
        offset = int(sys.argv[3])
        chromosome_name = sys.argv[4]
        start_region1 = int(sys.argv[5])
        end_region1 = int(sys.argv[6])
        start_region2 = int(sys.argv[7])
        end_region2 = int(sys.argv[8])


        reads_dict = getReads(readsFile, readSize, chromosome_name, offset)

        region1_total = getRegionReads(reads_dict, start_region1, end_region1)
        region2_total = getRegionReads(reads_dict, start_region2, end_region2)

        #next, calculate the ratio of region1 reads to region2
        print(f"region 1 total reads: {region1_total}.\nregion 2 total reads: {region2_total}")
        PA_score = region1_total / region2_total
        print(f"{PA_score} is the ratio between gene region 1 from {start_region1} to {end_region1} and gene region 2 from {start_region2} to {end_region2}")
