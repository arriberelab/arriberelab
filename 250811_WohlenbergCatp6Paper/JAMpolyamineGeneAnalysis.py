"""
Chloe Wohlenberg 3/23/24
Python3

Script to calculate a score representing the prediction of low vs high polyamine concentrations based on ribosome occupancy in different regions of polyamine synthesis genes

Input:
readsFile - path to library jam-formatted file
readSizes - as a list of comma-delimited sizes. Will only consider reads of these sizes. By default, align by the 3'end.
offset - ribosome offset size needed to adjust location of to desired ribosome site (in nts) translating.
gene_name - Wormbase gene sequence name in which regions come from
transcript bounds of (two) comparable regions for gene of choice (ex. ORF1 vs ORF2)
	start_region1 - start position of region 1 of gene
	end_region1 - end position of region 1 of gene
	start_region2 - start position of region 2 of gene
	end_region2 - end position of region 2 of gene

Output: numerical ratio score for gene of interest in library of interest

Input: path to library jam file, genomic bounds of (two!) comparable regions for gene of choice (ex. ORF1 vs ORF2), desired FP size, offset size (I think we want to be in Psite)
"""
import sys
import collections

#first we take the jam file and make a dictionary based on desired read size and offset. This is taken from metaGeneForaBunchofPositions3.py
def getReads(readsFile, readSizes, offset):
	"""
	Will convert a jam file, restricting to a list of readSizes to a dict of the format:
	{txtName:{position:[ct, readid_list]}} where position is relATG and of the 3'end of the read
	"""
	print('Using offset such that we record the modified positions of the desired offset.')
	#aa=collections.defaultdict(lambda:collections.defaultdict(int))
	aa = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, []]))
	##
	with open(readsFile,'r') as f:
		f.readline()#skips first line
		for line in f:
			line=line.strip().split('\t')
			##
			readLen=len(line[6])
			if readLen == readSizes:
				##
				gene=line[8]
				txts=line[9]
				readid=line[0]
				##
				strand=gene.split(':')[1]
				if strand=='S':
					txts=txts.split('|')
					totalTxts=len(txts)
					for txt in txts:
						txtName,posRelStart,blah=txt.split(':')
						##
						#aa[txtName][int(posRelStart)+offset]+=1#/totalTxts
						pos = int(posRelStart) + offset
						aa[txtName][pos][0] += 1  # increment count
						aa[txtName][pos][1].append(readid)  # store readid
	##
	return aa

#now we want to sum the total number of reads in region1 vs region2. This is also modified from metaGeneForaBunchofPositions3.py
def getRegionReads(reads_dict, gene_name,  start_region, end_region):
	"""
	reads_dict={positions:ct}, position is some number. Will get the read
	counts for all positions w/in a window of position. Will then calculate
	the total number of reads. If this is lower than some number (perhaps
	zero), will return a dict of zeros.
	"""
	print('Summing total reads per region')
	aa={}
	read_ids = []
	#for gene in reads_dict.keys():
	#	if gene == gene_name:
	#		print(gene)
	#if start_region < 0
	for k in range(start_region,end_region+1):
		if k in reads_dict[gene_name]:
			aa[k]=reads_dict[gene_name][k][0]
			#also collect read ids here
			read_ids.extend(reads_dict[gene_name][k][1])
		else:
			aa[k]=0
	##
	total = sum(aa.values())
	##
	if total > 0:
		for k, v in aa.items():
			aa[k] = v / total
	##
	return aa, total, read_ids


#where functions are called
if __name__ == "__main__":
	readsFile = sys.argv[1]
	readSize = int(sys.argv[2])
	offset = int(sys.argv[3])
	gene_name = sys.argv[4]
	start_region1 = int(sys.argv[5])
	end_region1 = int(sys.argv[6])
	start_region2 = int(sys.argv[7])
	end_region2 = int(sys.argv[8])

	reads_dict = getReads(readsFile, readSize, offset)

	#for gene in reads_dict:
		#if gene.startswith("F47G4"):
			#print(gene)

	region1_reads,region1_total,region1_read_ids = getRegionReads(reads_dict, gene_name, start_region1, end_region1)
	region2_reads,region2_total,region2_read_ids = getRegionReads(reads_dict, gene_name, start_region2, end_region2)

	print(f"Region 1 read IDs:")
	for read_id in region1_read_ids:
		print(read_id)

	#next, calculate the ratio of region1 reads to region2
	print(f"region 1 total reads: {region1_total}.\nregion 2 total reads: {region2_total}")
	PA_score = region1_total / region2_total
	print(f"{PA_score} is the ratio between gene region 1 from {start_region1} to {end_region1} and gene region 2 from {start_region2} to {end_region2}")
