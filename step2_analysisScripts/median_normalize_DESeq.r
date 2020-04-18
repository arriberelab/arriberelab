#Joshua Arribere Oct 7, 2013
#
#Script to Median normalize a table of gene cts
#
#Input: file.geneCt - output of medianNormalizedWithDESeq.writeGeneCtFile script
#   file.conditions - the conditions that correspond to each of the libraries. Strings.
#   file.libType - probably just tab-delimited 'single-end'
#
#Output: median normalized output of DESeq
#
#EDIT: Added arg[5] option, where it will run differential expression analysis if True
#EDIT: Dec 10, 2017 - JOSH updated to use DESeq2
#
#run as Rscript median_normalize_DESeq2.r file.geneCt file.conditions file.libType outFile

#This enables command line arguments
args<-commandArgs(TRUE)

#This calls the DESeq library
library("DESeq")

#inputs the gene count file into variable countTable
countTable=read.table( args[1], header=TRUE, row.names=1 )#this line works OK

#same with conditions file and libTypes
conditions <- strsplit(paste(readLines(args[2])),'\t')
libTypes <- strsplit(paste(readLines(args[3])),'\t')

#Next line is likely unnecessary and I comment it out.
#studyDesign=data.frame(row.names=colnames(countTable),condition=conditions,libType=libTypes)


condition=factor(conditions[[1]])
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )

#finally, write the output
write.table(counts( cds, normalized=TRUE ), args[4] , sep="\t")

#edit March 4, 2014
if(length(args)>4) {
    cds = estimateDispersions( cds )
    uniqueConditions <- unique(conditions[[1]])
    res = nbinomTest( cds, uniqueConditions[1], uniqueConditions[2])
    write.csv( res, file=paste(args[4],"diffExpression",sep="_"))
}
