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
#
#run as Rscript median_normalize_DESeq.r file.geneCt file.conditions file.libType outFile
#This enables command line arguments
args<-commandArgs(TRUE)
#variable names for each input
geneCtFileName <- args[1]
conditionsFileName <- args[2]
libTypeFileName <- args[3]
outFileName <- args[4]

#This calls the DESeq library
library("DESeq")

#inputs the gene count file into variable countData
countData <- read.table( geneCtFileName, header=TRUE, row.names=1 )

#same with conditions file and libTypes
conditions <- strsplit(paste(readLines(conditionsFileName)),'\t')
libTypes <- strsplit(paste(readLines(libTypeFileName)),'\t')

#Next line is likely unnecessary and I comment it out.
#studyDesign=data.frame(row.names=colnames(countData),condition=conditions,libType=libTypes)


condition=factor(conditions[[1]])
cds = newCountDataSet( countData, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )

#finally, write the output
write.table(counts( cds, normalized=TRUE ), outFileName , sep="\t")

#edit March 4, 2014
if(length(args)>4) {
    uniqueConditions <- unique(conditions[[1]])
    for (i in 2:length(uniqueConditions)) {
        cds <- newCountDataSet( countData, as.factor(c(uniqueConditions[1],uniqueConditions[i])) )
        cds = estimateSizeFactors( cds )
        sizeFactors( cds )
        res = nbinomTest( cds, uniqueConditions[1], uniqueConditions[i])
        resSig <- res[which(res$padj < 0.05),]
        write.csv( resSig, file=paste(outFileName + "_1_" + i,"diffExpression",sep="_"))
    }
}

