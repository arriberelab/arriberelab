##Joshua Arribere May 5, 2022
##
##Script to run DESeq2 given:
##inputCtFile
##conditionsFile
##treatedNum and untreated
##outPrefix
args <- commandArgs(TRUE)
##
cts <- as.matrix(read.csv(args[1],sep="\t",row.names=1))
##
coldata <- read.csv(args[2], sep="\t", header=FALSE, row.names=1)
##
colnames(coldata) <- c("condition", "type")
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
##
library("DESeq2")
##
dds <- DESeqDataSetFromMatrix(countData = cts, 
    colData = coldata, design = ~ condition)
##
dds <- DESeq(dds)
##
res <- results(dds, contrast=c("condition",args[3],args[4]))
##
write.csv(as.data.frame(res),file=args[5])
