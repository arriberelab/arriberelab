#########################
# Annie Courney 4/28/2020 
#Script to plot genecounts (output from medianNormalizer), highlighting specific genes in a csv file
#########################

library("plotly")
library("tidyverse")
library("htmlwidgets")
library("ggplot2")
library("gghighlight")

#set theme to get gghighlight to work
theme_set(theme_bw())

#assign csv files to variables
DESeqGeneCts <- read.csv("C:/Users/aecou/github/arribere_lab/working/smg6RnaSeq/200428_SJA246-248_S.DESeqgeneCts.csv")
GeneNames <- read.csv("C:/Users/aecou/github/arribere_lab/working/smg6RnaSeq/GeneNames.csv")

#make a list of important genes
importantGenes <- which(DESeqGeneCts[,1] %in% GeneNames)

#assign the ggplot to variable genePlot
genePlot <- ggplot(data=DESeqGeneCts, aes(x=X200210_SJA246Nugen, y=X200210_SJA248Nugen)) + geom_point(alpha=0.7, colour = "#51A0D5") + gghighlight(Gene_Name %in% importantGenes,
                                                                                                                                                   unhighlighted_params = list("darkred", 0.4),
                                                                                                                                                   use_direct_label = TRUE,
                                                                                                                                                   label_key = name,
                                                                                                                                                   label_params = list(size = 5)) + labs(x="wt", y="smg-6") + coord_trans(x="log10", y="log10")

#plot it
ggplotly(genePlot)


#attach(DESeqGeneCts)
#attach(GeneNames)
#par(mfrow=c(1,2), mar=c(4.1, 3.79, 3.1, .1))
#plot(x=DESeqGeneCts$X200210_SJA246Nugen, y=DESeqGeneCts$X200210_SJA248Nugen, xlab="wt", ylab="smg-6(srf2102; D1070A)", log="xy")
#plot(x=DESeqGeneCts$X200210_SJA247Nugen, y=DESeqGeneCts$X200210_SJA248Nugen, xlab="smg-1(e1228)", ylab=" ", log="xy")
#identify(x=DESeqGeneCts$X200210_SJA246Nugen, y=DESeqGeneCts$X200210_SJA248Nugen, labels=DESeqGeneCts$Gene_Name)
#identify(x=DESeqGeneCts$X200210_SJA247Nugen, y=DESeqGeneCts$X200210_SJA248Nugen, labels=DESeqGeneCts$Gene_Name)

