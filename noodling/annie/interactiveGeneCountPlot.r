library("plotly")
library("tidyverse")
library("htmlwidgets")
library("ggplot2")

DESeqGeneCts <- read.csv("C:/Users/aecou/github/arribere_lab/working/smg6RnaSeq/200428_SJA246-248_S.DESeqgeneCts.csv")
GeneNames <- read.csv("C:/Users/aecou/github/arribere_lab/working/smg6RnaSeq/GeneNames.csv")

genePlot <- ggplot(data=DESeqGeneCts, aes(x=X200210_SJA246Nugen, y=X200210_SJA248Nugen, text=paste("Gene: ", GeneNames))) + geom_point(alpha=0.7, colour = "#51A0D5") + labs(x="wt", y="smg-6") + theme_classic()

ggplotly(genePlot, tooltip=text)


#attach(DESeqGeneCts)
#attach(GeneNames)
#par(mfrow=c(1,2), mar=c(4.1, 3.79, 3.1, .1))
#plot(x=DESeqGeneCts$X200210_SJA246Nugen, y=DESeqGeneCts$X200210_SJA248Nugen, xlab="wt", ylab="smg-6(srf2102; D1070A)", log="xy")
#plot(x=DESeqGeneCts$X200210_SJA247Nugen, y=DESeqGeneCts$X200210_SJA248Nugen, xlab="smg-1(e1228)", ylab=" ", log="xy")
#identify(x=DESeqGeneCts$X200210_SJA246Nugen, y=DESeqGeneCts$X200210_SJA248Nugen, labels=DESeqGeneCts$Gene_Name)
#identify(x=DESeqGeneCts$X200210_SJA247Nugen, y=DESeqGeneCts$X200210_SJA248Nugen, labels=DESeqGeneCts$Gene_Name)

