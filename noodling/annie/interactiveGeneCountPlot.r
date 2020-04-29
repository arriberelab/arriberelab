#########################
# Annie Courney 4/28/2020 
# Script to plot genecounts (output from medianNormalizer), highlighting specific genes in a line-separated file
# where the first line should be the a column name
#########################

library("plotly")
library("tidyverse")
library("htmlwidgets")
library("ggplot2")
library("gghighlight")

# set theme to get gghighlight to work
theme_set(theme_bw())

# assign csv files to variables
DESeqGeneCts <- read.csv("C:/Users/aecou/github/arribere_lab/working/smg6RnaSeq/200428_SJA246-248_S.DESeqgeneCts.csv")
GeneNames <- read.csv("C:/Users/aecou/github/arribere_lab/working/smg6RnaSeq/GeneNames.csv")

# assign the ggplot to variable genePlot
genePlot <- ggplot(data=DESeqGeneCts, aes(x=X200210_SJA246Nugen, y=X200210_SJA248Nugen)) + geom_point(aes(colour=Gene_Name)) + gghighlight(DESeqGeneCts[,1] %in% GeneNames[,1]) + labs(x="wt", y="smg-6") + coord_trans(x="log10", y="log10")

# plot genePlot
ggplotly(genePlot)

