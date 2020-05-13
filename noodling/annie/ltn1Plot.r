#read in data
r293Data <- read.csv("C:/Users/aecou/github/arribere_lab/working/ltn1/r293Data.csv")
r274Data <- read.csv("C:/Users/aecou/github/arribere_lab/working/ltn1/r274Data.csv")
tra2Data <- read.csv("C:/Users/aecou/github/arribere_lab/working/ltn1/tra2Data.csv")

#import libraries
library("ggplot2")
library("gridExtra")
library("tidyverse")

#put all the plots into a list first to arrange them later
ltn1_plots = list()

grid.arrange(ltn1_plots$sc1, ltn1_plots$sc2, ncol = 3)

#Note: do NOT use geom_point and geom_jitter or data points will be doubled
ltn1_plots$r293Plot <- ggplot(data=r293Data, aes( x=Mutation, y=BPM, group=Mutation, colour=Mutation  )) + geom_jitter(width = 0.3) + labs( title="A", y="Beats Per Minute of Animals in EN50") +  theme(legend.position = "none", axis.text.x = element_text(color = "grey20", size = 8, angle = 45, hjust = .5, vjust = .5, face = "plain"), axis.title.y = element_text(color = "grey20", size = 7, angle = 90, hjust = .5, vjust = .5, face = "plain"))
ltn1_plots$r274Plot <- ggplot(data=r274Data, aes( x=Mutation, y=Progeny, group=Mutation, colour=Mutation )) + geom_jitter(width = 0.3, height = 0) + labs( title="B", y="Progeny Per Animal" ) +  theme(legend.position = "none", axis.text.x = element_text(color = "grey20", size = 8, angle = 45, hjust = .5, vjust = .5, face = "plain"), axis.title.y = element_text(color = "grey20", size = 7, angle = 90, hjust = .5, vjust = .5, face = "plain"))
ltn1_plots$tra2Plot <- ggplot(data=tra2Data, aes(x=Mutation, y=Progeny, group=Mutation, colour=Mutation ))  + geom_jitter(width = 0.3, height = 0)  + labs( title="C", y="Progeny Per Animal" ) +  theme(legend.position = "none", axis.text.x = element_text(color = "grey20", size = 8, angle = 45, hjust = .5, vjust = .5, face = "plain"), axis.title.y = element_text(color = "grey20", size = 7, angle = 90, hjust = .5, vjust = .5, face = "plain"))

grid.arrange(ltn1_plots$r293Plot, ltn1_plots$r274Plot, ltn1_plots$tra2Plot, ncol = 4)

