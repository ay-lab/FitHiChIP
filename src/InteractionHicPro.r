#!/usr/bin/env Rscript

#===========================================================
# R script for creating an interaction matrix (chromosome intervals + contact count)
# from a given HiC-pro output matrix (generated using the HiC-pro utility)
# currently cis interactions are returned

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================

library(tools)
library(data.table)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(TRUE)

# file storing the indices of individual bins of chromosomes
IntervalFile <- args[1]

# file storing the interactions between bins (matrix)
MatrixFile <- args[2]

# Initial file of interactions 
# (where the intervals + contact count will be written)
# both cis and trans interactions (with no distance constraints)
InitialInteractionFile <- args[3]

IntervalMat <- data.table::fread(IntervalFile, header=F, sep="\t", stringsAsFactors=F)	
colnames(IntervalMat) <- c("chr1","s1","e1","idx")

InpInteraction <- data.table::fread(MatrixFile, header=F, sep="\t", stringsAsFactors=F)
colnames(InpInteraction) <- c("idx1","idx2","cc")

# merge with respect to the index of 1st chromosome interval
# the output data frame has the 1st chromosome intervals in the last three columns
df1 <- merge(x=InpInteraction, y=IntervalMat, by.x=colnames(InpInteraction)[1], by.y=colnames(IntervalMat)[4])
colnames(df1) <- c(colnames(InpInteraction), "chr1","s1","e1")
# write.table(df1, 'merge1.bed', row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

# merge with respect to the index of 2nd chromosome interval
# the output data frame has the 2nd chromosome intervals in the last three columns
df2 <- merge(x=df1, y=IntervalMat, by.x=colnames(InpInteraction)[2], by.y=colnames(IntervalMat)[4])
colnames(df2) <- c(colnames(df1), "chr2","s2","e2")
# write.table(df2, 'merge2.bed', row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

# select the 4th to last column (signifying the pair of chromosome intervals)
# and also the third column (contact count)
df3 <- df2[,c(4:ncol(df2),3)]
colnames(df3) <- c(colnames(df2)[4:ncol(df2)], "cc")
write.table(df3, InitialInteractionFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

