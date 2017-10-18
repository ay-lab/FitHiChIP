#!/usr/bin/env Rscript

#===========================================================
# R script for combining the mappability,GC content on the pre-specified RE fragment file

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================

args <- commandArgs(TRUE)

REFragFile <- args[1]
Temp_Mapp_File <- args[2]
Temp_GC_File <- args[3]
REFragMappGCFile <- args[4]

REFragData <- read.table(REFragFile, header=F)
Temp_Mapp_Data <- read.table(Temp_Mapp_File, header=F)
Temp_GC_Data <- read.table(Temp_GC_File, header=F)

write.table(cbind(REFragData[1:3], Temp_Mapp_Data, Temp_GC_Data, REFragData[4:ncol(REFragData)]), REFragMappGCFile, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

