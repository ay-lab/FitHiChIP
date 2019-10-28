#!/usr/bin/env Rscript

#======================
# script to sort the input locus pairs according to their genomic distance (ascending order)
# required for spline fitting
#======================

args <- commandArgs(TRUE)

InpFile <- args[1]
OutFile <- args[2]

tempDistFile <- paste0(dirname(OutFile), '/temp_Genomic_Distance.bed')
system(paste0("awk -F[\'\t\'] \'function abs(v) {return v < 0 ? -v : v} {if (NR>1) {print abs($5-$2)}}\' ", InpFile, " | sort -k1,1n | uniq > ", tempDistFile))

# copy the first line in the output file
system(paste("awk \'{if (NR==1) {print $0}}\' ", InpFile, " > ", OutFile))

DistData <- read.table(tempDistFile, header=F, sep="\t", stringsAsFactors=F)
DistValues <- sort(as.integer(DistData[,1]))

# temp_Dist_LocusPairFile <- paste0(dirname(OutFile), '/temp_Dist_LocusPair.bed')
for (i in 1:length(DistValues)) {
	currDist <- DistValues[i]
	cat(sprintf("\n Processing distance value for sorting locus pairs based on interaction distance: %s ", currDist))

	system(paste0("awk -F[\'\t\'] \'function abs(v) {return v < 0 ? -v : v} {if ((NR>1) && (abs($5-$2) == ", currDist, ")) {print $0}}\' ", InpFile, " | sort -k7,7nr >> ", OutFile))

}	# end distance value loop

system(paste("rm", tempDistFile))
# system(paste("rm", temp_Dist_LocusPairFile))

