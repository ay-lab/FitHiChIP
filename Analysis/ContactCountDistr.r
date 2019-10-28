#!/usr/bin/env Rscript

#===========================================================
# R script for finding the significant contact count distribution for individual peaks

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript ContactCountDistr.r $PeakFile $InteractionFile $PlotFile $OutText
#===========================================================

library(data.table)

options(scipen = 999)
options(datatable.fread.datatable=FALSE)

args <- commandArgs(TRUE)

PeakFile <- args[1]
InteractionFile <- args[2]
PlotFile <- args[3]
OutText <- args[4]

# peaks in the input peak file (PD = PeakData)
# PD <- read.table(PeakFile, header=FALSE)
PD <- data.table::fread(PeakFile, header=FALSE)

# peaks from the interaction file (PI = peaks from interactions)
# PI <- read.table(InteractionFile, header=FALSE)
PI <- data.table::fread(InteractionFile, header=FALSE)

# significant contact count (filtered according to Q value) for individual peaks
count <- c()

for (i in (1:nrow(PD))) {
	chr1 <- as.character(PD[i,1])
	s1 <- as.integer(PD[i,2])
	e1 <- as.integer(PD[i,3])

	# cat(sprintf("\n\n ===>>> checking chr1: %s s1: %s e1: %s \n ", chr1, s1, e1))

	# overlap between peaks in the peak file and in the interaction file
	# comparison with respect to three columns
	# note: use of '&' operator
	c <- 0

	# set of indices matching with chr1
	idx <- which(PI[,1] == chr1)

	# cat(sprintf("\n idx: %s \n length of idx: %s \n ", idx, length(idx)))

	if (length(idx) > 0) {
		for (j_idx in (1:length(idx))) {
			j <- idx[j_idx]

			chr2 <- as.character(PI[j,1])
			s2 <- as.integer(PI[j,2])
			e2 <- as.integer(PI[j,3])

			# cat(sprintf("\n --- individual  chr1: %s  s1: %s  e1: %s chr2: %s s2: %s e2: %s \n ", chr1, s1, e1, chr2, s2, e2)) 		

			if ((!is.na(s1)) & (!is.na(e1)) & (!is.na(s2)) & (!is.na(e2))) {
				if ((s2 <= s1) & (e2 >= e1)) {
					c <- c + 1
				} else if ((s2 >= s1) & (e2 <= e1)) {
					c <- c + 1
				} else if ((s2 <= s1) & (e2 >= s1)) {
					c <- c + 1
				} else if ((s2 <= e1) & (e2 >= e1)) {
					c <- c + 1
				}			
			}
		}		
	}
	count[i] <- c
}

# dump the original peaks vs significant contact count in a text file
write.table(cbind(PD[,1:3], count), OutText, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)

# first sort the count vector and then draw the frequency distribution
pdf(PlotFile, width=14, height=10)
plot(sort(count, decreasing=TRUE), type = "o", cex=0.5, col="red", xlab="Peak counter", ylab="Contact count (filtered)")
title(sub('\\.pdf$', '', basename(PlotFile)))
dev.off()


