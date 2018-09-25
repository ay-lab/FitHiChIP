#!/usr/bin/env Rscript

#===========================================================
# R script for differential analysis
# on FitHiChIP generated loops
# using EdgeR 
# only 3D differences, excluding 1D differences, are derived

# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
#===========================================================

# package to compute the overlap among intervals
suppressMessages(library(GenomicRanges))
library(optparse)
# library(DESeq2)
library(edgeR)
library(dplyr)
# library(UtilRPckg)
library(ggplot2)
# library(ggsignif)
library(tools)

# ggplot parameters
FontSize=18
PlotWidth=10
PlotHeight=6

#==========================
# getting number of lines of a file
#==========================
GetNumLines <- function(inpfile) {
	nline <- as.integer(system(paste("cat", inpfile, "| wc -l"), intern = TRUE))
	return(nline)
}

#=================
# computing overlap between peaks and loops
#=================
Overlap1Dwith2D <- function(Peakdata, Loopdata, boundary=1, offset=0, uniqov=TRUE) {

	ov1 <- as.data.frame(findOverlaps(GRanges(Peakdata[,1], IRanges(Peakdata[,2]+boundary-offset, Peakdata[,3]-boundary+offset)),GRanges(Loopdata[,1], IRanges(Loopdata[,2]+boundary-offset, Loopdata[,3]-boundary+offset))))

	ov2 <- as.data.frame(findOverlaps(GRanges(Peakdata[,1], IRanges(Peakdata[,2]+boundary-offset, Peakdata[,3]-boundary+offset)),GRanges(Loopdata[,4], IRanges(Loopdata[,5]+boundary-offset, Loopdata[,6]-boundary+offset))))

	if (uniqov == TRUE) {
		ov_A_B1 <- unique(ov1[,1])
		ov_A_B2 <- unique(ov2[,1])
		ov_B1_A <- unique(ov1[,2])
		ov_B2_A <- unique(ov2[,2])		
	} else {
		ov_A_B1 <- ov1[,1]
		ov_A_B2 <- ov2[,1]
		ov_B1_A <- ov1[,2]
		ov_B2_A <- ov2[,2]
	}

	Nonov_A_B1 <- setdiff(seq(1, nrow(Peakdata)), ov_A_B1)
	Nonov_A_B2 <- setdiff(seq(1, nrow(Peakdata)), ov_A_B2)
	Nonov_B1_A <- setdiff(seq(1, nrow(Loopdata)), ov_B1_A)
	Nonov_B2_A <- setdiff(seq(1, nrow(Loopdata)), ov_B2_A)

	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B1 = ov_A_B1, A_AND_B2 = ov_A_B2, B1_AND_A = ov_B1_A, B2_AND_A = ov_B2_A, A_MINUS_B1 = Nonov_A_B1, A_MINUS_B2 = Nonov_A_B2, B1_MINUS_A = Nonov_B1_A, B2_MINUS_A = Nonov_B2_A)
	return(newList)

}

#=================
# computing peak overlap
#=================
Overlap1D <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE) {

	ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))

	if (uniqov == TRUE) {
		ov_idx_file1 <- unique(ov1[,1])
		ov_idx_file2 <- unique(ov1[,2])		
	} else {
		ov_idx_file1 <- ov1[,1]
		ov_idx_file2 <- ov1[,2]
	}
	nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
	nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	return(newList)

}

#=================
# computing loop overlap
#=================
OverlapLoop <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE) {

	ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))

	ov2 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,4], IRanges(Inpdata1[,5]+boundary-offset, Inpdata1[,6]-boundary+offset)),GRanges(Inpdata2[,4], IRanges(Inpdata2[,5]+boundary-offset, Inpdata2[,6]-boundary+offset))))

	overlap_uniq_mat <- ov1[unique(which(paste(ov1[,1], ov1[,2], sep=".") %in% paste(ov2[,1], ov2[,2], sep="."))),]

	if (uniqov == TRUE) {
		ov_idx_file1 <- unique(overlap_uniq_mat[,1])
		ov_idx_file2 <- unique(overlap_uniq_mat[,2])		
	} else {
		ov_idx_file1 <- overlap_uniq_mat[,1]
		ov_idx_file2 <- overlap_uniq_mat[,2]
	}
	nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
	nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2, A_AND_B.df = Inpdata1[ov_idx_file1, ], B_AND_A.df = Inpdata2[ov_idx_file2, ], A_MINUS_B.df = Inpdata1[nonov_idx_file1, ], B_MINUS_A.df = Inpdata2[nonov_idx_file2, ])
	return(newList)
}

#==========================
# extracting chromosome data 
#==========================

ExtractChrData <- function(InpFile, chrName, OutFile=NULL, header=TRUE, dist=c(-1,-1), mid=FALSE) {
	if (is.null(OutFile)) {
		OutFile <- paste0(dirname(InpFile), "/temp_Chr_data.bed")
	}

	# process the distance thresholds
	# and insert in two variables
	if ((dist[1] > 0) & (dist[2] > 0) & (dist[2] > dist[1])) {
		distthrlow <- dist[1]
		distthrhigh <- dist[2]		
	} else {
		distthrlow <- -1
		distthrhigh <- -1
	}

	# condition based on using gzipped input file
	# or plain text file

	if (file_ext(InpFile) == "gz") {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	} else {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	}

}	# end function


#======================================================
# this function returns loop indices
# which are significant in the specified number of replicates
#======================================================
GetSpecificSigIdx <- function(IntDF, CCCols, QValCols, nRep) {
	data <- matrix(0, nrow=nrow(IntDF), ncol=length(CCCols))
	for (i in (1:length(CCCols))) {
		# get the loop indices such that 
		# current replicate is significant
		curr_cc_col <- CCCols[i]
		curr_q_col <- QValCols[i]
		Idx <- which((IntDF[, curr_cc_col] > 0) & (IntDF[, curr_q_col] < 0.01))
		# set the corresponding matrix entries to 1
		if (length(Idx) > 0) {
			data[Idx, i] <- 1 	
		}
	}

	# now return those loop indices (rows)
	# such that the columns having 1 equals to the specified "nRep"
	s <- rowSums(data)
	FinalIdx <- which(s == nRep)
	return(FinalIdx)

}

#======================================================
# this function returns loop indices
# such that none of the replicates are significant
#======================================================
GetNoneSigIdx <- function(IntDF, CCCols, QValCols) {
	for (i in (1:length(CCCols))) {
		# get the loop indices such that 
		# current replicate is significant
		curr_cc_col <- CCCols[i]
		curr_q_col <- QValCols[i]
		Idx <- which((IntDF[, curr_cc_col] == 0) | (IntDF[, curr_q_col] > 0.1))
		if (i == 1) {
			FinalIdx <- Idx
		} else {
			FinalIdx <- intersect(FinalIdx, Idx)
		}
	}
	return(FinalIdx)
}

#======================================================
# this function fills the "TargetVec" 
# according to the overlap of "PeakDF" and "IntDF"
# "targetcols" is a vector of two elements
# first element is the column index storing the binary (1/0) overlap indication
# second element is the column index storing the average peak strength
#======================================================
FillPeakVector <- function(TargetVec, PeakDF, IntDF, targetcols) {

	# get overlap of peak data and the interaction data
	ov1 <- as.data.frame(findOverlaps(GRanges(IntDF[,1], IRanges(IntDF[,2], IntDF[,3])), GRanges(PeakDF[,1], IRanges(PeakDF[,2], PeakDF[,3]))))

	LoopIdx <- ov1[,1]
	PeakIdx <- ov1[,2]

	# create a data frame such that its first column is the loop index
	# and the second column is the peak strength
	# of the corresponding overlapping peak
	Z <- cbind.data.frame(LoopIdx, PeakDF[PeakIdx, ncol(PeakDF)])
	colnames(Z) <- c('I', 'P')

	# now for each loop, aggregate the peak strength of the 
	# overlapping peaks
	ZA <- aggregate(P~I, data=Z, FUN=mean)

	# ZA has two columns
	# first column: index of loop
	# second column: mean peak strength
	# in the TargetVec, fill the columns indicated by "targetcols"
	# with the peak overlap and mean peak strength values
	TargetVec[ZA[,1], targetcols[1]] <- 1
	TargetVec[ZA[,1], targetcols[2]] <- ZA[,2]

	return(TargetVec)

}	# end function

#======================================================
# this function merges all the FitHiChIP loops
# (first 6 fields)
# of all the input files (present in the array of file names LoopList)
# in the final file "UnionLoopFile"
#======================================================
MergeLoops <- function(LoopList, UnionLoopFile, UnionLoopTempFile, FDRThr) {
	for (i in (1:length(LoopList))) {
		inpfile <- LoopList[i]
		if (i == 1) {
			if (file_ext(inpfile) == "gz") {
				system(paste0("zcat ", inpfile, " | awk \'((NR>1) && ($NF < ", FDRThr, "))\' - | cut -f1-6 > ", UnionLoopTempFile))
			} else {
				system(paste0("awk \'((NR>1) && ($NF < ", FDRThr, "))\' ", inpfile, " | cut -f1-6 > ", UnionLoopTempFile))
			}
		} else {
			if (file_ext(inpfile) == "gz") {
				system(paste0("zcat ", inpfile, " | awk \'((NR>1) && ($NF < ", FDRThr, "))\' - | cut -f1-6 >> ", UnionLoopTempFile))
			} else {
				system(paste0("awk \'((NR>1) && ($NF < ", FDRThr, "))\' ", inpfile, " | cut -f1-6 >> ", UnionLoopTempFile))
			}
		}
	}
	# sort the file, remove duplicate loops
	system(paste("sort -k1,1 -k2,2n -k5,5n", UnionLoopTempFile, "| awk -F\"\t\" \'!seen[$1, $2, $5]++\' - >", UnionLoopFile))
	cat(sprintf("\n ***** Constructed sorted interactions ****** \n"))
}	# end function


#================================
# this function annotates individual loops in the merged union set of loops
# according to the feature vectors of individual loops 
# (FitHiChIP significance loops - significant or )
# parameters:
# 1) UnionLoopFile: file containing merged significant loops from all input replicates and all categories
# 2) AllLoopList: FitHiChIP significance files (significant + non-significant interactions)
# 3) AllPeakFileList: list of peak files for two categories
# 4) AllChIPCovFileList: list of files containing ChIP-seq coverage for the target bin size (e.g. 5 Kb) for both categories
# 5) BiasFileList: list of files containing the HiChIP specific bias values for the given bins
# 6) ChrList_NameNum: list of chromosome names
# 7) UseRawCC: if 1, raw contact count is used. Else, ratio of raw and the predicted contact count is used.
# 8) AllRepLabels: labels of individual categories
# 9) CategoryList: two categories experimented
#================================
# FillFeatureValues <- function(UnionLoopFile, AllLoopList, AllPeakFileList, AllChIPCovFileList, BiasFileList, ChrList_NameNum, UseRawCC, AllRepLabels, CategoryList) {

FillFeatureValues <- function(UnionLoopFile, AllLoopList, AllPeakFileList, ChrList_NameNum, UseRawCC, AllRepLabels, CategoryList) {

	for (chr_idx in (1:length(ChrList_NameNum))) {
		chrName <- ChrList_NameNum[chr_idx]
		
		# extract the merged loops with respect to current chromosome
		UnionLoopTempFile1 <- paste0(dirname(UnionLoopFile), '/UnionLoopTempFile.bed')
		ExtractChrData(UnionLoopFile, chrName, UnionLoopTempFile1, header=FALSE)
		
		# check if there is no loop for the current chromosome
		nreadCurr <- GetNumLines(UnionLoopTempFile1)
		if (nreadCurr == 0) {
			next
		}

		# read the loops for the current chromosome
		MergedIntTempData <- read.table(UnionLoopTempFile1, header=F, sep="\t", stringsAsFactors=F)

		# list of vectors for storing the feature values of individual loops
		# with respect to current chromosome
		ContactCount_Categ <- list()
		QVal_Categ <- list()
		Bias1_Categ <- list()
		Bias2_Categ <- list()

		# process individual files of FitHiChIP significance 
		# for all loops
		for (i in (1:length(AllLoopList))) {

			# allocate the following vectors
			ccvec <- rep(0, nrow(MergedIntTempData))
			qvec <- rep(1, nrow(MergedIntTempData))
			bias1vec <- rep(0, nrow(MergedIntTempData))
			bias2vec <- rep(0, nrow(MergedIntTempData))

			# input file containing FitHiChIP significance for all loops
			# i.e. both significant and insignificant loops
			inpfile <- AllLoopList[i]

			# extract the input data for the current chromosome
			InpTempFile <- paste0(dirname(UnionLoopFile), '/InpTempFile.bed')
			ExtractChrData(inpfile, chrName, InpTempFile, header=TRUE)

			# check the number of loops
			nreadInp <- GetNumLines(InpTempFile)
			if (nreadInp > 0) {
				# read the input loops for the current chromosome
				InpTempData <- read.table(InpTempFile, header=F, sep="\t", stringsAsFactors=F)
				cat(sprintf("\n ***** Computing overlap of merged file with the input file: %s for the chromosome : %s ***** \n", i, chrName))

				# find the overlap between the merged set of loops
				# and the input loops
				CurrOv <- OverlapLoop(MergedIntTempData, InpTempData)

				# assign the contact count vector
				if (opt$UseRawCC == 1) {
					# here use the absolute contact count
					ccvec[CurrOv$A_AND_B] <- InpTempData[CurrOv$B_AND_A, 7]	
				} else {
					# here use the ratio of observed contact count 
					# and the expected contact count
					# and multiply with the observed contact count
					# and map to the nearest integer
					ccvec[CurrOv$A_AND_B] <- round(InpTempData[CurrOv$B_AND_A, 7] * (InpTempData[CurrOv$B_AND_A, 7] / InpTempData[CurrOv$B_AND_A, 21]))
				}

				# assign the q-value vector
				qvec[CurrOv$A_AND_B] <- InpTempData[CurrOv$B_AND_A, 25]
				
				# assign the bias 1 vector
				bias1vec[CurrOv$A_AND_B] <- InpTempData[CurrOv$B_AND_A, 10]

				# assign the bias 2 vector
				bias2vec[CurrOv$A_AND_B] <- InpTempData[CurrOv$B_AND_A, 16]
			
			}	# end number of reads condition

			# assign the contact count vector to the list of CC
			ContactCount_Categ[[i]] <- ccvec
			QVal_Categ[[i]] <- qvec
			Bias1_Categ[[i]] <- bias1vec
			Bias2_Categ[[i]] <- bias2vec
			cat(sprintf("\n ***** Finished Computing overlap of merged file with the input file: %s for the chromosome : %s ***** \n", i, chrName))

		}	# end loop FitHiChIP significance files

		# now adjust the data frame of the merged loops by appending the feature vectors
		# note that this is for a single chromosome
		for (i in (1:length(AllRepLabels))) {
			MergedIntTempData <- cbind.data.frame(MergedIntTempData, ContactCount_Categ[[i]], QVal_Categ[[i]], Bias1_Categ[[i]], Bias2_Categ[[i]])
		}

		# now append the contact count and q-value information
		# for the current chromosome and for all the input data
		# in the final data frame
		if (chr_idx == 1) {
			Merged_IntData <- MergedIntTempData
		} else {
			Merged_IntData <- rbind.data.frame(Merged_IntData, MergedIntTempData)
		}
	
	} 	# end chromosome index loop

	#==========================
	# now assign the peak strength information
	# and also the "isPeak" information
	# i.e. whether the interacting segment is peak or not

	for (i in (1:length(AllPeakFileList))) {
		currpeakfile <- AllPeakFileList[i]
		cat(sprintf("\n Merging peak strength: processing the peak file number :  name :  %s ", i, currpeakfile))

		# sourya
		# compute according to the number of features used for each cell type
		curr_col_end <- 6 + (i-1) * 4
		cat(sprintf("\n curr_col_end : %s ", curr_col_end))
		
		if (file_ext(currpeakfile) == "gz") {
			PeakData1 <- read.table(gzfile(currpeakfile), header=F, sep="\t", stringsAsFactors=F)
		} else {
			PeakData1 <- read.table(currpeakfile, header=F, sep="\t", stringsAsFactors=F)
		}
		PeakData1 <- PeakData1[, c(1,2,3,9)]

		# extract the loops
		temp_IntData <- Merged_IntData[,1:6]

		# form a vector of size = number of loops
		# corresponds to peak file of one specific cell type
		# first column: 0/1 - overlap of peak with the first interacting segment
		# second column: 0/1 - overlap of peak with the second interacting segment
		# third column: peak strength (overlapping peak) for the first interacting segment
		# fourth column: peak strength (overlapping peak) for the second interacting segment
		Peak1_Mat <- matrix(0, nrow=nrow(temp_IntData), ncol=4)

		# get overlap of first peak file with the first part of the interactions
		Peak1_Mat <- FillPeakVector(Peak1_Mat, PeakData1, temp_IntData[,1:3], c(1,3))

		# get overlap of first peak file with the second part of the interactions
		Peak1_Mat <- FillPeakVector(Peak1_Mat, PeakData1, temp_IntData[,4:6], c(2,4))

		# append the peak strength information for the current data
		OutDF <- cbind.data.frame(Merged_IntData[,1:curr_col_end], Peak1_Mat[,1:4], Merged_IntData[, (curr_col_end+1):ncol(Merged_IntData)])	

		# re-assign the data frame
		Merged_IntData <- OutDF

	}	# end peak file processing loop
	
	# #==========================
	# # now assign the ChIP-seq coverage information
	# # for the specified interacting bins
	# # and for both of the interacting categories
	# for (i in (1:length(AllChIPCovFileList))) {
	# 	currcovfile <- AllChIPCovFileList[i]
	# 	cat(sprintf("\n Merging ChIP-seq cooverage information: processing the coverage file number :  name :  %s ", i, currcovfile))	

	# 	# sourya - get the column number
	# 	# after which all of these information is to be entered
	# 	curr_col_end <- 6 + length(AllPeakFileList) * 4 + (i-1) * 2
	# 	cat(sprintf("\n curr_col_end : %s ", curr_col_end))

	# 	CovData1 <- read.table(currcovfile, header=F, sep="\t", stringsAsFactors=F)

	# 	ov1 <- Overlap1D(Merged_IntData[,1:3], CovData1[,1:3], boundary=1, offset=0, uniqov=FALSE)
	# 	ov2 <- Overlap1D(Merged_IntData[,4:6], CovData1[,1:3], boundary=1, offset=0, uniqov=FALSE)
	# 	coverage_Seg1 <- CovData1[ov1$B_AND_A, 4]
	# 	coverage_Seg2 <- CovData1[ov2$B_AND_A, 4]

	# 	# append the ChIP-seq coverage information for the current data
	# 	OutDF <- cbind.data.frame(Merged_IntData[,1:curr_col_end], coverage_Seg1, coverage_Seg2, Merged_IntData[, (curr_col_end+1):ncol(Merged_IntData)])

	# 	# re-assign the data frame
	# 	Merged_IntData <- OutDF

	# }	# end coverage file processing loop

	# #==========================
	# # also check the input coverage bias files
	# # for individual replicates of individual categories
	# # and insert the coverage bias values for the individual interacting segments
	# for (i in (1:length(BiasFileList))) {
	# 	currbiasfile <- BiasFileList[i]
	# 	cat(sprintf("\n Merging bias values: processing the bias file number : %s  name : %s ", i, currbiasfile))

	# 	BiasData <- read.table(currbiasfile, header=T, sep="\t", stringsAsFactors=F)
	# 	ov1 <- Overlap1D(Merged_IntData[,1:3], BiasData[,1:3], boundary=1, offset=0, uniqov=FALSE)
	# 	ov2 <- Overlap1D(Merged_IntData[,4:6], BiasData[,1:3], boundary=1, offset=0, uniqov=FALSE)
	# 	biacol1 <- 6 + 4 * length(AllPeakFileList) + 2 * length(AllChIPCovFileList) + (i-1) * 4 + 3
	# 	biacol2 <- biacol1 + 1
	# 	Merged_IntData[ov1$A_AND_B, biacol1] <- BiasData[ov1$B_AND_A, 6]
	# 	Merged_IntData[ov2$A_AND_B, biacol2] <- BiasData[ov2$B_AND_A, 6]

	# }	# end bias file processing loop
	# #==========================

	# assign the column names 
	appendnamevec <- c()
	for (i in (1:length(AllRepLabels))) {
		for (v in c('_CC', '_QVal', '_Bias1', '_Bias2')) {
			appendnamevec <- c(appendnamevec, paste0(AllRepLabels[i], v))
		}
	}

	# column names of the final list of interactions
	
	# namesvec <- c("chr1", "start1", "end1", "chr2", "start2", "end2", paste0(CategoryList[1], c('_isPeak1', '_isPeak2', '_PS1', '_PS2')), paste0(CategoryList[2], c('_isPeak1', '_isPeak2', '_PS1', '_PS2')), paste0(CategoryList[1], c('_ChIPCov1', '_ChIPCov2')), paste0(CategoryList[2], c('_ChIPCov1', '_ChIPCov2')), appendnamevec)
	namesvec <- c("chr1", "start1", "end1", "chr2", "start2", "end2", paste0(CategoryList[1], c('_isPeak1', '_isPeak2', '_PS1', '_PS2')), paste0(CategoryList[2], c('_isPeak1', '_isPeak2', '_PS1', '_PS2')), appendnamevec)

	colnames(Merged_IntData) <- namesvec
	write.table(Merged_IntData, UnionLoopFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)	

}	# end function

#======================================================
# this function applies EdgeR to the input ChIP-seq coverage data
# where for individual replicates, there is a input bed formatted data
# first three columns are the chromosome interval
# and the fourth column is the ChIP-seq coverage
# parameters:
# InpFileList: list of input files storing the ChIP seq segments and coverage
# HeaderInpList: header information of the input files
# ContactColList: column numbers for individual input files storing the count data
#======================================================
ApplyEdgeR_ChIPCoverage <- function(InpFileList, HeaderInpList, ContactColList, CategoryList, ReplicaCount, ResFile, bcv=0.4) {

	tempfile <- paste0(dirname(ResFile), '/temp_Peak_Union.bed')
	tempfile1 <- paste0(dirname(ResFile), '/temp_Peak_Union1.bed')

	for (i in (1:length(InpFileList))) {
		if (HeaderInpList[i] == 1) {
			if (i == 1) {
				if (file_ext(InpFileList[i]) == "gz") {
					system(paste("zcat", InpFileList[i], "| cut -f1-3 | awk \'(NR>1)\' - >", tempfile))
				} else {
					system(paste("cat", InpFileList[i], "| cut -f1-3 | awk \'(NR>1)\' - >", tempfile))
				}
			} else {
				if (file_ext(InpFileList[i]) == "gz") {
					system(paste("zcat", InpFileList[i], "| cut -f1-3 | awk \'(NR>1)\' - >>", tempfile))
				} else {
					system(paste("cat", InpFileList[i], "| cut -f1-3 | awk \'(NR>1)\' - >>", tempfile))
				}
			}
		} else {
			if (i == 1) {
				if (file_ext(InpFileList[i]) == "gz") {
					system(paste("zcat", InpFileList[i], "| cut -f1-3 >", tempfile))
				} else {
					system(paste("cat", InpFileList[i], "| cut -f1-3 >", tempfile))
				}
			} else {
				if (file_ext(InpFileList[i]) == "gz") {
					system(paste("zcat", InpFileList[i], "| cut -f1-3 >>", tempfile))
				} else {
					system(paste("cat", InpFileList[i], "| cut -f1-3 >>", tempfile))
				}
			}
		}
	}

	# sort the file, remove duplicate loops
	system(paste("sort -k1,1 -k2,2n", tempfile, "| awk -F\"\t\" \'!seen[$1, $2, $3]++\' - >", tempfile1))

	OutDF <- read.table(tempfile1, header=T, sep="\t", stringsAsFactors=F)

	# construct the coverage values for all the input data
	for (i in (1:length(InpFileList))) {
		if (HeaderInpList[i] == 1) {
			if (file_ext(InpFileList[i]) == "gz") {
				x <- read.table(gzfile(InpFileList[i]), header=T, sep="\t", stringsAsFactors=F)
			} else {
				x <- read.table(InpFileList[i], header=T, sep="\t", stringsAsFactors=F)
			}
		} else {
			if (file_ext(InpFileList[i]) == "gz") {
				x <- read.table(gzfile(InpFileList[i]), header=F, sep="\t", stringsAsFactors=F)
			} else {
				x <- read.table(InpFileList[i], header=F, sep="\t", stringsAsFactors=F)
			}
		}
		y <- rep(0, nrow(OutDF))
		ov <- Overlap1D(OutDF[,1:3], x[,1:3], boundary=1, offset=0, uniqov=FALSE)
		y[ov$A_AND_B] <- x[ov$B_AND_A, ContactColList[i]]
		OutDF <- cbind.data.frame(OutDF, y)
	}

	# construct the final data frame to be used as a count matrix
	CountData <- cbind.data.frame(seq(1,nrow(OutDF)), OutDF[,1], (OutDF[,2] + OutDF[,3])/2, OutDF[,4:ncol(OutDF)])
	colnames(CountData) <- c("Idx", "chr", "mid", paste0(CategoryList[1], "_", seq(1,ReplicaCount[1])), paste0(CategoryList[2], "_", seq(1,ReplicaCount[2])))

	CountDataFile <- paste0(dirname(ResFile), '/Dumped_CountData.bed')
	write.table(CountData, CountDataFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

	# read the raw count data
	# and assign the replicate information
	# this vector stores the category and replicate distribution
	# to be used for edgeR execution
	GroupDistrVec <- c(rep(CategoryList[1], ReplicaCount[1]), rep(CategoryList[2], ReplicaCount[2]))
	cat(sprintf("\n ===>>> GroupDistrVec : %s ", paste0(GroupDistrVec, sep="\t")))

	# create the EdgeR count data structure
	y <- DGEList(counts=CountData[,4:ncol(CountData)], group=GroupDistrVec)

	if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
		y <- estimateDisp(y)
		et <- exactTest(y, dispersion="trended", pair=c(CategoryList[1], CategoryList[2]))
	} else {
		et <- exactTest(y, dispersion=bcv^2)
	}

	# estimate the q-value
	Qval <- p.adjust(et$table$PValue, method = "BH")

	# create a data frame containing et
	# and also the adjusted FDR
	EdgeRRes <- data.frame(CountData, logFC  = et$table$logFC, logCPM = et$table$logCPM, PValue = et$table$PValue, QValue=Qval)

	# the above structure does not contain the chromosome full interval
	# but rather the midpoints
	# create a data structure which will contain the original UnionLoopData
	# and the EdgeR results
	EdgeRRes_NEW <- cbind.data.frame(OutDF[,1:3], EdgeRRes[, ((ncol(EdgeRRes) - 3):ncol(EdgeRRes))])
	colnames(EdgeRRes_NEW) <- c("chr", "start", "end", colnames(EdgeRRes)[((ncol(EdgeRRes) - 3):ncol(EdgeRRes))])
	write.table(EdgeRRes_NEW, ResFile, row.names=F, col.names=T, quote=F, sep="\t")	

	if (file.exists(tempfile)) {
		system(paste("rm", tempfile))
	}

	if (file.exists(tempfile1)) {
		system(paste("rm", tempfile1))
	}	

}	# end function

#======================================================
# this function applies EdgeR to the input loops
#======================================================
ApplyEdgeR <- function(UnionLoopFile, MainDir, CountCol, CategoryList, ReplicaCount, ResFile, bcv=0.4) {

	# read the input set of loops
	# with the feature values as well
	UnionLoopData <- read.table(UnionLoopFile, header=T, sep="\t", stringsAsFactors=F)

	# first write a text file
	# which will contain the interacting segments 
	# along with the index of those interactions
	IntervalTextFile <- paste0(MainDir, '/Interacting_segments.bed')
	system(paste("awk \'{if (NR>1) {print $1\"\t\"$2\"\t\"$3; print $4\"\t\"$5\"\t\"$6}}\'",  UnionLoopFile, "| sort -k1,1 -k2,2n | uniq | awk \'{print NR\"\t\"$0}\' - >", IntervalTextFile))

	IntervalData <- read.table(IntervalTextFile, header=F, sep="\t", stringsAsFactors=F)
	CountDataFile <- paste0(MainDir, '/Dumped_CountData.bed')

	# merge the first three columns of the UnionLoopData
	# with the columns 2-4 of the IntervalData
	# to map the index of corresponding chromosomal segment
	colnames(IntervalData) <- c("Idx1", "chr1", "start1", "end1")
	df1 <- inner_join(IntervalData, UnionLoopData, by.x=colnames(IntervalData)[2:4], by.y=colnames(UnionLoopData)[1:3])

	# merge the columns 4-6 of the UnionLoopData
	# with the columns 2-4 of the IntervalData
	# to map the index of corresponding chromosomal segment
	colnames(IntervalData) <- c("Idx2", "chr2", "start2", "end2")
	df2 <- inner_join(IntervalData, UnionLoopData, by.x=colnames(IntervalData)[2:4], by.y=colnames(UnionLoopData)[4:6])

	# only consider individual replicates
	# for creating the count matrix
	# consider the contact count values
	# also, edgeR requires the indices and data 
	# for individual chromosomal segments
	CountData <- cbind.data.frame(df1[,1], df2[,1], UnionLoopData[,1], (UnionLoopData[,2] + UnionLoopData[,3])/2, UnionLoopData[,4], (UnionLoopData[,5] + UnionLoopData[,6])/2, UnionLoopData[, CountCol])
	colnames(CountData) <- c(colnames(df1)[1], colnames(df2)[1], colnames(UnionLoopData)[1], "mid1", colnames(UnionLoopData)[4], "mid2", colnames(UnionLoopData)[CountCol])
	cn <- colnames(CountData)

	# dump the count data
	write.table(CountData, CountDataFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

	# read the raw count data
	# and assign the replicate information

	# this vector stores the category and replicate distribution
	# to be used for edgeR execution
	GroupDistrVec <- c(rep(CategoryList[1], ReplicaCount[1]), rep(CategoryList[2], ReplicaCount[2]))
	cat(sprintf("\n ===>>> GroupDistrVec : %s ", paste0(GroupDistrVec, sep="\t")))

	# create the EdgeR count data structure
	y <- DGEList(counts=CountData[,7:ncol(CountData)], group=GroupDistrVec)

	if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
		y <- estimateDisp(y)
		et <- exactTest(y, dispersion="trended", pair=c(CategoryList[1], CategoryList[2]))
	} else {
		et <- exactTest(y, dispersion=bcv^2)
	}

	# estimate the q-value
	Qval <- p.adjust(et$table$PValue, method = "BH")

	# create a data frame containing et
	# and also the adjusted FDR
	EdgeRRes <- data.frame(CountData, logFC  = et$table$logFC, logCPM = et$table$logCPM, PValue = et$table$PValue, QValue=Qval)
	# write.table(EdgeRRes, ResFile, row.names=F, col.names=T, quote=F, sep="\t")

	# the above structure does not contain the chromosome full interval
	# but rather the midpoints
	# create a data structure which will contain the original UnionLoopData
	# and the EdgeR results
	EdgeRRes_NEW <- cbind.data.frame(UnionLoopData, EdgeRRes[, ((ncol(EdgeRRes) - 3):ncol(EdgeRRes))])
	colnames(EdgeRRes_NEW) <- c(colnames(UnionLoopData), colnames(EdgeRRes)[((ncol(EdgeRRes) - 3):ncol(EdgeRRes))])
	write.table(EdgeRRes_NEW, ResFile, row.names=F, col.names=T, quote=F, sep="\t")

}	# end EdgeR function

# #==============================
# # this plot combines the replicate specific boxplots 
# # produced in the function "WritePromCentricLoops"
# # to generate a box plot of the given categories
# # for P-P interactions and P-E interactions
# #==============================
# CombineExclReplicateLoopsBoxPlot <- function(P_P_CurrDir, P_E_CurrDir, CategoryList, ReplicaCountVecCat1, ReplicaCountVecCat2) {

# 	if ((dir.exists(P_P_CurrDir) == FALSE) | (dir.exists(P_E_CurrDir) == FALSE)) {
# 		return()
# 	}

# 	# data structures for plotting the boxplot (consolidated for P-P and P-E loops)
# 	Sets <- c()
# 	Categories <- c()
# 	DataVec <- c()
	
# 	# stores the statistical significance of gene expression ratio values
# 	PVal_Sig_List <- c(1,1,1,1)
# 	# stores the significance symbols (< 0.01: **, < 0.05: *)
# 	PVal_Sig_Symbol_List <- c(rep(paste0(''), 4))
# 	# label vector of p values
# 	PVal_Labels <- c()

# 	# labels storing count of 
# 	# different categories
# 	Cnt_Labels <- c()

# 	cat(sprintf("\n P_P_CurrDir : %s ", P_P_CurrDir))

# 	# also writing some summary statistics
# 	fp_out <- file(paste0(P_P_CurrDir, '/CONSOLIDATED_P_P_and_P_E_Significance_Tests.log'), "w")
# 	outtext <- paste0("\n Category 1: ", CategoryList[1], "\n Category 2: ", CategoryList[2])
# 	writeLines(outtext, con=fp_out, sep="\n")
	
# 	# P-P - exclusive to category 1
# 	# gene expr ratio values
# 	gene_expr_ratiodata <- c()
# 	for (i in (1:length(ReplicaCountVecCat1))) {
# 		cat1_rc <- ReplicaCountVecCat1[i]
# 		cat2_rc <- 0
# 		currfile <- paste0(P_P_CurrDir, '/Loops_', CategoryList[1], '_', cat1_rc, '_', CategoryList[2], '_', cat2_rc, '.bed')
# 		if (file.exists(currfile) == TRUE) {
# 			currdata <- read.table(currfile, header=T, sep="\t", stringsAsFactors=F)
# 			# create ratio of gene expression
# 			currratiodata <- log2((currdata[, (ncol(currdata) - 1)] + 1) / (currdata[, ncol(currdata)] + 1))
# 			# append the ratio values in the designated vector
# 			gene_expr_ratiodata <- c(gene_expr_ratiodata, currratiodata)			
# 		}
# 	}

# 	if (length(unique(gene_expr_ratiodata)) > 1) {
		
# 		# significance of the gene expression ratio values
# 		# for P-P exclusive loops in category 1
# 		GE_One_t_test_df <- t.test(gene_expr_ratiodata, mu=0.5, alternative = c("greater"))	
		
# 		# dump the results
# 		outtext <- paste0("\n P-P Loops exclusive to ", CategoryList[1], "\n number of entries: ", length(gene_expr_ratiodata), "\n Ratio expression one sample test p value: ", GE_One_t_test_df$p.value)
# 		writeLines(outtext, con=fp_out, sep="\n")

# 		# append in the final data structures
# 		Sets <- c(Sets, rep(paste0("Excl_", CategoryList[1]), length(gene_expr_ratiodata)))
# 		Categories <- c(Categories, rep(paste0("P_P"), length(gene_expr_ratiodata)))
# 		DataVec <- c(DataVec, gene_expr_ratiodata)

# 		# add the gene expression P value significance
# 		PVal_Sig_List[1] <- GE_One_t_test_df$p.value
# 		if (GE_One_t_test_df$p.value < 0.01) {
# 			PVal_Sig_Symbol_List[1] <- '**'
# 		} else if (GE_One_t_test_df$p.value < 0.05) {
# 			PVal_Sig_Symbol_List[1] <- '*'
# 		}
# 		PVal_Labels <- c(PVal_Labels, rep(PVal_Sig_Symbol_List[1], length(gene_expr_ratiodata)))

# 		# adjust the count of samples for this category
# 		Cnt_Labels <- c(Cnt_Labels, rep(paste0("(", length(gene_expr_ratiodata), ")"), length(gene_expr_ratiodata)))

# 	}

# 	# P-P - exclusive to category 2
# 	# gene expr ratio values
# 	gene_expr_ratiodata <- c()
# 	for (i in (1:length(ReplicaCountVecCat2))) {
# 		cat2_rc <- ReplicaCountVecCat2[i]
# 		cat1_rc <- 0
# 		currfile <- paste0(P_P_CurrDir, '/Loops_', CategoryList[1], '_', cat1_rc, '_', CategoryList[2], '_', cat2_rc, '.bed')
# 		if (file.exists(currfile) == TRUE) {
# 			currdata <- read.table(currfile, header=T, sep="\t", stringsAsFactors=F)
# 			# create ratio of gene expression
# 			currratiodata <- log2((currdata[, (ncol(currdata) - 1)] + 1) / (currdata[, ncol(currdata)] + 1))
# 			# append the ratio values in the designated vector
# 			gene_expr_ratiodata <- c(gene_expr_ratiodata, currratiodata)
# 		}
# 	}

# 	if (length(unique(gene_expr_ratiodata)) > 1) {
# 		# significance of the gene expression ratio values
# 		# for P-P exclusive loops in category 2
# 		GE_One_t_test_df <- t.test(gene_expr_ratiodata, mu=-0.5, alternative = c("less"))	
		
# 		# dump the results
# 		outtext <- paste0("\n P-P Loops exclusive to ", CategoryList[2], "\n number of entries: ", length(gene_expr_ratiodata), "\n Ratio expression one sample test p value: ", GE_One_t_test_df$p.value)
# 		writeLines(outtext, con=fp_out, sep="\n")	

# 		# append in the final data structures
# 		Sets <- c(Sets, rep(paste0("Excl_", CategoryList[2]), length(gene_expr_ratiodata)))
# 		Categories <- c(Categories, rep(paste0("P_P"), length(gene_expr_ratiodata)))
# 		DataVec <- c(DataVec, gene_expr_ratiodata)	

# 		# add the gene expression P value significance
# 		PVal_Sig_List[2] <- GE_One_t_test_df$p.value
# 		if (GE_One_t_test_df$p.value < 0.01) {
# 			PVal_Sig_Symbol_List[2] <- '**'
# 		} else if (GE_One_t_test_df$p.value < 0.05) {
# 			PVal_Sig_Symbol_List[2] <- '*'
# 		}
# 		PVal_Labels <- c(PVal_Labels, rep(PVal_Sig_Symbol_List[2], length(gene_expr_ratiodata)))

# 		# adjust the count of samples for this category
# 		Cnt_Labels <- c(Cnt_Labels, rep(paste0("(", length(gene_expr_ratiodata), ")"), length(gene_expr_ratiodata)))

# 	}


# 	# P-E - exclusive to category 1
# 	# gene expr ratio values
# 	gene_expr_ratiodata <- c()
# 	for (i in (1:length(ReplicaCountVecCat1))) {
# 		cat1_rc <- ReplicaCountVecCat1[i]
# 		cat2_rc <- 0
# 		currfile <- paste0(P_E_CurrDir, '/Loops_', CategoryList[1], '_', cat1_rc, '_', CategoryList[2], '_', cat2_rc, '.bed')
# 		if (file.exists(currfile) == TRUE) {
# 			currdata <- read.table(currfile, header=T, sep="\t", stringsAsFactors=F)
# 			# create ratio of gene expression
# 			currratiodata <- log2((currdata[, (ncol(currdata) - 1)] + 1) / (currdata[, ncol(currdata)] + 1))
# 			# append the ratio values in the designated vector
# 			gene_expr_ratiodata <- c(gene_expr_ratiodata, currratiodata)
# 		}
# 	}

# 	if (length(unique(gene_expr_ratiodata)) > 1) {
# 		# significance of the gene expression ratio values
# 		# for P-E exclusive loops in category 1
# 		GE_One_t_test_df <- t.test(gene_expr_ratiodata, mu=0.5, alternative = c("greater"))	

# 		# dump the results
# 		outtext <- paste0("\n P-E Loops exclusive to ", CategoryList[1], "\n number of entries: ", length(gene_expr_ratiodata), "\n Ratio expression one sample test p value: ", GE_One_t_test_df$p.value)
# 		writeLines(outtext, con=fp_out, sep="\n")		

# 		# append in the final data structures
# 		Sets <- c(Sets, rep(paste0("Excl_", CategoryList[1]), length(gene_expr_ratiodata)))
# 		Categories <- c(Categories, rep(paste0("P_E"), length(gene_expr_ratiodata)))
# 		DataVec <- c(DataVec, gene_expr_ratiodata)	

# 		# add the gene expression P value significance
# 		PVal_Sig_List[3] <- GE_One_t_test_df$p.value
# 		if (GE_One_t_test_df$p.value < 0.01) {
# 			PVal_Sig_Symbol_List[3] <- '**'
# 		} else if (GE_One_t_test_df$p.value < 0.05) {
# 			PVal_Sig_Symbol_List[3] <- '*'
# 		}
# 		PVal_Labels <- c(PVal_Labels, rep(PVal_Sig_Symbol_List[3], length(gene_expr_ratiodata)))

# 		# adjust the count of samples for this category
# 		Cnt_Labels <- c(Cnt_Labels, rep(paste0("(", length(gene_expr_ratiodata), ")"), length(gene_expr_ratiodata)))

# 	}


# 	# P-E - exclusive to category 2
# 	# gene expr ratio values
# 	gene_expr_ratiodata <- c()
# 	for (i in (1:length(ReplicaCountVecCat2))) {
# 		cat2_rc <- ReplicaCountVecCat2[i]
# 		cat1_rc <- 0
# 		currfile <- paste0(P_E_CurrDir, '/Loops_', CategoryList[1], '_', cat1_rc, '_', CategoryList[2], '_', cat2_rc, '.bed')
# 		if (file.exists(currfile) == TRUE) {
# 			currdata <- read.table(currfile, header=T, sep="\t", stringsAsFactors=F)
# 			# create ratio of gene expression
# 			currratiodata <- log2((currdata[, (ncol(currdata) - 1)] + 1) / (currdata[, ncol(currdata)] + 1))
# 			# append the ratio values in the designated vector
# 			gene_expr_ratiodata <- c(gene_expr_ratiodata, currratiodata)
# 		}
# 	}

# 	if (length(unique(gene_expr_ratiodata)) > 1) {
# 		# significance of the gene expression ratio values
# 		# for P-E exclusive loops in category 2
# 		GE_One_t_test_df <- t.test(gene_expr_ratiodata, mu=-0.5, alternative = c("less"))	
		
# 		# dump the results
# 		outtext <- paste0("\n P-E Loops exclusive to ", CategoryList[2], "\n number of entries: ", length(gene_expr_ratiodata), "\n Ratio expression one sample test p value: ", GE_One_t_test_df$p.value)
# 		writeLines(outtext, con=fp_out, sep="\n")	

# 		# append in the final data structures
# 		Sets <- c(Sets, rep(paste0("Excl_", CategoryList[2]), length(gene_expr_ratiodata)))
# 		Categories <- c(Categories, rep(paste0("P_E"), length(gene_expr_ratiodata)))
# 		DataVec <- c(DataVec, gene_expr_ratiodata)

# 		# add the gene expression P value significance
# 		PVal_Sig_List[4] <- GE_One_t_test_df$p.value
# 		if (GE_One_t_test_df$p.value < 0.01) {
# 			PVal_Sig_Symbol_List[4] <- '**'
# 		} else if (GE_One_t_test_df$p.value < 0.05) {
# 			PVal_Sig_Symbol_List[4] <- '*'
# 		}
# 		PVal_Labels <- c(PVal_Labels, rep(PVal_Sig_Symbol_List[4], length(gene_expr_ratiodata)))

# 		# adjust the count of samples for this category
# 		Cnt_Labels <- c(Cnt_Labels, rep(paste0("(", length(gene_expr_ratiodata), ")"), length(gene_expr_ratiodata)))

# 	}

# 	# close the text file
# 	close(fp_out)

# 	# plot the ratio of gene expressions
# 	if (length(DataVec) > 0) {

# 		plotfile2 <- paste0(P_P_CurrDir, '/CONSOLIDATED_P_P_and_P_E_ratio_gene_expr_boxplot.pdf')

# 		CompleteData <- data.frame(y=DataVec, x=Categories, fill=Sets)
# 		cutoff <- data.frame(yintercept=0, cutoff=factor(0))
# 		ylim1 = boxplot.stats(CompleteData$y)$stats

# 		currplot <- ggplot(CompleteData, aes(x=Categories, y=DataVec, fill=Sets)) + geom_boxplot(width=0.2, outlier.shape = NA) + scale_y_continuous(name = paste0("log2(", CategoryList[1], " / ", CategoryList[2], ") Expression"), limits = c(-5, 5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16), legend.position="bottom") + geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff, show.legend=FALSE) + geom_text(aes(label = PVal_Labels, y=(ylim1[4] + 1), group = Sets), position = position_dodge(0.9), size=8) + geom_text(aes(label = Cnt_Labels, y=(ylim1[1] - 1), group = Sets), position = position_dodge(0.9), size=3.5)
		
# 		# currplot + ggtitle(paste0("Ratio gene expression - ", CategoryList[1], "/", CategoryList[2]))
# 		ggsave(plotfile2, width=5, height=5)

# 	}

# }	# end function


#==============================
# this function writes the promoter centric loops union, and exclusive to either categories
# also it plots the gene expressions for the differential loops
# exclusive to one of the two categories
#==============================
WritePromCentricLoops <- function(CurrDF, LoopFile, CategoryList, ReplicaCountVecCat1, ReplicaCountVecCat2, ContactCountCol1, ContactCountCol2, QvalCol1, QvalCol2, GeneExprDataInp=TRUE) {

	GroupLabels_1 <- c()
	DataLabels_1 <- c()
	DataVec_1 <- c()
	pval_1 <- c()
	pval_2 <- c()

	outRatioDF_Filled <- FALSE

	# also writing some summary statistics
	if (GeneExprDataInp == TRUE) {
		fp_out <- file(paste0(dirname(LoopFile), '/Significance_Tests.log'), "w")
		outtext <- paste0("\n Category 1: ", CategoryList[1], "\n Category 2: ", CategoryList[2])
		writeLines(outtext, con=fp_out, sep="\n")		
	}

	# write the union of loops for both categories
	write.table(CurrDF, LoopFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

	# form a master list of vectors
	# containing the replicate count for both categories
	RepCountCat1List <- c(ReplicaCountVecCat1, rep(0, length(ReplicaCountVecCat2)))
	RepCountCat2List <- c(rep(0, length(ReplicaCountVecCat1)), ReplicaCountVecCat2)

	# now loop through the category specific replicate count information
	for (i in (1:length(RepCountCat1List))) {
		cat1_rep_cnt <- RepCountCat1List[i]
		cat2_rep_cnt <- RepCountCat2List[i]
		targetfile <- paste0(sub('\\.bed$', '', LoopFile), '_', CategoryList[1], '_', cat1_rep_cnt, '_', CategoryList[2], '_', cat2_rep_cnt, '.bed')
		# cat(sprintf("\n target file: %s ", targetfile))
		
		# get the loop indices satisfying the 
		# significant replicate count for these two categories
		Cat1_Sig_Idx <- GetSpecificSigIdx(CurrDF, c(ContactCountCol1), c(QvalCol1), cat1_rep_cnt)
		Cat2_Sig_Idx <- GetSpecificSigIdx(CurrDF, c(ContactCountCol2), c(QvalCol2), cat2_rep_cnt)
		Idx_Sig <- intersect(Cat1_Sig_Idx, Cat2_Sig_Idx)
		
		if (length(Idx_Sig) > 0) {

			# first dump the loops which are significant in the mentioned
			# replciate counts of either categories
			OutData <- CurrDF[Idx_Sig, ]
			write.table(OutData, targetfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

			if (GeneExprDataInp == TRUE) {

				# check if t-test is possible
				tempDF <- cbind.data.frame(OutData[, (ncol(OutData) - 1)], OutData[, ncol(OutData)])
				tempDF1 <- unique(tempDF)
				
				if (nrow(tempDF1) > 1) {
					# check the raw gene expression values
					# for both categories
					# employ paired t test
					rawdata1 <- OutData[, (ncol(OutData) - 1)]
					rawdata2 <- OutData[, ncol(OutData)]
					GE_t_test_df <- t.test(rawdata1, rawdata2, paired=TRUE, conf.level=0.95)
					pval_1 <- c(pval_1, GE_t_test_df$p.value)
					
					# perform one sample t test for checking the log ratio cat1 / cat2
					# for exclusive interactions with respect to the first category
					# alternative hypothesis should be "greater"
					# for exclusive interactions with respect to the second category
					# alternative hypothesis should be "less"
					ratiodata <- log2((rawdata1 + 1) / (rawdata2 + 1))
					if ((cat1_rep_cnt > 0) & (cat2_rep_cnt == 0)) {
						GE_One_t_test_df <- t.test(ratiodata, mu=0.5, alternative = c("greater"))	
					} else if ((cat1_rep_cnt == 0) & (cat2_rep_cnt > 0)) {
						GE_One_t_test_df <- t.test(ratiodata, mu=-0.5, alternative = c("less"))
					} else {
						GE_One_t_test_df <- t.test(ratiodata, mu=0, alternative = c("two.sided"))
					}
					
					pval_2 <- c(pval_2, GE_One_t_test_df$p.value)

					# data frame for first plot
					GroupLabels_1 <- c(GroupLabels_1, rep(CategoryList[1], length(Idx_Sig)), rep(CategoryList[2], length(Idx_Sig)))
					DataLabels_1 <- c(DataLabels_1, rep(paste0(CategoryList[1], "_", cat1_rep_cnt, "\n", CategoryList[2], '_', cat2_rep_cnt, "\n", GE_t_test_df$p.value), (2 * length(Idx_Sig))))
					DataVec_1 <- c(DataVec_1, rawdata1, rawdata2)

					# data frame for the ratio plot
					currRatioDF <- data.frame(group=paste0(CategoryList[1], "_", cat1_rep_cnt, "\n", CategoryList[2], '_', cat2_rep_cnt, "\n(", length(Idx_Sig), ")\n", GE_One_t_test_df$p.value), value=ratiodata)
					if (outRatioDF_Filled == FALSE) {
						outRatioDF <- currRatioDF
						outRatioDF_Filled <- TRUE
					} else {
						outRatioDF <- rbind.data.frame(outRatioDF, currRatioDF)
					}

					# dump the results
					outtext <- paste0("\n output loop file: ", targetfile, "\n category 1 - no of significant replicates : ", cat1_rep_cnt, "\n category 2 - no of significant replicates : ", cat2_rep_cnt, "\n number of genes (to bins): ", length(Idx_Sig), "\n raw expression t test p val: ", GE_t_test_df$p.value, "\n Ratio expression wilcoxon one sample signed test p value: ", GE_One_t_test_df$p.value)
					writeLines(outtext, con=fp_out, sep="\n")

				}	# end nrow condition if

			}	# end gene expression condition

		}	# end idx_sig condition
	}	# end replicata condition

	if (GeneExprDataInp == TRUE) {

		# now dump the complete set of loops (without any exclusive replicates)
		OutData <- CurrDF
		# check if t-test is possible
		tempDF <- cbind.data.frame(OutData[, (ncol(OutData) - 1)], OutData[, ncol(OutData)])
		tempDF1 <- unique(tempDF)
		
		if (nrow(tempDF1) > 1) {
			# check the raw gene expression values
			# for both categories
			# employ paired t test
			rawdata1 <- OutData[, (ncol(OutData) - 1)]
			rawdata2 <- OutData[, ncol(OutData)]
			GE_t_test_df <- t.test(rawdata1, rawdata2, paired=TRUE, conf.level=0.95)
			pval_1 <- c(pval_1, GE_t_test_df$p.value)
			
			# perform one sample t test for checking the log ratio cat1 / cat2
			ratiodata <- log2((rawdata1 + 1) / (rawdata2 + 1))
			GE_One_t_test_df <- t.test(ratiodata, mu=0, alternative = c("two.sided"))
			
			pval_2 <- c(pval_2, GE_One_t_test_df$p.value)

			# data frame for first plot
			GroupLabels_1 <- c(GroupLabels_1, rep(CategoryList[1], nrow(OutData)), rep(CategoryList[2], nrow(OutData)))
			DataLabels_1 <- c(DataLabels_1, rep(paste0("ALL samples\n", GE_t_test_df$p.value), (2 * nrow(OutData))))
			DataVec_1 <- c(DataVec_1, rawdata1, rawdata2)

			# data frame for the ratio plot
			currRatioDF <- data.frame(group=paste0("ALL samples\n(", nrow(OutData), ")\n", GE_One_t_test_df$p.value), value=ratiodata)
			if (outRatioDF_Filled == FALSE) {
				outRatioDF <- currRatioDF
				outRatioDF_Filled <- TRUE
			} else {
				outRatioDF <- rbind.data.frame(outRatioDF, currRatioDF)
			}

			# dump the results
			outtext <- paste0("Considering ALL SAMPLES \n number of genes (to bins): ", nrow(OutData), "\n raw expression t test p val: ", GE_t_test_df$p.value, "\n Ratio expression wilcoxon one sample signed test p value: ", GE_One_t_test_df$p.value)
			writeLines(outtext, con=fp_out, sep="\n")
		}

		# the final category is to first get the loops exclusively to one category
		# (one or more replicates) from both sides
		# and then get equal number of samples from them
		# to plot the ratio of gene expression
		# loop indices exclusive to first category (one or more replicates)
		Cat1_Excl_Cat2_Zero_Idx <- GetSpecificSigIdx(CurrDF, c(ContactCountCol2), c(QvalCol2), 0)
		# loop indices exclusive to second category (one or more replicates)
		Cat2_Excl_Cat1_Zero_Idx <- GetSpecificSigIdx(CurrDF, c(ContactCountCol1), c(QvalCol1), 0)
		# number of sample to randomly pick (75%)
		nsample <- as.integer((min(length(Cat1_Excl_Cat2_Zero_Idx), length(Cat2_Excl_Cat1_Zero_Idx)) * 3) / 4)

		cat(sprintf("\n length Cat1_Excl_Cat2_Zero_Idx : %s  \n length Cat2_Excl_Cat1_Zero_Idx : %s \n Sampling equal plus exclusive : number of samples for both categories: %s ", length(Cat1_Excl_Cat2_Zero_Idx), length(Cat2_Excl_Cat1_Zero_Idx), nsample))

		Cat1_Excl_Sampled_Idx <- sample(Cat1_Excl_Cat2_Zero_Idx, nsample, replace=FALSE)
		Cat2_Excl_Sampled_Idx <- sample(Cat2_Excl_Cat1_Zero_Idx, nsample, replace=FALSE)
		
		OutData <- CurrDF[union(Cat1_Excl_Sampled_Idx, Cat2_Excl_Sampled_Idx), ]
		tempDF <- cbind.data.frame(OutData[, (ncol(OutData) - 1)], OutData[, ncol(OutData)])
		tempDF1 <- unique(tempDF)

		if (nrow(tempDF1) > 1) {
			# check the raw gene expression values
			# for both categories
			# employ paired t test
			rawdata1 <- OutData[, (ncol(OutData) - 1)]
			rawdata2 <- OutData[, ncol(OutData)]
			GE_t_test_df <- t.test(rawdata1, rawdata2, paired=TRUE, conf.level=0.95)
			pval_1 <- c(pval_1, GE_t_test_df$p.value)
			
			# perform one sample t test for checking the log ratio cat1 / cat2
			ratiodata <- log2((rawdata1 + 1) / (rawdata2 + 1))
			GE_One_t_test_df <- t.test(ratiodata, mu=0, alternative = c("two.sided"))
			
			pval_2 <- c(pval_2, GE_One_t_test_df$p.value)

			# data frame for first plot
			GroupLabels_1 <- c(GroupLabels_1, rep(CategoryList[1], nrow(OutData)), rep(CategoryList[2], nrow(OutData)))
			DataLabels_1 <- c(DataLabels_1, rep(paste0("Equal +\nExclusive\nsamples\n", GE_t_test_df$p.value), (2 * nrow(OutData))))
			DataVec_1 <- c(DataVec_1, rawdata1, rawdata2)

			# data frame for the ratio plot
			currRatioDF <- data.frame(group=paste0("Equal +\nExclusive\nsamples\n(", nrow(OutData), ")\n", GE_One_t_test_df$p.value), value=ratiodata)
			if (outRatioDF_Filled == FALSE) {
				outRatioDF <- currRatioDF
				outRatioDF_Filled <- TRUE
			} else {
				outRatioDF <- rbind.data.frame(outRatioDF, currRatioDF)
			}

			# dump the results
			outtext <- paste0("Considering Equal + Exclusive samples \n number of genes (to bins): ", nrow(OutData), "\n raw expression t test p val: ", GE_t_test_df$p.value, "\n Ratio expression wilcoxon one sample signed test p value: ", GE_One_t_test_df$p.value)
			writeLines(outtext, con=fp_out, sep="\n")
		}

		if (outRatioDF_Filled == TRUE) {

			# plot the absolute gene expressions
			plotfile1 <- paste0(dirname(LoopFile), '/absolute_gene_expr_boxplot.pdf')		
			CompleteData <- data.frame(y=DataVec_1, x=DataLabels_1, fill=GroupLabels_1)
			ylim1 = boxplot.stats(CompleteData$y)$stats
			ggplot(CompleteData, aes(x=DataLabels_1, y=DataVec_1, fill=GroupLabels_1)) + geom_boxplot(width=0.2, outlier.shape = NA) + scale_y_continuous(name = "Gene Expression", limits = c(ylim1[1], ylim1[4] * 1.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16), axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom") #+ guides(fill=FALSE) + geom_signif(comparisons = list(c(CategoryList[1], CategoryList[2])), map_signif_level=FALSE)
			ggsave(plotfile1, width=10, height=7)

			# plot the ratio of gene expressions
			plotfile2 <- paste0(dirname(LoopFile), '/ratio_gene_expr_boxplot.pdf')
			cutoff <- data.frame(yintercept=0, cutoff=factor(0))
			currplot <- ggplot(outRatioDF, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.2, outlier.shape = NA) + scale_y_continuous(name = "(log2) expression ratio") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16), axis.text.x=element_text(angle=45, hjust=1)) + guides(fill=FALSE) + geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff, show.legend=FALSE)
			currplot + ggtitle(paste0("Ratio gene expression - ", CategoryList[1], "/", CategoryList[2]))
			ggsave(plotfile2, width=10, height=6)

		}

		# close the text file
		close(fp_out)	

	}	# end gene expression condition

}	# end WritePromCentricLoops function


#==============================
# this function processes only those loops (promoter specific) which belong to 
# genes which do not have exclusive loops in both of the input categories
#==============================
ProcessDiffLoopsExclSingleCategory <- function(RefGeneData, InpData, InpDir, CategoryList, ReplicaCount, CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp=TRUE) {

	Idx_ExclSingleCat <- which(InpData[,7] %in% RefGeneData)
	if (length(Idx_ExclSingleCat) > 0) {
		InpData_ExclSingleCat <- InpData[Idx_ExclSingleCat, ]
		InpDir_ExclSingleCat <- paste0(InpDir, '/Genes_Not_Excl_Both_Categories')
		system(paste("mkdir -p", InpDir_ExclSingleCat))
		LoopFile <- paste0(InpDir_ExclSingleCat, '/Loops.bed')
		WritePromCentricLoops(InpData_ExclSingleCat, LoopFile, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]), CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp)
	}
}

# #==============================
# # this function processes only those loops (promoter specific) which belong to 
# # genes uniquely present in the list of differential loops
# #==============================
# ProcessUniqueGeneDiffLoops <- function(CntGene_DiffLoopData_uniq, InpData, InpDir, CategoryList, ReplicaCount, CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2) {

# 	Idx_UniqGene <- which(InpData[,7] %in% CntGene_DiffLoopData_uniq[,1])
# 	if (length(Idx_UniqGene) > 0) {
# 		InpData_UniqGene <- InpData[Idx_UniqGene, ]
# 		InpDir_UniqGene <- paste0(InpDir, '/Genes_Count_1')
# 		system(paste("mkdir -p", InpDir_UniqGene))
# 		LoopFile <- paste0(InpDir_UniqGene, '/Loops.bed')
# 		WritePromCentricLoops(InpData_UniqGene, LoopFile, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]), CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2)
# 	}
# }

#======================
# function to get the promoter and enhancers 
# from one interacting segment (out of two interacting segments)
# in a FitHiChIP loop
#======================
GetPromEnh_OneSegment <- function(Segdata, TSSdata, ChIPSeqData, Ov_Offset=5000) {

	# overlap of the interacting segment
	# with the TSS data
	# offset of 5 Kb is used for promoters
	OvTSS_Seg1 <- Overlap1D(Segdata, TSSdata, boundary=0, offset=Ov_Offset, uniqov=FALSE)
	
	# index of promoters
	PromIdx_Seg1 <- OvTSS_Seg1$A_AND_B

	# overlap of first part of interacting segment
	# with the reference ChIP seq file
	# to get the enhancers index
	OvPeak_Seg1 <- Overlap1D(Segdata, ChIPSeqData, boundary=0, offset=0, uniqov=FALSE)
	EnhIdx_Seg1 <- setdiff(OvPeak_Seg1$A_AND_B, PromIdx_Seg1)

	# cat(sprintf("\n -->> exiting the function GetPromEnh_OneSegment \n"))

	newList <- list(Prom = PromIdx_Seg1, Enh = EnhIdx_Seg1)
	return(newList)
}

#=========================
# this function annotates FitHiChIP or similar types of loop files
# using promoter or enhancer annotations
# by using the input GTF annotation file
# and the reference ChIP seq files
# parameters:
# InpLoopFile: input FitHiChIP loop file
# InpTSSFile: file containing gene (TSS) information
# ListChIPSeqFiles: list of reference ChIP-seq information for different categories
# AnnotatedLoopFile: output loop file, similar to the FitHiChIP loop file format, 
# 					 but with two additional columns: labels of each bins (P / E / O)
# Ov_Offset: default 5000 means that at most 5 kb gap from each segment is to be provided
#=========================
Annotate_Loops_P_E <- function(InpLoopFile, InpTSSFile, ListChIPSeqFiles, AnnotatedLoopFile, CategoryList, Ov_Offset=5000) {

	# read the FitHiChIP loops
	if (file_ext(InpLoopFile) == "gz") {
		FitHiChIPLoopData <- read.table(gzfile(InpLoopFile), header=T, sep="\t", stringsAsFactors=F)
	} else {
		FitHiChIPLoopData <- read.table(InpLoopFile, header=T, sep="\t", stringsAsFactors=F)
	}

	# read the GTF file containing TSS annotation
	# first three columns are assumed to have the chromosome segment information
	if (file_ext(InpTSSFile) == "gz") {
		InpGTFData <- read.table(gzfile(InpTSSFile), header=T, sep="\t", stringsAsFactors=F)
	} else {
		InpGTFData <- read.table(InpTSSFile, header=T, sep="\t", stringsAsFactors=F)
	}

	for (i in (1:length(ListChIPSeqFiles))) {
		# current peak file
		InpChIPSeqFile <- ListChIPSeqFiles[i]
		# column number after which the label 
		# for the current category will be inserted
		currcol <- 6 + (i-1) * 2
		# read the input ChIP-seq file 
		InpPeakData <- read.table(InpChIPSeqFile, header=F, sep="\t", stringsAsFactors=F)
		# initialize labels of interacting segments
		# default: "O" (others)
		LabelVec_Seg1 <- rep("O", nrow(FitHiChIPLoopData))
		LabelVec_Seg2 <- rep("O", nrow(FitHiChIPLoopData))

		# overlap of first part of interacting segment
		# with the TSS data
		# and the reference ChIP seq file
		# to get promoters and enhancers
		IdxSet1 <- GetPromEnh_OneSegment(FitHiChIPLoopData[,1:3], InpGTFData[,1:3], InpPeakData[,1:3], Ov_Offset)
		PromIdx_Seg1 <- IdxSet1$Prom
		EnhIdx_Seg1 <- IdxSet1$Enh

		# assign the labels for promoters and enhancers
		# for the first part of interacting segment
		LabelVec_Seg1[PromIdx_Seg1] <- "P"
		LabelVec_Seg1[EnhIdx_Seg1] <- "E"
	
		# overlap of second part of interacting segment
		# with the TSS data
		# and the reference ChIP seq file
		# to get promoters and enhancers
		IdxSet2 <- GetPromEnh_OneSegment(FitHiChIPLoopData[,4:6], InpGTFData[,1:3], InpPeakData[,1:3], Ov_Offset)
		PromIdx_Seg2 <- IdxSet2$Prom
		EnhIdx_Seg2 <- IdxSet2$Enh

		# assign the labels for promoters and enhancers
		# for the first part of interacting segment
		LabelVec_Seg2[PromIdx_Seg2] <- "P"
		LabelVec_Seg2[EnhIdx_Seg2] <- "E"

		# now insert the labels in the FitHiChIP loops
		# and also adjust the column vectors
		if (i == 1) {
			AnnotatedLoopData <- cbind.data.frame(FitHiChIPLoopData[,1:currcol], LabelVec_Seg1, LabelVec_Seg2, FitHiChIPLoopData[,(currcol+1):ncol(FitHiChIPLoopData)])
		} else {
			AnnotatedLoopData <- cbind.data.frame(AnnotatedLoopData[,1:currcol], LabelVec_Seg1, LabelVec_Seg2, AnnotatedLoopData[,(currcol+1):ncol(AnnotatedLoopData)])
		}
	}	# end ChIP-seq file traversal loop

	# set the column names of the output loops
	CN <- colnames(FitHiChIPLoopData)[1:6]
	for (i in (1:length(ListChIPSeqFiles))) {
		CN <- c(CN, paste0(CategoryList[i], "_Label1"), paste0(CategoryList[i], "_Label2"))
	}
	CN <- c(CN, colnames(FitHiChIPLoopData)[7:ncol(FitHiChIPLoopData)])
	colnames(AnnotatedLoopData) <- CN

	# now write the final annotated loop in the given file
	write.table(AnnotatedLoopData, AnnotatedLoopFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

}	# end function

#===================================
# this function analyzes the input FitHiChIP loops
# and finds out the promoter centric loops
# then it prints the promoter (gene) information
# along with the other interacting segments
# parameters:
# InpGeneExprFile: file containing gene expression for the current input data / cell type
# ColGeneName:  column number within the "InpGeneExprFile", storing the gene names
# ColGeneExpr: 	column number within the "InpGeneExprFile", storing the gene expression values
# Ov_Offset: offset used to compute overlap between TSS and the interacting segments
# 			 to declare if it is a promoter
# 			 default value of this offset is 5 Kb
# OutPromoterLoopFile: output file storing the results
#===================================
Get_Promoter_Loops_FitHiChIP <- function(AnnotatedLoopFile, InpTSSFile, OutPromoterLoopFile, CategoryList, Ov_Offset=5000, GeneExprFileList=c(), GeneExprFileColGeneNameList=c(), GeneExprFileColGeneExprList=c()) {

	# read the FitHiChIP loops
	if (file_ext(AnnotatedLoopFile) == "gz") {
		AnnotatedLoopData <- read.table(gzfile(AnnotatedLoopFile), header=T, sep="\t", stringsAsFactors=F)
	} else {
		AnnotatedLoopData <- read.table(AnnotatedLoopFile, header=T, sep="\t", stringsAsFactors=F)
	}
	CN <- colnames(AnnotatedLoopData)

	# get the columns (indices) containing feature vectors of first segment
	Col_idx_Feat1 <- which(grepl('1$', CN))

	# get the columns (indices) containing feature vectors of second interacting segment
	Col_idx_Feat2 <- which(grepl('2$', CN))

	# column names which are generalized (from being specific to interacting segment)
	CN1 <- CN[Col_idx_Feat1]
	CN1 <- gsub('.{1}$', '', CN1)

	# get the columns not specific to any of the interacting segments
	# (apart from the first 10 columns)
	Col_Feat_Gen <- setdiff(setdiff(seq(1, ncol(AnnotatedLoopData)), Col_idx_Feat1), Col_idx_Feat2)

	# read the GTF file containing TSS annotation
	if (file_ext(InpTSSFile) == "gz") {
		InpGTFData <- read.table(gzfile(InpTSSFile), header=T, sep="\t", stringsAsFactors=F)
	} else {
		InpGTFData <- read.table(InpTSSFile, header=T, sep="\t", stringsAsFactors=F)
	}

	# first find out the promoter specific loops
	# from the annotated interactions

	# loop indices whose first part is promoter
	Loop_Idx_Prom1 <- which((AnnotatedLoopData[, 7] == "P") & (AnnotatedLoopData[, 9] == "P"))

	# loop indices whose second part is promoter
	Loop_Idx_Prom2 <- which((AnnotatedLoopData[, 8] == "P") & (AnnotatedLoopData[, 10] == "P"))

	#====================
	# now find the promoters and corresponding interacting segments
	#====================

	# initialize - important for error cases - sourya
	GTF_Idx_1 <- c()
	Loop_Idx_1 <- c()

	if (length(Loop_Idx_Prom1) > 0) {
		AnnotatedLoopData_Curr <- AnnotatedLoopData[Loop_Idx_Prom1, ]
		CurrOv1 <- Overlap1D(InpGTFData[,1:3], AnnotatedLoopData_Curr[,1:3], boundary=0, offset=Ov_Offset, uniqov=FALSE)
		# cat(sprintf("\n -->> Computed overlap of GTF with the first interacting segment"))
		GTF_Idx_1 <- CurrOv1$A_AND_B
		Loop_Idx_1 <- CurrOv1$B_AND_A
		if ((length(GTF_Idx_1) > 0) & (length(Loop_Idx_1) > 0)) {
			# construct a data frame
			# by using the GTF gene information
			# only the interacting segment of the first chromosome (where gene promoter overlaps)
			# and the second interacting segment of the overlapping loops
			OutDF1 <- cbind.data.frame(InpGTFData[GTF_Idx_1, c(1:6, 8)], AnnotatedLoopData_Curr[Loop_Idx_1, 1:3], AnnotatedLoopData_Curr[Loop_Idx_1, c(Col_idx_Feat2, Col_Feat_Gen)])
			colnames(OutDF1) <- c(colnames(InpGTFData)[c(1:6, 8)], "chrP", "startP", "endP", CN1, CN[Col_Feat_Gen])
		}	# end condition GTF indices
	}
	
	# initialize - important for error cases - sourya
	GTF_Idx_2 <- c()
	Loop_Idx_2 <- c()

	if (length(Loop_Idx_Prom2) > 0) {
		AnnotatedLoopData_Curr <- AnnotatedLoopData[Loop_Idx_Prom2, ]
		CurrOv2 <- Overlap1D(InpGTFData[,1:3], AnnotatedLoopData_Curr[,4:6], boundary=0, offset=Ov_Offset, uniqov=FALSE)
		# cat(sprintf("\n -->> Computed overlap of GTF with the second interacting segment"))
		GTF_Idx_2 <- CurrOv2$A_AND_B
		Loop_Idx_2 <- CurrOv2$B_AND_A
		if ((length(GTF_Idx_2) > 0) & (length(Loop_Idx_2) > 0)) {
			# construct a data frame
			# by using the GTF gene information
			# only the interacting segment of the second chromosome (where gene promoter overlaps)
			# and the first interacting segment of the overlapping loops
			OutDF2 <- cbind.data.frame(InpGTFData[GTF_Idx_2, c(1:6, 8)], AnnotatedLoopData_Curr[Loop_Idx_2, 4:6], AnnotatedLoopData_Curr[Loop_Idx_2, c(Col_idx_Feat1, Col_Feat_Gen)])
			colnames(OutDF2) <- c(colnames(InpGTFData)[c(1:6, 8)], "chrP", "startP", "endP", CN1, CN[Col_Feat_Gen])
		}	# end condition GTF indices
	}

	# combine these two data frames
	if ((length(GTF_Idx_1) > 0) & (length(Loop_Idx_1) > 0) & (length(GTF_Idx_2) > 0) & (length(Loop_Idx_2) > 0)) {
		MasterSheetDF <- rbind.data.frame(OutDF1, OutDF2)	
	} else if ((length(GTF_Idx_1) > 0) & (length(Loop_Idx_1) > 0)) {
		MasterSheetDF <- OutDF1
	} else if ((length(GTF_Idx_2) > 0) & (length(Loop_Idx_2) > 0)) {
		MasterSheetDF <- OutDF2
	}

	# sort this master data frame
	# with respect to the genes
	if (((length(GTF_Idx_1) > 0) & (length(Loop_Idx_1) > 0)) | ((length(GTF_Idx_2) > 0) & (length(Loop_Idx_2) > 0))) {
		# cat(sprintf("\n -->> before sorting --- MasterSheetDF dimension -- nrow: %s  ncol : %s ", nrow(MasterSheetDF), ncol(MasterSheetDF)))

		# sort the data frame with respect to 
		# columns 1 (chr), 2 (number), 11 (chr), 12 (number)
		MasterSheetDF <- MasterSheetDF[ order( MasterSheetDF[,1], MasterSheetDF[,2], MasterSheetDF[,11], MasterSheetDF[,12] ), ]
		# cat(sprintf("\n -->> After sorting --- MasterSheetDF dimension -- nrow: %s  ncol : %s ", nrow(MasterSheetDF), ncol(MasterSheetDF)))

		if ((length(GeneExprFileList) > 0) & (length(GeneExprFileColGeneNameList) > 0) & (length(GeneExprFileColGeneExprList) > 0)) {

			# now also merge the gene expression data
			# for individual promoters
			# and for individual cell types
			for (i in (1:length(GeneExprFileList))) {
				# current gene expression files
				# and the associated column information
				InpGeneExprFile <- GeneExprFileList[i]
				ColGeneName <- GeneExprFileColGeneNameList[i]
				ColGeneExpr <- GeneExprFileColGeneExprList[i]

				if (file_ext(InpGeneExprFile) == "gz") {
					GeneExprData <- read.table(gzfile(InpGeneExprFile), header=T, sep="\t", stringsAsFactors=FALSE)
				} else {
					GeneExprData <- read.table(InpGeneExprFile, header=T, sep="\t", stringsAsFactors=FALSE)
				}

				# first initialize a vector which will contain the expression value
				# for all the interactions
				ExprVal <- matrix(0, nrow=nrow(MasterSheetDF), ncol=1)
				colnames(ExprVal) <- paste0("GeneExpr_", CategoryList[i])

				# check the gene names belonging to the 
				# MasterSheetDF (7th column)
				# with the gene expression data
				m <- match(MasterSheetDF[,7], GeneExprData[, ColGeneName])

				# get the indices (rows) of the master sheet
				# whose gene names match with the genes 
				# provided in the gene expression data 
				# (removing the NA cases from the match output)
				idx_MasterSheet_Ov_GeneExpr <- which(!is.na(m))

				# get the indices in the gene expression data
				# matching with the master data frame
				idx_GeneExpr_Ov_MasterSheet <- m[!is.na(m)]

				# copy the gene expression value
				ExprVal[idx_MasterSheet_Ov_GeneExpr] <- GeneExprData[idx_GeneExpr_Ov_MasterSheet, ColGeneExpr]

				# append the gene expression data in the last column 
				# of the master sheet data frame
				MasterSheetDF <- cbind.data.frame(MasterSheetDF, ExprVal)
			
			}	# end loop gene expression data

		}	# end condition on availability of gene expression file

		# write the overlapping interactions of the TSS
		# to output file		
		write.table(MasterSheetDF, OutPromoterLoopFile, row.names=FALSE, col.names=TRUE, sep = "\t", quote=FALSE, append=FALSE)	
		# cat(sprintf("\n -->> Written the output master sheet file"))

	}	# end condition 

}	# end function

#===========================================================
option_list = list(
	# make_option(c("--LoopList"), type="character", default=NULL, help="Comma or colon separated list of interaction / loop files containing FitHiChIP SIGNIFICANT loops. Mandatory parameter."),
	make_option(c("--AllLoopList"), type="character", default=NULL, help="Comma or colon separated list of all possible interactions, along with their FItHiChIP significance value (would include both significant and insignificant loops as well). Mandatory parameter."),
	make_option(c("--FDRThrLoop"), type="numeric", default=0.01, help="FDR threshold used for FitHiChIP loops. Default = 0.01 (same as FitHiChIP implementation)."),

	make_option(c("--OutDir"), type="character", default=NULL, help="Base Output directory. Mandatory parameter."),

  	make_option(c("--UseRawCC"), type="integer", action="store", default=0, help="If 1, uses the raw contact count for replicate analysis. If 0, uses the observed contact count (c), multiplied by the ratio (c/e) (upto nearest integer) where e is the expected contact count of this locus pair, according to the spline fit and bias regression model of FitHiChIP (listed in the field exp_cc_bias of FitHiChIP loop file). Default value of this parameter is 0."),

	# make_option(c("--BiasFileList"), type="character", default=NULL, help="Comma or colon separated list of bias files generated by FitHiChIP loops. Mandatory parameter."),	

	make_option(c("--PeakFileCat1"), type="character", default=NULL, help="Peak file used for the samples in the first category (to infer peak specific HiChIP loops). Mandatory parameter."),
	make_option(c("--PeakFileCat2"), type="character", default=NULL, help="Peak file used for the samples in the second category (to infer peak specific HiChIP loops). Mandatory parameter."),	

	make_option(c("--CategoryList"), type="character", default=NULL, help="Comma or colon separated list of the two main categories (whose replicates are present). Default: Category1, Category2."),
	
	make_option(c("--ReplicaCount"), type="character", default=NULL, help="Comma or colon separated list of the count of replicates for individual categories. Default: 1,1 (means that we are considering one replicate per sample)."),
	make_option(c("--ReplicaLabels1"), type="character", default=NULL, help="Comma or colon separated list of the label of replicates for the first category. Default: R1,R2, etc."),
	make_option(c("--ReplicaLabels2"), type="character", default=NULL, help="Comma or colon separated list of the label of replicates for the second category. Default: R1,R2, etc."),

	# make_option(c("--ChIPCovFileCat1"), type="character", default=NULL, help="File storing the ChIP coverage for the targeted bin size, corresponding to the first category. Mandatory parameter."),
	# make_option(c("--ChIPCovFileCat2"), type="character", default=NULL, help="File storing the ChIP coverage for the targeted bin size, corresponding to the second category. Mandatory parameter."),

	make_option(c("--BinCoverageList"), type="character", default=NULL, help="List of files storing the ChIP-seq coverage for the specified bin size. Mandatory parameter."),
	
	make_option(c("--InpTSSFile"), type="character", default=NULL, help="TSS containing file for the reference genome. Mandatory parameter."),

	make_option(c("--GeneExprFileList"), type="character", default=NULL, help="Comma or colon separated list of gene expression containing files, for the cell lines checked."),
	
	make_option(c("--GeneNameColList"), type="character", default=NULL, help="Comma or colon separated list of columns containing the gene names, for the above mentioned gene expression files."),
	
	make_option(c("--ExprValColList"), type="character", default=NULL, help="Comma or colon separated list of columns containing the gene expression values, for the above mentioned gene expression files."),

  	# make_option(c("--UseDESeq"), type="integer", action="store", default=0, help="If 1, DESeq2 is used for differential analysis. Else, EdgeR is used. Default value of this parameter is 0 (means EdgeR is used)."),	

  	make_option(c("--FoldChangeThr"), type="integer", action="store", default=3, help="DESeq / EdgeR fold change threshold - log2 of this value is used. Default = 3, means that log2(3) is used as the fold change threshold."),	
  	
  	make_option(c("--FDRThr"), type="numeric", default=0.05, help="FDR threshold for DESeq / EdgeR. Default is 0.05, means that loops with FDR < 0.05, and fold change >= log2(FoldChangeThr) would be considered as differential."),

  	make_option(c("--bcv"), type="numeric", default=0.4, help="If EdgeR is used with single samples (replica count = 1 for any of the categories), this value is the square-root-dispersion. For datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates. For details, see the edgeR manual. By default, the value is set as 0.4.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#====================================
# process the input arguments

USEDESEQ=0	#opt$UseDESeq

# significant loops by FitHiChIP
# if (is.null(opt$LoopList)) {
# 	print_help(parser)
# 	stop("Significant interaction / loop files of FitHiChIP are not provided - check the option --LoopList \n", call.=FALSE)
# } else {
# 	LoopList <- as.character(unlist(strsplit(opt$LoopList,"[,:]")))
# }


# all loops by FitHiChIP
if (is.null(opt$AllLoopList)) {
	print_help(parser)
	stop("Files containing all (significant or not) loops along with their FitHiChIP significance are not provided - check the option --AllLoopList \n", call.=FALSE)
} else {
	AllLoopList <- as.character(unlist(strsplit(opt$AllLoopList,"[,:]")))
	# if (length(AllLoopList) != length(LoopList)) {
	# 	stop("Number of Files in the --AllLoopList parameter is not same as the number of files in the --LoopList parameter !! quit !!! \n", call.=FALSE)
	# }
}

if (is.null(opt$OutDir)) {
	print_help(parser)
	stop("Output directory is not provided - check the option --OutDir \n", call.=FALSE)
}

if (is.null(opt$CategoryList)) {
	CategoryList <- paste0('Category', seq(1:2))
} else {
	CategoryList <- as.character(unlist(strsplit(opt$CategoryList,"[,:]")))
}

if (is.null(opt$ReplicaCount)) {
	ReplicaCount <- c(1,1)
} else {
	ReplicaCount <- as.integer(unlist(strsplit(opt$ReplicaCount,"[,:]")))
}

if (sum(ReplicaCount) != length(AllLoopList)) {
	stop("Count of the replicates does not match with the number of input loops files - check the option --ReplicaCount \n", call.=FALSE)
}

if (is.null(opt$ReplicaLabels1)) {
	ReplicaLabels1 <- paste0('R', seq(1:ReplicaCount[1]))
} else {
	ReplicaLabels1 <- as.character(unlist(strsplit(opt$ReplicaLabels1,"[,:]")))
	if (length(ReplicaLabels1) != ReplicaCount[1]) {
		stop("Number of labels for the replicates in the first category does not match with the number of replicates provided for the first category - check the options --ReplicaCount and --ReplicaLabels1 \n", call.=FALSE)
	}
}

if (is.null(opt$ReplicaLabels2)) {
	ReplicaLabels2 <- paste0('R', seq(1:ReplicaCount[2]))
} else {
	ReplicaLabels2 <- as.character(unlist(strsplit(opt$ReplicaLabels2,"[,:]")))
	if (length(ReplicaLabels2) != ReplicaCount[2]) {
		stop("Number of labels for the replicates in the first category does not match with the number of replicates provided for the first category - check the options --ReplicaCount and --ReplicaLabels2 \n", call.=FALSE)
	}
}

# parse the list of input files 
# storing the ChIP-seq coverage for the specified bin size
if (is.null(opt$BinCoverageList)) {
	print_help(parser)
	stop("File list containing the ChIP-seq coverage for the specified bin size is not provided - check the option --BinCoverageList \n", call.=FALSE)
} else {
	BinCoverageList <- as.character(unlist(strsplit(opt$BinCoverageList,"[,:]")))
}

# # parse the bias file list
# if (is.null(opt$BiasFileList)) {
# 	print_help(parser)
# 	stop("Files containing bias values are not provided - check the option --BiasFileList \n", call.=FALSE)
# } else {
# 	BiasFileList <- as.character(unlist(strsplit(opt$BiasFileList,"[,:]")))
# 	if (length(BiasFileList) != length(AllLoopList)) {
# 		stop("Number of Files in the --BiasFileList parameter is not same as the number of files in the --AllLoopList parameter !! quit !!! \n", call.=FALSE)
# 	}
# }

if (is.null(opt$PeakFileCat1)) {
	print_help(parser)
	stop("Peak file for the first category is not provided - check the option --PeakFileCat1 \n", call.=FALSE)
}

if (is.null(opt$PeakFileCat2)) {
	print_help(parser)
	stop("Peak file for the second category is not provided - check the option --PeakFileCat2 \n", call.=FALSE)
}

# if (is.null(opt$ChIPCovFileCat1)) {
# 	print_help(parser)
# 	stop("ChIP-seq coverage file (with respect to the given bin size) for the first category is not provided - check the option --ChIPCovFileCat1 \n", call.=FALSE)
# }

# if (is.null(opt$ChIPCovFileCat2)) {
# 	print_help(parser)
# 	stop("ChIP-seq coverage file (with respect to the given bin size) for the second category is not provided - check the option --ChIPCovFileCat2 \n", call.=FALSE)
# }

if (is.null(opt$InpTSSFile)) {
	print_help(parser)
	stop("TSS containing file is not provided - check the option --InpTSSFile \n", call.=FALSE)
}

#===============================
# gene expression related processing

# parse the gene expression file list
if (is.null(opt$GeneExprFileList)) {
# 	print_help(parser)
# 	stop("Files containing gene expression values are not provided  - check the option --GeneExprFileList \n", call.=FALSE)
	# boolean variable indicating whether gene expression data is provided or not
	GeneExprDataInp <- FALSE
} else {
	GeneExprFileList <- as.character(unlist(strsplit(opt$GeneExprFileList,"[,:]")))
	# boolean variable indicating whether gene expression data is provided or not
	GeneExprDataInp <- TRUE

	# now parse the list of column indices (integers)
	# containing the name of genes in the gene expression files
	if (is.null(opt$GeneNameColList)) {
		print_help(parser)
		stop("Gene expression files are provided but the Columns containing gene names within gene expression files are not provided  - check the option --GeneNameColList \n", call.=FALSE)
	} else {
		GeneNameColList <- as.integer(unlist(strsplit(opt$GeneNameColList,"[,:]")))
		if (length(GeneNameColList) < length(GeneExprFileList)) {
			GeneNameColList <- c(GeneNameColList, rep(GeneNameColList[length(GeneNameColList)], (length(GeneExprFileList) - length(GeneNameColList))))
		}		
	}

	# now parse the list of column indices (integers)
	# containing the expression of genes in the gene expression files
	if (is.null(opt$ExprValColList)) {
		print_help(parser)
		stop("Gene expression files are provided but the Columns containing gene expression values are not provided  - check the option --ExprValColList \n", call.=FALSE)
	} else {
		ExprValColList <- as.integer(unlist(strsplit(opt$ExprValColList,"[,:]")))
		if (length(ExprValColList) < length(GeneExprFileList)) {
			ExprValColList <- c(ExprValColList, rep(ExprValColList[length(ExprValColList)], (length(GeneExprFileList) - length(ExprValColList))))
		}
	}
}

#===============================

# # check the replica count for both categories
# # if any of the groups involve no replicate 
# # and DESeq is used 
# # then stop 
# if (((ReplicaCount[1] == 1) | (ReplicaCount[2] == 1)) & (USEDESEQ == 1)) {
# 	stop("At least one of the input categories do not involve multiple replicates - so DESeq cannot be used for the differential analysis - exit !!! \n", call.=FALSE)
# }


# DESEq / EdgeR thresholds
FOLD_Change_Thr <- as.integer(opt$FoldChangeThr)
FDR_Th_DESeq <- as.numeric(opt$FDRThr)

# combining all the replica labels
if (is.null(opt$ReplicaLabels1)) {
	AllRepLabels <- c(paste0(CategoryList[1], '_', ReplicaLabels1))
} else {
	AllRepLabels <- c(ReplicaLabels1)
}

if (is.null(opt$ReplicaLabels2)) {
	AllRepLabels <- c(AllRepLabels, c(paste0(CategoryList[2], '_', ReplicaLabels2)))
} else {
	AllRepLabels <- c(AllRepLabels, c(ReplicaLabels2))
}

cat(sprintf("\n ===>>> CategoryList : %s ", paste0(CategoryList, sep="\t")))
cat(sprintf("\n ===>>> ReplicaCount : %s ", paste0(ReplicaCount, sep="\t")))
cat(sprintf("\n ===>>> ReplicaLabels1 : %s ", paste0(ReplicaLabels1, sep="\t")))
cat(sprintf("\n ===>>> ReplicaLabels2 : %s ", paste0(ReplicaLabels2, sep="\t")))
cat(sprintf("\n ===>>> AllRepLabels : %s ", paste0(AllRepLabels, sep="\t")))


#====================================
# set of chromosomes which would be individually analyzed for testing overlap
ChrList_NameNum <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")

# number of features for individual loops
# contact count, q-val, bias1, bias2 = 4,
numFeat <- 4

# directory designated for differential analysis
if (opt$UseRawCC == 0) {
	OutTableDir <- paste0(opt$OutDir, '/ExpCC')
} else {
	OutTableDir <- paste0(opt$OutDir, '/RawCC')
}

# if (USEDESEQ == 0) {
MainDir <- paste0(OutTableDir, '/EdgeR')
# } else {
# 	MainDir <- paste0(OutTableDir, '/DESeq')
# }
system(paste("mkdir -p", MainDir))

#================================
# at first create a file which contains the union of loops
# for all the replicates of the two categories
# Note: these loops are significant in at least one of the input files
#================================
UnionLoopFile <- paste0(OutTableDir, '/Union_', CategoryList[1], '_', CategoryList[2], '_Loops.bed')
UnionLoopTempFile <- paste0(OutTableDir, '/Union_', CategoryList[1], '_', CategoryList[2], '_Loops_Temp.bed')

if (1) { #(file.exists(UnionLoopFile) == FALSE) {

	# merge all the loops for individual files within a file
	MergeLoops(AllLoopList, UnionLoopFile, UnionLoopTempFile, opt$FDRThrLoop)

	# fill the following information
	# for individual loops from the union set
	# contact count, bias, significance, isPeak
	# with respect to individual FitHiChIP significance files 
	# (containing significant + non-significant interactions)
	
	# FillFeatureValues(UnionLoopFile, AllLoopList, c(opt$PeakFileCat1, opt$PeakFileCat2), c(opt$ChIPCovFileCat1, opt$ChIPCovFileCat2), BiasFileList, ChrList_NameNum, opt$UseRawCC, AllRepLabels, CategoryList)

	# removing bias and CHIP-seq coverage statistics
	FillFeatureValues(UnionLoopFile, AllLoopList, c(opt$PeakFileCat1, opt$PeakFileCat2), ChrList_NameNum, opt$UseRawCC, AllRepLabels, CategoryList)	

}	# end file exist condition

#====================================
# now apply EdgeR to get the differential loops between two categories
# with respect to the loop file "EdgeR_CC_Significance.bed"
#===================================

# columns containing counts
# first 14 columns contain the interacting segments (6 columns)
# peak information for both segments (8 columns)
CountCol <- as.vector(14 + (seq(1:length(AllRepLabels)) - 1) * numFeat + 1)
cat(sprintf("\n ===>>> CountCol : %s ", paste0(CountCol, sep="\t")))

# if (USEDESEQ == 0) {
# file name prefix string
PrefixStr <- 'EdgeR'
# } else {
# 	# file name prefix string
# 	PrefixStr <- 'DESeq'
# }

# file names storing different loops (significant / all)
ResFile <- paste0(MainDir, '/', PrefixStr, '_CC_Significance.bed')
ResFile_Expand_FINAL_Sig <- paste0(MainDir, '/', PrefixStr, '_CC_Significance_Expanded_FINAL_Sig.bed')

# depending on either DESeq or EdgeR 
# the "ResFile" contains the significance output
if (1) { #(file.exists(ResFile) == FALSE) {
	# if (USEDESEQ == 0) {
	# apply EdgeR function
	ApplyEdgeR(UnionLoopFile, MainDir, CountCol, CategoryList, ReplicaCount, ResFile, opt$bcv)
	# } else {
	# 	# apply DEseq on the generated dataset
	# 	ApplyDESeq(UnionLoopFile, MainDir, CountCol, CategoryList, ReplicaCount, ResFile, FOLD_Change_Thr)
	# }	# end DESeq / EdgeR condition
}

# read the EdgeR / DESeq output
FinalOutDF <- read.table(ResFile, header=T, sep="\t", stringsAsFactors=F)

# also extract the DEseq significant loops
# based on FDR and log2 fold change threshold
# if (USEDESEQ == 0) {
# edgeR output
SigIdx <- which((FinalOutDF[,ncol(FinalOutDF)] < FDR_Th_DESeq) & (abs(FinalOutDF[,(ncol(FinalOutDF) - 3)]) >= log2(FOLD_Change_Thr)))
# } else {
# 	# DEseq output
# 	SigIdx <- which((FinalOutDF[,ncol(FinalOutDF)] < FDR_Th_DESeq) & (abs(FinalOutDF[,(ncol(FinalOutDF) - 4)]) >= log2(FOLD_Change_Thr)))
# }
write.table(FinalOutDF[SigIdx, ], ResFile_Expand_FINAL_Sig, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

MainDF <- read.table(ResFile_Expand_FINAL_Sig, header=T, sep="\t", stringsAsFactors=F)

CurrMainBaseDir <- MainDir 	#paste0(MainDir, '/Diffloop_HiChIP')
system(paste("mkdir -p", CurrMainBaseDir))

#====================================
# also apply EdgeR to get the differential 1D ChIP-seq bins
#====================================

PrefixStr <- 'EdgeR'
EdgeRCovFile <- paste0(CurrMainBaseDir, '/', PrefixStr, '_ChIPCoverage_Significance.bed')
EdgeRCovSigFile <- paste0(CurrMainBaseDir, '/', PrefixStr, '_ChIPCoverage_Significance_FINAL_Sig.bed')

# apply EdgeR on the bins: ChIP-seq coverage is the count measure
if (1) { #(file.exists(EdgeRCovFile) == FALSE) {
	ApplyEdgeR_ChIPCoverage(BinCoverageList, rep(1,length(BinCoverageList)), rep(4,length(BinCoverageList)), CategoryList, ReplicaCount, EdgeRCovFile, bcv=0.4)
}

# get the significant differential bins
# fold change threshold of 2 is employ
EdgeRCovOutDF <- read.table(EdgeRCovFile, header=T, sep="\t", stringsAsFactors=F)

# using the input parameter fold change
EdgeRCovSigIdx <- which((EdgeRCovOutDF[,ncol(EdgeRCovOutDF)] < FDR_Th_DESeq) & (abs(EdgeRCovOutDF[,(ncol(EdgeRCovOutDF) - 3)]) >= log2(FOLD_Change_Thr)))

cat(sprintf("\n length EdgeRCovSigIdx : %s  ", length(EdgeRCovSigIdx)))

# store the differential 1D segments 
EdgeRCov_Differential_Bin <- EdgeRCovOutDF[EdgeRCovSigIdx, ]
write.table(EdgeRCov_Differential_Bin, EdgeRCovSigFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

#===========================================
# processing the differential loops which are different in 3D
# but not different in 1D
#===========================================

# with respect to the promoter (gene) specific loop file (Annot_DiffLoops_Promoter_Specific.bed)
# first 19 columns are for different interacting segments
# and the ChIP-seq related features like peak strength, etc.

# # columns containing bias values for the first category
# biascol_cat1 <- 21 + seq(1, ReplicaCount[1])
# # columns containing bias values for the second category
# biascol_cat2 <- 21 + ReplicaCount[1] + seq(1, ReplicaCount[2])
# columns containing contact count for different replicates of the first category
CCColList_Cat1 <- 19 + ((seq(1, ReplicaCount[1]) - 1) * 2 + 1)
# columns containing contact count for different replicates of the second category
CCColList_Cat2 <- 19 + 2 * ReplicaCount[1] + ((seq(1, ReplicaCount[2]) - 1) * 2 + 1)
# columns containing q-val for different replicates of the first category
QValColList_Cat1 <- CCColList_Cat1 + 1
# columns containing q-val for different replicates of the second category
QValColList_Cat2 <- CCColList_Cat2 + 1

# cat(sprintf("\n\n biascol_cat1 : %s  ", paste(biascol_cat1, collapse=" ")))
# cat(sprintf("\n\n biascol_cat2 : %s  ", paste(biascol_cat2, collapse=" ")))
cat(sprintf("\n\n CCColList_Cat1 : %s  ", paste(CCColList_Cat1, collapse=" ")))
cat(sprintf("\n\n CCColList_Cat2 : %s  ", paste(CCColList_Cat2, collapse=" ")))
cat(sprintf("\n\n QValColList_Cat1 : %s  ", paste(QValColList_Cat1, collapse=" ")))
cat(sprintf("\n\n QValColList_Cat2 : %s  ", paste(QValColList_Cat2, collapse=" ")))

#====================
# first annotate the differential loops obtained from EdgeR
# allow 5 Kb offset from the promoter
Annotated_Differential_LoopFile <- paste0(CurrMainBaseDir, '/Annot_DiffLoops.bed')
Annotate_Loops_P_E(ResFile_Expand_FINAL_Sig, opt$InpTSSFile, c(opt$PeakFileCat1, opt$PeakFileCat2), Annotated_Differential_LoopFile, CategoryList, 5000)

# read the annotated loops
Annot_DiffLoopData <- read.table(Annotated_Differential_LoopFile, header=T, sep="\t", stringsAsFactors=F)

# get the promoter specific loops
# and also list its expressions for the input categories
Promoter_Specific_DiffLoop_File <- paste0(CurrMainBaseDir, '/Annot_DiffLoops_Promoter_Specific.bed')
if (GeneExprDataInp == TRUE) {
	Get_Promoter_Loops_FitHiChIP(Annotated_Differential_LoopFile, opt$InpTSSFile, Promoter_Specific_DiffLoop_File, CategoryList, 5000, GeneExprFileList, GeneNameColList, ExprValColList)
} else {
	Get_Promoter_Loops_FitHiChIP(Annotated_Differential_LoopFile, opt$InpTSSFile, Promoter_Specific_DiffLoop_File, CategoryList, 5000)
}

# read these promoter specific loops 
Promoter_Specific_DiffLoop_Data <- read.table(Promoter_Specific_DiffLoop_File, header=T, sep="\t", stringsAsFactors=F)

# debug - sourya
if (0) {
	LoopFile <- paste0(CurrMainBaseDir, '/Loops.bed')
	WritePromCentricLoops(Promoter_Specific_DiffLoop_Data, LoopFile, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]), CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp)	
}

#*********************
# derived file - 1
# now get the count of genes associated with this data
# two column data frame: genes and their counts
CntGene_DiffLoopData <- as.data.frame(table(Promoter_Specific_DiffLoop_Data[,7]))

# get the genes having count = 1 (unique genes) in this data
# two column data frame: genes and their counts (1)
CntGene_DiffLoopData_uniq <- CntGene_DiffLoopData[which(CntGene_DiffLoopData[,2]==1), ]

#*********************
# derived file - 2
# get the genes such that they do not have exclusive loops in two categories

# loop indices having exclusive loops in the first category
Idx_Promoter_Specific_DiffLoop_ExclCat1 <- GetNoneSigIdx(Promoter_Specific_DiffLoop_Data, CCColList_Cat2, QValColList_Cat2)

# genes having exclusive loops in the first category
Genes_ExclCat1 <- c()
if (length(Idx_Promoter_Specific_DiffLoop_ExclCat1) > 0) {
	Genes_ExclCat1 <- unique(Promoter_Specific_DiffLoop_Data[Idx_Promoter_Specific_DiffLoop_ExclCat1, 7])
}

# loop indices having exclusive loops in the second category
Idx_Promoter_Specific_DiffLoop_ExclCat2 <- GetNoneSigIdx(Promoter_Specific_DiffLoop_Data, CCColList_Cat1, QValColList_Cat1)

# genes having exclusive loops in the second category
Genes_ExclCat2 <- c()
if (length(Idx_Promoter_Specific_DiffLoop_ExclCat2) > 0) {
	Genes_ExclCat2 <- unique(Promoter_Specific_DiffLoop_Data[Idx_Promoter_Specific_DiffLoop_ExclCat2, 7])
}

# genes having exclusive loops in both categories
genes_Excl_bothCat <- intersect(Genes_ExclCat1, Genes_ExclCat2)

# genes having exclusive loops in either of these two categories (but not both)
# or not having any exclusive loops
Genes_Excl_SingleCat <- setdiff(unique(Promoter_Specific_DiffLoop_Data[,7]), genes_Excl_bothCat)

# now check the annotated loop file
# divide interactions into three parts: 
# P-P, P-E
P_P_Dir_Main <- paste0(CurrMainBaseDir, '/Annot_P_P')
P_E_Dir_Main <- paste0(CurrMainBaseDir, '/Annot_P_E')
system(paste("mkdir -p", P_P_Dir_Main))
system(paste("mkdir -p", P_E_Dir_Main))

#--------------------------------
# first get the P-P loops 
# and dump 
P_P_DiffLoop_File <- paste0(P_P_Dir_Main, '/DiffLoops_Annot_ALL.bed')
P_P_Idx <- which((Annot_DiffLoopData[,7]=="P") & (Annot_DiffLoopData[,8]=="P") & (Annot_DiffLoopData[,9]=="P") & (Annot_DiffLoopData[,10]=="P"))

cat(sprintf("\n\n **** length P_P_Idx : %s  ", length(P_P_Idx)))

if (length(P_P_Idx) > 0) {

	P_P_DiffLoop <- Annot_DiffLoopData[P_P_Idx, ]
	write.table(P_P_DiffLoop, P_P_DiffLoop_File, row.names=FALSE, col.names=TRUE, sep = "\t", quote=FALSE, append=FALSE)

	# get the promoter specific loops
	# and also list its expressions for the input categories
	Promoter_Specific_P_P_DiffLoop_File <- paste0(P_P_Dir_Main, '/DiffLoops_Annot_ALL_Promoter_Specific.bed')
	if (GeneExprDataInp == TRUE) {
		Get_Promoter_Loops_FitHiChIP(P_P_DiffLoop_File, opt$InpTSSFile, Promoter_Specific_P_P_DiffLoop_File, CategoryList, 5000, GeneExprFileList, GeneNameColList, ExprValColList)
	} else {
		Get_Promoter_Loops_FitHiChIP(P_P_DiffLoop_File, opt$InpTSSFile, Promoter_Specific_P_P_DiffLoop_File, CategoryList, 5000)
	}

	if (0) {
		# debug - sourya	
		CurrData <- read.table(Promoter_Specific_P_P_DiffLoop_File, header=T, sep="\t", stringsAsFactors=F)
		LoopFile <- paste0(P_P_Dir_Main, '/Loops.bed')
		WritePromCentricLoops(CurrData, LoopFile, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]), CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp)		
	}

	P_P_Dir <- paste0(P_P_Dir_Main, '/Differential_3D_not_1D')
	system(paste("mkdir -p", P_P_Dir))

	# within the P-P loops, get the loops
	# which are differential in 3D 
	# but not differential in 1D
	# but get the overlap between the differential 1D bins
	# and both sides of the P-P loops
	ov <- Overlap1Dwith2D(EdgeRCov_Differential_Bin[,1:3], P_P_DiffLoop[,1:6])
	cat(sprintf("\n length ov$B1_AND_A : %s  length ov$B2_AND_A : %s  ", length(ov$B1_AND_A), length(ov$B2_AND_A)))

	# loop indices which contain differential ChIP-seq bins 
	Loop_DiffChIPIdx <- union(ov$B1_AND_A, ov$B2_AND_A)
	cat(sprintf("\n length Loop_DiffChIPIdx : %s  ", length(Loop_DiffChIPIdx)))

	# also get the loops where peak status (i.e. isPeak field) is same 
	# for both categories and both ends
	# i.e. for both the bins, isPeak fields are same
	IsPeakStatus_Same_Idx <- which((P_P_DiffLoop[,11] == P_P_DiffLoop[,15]) & (P_P_DiffLoop[,12] == P_P_DiffLoop[,16]))

	# impose further conditions on ChIP coverage
	# but this time the coverage is computed with respect to the 
	# merged ChIP-seq alignment of individual cell types 
	# (by combining all the replicates) 
	# for (cov_thr in c(0, 0.25, 0.5)) {
	for (cov_thr in c(0)) {
		# if (cov_thr == 0) {
		# 	P_P_SubDir <- paste0(P_P_Dir, '/No_Coverage_Condition')
		# } else {
		# 	P_P_SubDir <- paste0(P_P_Dir, '/Coverage_pct_', cov_thr)
		# }
		P_P_SubDir <- P_P_Dir
		system(paste("mkdir -p", P_P_SubDir))

		# cat(sprintf("\n\n ******** cov_thr : %s ", cov_thr))

		# get the P-P loop indices 
		# whose any segment does not overlap with the differential bins
		# and optionally, both segments maintain coverage constraint (on the merged ChIP-seq alignment)
		# also the loop indices should belong to the set "IsPeakStatus_Same_Idx"
		# if (cov_thr == 0) {
		Filt_P_P_Idx <- intersect(setdiff(seq(1,nrow(P_P_DiffLoop)), Loop_DiffChIPIdx), IsPeakStatus_Same_Idx)
		cat(sprintf("\n length Filt_P_P_Idx : %s  ", length(Filt_P_P_Idx)))
		# } else {
		# 	Cov_Fail_Idx <- which(((abs(P_P_DiffLoop[,19] - P_P_DiffLoop[,21]) / P_P_DiffLoop[,19]) > cov_thr) | ((abs(P_P_DiffLoop[,20] - P_P_DiffLoop[,22]) / P_P_DiffLoop[,20]) > cov_thr))
		# 	Filt_P_P_Idx <- intersect(setdiff(seq(1,nrow(P_P_DiffLoop)), union(Loop_DiffChIPIdx, Cov_Fail_Idx)), IsPeakStatus_Same_Idx)
		# 	cat(sprintf("\n length Cov_Fail_Idx : %s  ", length(Cov_Fail_Idx)))
		# 	cat(sprintf("\n length Filt_P_P_Idx : %s  ", length(Filt_P_P_Idx)))
		# }
		
		if (length(Filt_P_P_Idx) > 0) {

			# write the filtered loops (including annotations)
			P_P_DiffLoop_FiltFile <- paste0(P_P_SubDir, '/DiffLoops_Annot_Filt_ChIP1D.bed')
			write.table(P_P_DiffLoop[Filt_P_P_Idx, ], P_P_DiffLoop_FiltFile, row.names=FALSE, col.names=TRUE, sep = "\t", quote=FALSE, append=FALSE)

			# now from the filtered P-P loops 
			# get the promoter specific loops
			# and also list its expressions for the input categories
			Promoter_Specific_P_P_DiffLoop_FiltFile <- paste0(P_P_SubDir, '/DiffLoops_Annot_Filt_ChIP1D_Promoter_Specific.bed')
			if (GeneExprDataInp == TRUE) {
				Get_Promoter_Loops_FitHiChIP(P_P_DiffLoop_FiltFile, opt$InpTSSFile, Promoter_Specific_P_P_DiffLoop_FiltFile, CategoryList, 5000, GeneExprFileList, GeneNameColList, ExprValColList)
			} else {
				Get_Promoter_Loops_FitHiChIP(P_P_DiffLoop_FiltFile, opt$InpTSSFile, Promoter_Specific_P_P_DiffLoop_FiltFile, CategoryList, 5000)
			}

			InpData <- read.table(Promoter_Specific_P_P_DiffLoop_FiltFile, header=T, sep="\t", stringsAsFactors=F)

			# cat(sprintf("\n nrow InpData : %s  ", nrow(InpData)))

			if (0) {
				# case 1 - write all the genes and associated information 
				LoopFile <- paste0(P_P_SubDir, '/Loops.bed')
				WritePromCentricLoops(InpData, LoopFile, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]), CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp)				
			}

			# case 2 - now select only those genes 
			# such that they have only one loop
			# in the "Promoter_Specific_DiffLoop_Data" 
			# (containing gene information for all differential loops of EdgeR)
			
			# comment - sourya
			# ProcessUniqueGeneDiffLoops(CntGene_DiffLoopData_uniq, InpData, P_P_SubDir, CategoryList, ReplicaCount, CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2)

			# case 3 - now select only those genes 
			# such that they do not have exclusive loops in both categories
			ProcessDiffLoopsExclSingleCategory(Genes_Excl_SingleCat, InpData, P_P_SubDir, CategoryList, ReplicaCount, CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp)
		}

	}	# end coverage threshold loop

}	# end P-P loop condition


#--------------------------------
# now get the P-E loops
# and dump
P_E_DiffLoop_File <- paste0(P_E_Dir_Main, '/DiffLoops_Annot_ALL.bed')
P_E_Idx <- which(((Annot_DiffLoopData[,7]=="P") & (Annot_DiffLoopData[,8]=="E") & (Annot_DiffLoopData[,9]=="P") & (Annot_DiffLoopData[,10]=="E")) | ((Annot_DiffLoopData[,7]=="E") & (Annot_DiffLoopData[,8]=="P") & (Annot_DiffLoopData[,9]=="E") & (Annot_DiffLoopData[,10]=="P")))

cat(sprintf("\n\n ****** length P_E_Idx : %s  ", length(P_E_Idx)))

if (length(P_E_Idx) > 0) {

	P_E_DiffLoop <- Annot_DiffLoopData[P_E_Idx, ]
	write.table(P_E_DiffLoop, P_E_DiffLoop_File, row.names=FALSE, col.names=TRUE, sep = "\t", quote=FALSE, append=FALSE)

	# get the promoter specific loops
	# and also list its expressions for the input categories
	Promoter_Specific_P_E_DiffLoop_File <- paste0(P_E_Dir_Main, '/DiffLoops_Annot_ALL_Promoter_Specific.bed')
	if (GeneExprDataInp == TRUE) {
		Get_Promoter_Loops_FitHiChIP(P_E_DiffLoop_File, opt$InpTSSFile, Promoter_Specific_P_E_DiffLoop_File, CategoryList, 5000, GeneExprFileList, GeneNameColList, ExprValColList)
	} else {
		Get_Promoter_Loops_FitHiChIP(P_E_DiffLoop_File, opt$InpTSSFile, Promoter_Specific_P_E_DiffLoop_File, CategoryList, 5000)
	}

	if (0) {
		# debug - sourya
		CurrData <- read.table(Promoter_Specific_P_E_DiffLoop_File, header=T, sep="\t", stringsAsFactors=F)
		LoopFile <- paste0(P_E_Dir_Main, '/Loops.bed')
		WritePromCentricLoops(CurrData, LoopFile, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]), CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp)	
	}

	P_E_Dir <- paste0(P_E_Dir_Main, '/Differential_3D_not_1D')
	system(paste("mkdir -p", P_E_Dir))

	# within the P-E loops, get the loops
	# which are differential in 3D 
	# but not differential in 1D
	# but get the overlap between the differential 1D bins
	# and both sides of the P-E loops
	ov <- Overlap1Dwith2D(EdgeRCov_Differential_Bin[,1:3], P_E_DiffLoop[,1:6])
	cat(sprintf("\n length ov$B1_AND_A : %s  length ov$B2_AND_A : %s  ", length(ov$B1_AND_A), length(ov$B2_AND_A)))

	# loop indices which contain differential ChIP-seq bins 
	Loop_DiffChIPIdx <- union(ov$B1_AND_A, ov$B2_AND_A)	
	cat(sprintf("\n length Loop_DiffChIPIdx : %s  ", length(Loop_DiffChIPIdx)))

	# also get the loops where peak status (i.e. isPeak field) is same 
	# for both categories and both ends
	# i.e. for both the bins, isPeak fields are same
	IsPeakStatus_Same_Idx <- which((P_E_DiffLoop[,11] == P_E_DiffLoop[,15]) & (P_E_DiffLoop[,12] == P_E_DiffLoop[,16]))	

	# impose further conditions on ChIP coverage
	# but this time the coverage is computed with respect to the 
	# merged ChIP-seq alignment of individual cell types 
	# (by combining all the replicates) 
	# for (cov_thr in c(0, 0.25, 0.5)) {
	for (cov_thr in c(0)) {
		# if (cov_thr == 0) {
		# 	P_E_SubDir <- paste0(P_E_Dir, '/No_Coverage_Condition')
		# } else {
		# 	P_E_SubDir <- paste0(P_E_Dir, '/Coverage_pct_', cov_thr)
		# }
		P_E_SubDir <- P_E_Dir
		system(paste("mkdir -p", P_E_SubDir))

		# cat(sprintf("\n\n ******** cov_thr : %s ", cov_thr))

		# get the P-E loop indices 
		# whose any segment does not overlap with the differential bins
		# and optionally, both segments maintain coverage constraint (on the merged ChIP-seq alignment)
		# also the loop indices should belong to the set "IsPeakStatus_Same_Idx"
		# if (cov_thr == 0) {
		Filt_P_E_Idx <- intersect(setdiff(seq(1,nrow(P_E_DiffLoop)), Loop_DiffChIPIdx), IsPeakStatus_Same_Idx)
		cat(sprintf("\n length Filt_P_E_Idx : %s  ", length(Filt_P_E_Idx)))
		# } else {
		# 	Cov_Fail_Idx <- which(((abs(P_E_DiffLoop[,19] - P_E_DiffLoop[,21]) / P_E_DiffLoop[,19]) > cov_thr) | ((abs(P_E_DiffLoop[,20] - P_E_DiffLoop[,22]) / P_E_DiffLoop[,20]) > cov_thr))
		# 	Filt_P_E_Idx <- intersect(setdiff(seq(1,nrow(P_E_DiffLoop)), union(Loop_DiffChIPIdx, Cov_Fail_Idx)), IsPeakStatus_Same_Idx)
		# 	cat(sprintf("\n length Cov_Fail_Idx : %s  ", length(Cov_Fail_Idx)))
		# 	cat(sprintf("\n length Filt_P_E_Idx : %s  ", length(Filt_P_E_Idx)))
		# }

		if (length(Filt_P_E_Idx) > 0) {

			# write the filtered loops (including annotations)
			P_E_DiffLoop_FiltFile <- paste0(P_E_SubDir, '/DiffLoops_Annot_Filt_ChIP1D.bed')
			write.table(P_E_DiffLoop[Filt_P_E_Idx, ], P_E_DiffLoop_FiltFile, row.names=FALSE, col.names=TRUE, sep = "\t", quote=FALSE, append=FALSE)

			# now from the filtered P-E loops 
			# get the promoter specific loops
			# and also list its expressions for the input categories
			Promoter_Specific_P_E_DiffLoop_FiltFile <- paste0(P_E_SubDir, '/DiffLoops_Annot_Filt_ChIP1D_Promoter_Specific.bed')
			if (GeneExprDataInp == TRUE) {
				Get_Promoter_Loops_FitHiChIP(P_E_DiffLoop_FiltFile, opt$InpTSSFile, Promoter_Specific_P_E_DiffLoop_FiltFile, CategoryList, 5000, GeneExprFileList, GeneNameColList, ExprValColList)
			} else {
				Get_Promoter_Loops_FitHiChIP(P_E_DiffLoop_FiltFile, opt$InpTSSFile, Promoter_Specific_P_E_DiffLoop_FiltFile, CategoryList, 5000)
			}

			InpData <- read.table(Promoter_Specific_P_E_DiffLoop_FiltFile, header=T, sep="\t", stringsAsFactors=F)

			cat(sprintf("\n nrow InpData : %s  ", nrow(InpData)))

			if (0) {
				# case 1 - write all the genes and associated information 
				LoopFile <- paste0(P_E_SubDir, '/Loops.bed')
				WritePromCentricLoops(InpData, LoopFile, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]), CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp)
			}

			# case 2 - now select only those genes 
			# such that they have only one loop
			# in the "Promoter_Specific_DiffLoop_Data" 
			# (containing gene information for all differential loops of EdgeR)

			# comment - sourya
			# ProcessUniqueGeneDiffLoops(CntGene_DiffLoopData_uniq, InpData, P_E_SubDir, CategoryList, ReplicaCount, CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2)

			# case 3 - now select only those genes 
			# such that they do not have exclusive loops in both categories
			ProcessDiffLoopsExclSingleCategory(Genes_Excl_SingleCat, InpData, P_E_SubDir, CategoryList, ReplicaCount, CCColList_Cat1, CCColList_Cat2, QValColList_Cat1, QValColList_Cat2, GeneExprDataInp)
		}

	}	# end coverage threshold loop

}	# end P-E loop condition

# now combine the P-P loops and P-E loops
# which are exclusive to one of the two categories
# into a single plot
# Note: as there are multiple discriminating conditions (and subfolders)
# we compute this information for individual subfolders (i.e. individual condition level)

# if ((length(P_P_Idx) > 0) & (length(P_E_Idx) > 0)) {

# 	for (cov_thr in c(0, 0.25, 0.5)) {
# 		if (cov_thr == 0) {
# 			P_P_BaseDir <- paste0(P_P_Dir_Main, '/Differential_3D_not_1D/No_Coverage_Condition')
# 			P_E_BaseDir <- paste0(P_E_Dir_Main, '/Differential_3D_not_1D/No_Coverage_Condition')
# 		} else {
# 			P_P_BaseDir <- paste0(P_P_Dir_Main, '/Differential_3D_not_1D/Coverage_pct_', cov_thr)
# 			P_E_BaseDir <- paste0(P_E_Dir_Main, '/Differential_3D_not_1D/Coverage_pct_', cov_thr)
# 		}

# 		CombineExclReplicateLoopsBoxPlot(P_P_BaseDir, P_E_BaseDir, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]))		

# 		P_P_CurrDir <- paste0(P_P_BaseDir, '/Genes_Count_1')
# 		P_E_CurrDir <- paste0(P_E_BaseDir, '/Genes_Count_1')
# 		CombineExclReplicateLoopsBoxPlot(P_P_CurrDir, P_E_CurrDir, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]))

# 		P_P_CurrDir <- paste0(P_P_BaseDir, '/Genes_Not_Excl_Both_Categories')
# 		P_E_CurrDir <- paste0(P_E_BaseDir, '/Genes_Not_Excl_Both_Categories')
# 		CombineExclReplicateLoopsBoxPlot(P_P_CurrDir, P_E_CurrDir, CategoryList, seq(1, ReplicaCount[1]), seq(1, ReplicaCount[2]))

# 	}	# end coverage threshold loop

# }

