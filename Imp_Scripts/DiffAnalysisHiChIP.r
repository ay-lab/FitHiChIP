#!/usr/bin/env Rscript

#===========================================================
# R script for differential analysis between two categories of FitHiChIP loops
# each with multiple replicates

# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
#===========================================================

# package to compute the overlap among intervals
suppressMessages(library(GenomicRanges))
library(optparse)
library(edgeR)
library(plyr)
library(dplyr)
library(ggplot2)
library(tools)
library(data.table)

options(scipen = 999)
options(datatable.fread.datatable=FALSE)

# ggplot parameters
FONTSIZE=18
PLOTWIDTH=10
PLOTHEIGHT=6

CHRLIST_NAMENUM <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")

NUMFEAT <- 2
NUMCOL_SEGMENT_UNIONLOOP <- 12
TSS_OFFSET <- 5000
THR_DEFAULT_BCV <- 0.4

#======================================================
# this function returns a vector (equal to the size of number of loops)
# which counts for each loop the number of replicates showing significant interaction
# parameters:
# IntDF: loops from all replicates - DESEQ input structure like 
# CCCols: columns storing raw contact count
# QValCols: columns storing q-values
# FDR_thr_Loop: FDR threshold for FitHiChIP loop significance (default = 0.01)
#======================================================
GetCntNumRepVec <- function(IntDF, CCCols, QValCols, FDR_thr_Loop=0.01) {	
	CntVec <- matrix(0, nrow=nrow(IntDF), ncol=1)
	for (i in (1:length(CCCols))) {
		# get the loop indices such that 
		# current replicate is significant
		curr_cc_col <- CCCols[i]
		curr_q_col <- QValCols[i]
		Idx <- which((IntDF[, curr_cc_col] > 0) & (IntDF[, curr_q_col] <= FDR_thr_Loop))
		# increment the counter of those loops
		CntVec[Idx] <- CntVec[Idx] + 1
	}
	return(CntVec)
}

#=================================
# function to extract the number of lines
#=================================
GetNumLines <- function(inpfile) {
	nline <- as.integer(system(paste("cat", inpfile, "| wc -l"), intern = TRUE))
	return(nline)
}

#================================
# this function finds overlapping loops
#================================
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

#================================
# this function extracts the chromosome data
#================================
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


#=================================
# function to compute overlap of 1D bins
#=================================
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


#======================================================
# applies EdgeR to the input count matrix
# returns a set of differential loops / bins 
# parameters:
# ALLLoopData: matrix or master data frame containing information for all loops / bins 
# MainDir: main base directory for EdgeR outputs
# CountData: only the count matrix employed for EdgeR design
# CategoryList: two input categories
# ReplicaCount: list of two elements, corresponding to the count of replicates in each categories
# FDR_Th_DESeq: FDR threshold for differential analysis
# FOLD_Change_Thr: fold change threshold for differential analysis
# prefixstr: if 'Loops', computes differential loops (3D), else computes differential (Hi)ChIP (1D)
# RawCCCols: list of columns storing the raw contact count
# FDR_thr_Loop: FDR threshold of loops (FitHiChIP)
# bcv: default threshold for single replicates data
# NormUse: if 1, calculates normalization factors prior to EdgeR design
#======================================================
ApplyEdgeR <- function(ALLLoopData, MainDir, CountData, CategoryList, ReplicaCount, FDR_Th_DESeq, FOLD_Change_Thr, prefixstr, RawCCCols=c(), QCols=c(), FDR_thr_Loop=0.01, bcv=0.4, NormUse=1) {

	# directory to store the plots related to EdgeR
	PlotDir <- paste0(MainDir, '/Plots_EdgeR')
	system(paste("mkdir -p", PlotDir))

	# category and replicate information for EdgeR structure
	GroupDistrVec <- c(rep(CategoryList[1], ReplicaCount[1]), rep(CategoryList[2], ReplicaCount[2]))
	cat(sprintf("\n ===>>> GroupDistrVec : %s ", paste0(GroupDistrVec, sep="\t")))
	
	# create the EdgeR count data structure (7th column onwards)
	y <- DGEList(counts=CountData, group=GroupDistrVec)

	if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {

		# for data with multiple replicates, better to normalize the DGEList object
		if (NormUse == 1) {
			y <- calcNormFactors(y)			
		}

		# generate the design (according to the group distribution)
		design <- model.matrix(~GroupDistrVec) 

		# Estimate Common, Trended and Tagwise Negative Binomial dispersions 
		# by weighted likelihood empirical Bayes
		y <- estimateDisp(y, design)

		# approach 1 - identifying DE tags - qCML method
		# exact test (applicable for experiments with single factor) 
		# to generate matrix of pseudo counts
		et <- exactTest(y, dispersion="trended", pair=c(CategoryList[1], CategoryList[2]))	
		# get the DE genes - useful for plotting MA 
		et_tags <- topTags(et)			

	} else {
		# default exact test when no replicates are available
		et <- exactTest(y, dispersion=bcv^2)
		# get the DE genes - useful for plotting MA 
		et_tags <- topTags(et)			
	}

	#===============
	if ((ReplicaCount[1] > 1) & (ReplicaCount[2] > 1)) {
	
		# statistics
		# using plotMDS function to generate a plot in which distances between samples 
		# correspond to leading biological coefficient of variation (BCV) between those samples:
		MDSPlotFile <- paste0(PlotDir, '/EdgeR_MDS_Plot.pdf')
		pdf(MDSPlotFile, width=8, height=6)
		plotMDS(y)
		dev.off()

		# statistics
		# dispersion estimates can be viewed in a BCV plot
		BCVPlotFile <- paste0(PlotDir, '/EdgeR_BCV_Plot.pdf')
		pdf(BCVPlotFile, width=8, height=6)
		plotBCV(y)
		dev.off()

		# Plot log-fold change against log-counts per million, 
		# with DE genes highlighted
		# using qlf (output of glmQLFTest)
		et_MDPLotFile <- paste0(PlotDir, '/EdgeR_et_MD_Plot.pdf')
		pdf(et_MDPLotFile, width=8, height=6)
		plotMD(et, hl.cex = 0.1)
		abline(h=c(-1, 1), col="blue")
		dev.off()		

		# smear plot - MA for groups of samples
		# approach 1 - exact test
		SmearPlotFile <- paste0(PlotDir, '/EdgeR_Smear_Plot_et.pdf')
		pdf(SmearPlotFile, width=8, height=6)
		plotSmear(et, pair=c(CategoryList[1], CategoryList[2]), de.tags=et_tags)
		dev.off()

	}	# end replica count condition

	#===============
	# for each of the DE tag finding techniques (qCML or CR methods)
	# convert the p-values to the corresponding q-values, using BH correction
	#===============
	Qval_et <- p.adjust(et$table$PValue, method = "BH")		# approach 1 - exact test

	#============================
	# for individual loops, get the number of replicates for both categories
	# which exhibit significant loops
	#============================
	if (prefixstr == 'Loops') {				
		OutLoop_SigCount_Vec_Category1 <- GetCntNumRepVec(ALLLoopData, RawCCCols[1:ReplicaCount[1]], QCols[1:ReplicaCount[1]], FDR_thr_Loop)
		OutLoop_SigCount_Vec_Category2 <- GetCntNumRepVec(ALLLoopData, RawCCCols[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], QCols[(ReplicaCount[1]+1):(ReplicaCount[1]+ReplicaCount[2])], FDR_thr_Loop)
	}

	#==========================
	# merge the EdgeR results with the master sheet data
	# and also generate the significant loops
	#==========================

	# approach 1 - exact test
	if (prefixstr == 'Loops') {
		EdgeRRes_et <- cbind.data.frame(ALLLoopData, et$table$logFC, et$table$logCPM, et$table$PValue, Qval_et, OutLoop_SigCount_Vec_Category1, OutLoop_SigCount_Vec_Category2)
		colnames(EdgeRRes_et) <- c(colnames(ALLLoopData), 'logFC', 'logCPM', 'PValue', 'FDR', paste0(CategoryList[1], '_SigRepl'), paste0(CategoryList[2], '_SigRepl'))
	} else {
		EdgeRRes_et <- cbind.data.frame(ALLLoopData, et$table$logFC, et$table$logCPM, et$table$PValue, Qval_et)
		colnames(EdgeRRes_et) <- c(colnames(ALLLoopData), 'logFC', 'logCPM', 'PValue', 'FDR')
	}
	if (prefixstr == 'Loops') {
		# for differential loops, consider both fold change and FDR conditions
		SigIdx <- which((EdgeRRes_et$FDR <= FDR_Th_DESeq) & (abs(EdgeRRes_et$logFC) >= FOLD_Change_Thr))
	} else {
		# for 1D bins, consider only FDR condition
		SigIdx <- which(EdgeRRes_et$FDR <= FDR_Th_DESeq)
	}
	LoopFile <- paste0(MainDir, '/', prefixstr, '_EdgeR_Default.bed')
	SigLoopFile <- paste0(MainDir, '/', prefixstr, '_EdgeR_Default_SIG.bed')
	write.table(EdgeRRes_et, LoopFile, row.names=F, col.names=T, quote=F, sep="\t")
	write.table(EdgeRRes_et[SigIdx, ], SigLoopFile, row.names=F, col.names=T, quote=F, sep="\t")

	if (prefixstr != 'Loops') {
		# for 1D data, derive the insignificant (non-differential) bins
		NonSigIdx <- which(EdgeRRes_et$FDR >= 0.1)
		NonSigLoopFile <- paste0(MainDir, '/', prefixstr, '_EdgeR_Default_NonSIG.bed')
		write.table(EdgeRRes_et[NonSigIdx, ], NonSigLoopFile, row.names=F, col.names=T, quote=F, sep="\t")
	}

}	# end EdgeR function

#================================
# this function annotates individual loops in the merged union set of loops
# according to the feature vectors of individual loops 
# parameters:
# UnionLoopFile: file containing merged set of loops from all input replicates and all categories
# AllLoopList: Loops with FitHiChIP significance for different input replicates and categories
# ChIPCovFileList: ChIP-seq coverage file list (one per input category) for the bin size of loops
# AllRepLabels: labels of individual replicates of individual categories
# CategoryList: two categories experimented
#================================
FillFeatureValues <- function(UnionLoopFile, AllLoopList, BinSize, ChIPCovFileList, AllRepLabels, CategoryList) {
	
	# counter of chromosomes processed
	valid_chr_count <- 0

	# temporary loop + feature containing file used per iteration
	temp_final_UnionLoopFile <- paste0(dirname(UnionLoopFile), '/temp_final_mastersheet.bed')

	UnionLoopTempFile1 <- paste0(dirname(UnionLoopFile), '/temp_CurrChr_Merged_Loops.bed')

	for (chr_idx in (1:length(CHRLIST_NAMENUM))) {
		chrName <- CHRLIST_NAMENUM[chr_idx]
		cat(sprintf("\n ===>>> Within function FillFeatureValues --- processing chromosome : %s ", chrName))

		# first extract the loops involving current chromosome		
		ExtractChrData(UnionLoopFile, chrName, UnionLoopTempFile1, header=FALSE)
		nreadCurr <- GetNumLines(UnionLoopTempFile1)
		if (nreadCurr == 0) {
			# check if there is no loop for the current chromosome
			# in such case, continue
			next
		}
		
		# otherwise, read the loops for the current chromosome
		# MergedIntTempData <- read.table(UnionLoopTempFile1, header=F)
		MergedIntTempData <- data.table::fread(UnionLoopTempFile1, header=F)

		# also get the interacting bins (start position divided by the bin size)
		AllLoop_BinDF <- cbind.data.frame((MergedIntTempData[,2] / BinSize), (MergedIntTempData[,5] / BinSize))
		colnames(AllLoop_BinDF) <- c('B1', 'B2')

		#=====================
		# for the current chromosome, process loop files for different replicates and categories
		# containing FitHiChIP significance
		#=====================
		
		# list of vectors for storing the feature values of individual loops
		# with respect to current chromosome
		RawCC_Categ <- list()
		QVal_Categ <- list()

		# stores temporary FitHiChIP loops for the current chromosome
		InpTempFitHiChIPLoopFile <- paste0(dirname(UnionLoopFile), '/temp_CurrChr_InpFitHiChIPLoopFile.bed')

		for (i in (1:length(AllLoopList))) {
			# default initialization for the current chromosome
			# and current input FitHiChIP loop file
			rawccvec <- rep(0, nrow(MergedIntTempData))
			qvec <- rep(1, nrow(MergedIntTempData))

			# input file containing FitHiChIP significance for all loops
			inpfile <- AllLoopList[i]

			# extract the interacting bins
			# and the raw contact count and the q-values
			if (tools::file_ext(inpfile) == "gz") {
				system(paste0("zcat ", inpfile, " | awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print ($2/b)\"\t\"($5/b)\"\t\"$7\"\t\"$NF}}\' - > ", InpTempFitHiChIPLoopFile))
			} else {
				system(paste0("awk -v b=", BinSize, " \'{if ((NR>1) && ($1==\"", chrName, "\")) {print ($2/b)\"\t\"($5/b)\"\t\"$7\"\t\"$NF}}\' ", inpfile, " > ", InpTempFitHiChIPLoopFile))		
			}

			# # extract the input data for the current chromosome			
			# ExtractChrData(inpfile, chrName, InpTempFitHiChIPLoopFile, header=TRUE)

			# check the number of loops for the current chromosome 
			# and for the current sample
			nreadInp <- GetNumLines(InpTempFitHiChIPLoopFile)
			if (nreadInp > 0) {
				# InpTempData <- read.table(InpTempFitHiChIPLoopFile, header=F, sep="\t", stringsAsFactors=F)
				InpTempData <- data.table::fread(InpTempFitHiChIPLoopFile, header=F, sep="\t", stringsAsFactors=F)
				colnames(InpTempData) <- c('B1', 'B2', 'RawCC', 'qval')
				cat(sprintf("\n ***** Computing overlap of merged loops with the FitHiChIP loop file index: %s  name: %s for the chromosome : %s ***** \n", i, inpfile, chrName))

				# get the matching row indices with respect to the master (union) set of loops
				# i.e. which overlap with the current input file
				m <- plyr::match_df(AllLoop_BinDF, InpTempData)
				Match_RowIdx <- as.integer(rownames(m))
				cat(sprintf("\n ----> number of master loops : %s  number of overlapping loops in input file : %s ", nrow(AllLoop_BinDF), length(Match_RowIdx)))

				# also merge the data frames
				mergeDF <- dplyr::inner_join(AllLoop_BinDF, InpTempData)

				# now store the feature values
				rawccvec[Match_RowIdx] <- mergeDF$RawCC
				qvec[Match_RowIdx] <- mergeDF$qval

				# # find the overlap between the merged loops and input sample FitHiChIP loops
				# CurrOv <- OverlapLoop(MergedIntTempData, InpTempData)

				# # dump the feature vectors
				# rawccvec[CurrOv$A_AND_B] <- InpTempData[CurrOv$B_AND_A, 7]				
				# qvec[CurrOv$A_AND_B] <- InpTempData[CurrOv$B_AND_A, ncol(InpTempData)]
			
			}	# end number of reads condition

			# assign the feature vectors for the current input file
			# into the final list of feature arrays
			RawCC_Categ[[i]] <- rawccvec
			QVal_Categ[[i]] <- qvec

			cat(sprintf("\n ***** Assigned contact count and q-values for the FitHiChIP loop file index: %s for the chromosome : %s ***** \n", i, chrName))

		}	# end loop FitHiChIP significance files

		# delete the temporary file
		if (file.exists(InpTempFitHiChIPLoopFile) == TRUE) {
			system(paste("rm", InpTempFitHiChIPLoopFile))
		}

		#=====================
		# process the ChIP-seq coverage files given for two input categories
		#=====================
  		# list of vectors for storing the features of individual loops
  		# with respect to one interacting segment, and for the current chromosome
  		Seg1_ChIP_Coverage <- list()
		Seg2_ChIP_Coverage <- list()

		Seg1_ChIP_Label <- rep('HD', nrow(MergedIntTempData))
		Seg2_ChIP_Label <- rep('HD', nrow(MergedIntTempData))

		# stores temporary ChIP coverage for the current chromosome
		InpTempChIPCoverageFile <- paste0(dirname(UnionLoopFile), '/temp_CurrChr_InpChIPCoverageFile.bed')

		for (i in (1:length(ChIPCovFileList))) {
			seg1_chip_coverage_vec <- rep(0, nrow(MergedIntTempData))
			seg2_chip_coverage_vec <- rep(0, nrow(MergedIntTempData))

			currcovfile <- ChIPCovFileList[i]
			cat(sprintf("\n Merging reference ChIP-seq coverage values: file number : %s  file name : %s ", i, currcovfile))

			# # extract the coverage information only for the current chromosome
			# # check the first field of chromosome name
			# # and second and third fields as integers
			# system(paste0("awk \'(($1==\"", chrName, "\") && ($2 ~ /^[0-9]+$/) && ($3 ~ /^[0-9]+$/))\' ", currcovfile, " > ", InpTempChIPCoverageFile))

			# extract the coverage information only for the current chromosome
			# check the first field of chromosome name
			# and second and third fields as integers
			# extract the bin number and the coverage information (assume 4th field)
			system(paste0("awk -v b=", BinSize, " \'{if (($1==\"", chrName, "\") && ($2 ~ /^[0-9]+$/) && ($3 ~ /^[0-9]+$/)) {print ($2/b)\"\t\"$4}}\' ", currcovfile, " > ", InpTempChIPCoverageFile))

			# read the ChIP coverage for the current sample and for the current chromosome
			# ChIPCoverageData <- read.table(InpTempChIPCoverageFile, header=F, sep="\t", stringsAsFactors=F)
			ChIPCoverageData <- data.table::fread(InpTempChIPCoverageFile, header=F, sep="\t", stringsAsFactors=F)

			# *** Note: the below mentioned MATCH operation is only possible 
			# since we perform exact overlap of the segments
				
			# get the overlap of the first interacting bin of AllLoop_BinDF
			# with the bins in ChIPCoverageData
			m <- match(AllLoop_BinDF[,1], ChIPCoverageData[,1])
			idx_Loop <- which(!is.na(m))
			idx_Cov <- m[idx_Loop]
			# copy in the segment 1 ChIP coverage vector
			seg1_chip_coverage_vec[idx_Loop] <- ChIPCoverageData[idx_Cov, 2]
			if (i == 1) {
				Seg1_ChIP_Label[idx_Loop] <- ChIPCoverageData[idx_Cov, ncol(ChIPCoverageData)]
			}

			# get the overlap of the second interacting bin of AllLoop_BinDF
			# with the bins in ChIPCoverageData
			m <- match(AllLoop_BinDF[,2], ChIPCoverageData[,1])
			idx_Loop <- which(!is.na(m))
			idx_Cov <- m[idx_Loop]
			# copy in the segment 2 ChIP coverage vector
			seg2_chip_coverage_vec[idx_Loop] <- ChIPCoverageData[idx_Cov, 2]
			if (i == 1) {
				Seg2_ChIP_Label[idx_Loop] <- ChIPCoverageData[idx_Cov, ncol(ChIPCoverageData)]
			}

			# # overlap of the current set of merged loops (for the current chromosome)
			# # and the 1D segments (for the current chromosome)
			# ov1 <- Overlap1D(MergedIntTempData[,1:3], ChIPCoverageData[,1:3], boundary=1, offset=0, uniqov=FALSE)
			# ov2 <- Overlap1D(MergedIntTempData[,4:6], ChIPCoverageData[,1:3], boundary=1, offset=0, uniqov=FALSE)

			# # assign the ChIP coverage of overlapping bins 
			# # we assume that 4th field of ChIP coverage file stores the coverage information
			# seg1_chip_coverage_vec[ov1$A_AND_B] <- ChIPCoverageData[ov1$B_AND_A, 4]
			# seg2_chip_coverage_vec[ov2$A_AND_B] <- ChIPCoverageData[ov2$B_AND_A, 4]

			# assign the feature vectors for the current input file
			# into the final list of feature arrays
			Seg1_ChIP_Coverage[[i]] <- seg1_chip_coverage_vec
			Seg2_ChIP_Coverage[[i]] <- seg2_chip_coverage_vec

			# # assign the 1D label as well
			# if (i == 1) {
			# 	Seg1_ChIP_Label[ov1$A_AND_B] <- ChIPCoverageData[ov1$B_AND_A, ncol(ChIPCoverageData)]
			# 	Seg2_ChIP_Label[ov2$A_AND_B] <- ChIPCoverageData[ov2$B_AND_A, ncol(ChIPCoverageData)]
			# }

			cat(sprintf("\n ***** Assigned reference ChIP-seq coverage information for the file index: %s for the chromosome : %s ***** \n", i, chrName))
		}

		# delete the temporary file
		if (file.exists(InpTempChIPCoverageFile) == TRUE) {
			system(paste("rm", InpTempChIPCoverageFile))
		}

		#=======================
		# now for the current chromosome
		# merge the interacting regions along with the features accumulated
		#=======================
		# chromosome counter
		valid_chr_count <- valid_chr_count + 1

		namesvec <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

		# insert the ChIP-seq coverage information
		for (i in (1:2)) {
			MergedIntTempData <- cbind.data.frame(MergedIntTempData, Seg1_ChIP_Coverage[[i]], Seg2_ChIP_Coverage[[i]])
		}

		# insert the ChIP coverage differential labels (HD / LD / ND)
		MergedIntTempData <- cbind.data.frame(MergedIntTempData, Seg1_ChIP_Label, Seg2_ChIP_Label)

		# update the column labels as well
		namesvec <- c(namesvec, paste0(CategoryList[1], c('_ChIPCov1', '_ChIPCov2')), paste0(CategoryList[2], c('_ChIPCov1', '_ChIPCov2')), 'Bin1_Label', 'Bin2_Label')

		# insert for individual samples (replicates within categories)
		# different feature values
		# note that this is for a single chromosome
		for (i in (1:length(AllLoopList))) { 
			MergedIntTempData <- cbind.data.frame(MergedIntTempData, RawCC_Categ[[i]], QVal_Categ[[i]])
		}

		# column names
		appendnamevec <- c()
		for (i in (1:length(AllRepLabels))) {
			for (v in c('_RawCC', '_QVal')) {
				appendnamevec <- c(appendnamevec, paste0(AllRepLabels[i], v))
			}
		}
		namesvec <- c(namesvec, appendnamevec)

		# assign the column names into the data frame created for the current chromosome
		colnames(MergedIntTempData) <- namesvec

		# now write this data frame in the temporary master sheet file
		# sourya - should be modified - 
		# we need to insert a boolean variable to see if chromosome specific writing is enabled or not
		# if (chr_idx == 1) {
		if (valid_chr_count == 1) {
			write.table(MergedIntTempData, temp_final_UnionLoopFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		} else {
			write.table(MergedIntTempData, temp_final_UnionLoopFile, row.names=F, col.names=F, sep="\t", quote=F, append=T)
		}

	} 	# end chromosome index loop

	# rename the temporary master sheet file in the final file
	system(paste("mv", temp_final_UnionLoopFile, UnionLoopFile))

	if (file.exists(UnionLoopTempFile1) == TRUE) {
		system(paste("rm", UnionLoopTempFile1))
	}

}	# end function

#===========================================================
option_list = list(

	make_option(c("--AllLoopList"), type="character", default=NULL, help="Comma or colon separated list of FitHiChIP loops (without any significance thresholding) of all categories and all replicates. Check the file names $PREFIX$*.interactions.FitHiC.bed in the FitHiChIP output files. Loop files for the first category are to be mentioned first, followed by those of second category. For example, suppose there are two replicates for each category. Then, the file list would be like Cat1_Repl1_File:Cat1_Repl2_File:Cat2_Repl1_File:Cat2_Repl2_File where, Cat1 denotes input category 1, and Repl1 means replicate 1, and so on. Files can be gzipped as well. Mandatory parameter."),

	make_option(c("--ChrSizeFile"), type="character", default=NULL, help="File containing size of chromosomes for reference genome. Mandatory parameter."),

	make_option(c("--FDRThr"), type="numeric", default=0.01, help="FDR threshold used for determining FitHiChIP significant loops. Generally, 0.01 q-value threshold (default) is used."),	

	make_option(c("--CovThr"), type="integer", action="store", default=25, help="Threshold signifying the max allowed deviation of ChIP coverage between two categories, to consider those bins as ND (i.e. no difference). Default is 25, means 25% deviation of ChIP coverage is set as maximum, for a bin to be classified as ND. If user chooses 50, 50% maximum ChIP seq coverage deviation would be allowed."),

	make_option(c("--ChIPAlignFileList"), type="character", default=NULL, help="Comma or colon separated list of ChIP-seq files, either alignment (BAM) format, or in bedgraph (4 column) format. Default is NULL. User can 1) either provide two files, one for each input category, 2) or provide ChIP seq alignment / coverage files one for each input sample (i.e. for all input replicates). Mandatory parameter."),

	make_option(c("--OutDir"), type="character", default=NULL, help="Base Output directory. Mandatory parameter."),

	#==================
	make_option(c("--CategoryList"), type="character", default=NULL, help="Comma or colon separated list of strings, corresponding to the two main categories (whose replicates are present). Default: Category1, Category2."),

	make_option(c("--ReplicaCount"), type="character", default=NULL, help="Comma or colon separated list of the count of replicates for individual categories. Default: 1,1 (means that we are considering one replicate per sample)."),
	
	make_option(c("--ReplicaLabels1"), type="character", default=NULL, help="Comma or colon separated list of the label of replicates for the first category. Default: R1,R2, etc."),
	
	make_option(c("--ReplicaLabels2"), type="character", default=NULL, help="Comma or colon separated list of the label of replicates for the second category. Default: R1,R2, etc."),

	#==================
	make_option(c("--FoldChangeThr"), type="integer", action="store", default=2, help="EdgeR fold change threshold - log2 of this value is used. Default = 2, means that log2(2) = 1 is used as the fold change threshold in EdgeR output."),
	
  	make_option(c("--DiffFDRThr"), type="numeric", default=0.05, help="FDR threshold for EdgeR. Default is 0.05, means that loops with FDR < 0.05, and fold change >= log2(FoldChangeThr) would be considered as differential."),

  	make_option(c("--bcv"), type="numeric", default=0.4, help="If EdgeR is used with single samples (replica count = 1 for any of the categories), this value is the square-root-dispersion. For datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms, or 0.01 for technical replicates. For details, see the edgeR manual. By default, the value is set as 0.4.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#====================================
# process the input arguments
#====================================

# parse the FitHiChIP loop file list (containing loop significance, not constrained by FDR threshold)
if (is.null(opt$AllLoopList)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! FiHiChIP loop files (along with loop significance) are not provided - check the option --AllLoopList \n", call.=FALSE)
} else {
	AllLoopList <- as.character(unlist(strsplit(opt$AllLoopList,"[,:]")))
}

if (is.null(opt$ChrSizeFile)) {
	stop("ERROR !!!!!!! File containing chromosome sizes corresponding to the reference genome are not provided - check the option --ChrSizeFile \n", call.=FALSE)
}

if (!is.null(opt$ChIPAlignFileList)) {
	ChIPAlignFileList <- as.character(unlist(strsplit(opt$ChIPAlignFileList,"[,:]")))
	if ((length(ChIPAlignFileList) != 2) & (length(ChIPAlignFileList) != length(AllLoopList))) {
		stop("ERROR !!!!!!! Number of ChIP-seq alignment files in the --ChIPAlignFileList parameter should be either 2 (one for each category) or should be equal to the number of samples !! quit !!! \n", call.=FALSE)
	}
} else {
	print_help(opt_parser)
	stop("ERROR !!!!!!! Option --ChIPAlignFileList is empty. User should provide either 2 (one for each category) ChIP seq alignment files, or the number of files should be equal to the number of samples. Exit !!! \n", call.=FALSE)
}

if (is.null(opt$OutDir)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! Output directory is not provided - check the option --OutDir \n", call.=FALSE)
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
	stop("ERROR !!!!!!! Count of the replicates does not match with the number of input loops files - check the option --ReplicaCount \n", call.=FALSE)
}

if (is.null(opt$ReplicaLabels1)) {
	ReplicaLabels1 <- paste0('R', seq(1:ReplicaCount[1]))
} else {
	ReplicaLabels1 <- as.character(unlist(strsplit(opt$ReplicaLabels1,"[,:]")))
	if (length(ReplicaLabels1) != ReplicaCount[1]) {
		stop("ERROR !!!!!!! Number of labels for the replicates in the first category does not match with the number of replicates provided for the first category - check the options --ReplicaCount and --ReplicaLabels1 \n", call.=FALSE)
	}
}

if (is.null(opt$ReplicaLabels2)) {
	ReplicaLabels2 <- paste0('R', seq(1:ReplicaCount[2]))
} else {
	ReplicaLabels2 <- as.character(unlist(strsplit(opt$ReplicaLabels2,"[,:]")))
	if (length(ReplicaLabels2) != ReplicaCount[2]) {
		stop("ERROR !!!!!!! Number of labels for the replicates in the first category does not match with the number of replicates provided for the first category - check the options --ReplicaCount and --ReplicaLabels2 \n", call.=FALSE)
	}
}

if ((opt$CovThr <= 0) | (opt$CovThr > 100)) {
	stop("ERROR !!!!!!! The parameter CovThr should be a positive number between 1 to 100 - check !!! \n", call.=FALSE)
}

#======================
# assign other parameters
#======================
FDR_Th_FitHiChIP <- as.numeric(opt$FDRThr)	# FitHiChIP q-value threshold
FOLD_Change_Thr <- log2(opt$FoldChangeThr) 	# as.integer(opt$FoldChangeThr)
FDR_Th_DESeq <- as.numeric(opt$DiffFDRThr)
ChIP_Cov_Thr <- (opt$CovThr * 1.0) / 100	# max allowed deviation for ChIP seq coverage

# combining all the replica labels
AllRepLabels <- c(paste0(CategoryList[1], '_', ReplicaLabels1), paste0(CategoryList[2], '_', ReplicaLabels2))

# create the base output directory
system(paste("mkdir -p", opt$OutDir))

# count the number of ChIP-seq alignment files
if (length(ChIPAlignFileList) == 2) {
	ChIPAlignFileCountVec <- c(1,1)
} else {
	ChIPAlignFileCountVec <- ReplicaCount
}

#=========================
# get the bin size
#=========================
if (tools::file_ext(AllLoopList[1]) == "gz") {
	BinSize <- as.integer(system(paste("zcat", AllLoopList[1], " | awk \'{if (NR==2) {print ($3-$2)}}\' - "), intern = TRUE))
} else {
	BinSize <- as.integer(system(paste("awk \'{if (NR==2) {print ($3-$2)}}\'", AllLoopList[1]), intern = TRUE))	
}
cat(sprintf("\n **** Bin size of FitHiChIP loops : %s \n", BinSize))

#=========================
# read both input ChIP coverage files and scale + combine those coverage values
# consider only chromosomes common to both categories
#=========================
MergedChIPCovFile <- paste0(opt$OutDir, '/ChIP_Coverage_ALL.bed')

tempfile1 <- paste0(opt$OutDir, '/ChIP_Coverage_temp1.bed')
tempfile2 <- paste0(opt$OutDir, '/ChIP_Coverage_temp2.bed')

# first get the binned distribution of specified chromosome size file
TargetBinnedChrFile <- paste0(opt$OutDir, '/Binned_Chromosome_Size_RefGenome.bed')
system(paste("bedtools makewindows -g", opt$ChrSizeFile, "-w", BinSize, ">", TargetBinnedChrFile))

# process individual alignment files and first get the ChIP seq coverage for those files
# using the TargetBinnedChrFile
for (i in (1:length(ChIPAlignFileList))) {
	if (tools::file_ext(ChIPAlignFileList[i]) == "bam") {
		cat(sprintf("\n ===>> Merging ChIP coverage - processing ChIP-seq alignment file in BAM format : %s ====>> \n", ChIPAlignFileList[i]))
	} else if (tools::file_ext(ChIPAlignFileList[i]) == "gz") {
		cat(sprintf("\n ===>> Merging ChIP coverage - processing ChIP-seq alignment file in gzipped bedgraph (4 column) format : %s ====>> \n", ChIPAlignFileList[i]))
	} else {
		cat(sprintf("\n ===>> Merging ChIP coverage - processing ChIP-seq alignment file in bedgraph (4 column) format : %s ====>> \n", ChIPAlignFileList[i]))
	}

	if (tools::file_ext(ChIPAlignFileList[i]) == "gz") {
		# gzipped bedgraph format
		if (i == 1) {			
			system(paste("zcat", ChIPAlignFileList[i], "| bedtools coverage -a", TargetBinnedChrFile, "-b stdin -sorted -counts >", MergedChIPCovFile))
		} else {
			system(paste("zcat", ChIPAlignFileList[i], "| bedtools coverage -a", TargetBinnedChrFile, "-b stdin -sorted -counts | cut -f4 >", tempfile1))
			system(paste("paste", MergedChIPCovFile, tempfile1, ">", tempfile2))
			system(paste("mv", tempfile2, MergedChIPCovFile))			
		}
	} else {
		# either BAM file or plain bedgraph format
		if (i == 1) {			
			system(paste("bedtools coverage -a", TargetBinnedChrFile, "-b", ChIPAlignFileList[i], "-sorted -counts >", MergedChIPCovFile))
		} else {			
			system(paste("bedtools coverage -a", TargetBinnedChrFile, "-b", ChIPAlignFileList[i], "-sorted -counts | cut -f4 >", tempfile1))
			system(paste("paste", MergedChIPCovFile, tempfile1, ">", tempfile2))
			system(paste("mv", tempfile2, MergedChIPCovFile))
		}		
	}
}	# end input file loop

# now scale the ChIP coverages according to different categories
# Merged_ChIPCovData <- read.table(MergedChIPCovFile, header=F, sep="\t", stringsAsFactors=F)
Merged_ChIPCovData <- data.table::fread(MergedChIPCovFile, header=F, sep="\t", stringsAsFactors=F)

# ChIP coverage of the first category - row means operation
if (ChIPAlignFileCountVec[1] > 1) {
	ChIPCovCat1 <- rowMeans(Merged_ChIPCovData[, 4:(4+ChIPAlignFileCountVec[1]-1)])
} else {
	ChIPCovCat1 <- Merged_ChIPCovData[, 4]
}

# ChIP coverage of the second category - row means operation
if (ChIPAlignFileCountVec[2] > 1) {
	ChIPCovCat2 <- rowMeans(Merged_ChIPCovData[, (4+ChIPAlignFileCountVec[1]):(4+ChIPAlignFileCountVec[1]+ChIPAlignFileCountVec[2]-1)])
} else {
	ChIPCovCat2 <- Merged_ChIPCovData[, (4+ChIPAlignFileCountVec[1])]	
}

# scaling factor to be applied on the coverage values of first category
scaling_factor <- (sum(ChIPCovCat2) * 1.0) / sum(ChIPCovCat1)

cat(sprintf("\n ===>> Scaling ChIP coverage - scaling_factor : %s ====>> \n", scaling_factor))

# scale the ChIP-seq coverage of the first category
for (i in (1:ChIPAlignFileCountVec[1])) { 
	Merged_ChIPCovData[, (4+i-1)] <- as.integer(round(Merged_ChIPCovData[, (4+i-1)] * scaling_factor))
}

# now write the scaled coverage values
colnames(Merged_ChIPCovData) <- c('chr', 'start', 'end', paste0(CategoryList[1], '_ChIPCov_', seq(1,ChIPAlignFileCountVec[1])), paste0(CategoryList[2], '_ChIPCov_', seq(1,ChIPAlignFileCountVec[2])))
write.table(Merged_ChIPCovData, MergedChIPCovFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

# now remove the temporary files
if (file.exists(tempfile1) == TRUE) {
	system(paste("rm", tempfile1))
}
if (file.exists(tempfile2) == TRUE) {
	system(paste("rm", tempfile2))
}	

cat(sprintf("\n\n *** Performed scaling of ChIP coverage values for uniform ChIP coverage *** \n\n"))

#========================
# now apply EdgeR on the scaled ChIP coverage values of both categories
# to get the non-differential 1D bins
#========================

ChIPEdgeRDir <- paste0(opt$OutDir, '/ChIP_Coverage_EdgeR')
system(paste("mkdir -p", ChIPEdgeRDir))

# construct the count matrix for applying to EdgeR
# consists of bin indices, midpoints and their coverages (scaled) for both categories
CountData_1D <- cbind.data.frame(seq(1,nrow(Merged_ChIPCovData)), Merged_ChIPCovData[,1], (Merged_ChIPCovData[,2] + Merged_ChIPCovData[,3])/2, Merged_ChIPCovData[,4:ncol(Merged_ChIPCovData)])
colnames(CountData_1D) <- c("Idx", "chr", "mid", colnames(Merged_ChIPCovData)[4:ncol(Merged_ChIPCovData)])
CountDataColNames_1D <- colnames(CountData_1D)

# dump the count matrix 
CountDataFile_1D <- paste0(ChIPEdgeRDir, '/count_matrix_1D.bed')
write.table(CountData_1D, CountDataFile_1D, row.names=F, col.names=T, quote=F, sep="\t")

# apply the EdgeR routine on the 1D bins
# Note: here we supply the parameter ChIPAlignFileCountVec
# which contains the number of ChIP-seq alignment files provided
# for coverage estimation 
ApplyEdgeR(Merged_ChIPCovData, ChIPEdgeRDir, CountData_1D[, 4:ncol(CountData_1D)], CategoryList, ChIPAlignFileCountVec, FDR_Th_DESeq, FOLD_Change_Thr, 'ChIP_Coverage')

# file containing non-differential ChIP-seq bins (EdgeR)
ChIP_1D_EdgeR_NonSigFile <- paste0(ChIPEdgeRDir, '/ChIP_Coverage_EdgeR_Default_NonSIG.bed')
# ChIP_1D_EdgeR_NonSigData <- read.table(ChIP_1D_EdgeR_NonSigFile, header=T, sep="\t", stringsAsFactors=F)
ChIP_1D_EdgeR_NonSigData <- data.table::fread(ChIP_1D_EdgeR_NonSigFile, header=T, sep="\t", stringsAsFactors=F)

cat(sprintf("\n\n *** Performed EdgeR on ChIP coverage *** \n\n"))

# get the 1D bins having very low deviation in their ChIP-seq coverage
# subject to the input parameter threshold 
MaxCovVec <- pmax(Merged_ChIPCovData[,4], Merged_ChIPCovData[,5])
MinCovVec <- pmin(Merged_ChIPCovData[,4], Merged_ChIPCovData[,5])
BinIdxLowCovDiff <- which(MinCovVec > ((1 - ChIP_Cov_Thr) * MaxCovVec))

# label the ChIP bins as HD, LD or ND
# depending on their EdgeR output, and ChIP coverage comparison
ov <- Overlap1D(Merged_ChIPCovData[, 1:3], ChIP_1D_EdgeR_NonSigData[, 1:3])

ND_Label_Idx <- intersect(ov$A_AND_B, BinIdxLowCovDiff)
LD_Label_Idx <- setdiff(ov$A_AND_B, ND_Label_Idx)
HD_Label_Idx <- setdiff(setdiff(seq(1, nrow(Merged_ChIPCovData)), LD_Label_Idx), ND_Label_Idx)
cat(sprintf("\n Number of 1D bins: %s \n number of 1D bins classified as HD: %s \n number of 1D bins classified as LD: %s \n number of 1D bins classified as ND: %s ", nrow(Merged_ChIPCovData), length(HD_Label_Idx), length(LD_Label_Idx), length(ND_Label_Idx)))

# assign the 1D label information
ChIPLabelVec <- rep('HD', nrow(Merged_ChIPCovData))
ChIPLabelVec[ND_Label_Idx] <- 'ND'
ChIPLabelVec[LD_Label_Idx] <- 'LD'

# append the 1D label in the merged ChIP coverage data
# and write back to the merged ChIP coverage file
CN <- colnames(Merged_ChIPCovData)
Merged_ChIPCovData <- cbind.data.frame(Merged_ChIPCovData, ChIPLabelVec)
colnames(Merged_ChIPCovData) <- c(CN, 'Label')
write.table(Merged_ChIPCovData, MergedChIPCovFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

#==========================
# prepare the union of input loops
#==========================
UnionLoopFile <- paste0(opt$OutDir, '/MasterSheet_', CategoryList[1], '_', CategoryList[2], '_Loops.bed')

# if (file.exists(UnionLoopFile) == FALSE) {

	tempLoopFile1 <- paste0(opt$OutDir, '/tempLoop1.bed')
	 
	for (idx in (1:length(AllLoopList))) {
		CurrLoopFile <- AllLoopList[idx]		
		cat(sprintf("\n Creating master set of interactions - processing the loop file idx : %s  file name : %s ", idx, CurrLoopFile))
		if (tools::file_ext(CurrLoopFile) == "gz") {
			if (idx == 1) {
				system(paste0("zcat ", CurrLoopFile, " | awk \'(NR>1)\' - | cut -f1-6 > ", tempLoopFile1))
			} else {
				system(paste0("zcat ", CurrLoopFile, " | awk \'(NR>1)\' - | cut -f1-6 >> ", tempLoopFile1))
			}
		} else {
			if (idx == 1) {
				system(paste0("awk \'(NR>1)\' ", CurrLoopFile, " | cut -f1-6 > ", tempLoopFile1))
			} else {
				system(paste0("awk \'(NR>1)\' ", CurrLoopFile, " | cut -f1-6 >> ", tempLoopFile1))
			}			
		}
	}
	system(paste("sort -k1,1 -k2,2n -k5,5n", tempLoopFile1, " | uniq >", UnionLoopFile))
	system(paste("rm", tempLoopFile1))

	#===================
	# for individual loops, fill the features from input HiChIP loops
	# for different replicates for input categories
	# like contact count, bias, significance, isPeak
	# Note: here we create two new files for storing the scaled ChIP coverage values
	# and send it as arguments
	# when multiple ChIP-seq alignment files are provided as input
	# we use the mean of the ChIP-seq coverages
	#===================
	scaled_ChIPCovFile1 <- paste0(opt$OutDir, '/scaled_ChIP_Coverage_', CategoryList[1], '.bed')
	scaled_ChIPCovFile2 <- paste0(opt$OutDir, '/scaled_ChIP_Coverage_', CategoryList[2], '.bed')

	# mean scaled ChIP coverage of the first category
	if (ChIPAlignFileCountVec[1] > 1) {
		meanScaledChIPCovCat1 <- rowMeans(Merged_ChIPCovData[, 4:(4+ChIPAlignFileCountVec[1]-1)])
	} else {
		meanScaledChIPCovCat1 <- Merged_ChIPCovData[, 4]
	}
	# mean scaled ChIP coverage of the second category
	if (ChIPAlignFileCountVec[2] > 1) {
		meanScaledChIPCovCat2 <- rowMeans(Merged_ChIPCovData[, (4+ChIPAlignFileCountVec[1]):(4+ChIPAlignFileCountVec[1]+ChIPAlignFileCountVec[2]-1)])
	} else {
		meanScaledChIPCovCat2 <- Merged_ChIPCovData[, (4+ChIPAlignFileCountVec[1])]
	}

	# create two data frames using the scaled mean of ChIP coverage values
	covdf1 <- cbind.data.frame(Merged_ChIPCovData[, 1:3], meanScaledChIPCovCat1, Merged_ChIPCovData[, ncol(Merged_ChIPCovData)])
	colnames(covdf1) <- c(colnames(Merged_ChIPCovData)[1:3], paste0(CategoryList[1], '_meanChIPCov'), 'Label')

	covdf2 <- cbind.data.frame(Merged_ChIPCovData[, 1:3], meanScaledChIPCovCat2, Merged_ChIPCovData[, ncol(Merged_ChIPCovData)])
	colnames(covdf2) <- c(colnames(Merged_ChIPCovData)[1:3], paste0(CategoryList[2], '_meanChIPCov'), 'Label')

	write.table(covdf1, scaled_ChIPCovFile1, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	write.table(covdf2, scaled_ChIPCovFile2, row.names=F, col.names=T, sep="\t", quote=F, append=F)

	# now call the feature value function with the new set of ChIP coverage files
	FillFeatureValues(UnionLoopFile, AllLoopList, BinSize, c(scaled_ChIPCovFile1, scaled_ChIPCovFile2), AllRepLabels, CategoryList)

# }	# end if

cat(sprintf("\n\n *** Created master sheet of loops *** \n\n"))


#=========================
# paste the configuration options within this execution
#=========================
fp_out <- file(paste0(opt$OutDir, '/Input_Parameters_', gsub("-","_",gsub(" ","_",format(Sys.time(), "%F %H-%M"))), '.log'), "w")

outtext <- paste0("\n ********* Listing all the input parameters ****** ")
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> All FitHiChIP input loop files (without significant thresholding) : \n ", paste(AllLoopList, collapse="\n"))
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> Threshold of ChIP-seq coverage difference: ", opt$CovThr, " percent \n")
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> FDR threshold of FitHiChIP loops: ", FDR_Th_FitHiChIP)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> Files containing ChIP alignment information : \n ", paste(ChIPAlignFileList, collapse="\n"))
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> Output directory: ", opt$OutDir)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> Input categories : \t ", paste(CategoryList, collapse="\t"))
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> Count of replicates for first category : ", ReplicaCount[1], " and for the second category : ", ReplicaCount[2])
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> Labels associated with replicates of category 1 : \t ", paste(ReplicaLabels1, collapse="\t"))
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> Labels associated with replicates of category 2 : \t ", paste(ReplicaLabels2, collapse="\t"))
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> FDR threshold for differential analysis (DESeq / EdgeR) : ", FDR_Th_DESeq)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> Fold change threshold (logarithm base 2 of the input) for differential analysis (DESeq / EdgeR) : ", FOLD_Change_Thr)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n ===>>> bcv option (constant) for differential analysis (specific to EdgeR) : ", opt$bcv)
writeLines(outtext, con=fp_out, sep="\n")

# close the output text file
close(fp_out)

#===================
# read the master sheet for all loops  and their features
# then apply EdgeR based differential loop analysis
#===================
# MasterSheetData <- read.table(UnionLoopFile, header=T, sep="\t", stringsAsFactors=F)
MasterSheetData <- data.table::fread(UnionLoopFile, header=T, sep="\t", stringsAsFactors=F)

# re-initialize the column names
CN <- colnames(MasterSheetData)
colnames(MasterSheetData) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", CN[7:ncol(MasterSheetData)])

DiffLoopDir <- paste0(opt$OutDir, '/EdgeR_Loops_ALL')
system(paste("mkdir -p", DiffLoopDir))

# list of columns in the master sheet data storing the raw contact counts
RawCC_ColList <- NUMCOL_SEGMENT_UNIONLOOP + ((seq(1, length(AllRepLabels)) - 1) * NUMFEAT) + 1
QVal_ColList <- RawCC_ColList + 1

# file which would contain the EdgeR based differential loops
EdgeROutFile <- paste0(DiffLoopDir, '/Loops_EdgeR_Default_SIG.bed')

# raw count data (RawCC)
CountData <- MasterSheetData[, RawCC_ColList]

CountDataColNames <- colnames(CountData)

# first write a text file which will contain the interacting segments 
# along with the index of those interactions
# required for defining EdgeR count matrix
IntervalTextFile <- paste0(DiffLoopDir, '/Interacting_segments.bed')
system(paste("awk \'{if (NR>1) {print $1\"\t\"$2\"\t\"$3; print $4\"\t\"$5\"\t\"$6}}\'",  UnionLoopFile, "| sort -k1,1 -k2,2n | uniq | awk \'{print NR\"\t\"$0}\' - >", IntervalTextFile))

# file containing the count data structure for use in EdgeR
CountDataFile <- paste0(DiffLoopDir, '/count_matrix.bed')

# read the interacting bins
# IntervalData <- read.table(IntervalTextFile, header=F)
IntervalData <- data.table::fread(IntervalTextFile, header=F)
colnames(IntervalData) <- c("Idx1", "chr1", "start1", "end1")

# prepare a copy of this data frame
# used for joining operation
IntervalData2 <- data.frame(IntervalData)
colnames(IntervalData2) <- c("Idx2", "chr2", "start2", "end2")

# join the set of loops together with these bin index data frames
finaldf <- dplyr::inner_join(MasterSheetData, IntervalData) %>% dplyr::inner_join(IntervalData2)

# prepare the count data (for EdgeR structure)
# containing the indices as well
CountData <- cbind.data.frame(finaldf$Idx1, finaldf$Idx2, finaldf$chr1, ((finaldf$start1 + finaldf$end1) / 2), finaldf$chr2, ((finaldf$start2 + finaldf$end2) / 2), CountData)
colnames(CountData) <- c("Idx1", "Idx2", "chr1", "mid1", "chr2", "mid2", CountDataColNames)
CountDataColNames <- colnames(CountData)

# dump the count matrix 
write.table(CountData, CountDataFile, row.names=F, col.names=T, quote=F, sep="\t")

# apply EdgeR using the generated count data
# Note: we send the count matrix only
ApplyEdgeR(MasterSheetData, DiffLoopDir, CountData[, 7:ncol(CountData)], CategoryList, ReplicaCount, FDR_Th_DESeq, FOLD_Change_Thr, 'Loops', RawCC_ColList, QVal_ColList, FDR_Th_FitHiChIP)

cat(sprintf("\n\n *** Applied EdgeR for differential loop finding *** \n\n"))

#============================
# now find out the differential loops involving categories ND-ND, LD-ND, LD-LD, HD-LD/ND, and HD-HD
# for all the EdgeR based differential loops
#============================
categoryDir <- paste0(DiffLoopDir, '/ND_ND')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || (($11==\"ND\") && ($12==\"ND\"))\' ", EdgeROutFile, ">", allLoop_file))
system(paste("awk \'(NR==1) || ($(NF-5)<=0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || ($(NF-5)>0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[2], '.bed')))


categoryDir <- paste0(DiffLoopDir, '/LD_ND')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || ((($11==\"LD\") && ($12==\"ND\")) || (($11==\"ND\") && ($12==\"LD\")))\' ", EdgeROutFile, ">", allLoop_file))
system(paste("awk \'(NR==1) || ($(NF-5)<=0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || ($(NF-5)>0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[2], '.bed')))


categoryDir <- paste0(DiffLoopDir, '/LD_LD')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || (($11==\"LD\") && ($12==\"LD\"))\' ", EdgeROutFile, ">", allLoop_file))
system(paste("awk \'(NR==1) || ($(NF-5)<=0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || ($(NF-5)>0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[2], '.bed')))


categoryDir <- paste0(DiffLoopDir, '/HD_HD')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || (($11==\"HD\") && ($12==\"HD\"))\' ", EdgeROutFile, ">", allLoop_file))
system(paste("awk \'(NR==1) || ($(NF-5)<=0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || ($(NF-5)>0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[2], '.bed')))


categoryDir <- paste0(DiffLoopDir, '/HD_LD_or_ND')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || ((($11==\"HD\") && ($12!=\"HD\")) || (($11!=\"HD\") && ($12==\"HD\")))\' ", EdgeROutFile, ">", allLoop_file))
system(paste("awk \'(NR==1) || ($(NF-5)<=0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || ($(NF-5)>0)\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Up_', CategoryList[2], '.bed')))



#============================
# now find out the differential loops overlapping with FitHiChIP loops
# significant in at least one replicate
# involving categories ND-ND, LD-ND, LD-LD, HD-LD/ND, and HD-HD
# for all the EdgeR based differential loops
#============================
DiffLoopDir <- paste0(opt$OutDir, '/EdgeR_Loops_Ov_FitHiChIP_Sig_One_Repl')
system(paste("mkdir -p", DiffLoopDir))

OneSample_SigRepl_Loop_file <- paste0(categoryDir, '/DiffLoops_ALL.bed')
system(paste("awk \'((NR==1) || ($(NF-1)>0) || ($NF>0))\' ", EdgeROutFile, ">", OneSample_SigRepl_Loop_file))

categoryDir <- paste0(DiffLoopDir, '/ND_ND')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || (($11==\"ND\") && ($12==\"ND\"))\' ", OneSample_SigRepl_Loop_file, ">", allLoop_file))
system(paste("awk \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))


categoryDir <- paste0(DiffLoopDir, '/LD_ND')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || ((($11==\"LD\") && ($12==\"ND\")) || (($11==\"ND\") && ($12==\"LD\")))\' ", OneSample_SigRepl_Loop_file, ">", allLoop_file))
system(paste("awk \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))


categoryDir <- paste0(DiffLoopDir, '/LD_LD')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || (($11==\"LD\") && ($12==\"LD\"))\' ", OneSample_SigRepl_Loop_file, ">", allLoop_file))
system(paste("awk \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))


categoryDir <- paste0(DiffLoopDir, '/HD_HD')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || (($11==\"HD\") && ($12==\"HD\"))\' ", OneSample_SigRepl_Loop_file, ">", allLoop_file))
system(paste("awk \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))


categoryDir <- paste0(DiffLoopDir, '/HD_LD_or_ND')
system(paste("mkdir -p", categoryDir))
allLoop_file <- paste0(categoryDir, '/DiffLoops.bed')
system(paste("awk \'(NR==1) || ((($11==\"HD\") && ($12!=\"HD\")) || (($11!=\"HD\") && ($12==\"HD\")))\' ", OneSample_SigRepl_Loop_file, ">", allLoop_file))
system(paste("awk \'(NR==1) || (($(NF-1)>0) && ($NF==0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[1], '.bed')))
system(paste("awk \'(NR==1) || (($(NF-1)==0) && ($NF>0))\' ", allLoop_file, ">", paste0(categoryDir, '/DiffLoops_Excl_', CategoryList[2], '.bed')))

