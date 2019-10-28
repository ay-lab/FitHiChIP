#!/usr/bin/env Rscript

library(optparse)
library(data.table)
suppressMessages(library(GenomicRanges))

options(scipen = 999)
options(datatable.fread.datatable=FALSE)

#======================
# script to compute 1D bin specific HiChIP coverage
# and also their overlap with reference ChIP-seq / HiChIP peak files
#======================

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

#====================================================
option_list = list(
	make_option(c("--BinSize"), type="integer", action="store", default=5000, help="Bin size employed. Default 5000 (5 Kb)."),
	make_option(c("--ChrSizeFile"), type="character", default=NULL, help="File containing chromosome sizes corresponding to the reference genome."),
	make_option(c("--InpFile"), type="character", default=NULL, help="Input file with locus pairs and their contact count."),
	make_option(c("--Interval"), type="character", default=NULL, help="File to contain bin intervals."),
	make_option(c("--Matrix"), type="character", default=NULL, help="Input file storing matrix.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# create bin interval file
system(paste("mkdir -p", dirname(opt$Interval)))
TempBinnedChrFile <- paste0(dirname(opt$Interval), '/TempBinnedChrFile.bed')
system(paste("bedtools makewindows -g", opt$ChrSizeFile, "-w", opt$BinSize, ">", TempBinnedChrFile))
system(paste0("awk \'{print $0\"\t\"NR}\' ", TempBinnedChrFile, " > ", opt$Interval))

ChrSizeData <- read.table(opt$ChrSizeFile, header=F, sep="\t", stringsAsFactors=F)
ChrNames <- as.vector(ChrSizeData[,1])

CurrChrIntervalFile <- paste0(dirname(opt$Interval), '/TempBinnedInterval_CurrChr.bed')
CurrChrLoopFile <- paste0(dirname(opt$Interval), '/TempInpLoop_CurrChr.bed')

valid_chr_count <- 0

for (i in (1:length(ChrNames))) {
	CurrChr <- ChrNames[i]
	# export interval for the current chromosome
	system(paste0("awk \'{if ($1==\"", CurrChr, "\") {print $0}}\' ", opt$Interval, " > ", CurrChrIntervalFile))
	nBin <- as.integer(system(paste("cat", CurrChrIntervalFile, "| wc -l"), intern = TRUE))
	if (nBin == 0) {
		next
	}

	# export interactions for the current chromosome
	system(paste0("awk \'{if (($1==\"", CurrChr, "\") && ($4==\"", CurrChr, "\")) {print $0}}\' ", opt$InpFile, " > ", CurrChrLoopFile))

	nLoop <- as.integer(system(paste("cat", CurrChrLoopFile, "| wc -l"), intern = TRUE))
	cat(sprintf("\n Creating ICE compatible bin pair specific matrix file - processing chromosome : %s number of loops : %s ", CurrChr, nLoop))
	if (nLoop == 0) {
		next
	}

	# otherwise, process the loops
	valid_chr_count <- valid_chr_count + 1

	CurrChrIntervalData <- data.table::fread(CurrChrIntervalFile)
	CurrChrLoopData <- data.table::fread(CurrChrLoopFile)

	# define the output data frame
	outDF_CurrChr <- matrix(0, nrow=nrow(CurrChrLoopData), ncol=3)
	outDF_CurrChr[, 3] <- CurrChrLoopData[, ncol(CurrChrLoopData)]

	ov <- Overlap1D(CurrChrLoopData[,1:3], CurrChrIntervalData[,1:3], uniqov=FALSE)
	outDF_CurrChr[ov$A_AND_B, 1] <- CurrChrIntervalData[ov$B_AND_A, ncol(CurrChrIntervalData)]

	ov <- Overlap1D(CurrChrLoopData[,4:6], CurrChrIntervalData[,1:3], uniqov=FALSE)
	outDF_CurrChr[ov$A_AND_B, 2] <- CurrChrIntervalData[ov$B_AND_A, ncol(CurrChrIntervalData)]

	if (valid_chr_count == 1) {
		write.table(outDF_CurrChr, opt$Matrix, row.names=F, col.names=F, sep="\t", quote=F, append=F)
	} else {
		write.table(outDF_CurrChr, opt$Matrix, row.names=F, col.names=F, sep="\t", quote=F, append=T)
	}

}	# end chromosome loop

if (file.exists(TempBinnedChrFile) == TRUE) {
	system(paste("rm", TempBinnedChrFile))
}
if (file.exists(CurrChrIntervalFile) == TRUE) {
	system(paste("rm", CurrChrIntervalFile))
}
if (file.exists(CurrChrLoopFile) == TRUE) {
	system(paste("rm", CurrChrLoopFile))
}
