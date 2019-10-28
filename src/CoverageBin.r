#!/usr/bin/env Rscript

library(optparse)
library(GenomicRanges)
library(tools)
library(data.table)

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
  	make_option(c("--InpFile"), type="character", default=NULL, help="Input file with locus pairs and their contact count."),
  	make_option(c("--PeakFile"), type="character", default=NULL, help="Reference ChIP-seq or HiChIP peak file."),
  	make_option(c("--BinSize"), type="integer", action="store", default=5000, help="Bin size employed. Default 5000 (5 Kb)."),
  	make_option(c("--ChrSizeFile"), type="character", default=NULL, help="File containing chromosome sizes corresponding to the reference genome."),
  	make_option(c("--OutFile"), type="character", default=NULL, help="Output file storing 1D bin specific HiChIP coverage along with their overlap with input peak segments.")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

OutDir <- dirname(opt$OutFile)

# counter of valid chromosomes (containing one or more interactions)
valid_chr_count <- 0

# first use the reference chromosome size file and input bin size
# to extract all distinct 1D bins
temp_1D_binfile <- paste0(OutDir, '/temp_1D_ALLBins.bed')
system(paste("bedtools makewindows -g", opt$ChrSizeFile, "-w", opt$BinSize, ">", temp_1D_binfile))

# read the chromosomes of input ChrSizeFile
ChrData <- read.table(opt$ChrSizeFile, header=F, sep="\t", stringsAsFactors=F)
ChrNames <- sort(ChrData[,1])

# files for processing one chromosome
temp_1D_binfile_currChr <- paste0(OutDir, '/temp_1D_ALLBins_CurrChr.bed')

for (i in (1:length(ChrNames))) {
	currChr <- ChrNames[i]
	cat(sprintf("\n Processing chromosome for coverage computation: %s ", currChr))

	# get the bins for the current chromosome
	system(paste0("awk \'{if ($1==\"", currChr, "\") {print $0}}\' ", temp_1D_binfile, " > ", temp_1D_binfile_currChr))

	nline <- as.integer(system(paste("cat", temp_1D_binfile_currChr, "| wc -l"), intern = TRUE))
	if (nline == 0) {
		next
	}
	
	# now extract the HiChIP 1D coverage from the input locus pair and contact count file
	# and for the current chromosome
	temp_1D_Coveragefile <- paste0(OutDir, '/temp_1D_coverage.bed')
	system(paste0("awk -F\'[\t]\' -v OFS=\"\t\" \'{if ((NR>1) && ($1==\"", currChr, "\")) {a[$1\"\t\"$2\"\t\"$3]+=$7; a[$4\"\t\"$5\"\t\"$6]+=$7}} END {for (x in a) {print x,a[x]}}\' ", opt$InpFile, " | sort -k1,1 -k2,2n > ", temp_1D_Coveragefile))

	nline <- as.integer(system(paste("cat", temp_1D_Coveragefile, "| wc -l"), intern = TRUE))
	if (nline == 0) {
		next
	}

	# read the bins for the current chromosome
	temp_1D_bindata <- data.table::fread(temp_1D_binfile_currChr, header=F, sep="\t", stringsAsFactors=F)	
	
	# read the HiChIP coverage values for the current chromosome
	temp_1D_CoverageData <- data.table::fread(temp_1D_Coveragefile, header=F, sep="\t", stringsAsFactors=F)

	# append the coverage values of temp_1D_CoverageData 
	# into the bins in temp_1D_bindata
	coverageVec <- rep(0, nrow(temp_1D_bindata))
	ov <- Overlap1D(temp_1D_bindata, temp_1D_CoverageData, uniqov=FALSE)
	coverageVec[ov$A_AND_B] <- temp_1D_CoverageData[ov$B_AND_A, 4]

	# assign the peak information in the 1D bins
	peakvec <- rep(0, nrow(temp_1D_bindata))

	if (!is.null(opt$PeakFile)) {	
		# now read the peak file - supporting both gzipped and normal files
		temp_peakfile <- paste0(OutDir, '/temp_peaks.bed')
		if (tools::file_ext(opt$PeakFile) == "gz") {
			system(paste0("zcat ", opt$PeakFile, " | awk \'{if (($1==\"", currChr, "\") && ($2 ~ /^[0-9]+$/) && ($3 ~ /^[0-9]+$/)) {print $1\"\t\"$2\"\t\"$3}}\' - > ", temp_peakfile))
		} else {
			system(paste0("cat ", opt$PeakFile, " | awk \'{if (($1==\"", currChr, "\") && ($2 ~ /^[0-9]+$/) && ($3 ~ /^[0-9]+$/)) {print $1\"\t\"$2\"\t\"$3}}\' - > ", temp_peakfile))
		}
		temp_PeakData <- data.table::fread(temp_peakfile, header=F, sep="\t", stringsAsFactors=F)
		if (nrow(temp_PeakData) > 0) {
			ov1 <- Overlap1D(temp_1D_bindata, temp_PeakData, boundary=0, offset=0, uniqov=FALSE)
			peakvec[ov1$A_AND_B] <- 1		
		}
	}	# end peak file non-null condition

	# now construct the final data frame
	FinalDF <- cbind.data.frame(temp_1D_bindata, coverageVec, peakvec)
	colnames(FinalDF) <- c('Chr', 'Start', 'End', 'Coverage', 'IsPeak')

	# increment the valid chromosome counter
	valid_chr_count <- valid_chr_count + 1
	if (valid_chr_count == 1) {
		write.table(FinalDF, opt$OutFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	} else {
		write.table(FinalDF, opt$OutFile, row.names=F, col.names=F, sep="\t", quote=F, append=T)
	}

}	# end chromosome loop

# remove temporary files
system(paste("rm", temp_1D_binfile))
system(paste("rm", temp_1D_binfile_currChr))
system(paste("rm", temp_1D_Coveragefile))
if (!is.null(opt$PeakFile)) {
	system(paste("rm", temp_peakfile))
}

