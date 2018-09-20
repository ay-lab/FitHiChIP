#!/usr/bin/env Rscript

#===========================================================
# R script for producing bin specific coverage, ispeak, and bias information
# if the bias file is provided externally (ICE specific bias), it is appended
# otherwise, coverage specific bias is computed, and bias is separated for peaks and non-peaks

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================
library(optparse)

option_list = list(
  	make_option(c("--CoverageFile"), type="character", default=NULL, help="File containing the Coverage information for individual bins"),
  	make_option(c("--BiasFile"), type="character", default=NULL, help="If specified, contains the ICE bias for individual bins"),
	make_option(c("--OutFile"), type="character", default=NULL, help="Output file containing coverage + bias information")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read the genome coverage features
# Note: the coverage file has header information
# columns: chromosome interval, read coverage, is_peak
CoverageFeat <- read.table(opt$CoverageFile, header=T, sep="\t", stringsAsFactors=F)	
# colnames(CoverageFeat) <- c("chr1","s1","e1","coverage","isPeak")

# if bias file is not provided (no external ICE bias)
# then compute the bias from the coverage information
# separately process the peaks and non peaks
if (is.null(opt$BiasFile)) {

	# create a bias vector of the same size as the number of 
	# genomic intervals
	BiasVec <- matrix(0, nrow(CoverageFeat), 1)

	# find the indices corresponding to peak and non peaks
	PeakIdx <- which(CoverageFeat[,5] == 1)
	NonPeakIdx <- which(CoverageFeat[,5] == 0)

	# note the indices having peak information and non zero coverage
	PeakIDXNonZeroCoverage <- intersect(which(CoverageFeat[,4]> 0), which(CoverageFeat[,5]==1))

	# note the indices of non peak and non zero coverage
	nonPeakIDXNonZeroCoverage <- intersect(which(CoverageFeat[,4]> 0), which(CoverageFeat[,5]==0))

	# process the bias associated with peak segments
	if (length(PeakIDXNonZeroCoverage) > 0) {
		
		# coverage values (non zero) of the peak segments
		CoveragePeakVec <- CoverageFeat[PeakIDXNonZeroCoverage, 4]
		
		# mean value of these coverages
		meanCoverage <- mean(CoveragePeakVec)

		# copy the mean divided coverage values in the appropriate indices of the bias vector
		for (x in PeakIdx) {
			BiasVec[x] <- CoverageFeat[x, 4] / meanCoverage
		}
	}

	# process the bias associated with non-peak segments
	if (length(nonPeakIDXNonZeroCoverage) > 0) {

		# coverage values (non zero) of the non peak segments
		CoverageNonPeakVec <- CoverageFeat[nonPeakIDXNonZeroCoverage, 4]

		# mean value of these coverages
		meanCoverage <- mean(CoverageNonPeakVec)

		# copy the mean divided coverage values in the appropriate indices of the bias vector
		for (x in NonPeakIdx) {
			BiasVec[x] <- CoverageFeat[x, 4] / meanCoverage
		}
	}

	# now append the bias values with the original coverage statistics
	# and overwrite the existing coverage file
	OutCoverage.df <- cbind(CoverageFeat, BiasVec)
	colnames(OutCoverage.df) <- c(colnames(CoverageFeat), 'Bias')
	write.table(OutCoverage.df, opt$OutFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE) 

} else {

	# first read the external bias information (without any header)
	BiasFeat <- read.table(opt$BiasFile, header=F, sep="\t", stringsAsFactors=F)	
	colnames(BiasFeat) <- c("Chr","Start","End","Bias")

	# merge with respect to the 1st three fields of either data (chromosome interval)
	OutCoverage.df <- merge(x=CoverageFeat, y=BiasFeat, by.x=colnames(CoverageFeat)[1:3], by.y=colnames(BiasFeat)[1:3])
	colnames(OutCoverage.df) <- c(colnames(CoverageFeat), "Bias")
	
	# write the data frame in a temporary file
	# since the data frame may not be sorted
	temp_outfile <- paste0(opt$OutFile, '_temp')
	write.table(OutCoverage.df, temp_outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE) 

	# sort the output data frame with respect to first three columns
	# and write it in the final output file
	system(paste('sort -k1,1 -k2,2n', temp_outfile, '>', opt$OutFile))
	system(paste('rm', temp_outfile))

}


