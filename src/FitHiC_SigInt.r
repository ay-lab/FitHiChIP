#!/usr/bin/env Rscript

#===========================================================
# R code for FitHiChIP implementation 
# Uses spine fit (equal occupancy binning) as implemented in FitHiC paper Ay et. al. Genome Research 2014 
# then uses novel bias regression mechanism
# supports all of the following parameters:
# 1) zero contact count based locus pair use / only non zero contact based  locus pair use 
# 2) Using either all background or only peak to peak background
# 3) Specification of pre-filtering interactions subject to bias thresholds
# 4) Bias correction / no bias correction
# 		A) For bias correction, either multiply the probabilities or
# 		B) Use regression model
# 			B.1) Residual = 0, and Equal occupancy binning is either 0 or 1

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================

library(splines)
library(fdrtool)
library(parallel)	# library for parallel processing
library(optparse)
library(ggplot2)
library(data.table)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

#====================================================
# function to compute the number of possible interacting pairs (including zero contacts)
# given a particular genomic distance
# applicable for interactions involving peak segments
# parameters:
# 1) InpDistSet: set of input genomic distance values for any locus pairs
# 2) binsize
# 3) ChrSpecBin.df: chromosome specific bins and peak information
AllPossibleContacts <- function(InpDistSet, binsize, ChrSpecBin.df, IntType) {	

	# list of all chromosomes
	# ChrList_NameNum <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")
	ChrList_NameNum <- as.vector(unique(ChrSpecBin.df[,1]))
	# cat(sprintf("\n ===>>> within function AllPossibleContacts -- checking all possible locus pairs for these chromosomes : %s ", paste(ChrList_NameNum, collapse=" ")))

	# this vector stores number of locus pairs for individual genomic distance
	CountLocusPairDist <- rep(0, length(InpDistSet))

	# compute all possible locus pairs for individual chromosomes
	res_chrCntLocPair <- as.data.frame(parallel:::mclapply( 1:length(ChrList_NameNum) , mc.cores = ncore , function(chr_idx){

		# current chromosome name
		chrName <- ChrList_NameNum[chr_idx]
		if (0) {
			cat(sprintf("\n AllPossibleContacts - processing chromosome: %s ", chrName))	
		}		
		# vector stores the number of locus pairs for individual genomic distance 
		ChrCountLocusPairDist <- rep(0, length(InpDistSet))
		# extract the peak bins for this chromosome
		PeakBin <- ChrSpecBin.df[which((ChrSpecBin.df[,1] == chrName) & (ChrSpecBin.df[,3] == 1)), 2]
		# extract the non-peak bins for this chromosome
		NonPeakBin <- ChrSpecBin.df[which((ChrSpecBin.df[,1] == chrName) & (ChrSpecBin.df[,3] == 0)), 2]
		if (0) {
			cat(sprintf("\n Num PeakBin : %s Num NonPeakBin : %s ", length(PeakBin), length(NonPeakBin)))
		}

		# now find the locus pairs having the particular distance
		# and of the current interaction type
		Sort_PeakBin <- sort(PeakBin)
		Sort_NonPeakBin <- sort(NonPeakBin)
		if (length(Sort_PeakBin) > 1) {
			# peak to peak or peak to all interactions
			if ((IntType == 1) | (IntType == 3)) {
				# process individual distance values
				for (dist_idx in (1:length(InpDistSet))) {
					inpdist <- InpDistSet[dist_idx]
					nbins_thr <- (inpdist / binsize)
					# for peak to peak interactions, the other peak 
					# should have a distance of "nbins_thr" from the first peak (+ve distance)
					other_peak_bin <- Sort_PeakBin + nbins_thr
					# check which of these other end bins are actually peak bins
					# those correspond to the peak to peak interactions
					ChrCountLocusPairDist[dist_idx] <- ChrCountLocusPairDist[dist_idx] + length(intersect(other_peak_bin, Sort_PeakBin))
				}	# end distance loop
			}	# end peak to peak or peak to all interactions
		}	# end condition - presence of peak bins

		if ((length(Sort_PeakBin) > 0) & (length(Sort_NonPeakBin) > 0)) {
			# peak to non-peak or peak to all interactions
			if ((IntType == 2) | (IntType == 3)) {
				# process individual distance values
				for (dist_idx in (1:length(InpDistSet))) {
					inpdist <- InpDistSet[dist_idx]
					nbins_thr <- (inpdist / binsize)
					# for peak to peak interactions, the other peak 
					# should have a distance of "nbins_thr" from the first peak 
					# first explore the +ve distance
					other_peak_bin <- Sort_PeakBin + nbins_thr
					# check which of these other end bins are actually non-peak bins
					# those correspond to the peak to non-peak interactions
					ChrCountLocusPairDist[dist_idx] <- ChrCountLocusPairDist[dist_idx] + length(intersect(other_peak_bin, Sort_NonPeakBin))
					# then explore the -ve distance
					other_peak_bin <- Sort_PeakBin - nbins_thr
					# check which of these other end bins are actually non-peak bins
					# those correspond to the peak to non-peak interactions
					ChrCountLocusPairDist[dist_idx] <- ChrCountLocusPairDist[dist_idx] + length(intersect(other_peak_bin, Sort_NonPeakBin))
				}	# end distance loop
			}	# end peak to non peak or peak to all interactions
		}	# end condition - presence of peak and non-peak bins
		
		# we return the chromosome index 
		# and the distance specific peak and non peak count distribution
		# for this chromosome
		return(c(chr_idx,ChrCountLocusPairDist))
	} ))

	# now accumulate all the chromosome specific count (for individual distance)
	# in the final vector
	for (dist_idx in (1:length(InpDistSet))) {
		CountLocusPairDist[dist_idx] <- sum(as.integer(res_chrCntLocPair[dist_idx + 1, 1:ncol(res_chrCntLocPair)]))
	}

	return(CountLocusPairDist)
}

#====================================================
# function to compute the number of possible interacting pairs
# (including zero contacts)
# given a particular genomic distance
# applicable for all to all interactions
AllPossbleContactsGeneDist <- function(inpdist, binsize, MappableBinCount.df) {
	nbins_thr <- (inpdist / binsize)
	sum_val <- 0
	for (i in 1:nrow(MappableBinCount.df)) {
		curr_chr <- MappableBinCount.df[i,1]
		bincnt_curr_chr <- MappableBinCount.df[i,2]
		sum_val <- sum_val + (bincnt_curr_chr - nbins_thr)
	}
	return(sum_val)
}

#====================================================
option_list = list(
  	make_option(c("--InpFile"), type="character", default=NULL, help="File with interactions + normalization features among individual genomic bins, preferably sorted with respect to interaction distance"),
  	make_option(c("--headerInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input interaction file has the header. Default FALSE."),
	make_option(c("--OutFile"), type="character", default=NULL, help="Output file name storing interactions + probability, P value and Q value"),
	make_option(c("--CoverageFile"), type="character", default=NULL, help="Input file storing the chromosome bin and coverage values"),
	make_option(c("--BinSize"), type="integer", action="store", default=5000, help="Bin size employed. Default 5000 (5 Kb)."),
	
	make_option(c("--IntType"), type="integer", action="store", default=4, help="Can be 1 (peak to peak), 2 (peak to non peak), 3 (peak to all), or 4 (all to all). Default 4."),

	make_option(c("--UseNonzeroContacts"), type="integer", action="store", default=0, help="Can be 1 or 0. If 1, uses only non zero contacts for spline fitting. Else, uses all locus pairs (having zero contacts also). Default 0."),

	make_option(c("--P2P"), type="integer", action="store", default=0, help="Can be 0 or 1. If 1, spline fit is modeled using only peak to peak interactions. Default 0."),
	
  	make_option(c("--BiasCorr"), type="integer", action="store", default=0, help="If 0, raw contact count is used. Else, bias corrected contact count is used before applying FitHiC. Default 0."),
  	
  	make_option(c("--BiasType"), type="integer", action="store", default=1, help="If 1, coverage bias is used. If 2, ICE bias is used. Default 1."),  	

	make_option(c("--BiasLowThr"), type="numeric", default=0.2, help="Lower threshold of bias. Default 0.2. Used to model the bias specific regression, by filtering the training data with respect to this bias value. Optionally used for pre-filtering the interactions."),
	make_option(c("--BiasHighThr"), type="numeric", default=5, help="Higher threshold of bias. Default 5. Used to model the bias specific regression, by filtering the training data with respect to this bias value. Optionally used for pre-filtering the interactions."),  	
  	make_option(c("--nbins"), type="integer", action="store", default=200, help="Number of bins employed for FitHiC.", metavar="nbins"),
  	make_option(c("--Draw"), type="logical", action="store_true", default=FALSE, help="If TRUE, spline fit of FitHiC is plotted. Default FALSE."),
  	
  	make_option(c("--cccol"), type="integer", action="store", default=7, help="Column number of the file storing the contact count"),
  	make_option(c("--TimeFile"), type="character", default=NULL, help="If specified, denotes the file which will contain time profiling information"),

  	make_option(c("--BiasFilt"), type="integer", action="store", default=0, help="If 1, interactions are pre-filtered according to BiasLowThr and BiasHighThr, before processing. Default 0"),
	make_option(c("--MultBias"), type="integer", action="store", default=0, help="If 1, probability values are multiplied with the bias values. Default 0."),

	make_option(c("--Resid"), type="integer", action="store", default=0, help="If 1, and bias correction is enabled, models the residual contact count for regression. Default 0 (tested best performance)."),
	make_option(c("--EqOcc"), type="integer", action="store", default=0, help="If 1, and bias correction is enabled, models the regression with respect to equal occupancy binning. Else, equal distance based regression is modeled. Default 0 (tested best performance).")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$TimeFile)) {
	timeprof <- 0
} else {
	timeprof <- 1
}

# number of bins used in FitHiC
nbins <- as.integer(opt$nbins) 

# peak to peak background usage option
Peak2PeakBackg <- as.integer(opt$P2P)
cat(sprintf("\n Peak to peak background usage for spline fit: %s ", Peak2PeakBackg))

# input interaction type
IntType <- as.integer(opt$IntType)

# number of processors within the system
ncore <- detectCores()
cat(sprintf("\n Number of cores in the system: %s ", ncore))

# directory containing the input interaction file
inpdir <- dirname(opt$InpFile)

# directory to contain the spline fitted output interaction file
OutIntDir <- dirname(opt$OutFile)

# this directory will store the spline graph according to the model
outdir <- paste0(OutIntDir,'/Plots')
system(paste('mkdir -p', outdir))

if (timeprof == 1) {
	starttime <- Sys.time()
}

# list of parameters
if (0) {
	cat(sprintf("\n\n ******* list of parameters - \n BinSize : %s \n IntType : %s \n UseNonzeroContacts : %s \n P2P : %s \n BiasCorr : %s \n BiasType : %s \n BiasLowThr : %s \n BiasHighThr : %s \n nbins : %s \n BiasFilt : %s \n MultBias : %s \n Resid : %s \n EqOcc : %s \n\n", opt$BinSize, opt$IntType, opt$UseNonzeroContacts, opt$P2P, opt$BiasCorr, opt$BiasType, opt$BiasLowThr, opt$BiasHighThr, opt$nbins, opt$BiasFilt, opt$MultBias, opt$Resid, opt$EqOcc))
}

#========================================================
# applicable if the spline fitting is performed with all possible locus pairs (including zero contacts)
# input: coverage file 
# case 1: if all-to-all interactions are considered (IntType == 4), all bins are counted
# here the output: chromosome specific bin count - a table with 2 columns: chromosome, no_of_bins
# case 2: for other interaction types (involving peaks)
# output: chromosome specific bin and peak information
# also check for bias fltering, bias thresholds

if (opt$UseNonzeroContacts == 0) {
	MappableBinCountChrFile <- paste0(OutIntDir, '/MappableBinCountChr.bed')
	if (IntType == 4) {
		# all to all interactions
		if ((opt$BiasCorr == 1) & (opt$BiasFilt == 1)) {
			system(paste0("awk -F\'[\t]\' -v l=", opt$BiasLowThr, " -v h=", opt$BiasHighThr, " \'((NR>1) && ($6>=l) && ($6<=h))\' ", opt$CoverageFile, " | cut -f1 | uniq -c | awk \'{print $2\"\t\"$1}\' - > ", MappableBinCountChrFile))
		} else {
			system(paste0("awk -F\'[\t]\' \'(NR>1)\' ", opt$CoverageFile, " | cut -f1 | uniq -c | awk \'{print $2\"\t\"$1}\' - > ", MappableBinCountChrFile))
		}
		MappableBinCount.df <- data.table::fread(MappableBinCountChrFile, header=F, sep="\t", stringsAsFactors=F)
	} else {
		# dump individual chromosomes, their bins and corresponding peak information
		if ((opt$BiasCorr == 1) & (opt$BiasFilt == 1)) {
			system(paste0("awk -F\'[\t]\' -v b=", opt$BinSize, " -v l=", opt$BiasLowThr, " -v h=", opt$BiasHighThr, " \'{if ((NR>1) && ($6>=l) && ($6<=h)) {print $1\"\t\"($2/b)\"\t\"$5}}\' ", opt$CoverageFile, " > ", MappableBinCountChrFile))
		} else {
			system(paste0("awk -F\'[\t]\' -v b=", opt$BinSize, " \'{if (NR>1) {print $1\"\t\"($2/b)\"\t\"$5}}\' ", opt$CoverageFile, " > ", MappableBinCountChrFile))
		}
		# read the chromosome specific bin information
		ChrSpecBin.df <- data.table::fread(MappableBinCountChrFile, header=F, sep="\t", stringsAsFactors=F)
	}
}
#========================================================
# load the interaction data matrix (its header information is provided in opt$headerInp)
interaction.data <- data.table::fread(opt$InpFile, header=opt$headerInp, sep="\t", stringsAsFactors=F)
if (opt$headerInp == FALSE) {
	colnames(interaction.data) <- c("chr1", "s1", "e1", "chr2", "s2", "e2", "cc", "d1", "isPeak1", "Bias1", "mapp1", "gc1", "cut1", "d2", "isPeak2", "Bias2", "mapp2", "gc2", "cut2")
}

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time to load the interaction data file in the R structure: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}

cat(sprintf("\n ===>> Total Number of input interactions (locus pairs): %s ", nrow(interaction.data)))
# return if the number of interaction is 0
if (nrow(interaction.data) == 0) {
	cat(sprintf("\n ****** No of input interactions : 0 - exit !!!"))
	return()
}
if (ncol(interaction.data) < opt$cccol) {
	cat(sprintf("\n ****** There is no contact count column or the formatting of the input file has problems - exit !!!"))
	return()
}

# vector of initial (prior) probabilities for contact counts
Prob_Val <- c()
# probability of the observed contact count for a single locus pair, from the spline fit
Spline_Binom_Prob_CC <- c()
# P-value of the observed contact count for a single locus pair
# binomial distribution + spline based estimation
Spline_Binom_P_Val_CC <- c()
# for each bin, count the no of distinct locus pairs (having nonzero contacts)
no_distinct_loci <- c()
# for each bin, count the no of all possible locus pairs (including zero contacts)
no_possible_pair_loci <- c()
# for each bin, no of observed contact count
NumContact <- c()
# for each bin, stores the average contact count per interaction
avg_contact <- c()
# prior contact probability for a specific locus pair falling within a bin
prior_contact_prob <- c()
# average interaction distance for all locus pairs falling within a bin
avg_int_dist <- c()
# specifying the columns in the interaction data 
# where the bias information is provided
bias1.col <- 10
bias2.col <- 16
# few other vectors specifically employed for bias correction
Prob_BiasRegr <- c()
Exp_CC_BiasRegr <- c()

#======================================================
# if bias specific filtering is enabled, 
# discard the interactions where either loci has bias values outside specified thresholds
if ((opt$BiasCorr == 1) & (opt$BiasFilt == 1)) {
	if (ncol(interaction.data) >= bias2.col) {
		interaction.data <- interaction.data[which((interaction.data[,bias1.col] >= opt$BiasLowThr) & (interaction.data[,bias1.col] <= opt$BiasHighThr) & (interaction.data[,bias2.col] >= opt$BiasLowThr) & (interaction.data[,bias2.col] <= opt$BiasHighThr)), ]
		cat(sprintf("\n ===>> Bias specific filtering is enabled - Number of interactions where both loci satisfy bias criterion (within low and high limits): %s ", nrow(interaction.data)))
	}
}

#======================================================
# absolute genomic distance for an interaction instance
gene.dist <- abs(interaction.data[,2] - interaction.data[,5])

# no of interactions pairs
numPairs <- length(gene.dist)

# total number of contacts for all the interactions
TotContact <- sum(interaction.data[,opt$cccol])

if (timeprof == 1) {
	starttime <- Sys.time()
}

# check if the interaction file is not sorted - return
if (is.unsorted(gene.dist) == TRUE) {
	cat(sprintf("\n **** The interaction file is unsorted --- check - FitHiChIP quits !!!! ****"))
	return()
}
	
#======================================================
# form the training data out of the total number of interactions
# if --P2P is enabled, training data consists of only peak to peak interactions
# else training data consists of the complete interaction set
if ((Peak2PeakBackg == 1) & (ncol(interaction.data) >= 15)) {
	TrainingData <- interaction.data[which((interaction.data[,9] == 1) & (interaction.data[,15] == 1)), ]
	# append the interaction distance as the last column of the training data
	TrainingData <- cbind(TrainingData, abs(TrainingData[,2] - TrainingData[,5]))
	# now sort the data with respect to increasing interaction distance
	TrainingData <- TrainingData[order( TrainingData[,ncol(TrainingData)] ), ]
	# remove the last column
	TrainingData <- TrainingData[,1:(ncol(TrainingData)-1)]
} else {
	TrainingData <- interaction.data
}

numTrainingPairs <- nrow(TrainingData)
if (numTrainingPairs == 0) {
	cat(sprintf("\n ********* No training data for FitHiChIP spline fit - error - FitHiChIP quits !!! "))
	return()
}

cat(sprintf("\n ****** Total number of training interactions: %s ********", numTrainingPairs))

# interaction distance vector for the training data only
Trainingdata_GeneDist <- abs(TrainingData[,2] - TrainingData[,5])

# contact count sum for the training data
TotTrainingDataContact <- sum(TrainingData[,opt$cccol])

#=====================================================
# divide the genomic distance into b equal occupancy bins
# (approximately same sum of contact count, but variable number of contacts per bins)

# error condition - sourya
# for very small number of contact count, enforce only a small number of bins
if (TotTrainingDataContact < nbins) {
	nbins <- numTrainingPairs
}
ncontactbin <- floor(TotTrainingDataContact / nbins)
cat(sprintf("\n ****** Number of contacts per bin (allowed for equal occupancy binning): %s ******** \n", ncontactbin))

#=====================================================
# applicable if the spline fitting is performed with all possible locus pairs
# including zero contacts
if (opt$UseNonzeroContacts == 0) {

	time_start <- Sys.time()

	# if interaction type is not all to all
	# then pre-compute the number of locus pairs (including zero contacts)
	# for each of different distance values considered
	if (IntType < 4) {
		# unique genomic distance
		Unique_GeneDist_TrainingData <- unique(Trainingdata_GeneDist)

		if (0) {
			cat(sprintf("\n ******** Unique_GeneDist_TrainingData : %s ", paste(as.vector(Unique_GeneDist_TrainingData))))
			cat(sprintf("\n ********** length Unique_GeneDist_TrainingData: %s ", length(Unique_GeneDist_TrainingData)))		
		}
		# get the vector containing number of locus pairs
		# for each specific distance
		if ((Peak2PeakBackg == 1) & (ncol(interaction.data) >= 15)) {
			LocusPairCountGeneDist <- AllPossibleContacts(Unique_GeneDist_TrainingData, opt$BinSize, ChrSpecBin.df, 1)
		} else {
			LocusPairCountGeneDist <- AllPossibleContacts(Unique_GeneDist_TrainingData, opt$BinSize, ChrSpecBin.df, IntType)
		}
		if (0) {
			cat(sprintf("\n ******** LocusPairCountGeneDist : %s ", paste(as.vector(LocusPairCountGeneDist))))
		}
	}

	time_end <- Sys.time()
	if (0) {
		cat(sprintf("\n\n ==>>> Parallel processing --- time taken for AllPossibleContacts --- %s \n\n", (time_end - time_start)))
	}
}
# #=====================================================

#==============================================
# equal occupancy (contact count) bin
# spline fitting of FitHiC
# interaction distance vs average interaction probability
# Note: the spline fit is performed with respect to the training data
# for peak to peak background, the training data is a subset of the original interactions
#==============================================

# keeps track of total no of contact counts for a particular bin
cumContactCount <- 0
# total number of interacting pairs having nonzero contacts (interactions)
nelem <- 0
# end index of a particular bin
ei <- 0
# index of a particular equal occupancy bin
bin_idx <- 0

if (opt$UseNonzeroContacts == 0) {
	# total number of possible interacting pairs within a bin
	# including pairs of zero contact count
	nelem_all_possble <- 0
}

# keeps track of the no of different genomic distance values encountered within this bin
# and also the no of elements (locus pairs) having that particular genomic distance
ValGeneDist <- c()
nelemSameGeneDist <- c()

while (ei < numTrainingPairs) {
	# start index of this particular bin
	si <- ei + 1
	curr_gene_dist <- abs(TrainingData[si,2] - TrainingData[si,5])
	idx_list <- which(Trainingdata_GeneDist == curr_gene_dist)
	# number of elements (interacting segment pairs) with this particular genomic distance
	# and having nonzero contact count
	nelem <- nelem + length(idx_list)
	ei <- ei + length(idx_list)
	cumContactCount <- cumContactCount + sum(TrainingData[si:ei, opt$cccol])
	
	if (opt$UseNonzeroContacts == 0) {
		# number of possible contacts (including zero)
		# with respect to this particular genomic distance
		if (IntType == 4) {
			# all to all interactions
			all_possible_elem_curr_dist <- AllPossbleContactsGeneDist(curr_gene_dist, opt$BinSize, MappableBinCount.df)	
		} else {
			all_possible_elem_curr_dist <- LocusPairCountGeneDist[Unique_GeneDist_TrainingData == curr_gene_dist]		
		}
		nelem_all_possble <- nelem_all_possble + all_possible_elem_curr_dist
	}
	
	# insert the current gene distance value into a vector
	# which accounts for the different genomic distance values employed for this bin	
	ValGeneDist <- c(ValGeneDist, curr_gene_dist)
	
	# also account for the number of nonzero contacts for this particular genomic distance
	if (opt$UseNonzeroContacts == 0) {
		# considering all locus pairs (including zero contact)
		nelemSameGeneDist <- c(nelemSameGeneDist, all_possible_elem_curr_dist)
	} else {
		# considering only nonzero contact based locus pairs
		nelemSameGeneDist <- c(nelemSameGeneDist, length(idx_list))
	}

	if (cumContactCount >= ncontactbin) {
		# current cumulative contact exceeds the expected average contact per bin
		# and also, all the contacts with this specified genomic distance is covered
		# so, we fix this bin 
		bin_idx <- bin_idx + 1
		if (0) {
			cat(sprintf("\n processing bin_idx : %s ", bin_idx))
		}
		
		# number of interacting pairs having nonzero contact, for this particular bin
		no_distinct_loci[bin_idx] <- nelem
		if (0) {
			cat(sprintf("\n no_distinct_loci[bin_idx] : %s ", no_distinct_loci[bin_idx]))
		}
		
		if (opt$UseNonzeroContacts == 0) {
			# number of all possible interacting pairs (including zero contact) for this particular bin
			no_possible_pair_loci[bin_idx] <- nelem_all_possble
			if (0) {
				cat(sprintf("\n no_possible_pair_loci[bin_idx] : %s ", no_possible_pair_loci[bin_idx]))
			}
		}
		
		# total number of observed contact count for this particular bin
		NumContact[bin_idx] <- cumContactCount
		if (0) {
			cat(sprintf("\n NumContact[bin_idx] : %s ", NumContact[bin_idx]))
		}
		
		# average contact count for this particular bin
		# modification - sourya
		# in the previous implementation of FitHiC, 
		# average contact was computed using only the interacting pairs having nonzero contact
		# in the current implementation, 
		# all possible interacting pairs (including zero contacts) are also accounted
		if (opt$UseNonzeroContacts == 0) {
			avg_contact[bin_idx] <- NumContact[bin_idx] / no_possible_pair_loci[bin_idx]
			if (0) {
				cat(sprintf("\n no_possible_pair_loci[bin_idx] : %s ", no_possible_pair_loci[bin_idx]))
			}
		} else {
			avg_contact[bin_idx] <- NumContact[bin_idx] / no_distinct_loci[bin_idx]
			if (0) {
				cat(sprintf("\n no_distinct_loci[bin_idx] : %s ", no_distinct_loci[bin_idx]))
			}
		}
		if (0) {
			cat(sprintf("\n avg_contact[bin_idx] : %s ", avg_contact[bin_idx]))
		}

		# prior probability with respect to a particular bin
		# is computed by normalizing the average contact count 
		# (contact for individual interacting pair)
		# with the total contact count for all the interactions of all the bins
		prior_contact_prob[bin_idx] <- (avg_contact[bin_idx] / TotTrainingDataContact)
		if (0) {
			cat(sprintf("\n TotTrainingDataContact : %s ", TotTrainingDataContact))
			cat(sprintf("\n prior_contact_prob[bin_idx] : %s ", prior_contact_prob[bin_idx]))
		}

		# computing the average interaction distance for this particular bin		
		temp_sum <- 0
		for (j in (1:length(ValGeneDist))) {
			temp_sum <- temp_sum + (nelemSameGeneDist[j] * ValGeneDist[j])
		}

		if (opt$UseNonzeroContacts == 0) {
			# considering all locus pairs (including zero contact)
			avg_int_dist[bin_idx] <- as.integer(temp_sum / nelem_all_possble)
			if (0) {
				cat(sprintf("\n temp_sum : %s nelem_all_possble : %s avg_int_dist[bin_idx] : %s ", temp_sum, nelem_all_possble, avg_int_dist[bin_idx]))
			}
		} else {
			# considering only nonzero contact based locus pairs
			avg_int_dist[bin_idx] <- as.integer(temp_sum / nelem)
			if (0) {
				cat(sprintf("\n temp_sum : %s nelem : %s avg_int_dist[bin_idx] : %s ", temp_sum, nelem, avg_int_dist[bin_idx]))
			}
		}

		# reset the couners
		cumContactCount <- 0
	 	ValGeneDist <- c()
	 	nelemSameGeneDist <- c()
	 	nelem <- 0
	 	if (opt$UseNonzeroContacts == 0) {
	 		nelem_all_possble <- 0
	 	}
	}
}

# now dump the data in a text file
# bin specific statistics
OutBinfile <- paste0(OutIntDir, '/Bin_Info.log')
if (opt$UseNonzeroContacts == 0) {
	outBinInfoDF <- cbind.data.frame(avg_int_dist, no_distinct_loci, no_possible_pair_loci, NumContact, avg_contact, prior_contact_prob)
	colnames(outBinInfoDF) <- c("AvgIntDist", "NumLociNonZeroContact", "NumLoci_AllPossible", "NumContact","AvgContact","PriorProb")
	write.table(outBinInfoDF, OutBinfile, row.names=F, col.names=T, sep="\t", quote=F, append=F) 
} else {
	outBinInfoDF <- cbind.data.frame(avg_int_dist, no_distinct_loci, NumContact, avg_contact, prior_contact_prob)
	colnames(outBinInfoDF) <- c("AvgIntDist", "NumLociNonZeroContact", "NumContact","AvgContact","PriorProb")
	write.table(outBinInfoDF, OutBinfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
}

cat(sprintf("\n ****** WRITTEN BIN INFORMATION in file : %s ****\n", OutBinfile))

#===============================
# now prepare a spline fit where x axis is the 'avg_int_dist'
# and y axis is the 'prior_contact_prob'
# Note: the spline fit is performed with respect to the training data
# for peak to peak background, the training data is a subset of the original interactions
#===============================

# add - sourya
# check non-NA elements from "avg_int_dist" and "prior_contact_prob"
# and use them for spline fitting
non_NA_idx <- which((!is.na(avg_int_dist)) & (!is.na(prior_contact_prob)) & (is.finite(avg_int_dist)) & (is.finite(prior_contact_prob)))
if (length(non_NA_idx) <= 2) {
	cat(sprintf("\n\n ********* all NA entries in average interaction distance AND / OR prior contact probability -- spline fit is not possible - FitHiChIP quits !!!  **** \n\n "))
	return()
}

# the second spline works according to the cross validation principle
# fit2 <- smooth.spline(avg_int_dist, prior_contact_prob, cv=TRUE)
fit2 <- smooth.spline(avg_int_dist[non_NA_idx], prior_contact_prob[non_NA_idx], cv=TRUE)
cat(sprintf("\n ****** After fit2 smooth.spline ****\n")) 	
if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time for spline (CV): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}

# perform anti-tonic regression on the second spline
fit2.mr <- monoreg(fit2$x, fit2$y, type="antitonic")
cat(sprintf("\n ****** After fit2.mr monoreg ****\n")) 
if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time for antitonic regression on spline (CV): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}

#=========================
# now use this fitted spline
# to get predicted probability values for all distance values of the interactions
# in the unit of bin size employed in the current file)
# here the spline is adapted on the full set of interactions
#=========================
binsize <- abs(interaction.data[1,2] - interaction.data[1,5])
TestDistVal <- seq(min(gene.dist), max(gene.dist), binsize)

# probability obtained from the spline
pp.fit2 <- predict(fit2, TestDistVal)

# use this computed probability and the input distance to fit a spline again
fit2_new <- smooth.spline(pp.fit2$x, pp.fit2$y, cv=TRUE)
cat(sprintf("\n ****** After fit2_new smooth.spline ****\n")) 

# perform antitonic regression on this spline
fit2_new.mr <- monoreg(fit2_new$x, fit2_new$y, type="antitonic")
cat(sprintf("\n ****** After fit2_new.mr monoreg ****\n")) 

# plot the fitted spline and the smoothing regression
# if (opt$Draw) {
if (1) {
	plotfile1 <- paste0(outdir,'/EqOccBin_SplinePass1.png')
	# a <- data.frame(group = paste("Original points"), x = avg_int_dist, y = prior_contact_prob)
	a <- data.frame(group = paste("Original points"), x = avg_int_dist[non_NA_idx], y = prior_contact_prob[non_NA_idx])
	b <- data.frame(group = paste("Fitted spline"), x = fit2_new.mr$x, y = fit2_new.mr$yf)
	curr_plotA <- ggplot(rbind(a, b), aes(x=x, y=y, color=group, fill=group)) + geom_point(size=0.1) + geom_line(size=0.3) + xlab('Average interaction distance') + ylab('Prior contact probability') + scale_colour_manual(values = c("blue", "red"))
	curr_plotA + ggtitle("Spline fit")
	ggsave(plotfile1, plot = curr_plotA, width=5, height=5)
	cat(sprintf("\n ****** After plotting file : %s ****\n", plotfile1)) 
}

# plot the fitted spline and the smoothing regression (in log scale)
if (1) {	
	plotfile1 <- paste0(outdir,'/EqOccBin_SplinePass1_LOGScale.pdf')	
	# a <- data.frame(group = paste("Original"), x = log10(avg_int_dist), y = log10(prior_contact_prob))
	a <- data.frame(group = paste("Original"), x = log10(avg_int_dist[non_NA_idx]), y = log10(prior_contact_prob[non_NA_idx]))
	b <- data.frame(group = paste("spline (fit2) - ", as.integer(fit2$df), "df (cv)"), x = log10(fit2$x), y = log10(fit2$y))
	c <- data.frame(group = paste("MR on spline (fit2.mr)"), x = log10(fit2.mr$x), y = log10(fit2.mr$yf))
	d <- data.frame(group = paste("5 Kb uniform spline (fit2_new) - ", as.integer(fit2_new$df), "df (cv)"), x = log10(fit2_new$x), y = log10(fit2_new$y))
	e <- data.frame(group = paste("MR on 5 Kb uniform spline (fit2_new.mr)"), x = log10(fit2_new.mr$x), y = log10(fit2_new.mr$yf))
	curr_plotA <- ggplot(rbind(a, b, c, d, e), aes(x=x, y=y, color=group, fill=group)) + geom_point(size=0.1) + geom_line(size=0.3) + xlab('Average interaction distance (log)') + ylab('Prior contact probability (log)') + scale_colour_manual(values = c("black", "yellow", "blue", "green", "red"))
	curr_plotA + ggtitle("Smooth spline - antitonic regression - pass 1 - log scale")
	ggsave(plotfile1, plot = curr_plotA, width=8, height=6)
	cat(sprintf("\n ****** After plotting file : %s ****\n", plotfile1)) 
}

if (0) {
	# create a data frame to write the old spline and old MR regression results
	old_sample.df <- cbind.data.frame(fit2$x, fit2$y, fit2.mr$x, fit2.mr$y, fit2.mr$yf)
	colnames(old_sample.df) <- c('fit2$x', 'fit2$y', 'fit2.mr$x', 'fit2.mr$y', 'fit2.mr$yf')
	write.table(old_sample.df, paste0(outdir,'/EqOccBin_OldSpline_Regression.txt'), row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

	# create a data frame to write the old spline and old MR regression results
	new_sample.df <- cbind.data.frame(fit2_new$x, fit2_new$y, fit2_new.mr$x, fit2_new.mr$y, fit2_new.mr$yf)
	colnames(new_sample.df) <- c('fit2_new$x', 'fit2_new$y', 'fit2_new.mr$x', 'fit2_new.mr$y', 'fit2_new.mr$yf')
	write.table(new_sample.df, paste0(outdir,'/EqOccBin_NewSpline_Regression.txt'), row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)
}

#===========================================
# now predict the significance, P value etc. 
# (binomial distribution) of individual interactions
# check whether bias correction is enabled or not
#===========================================

# from the sorted gene distance values (of the full set of input interactions)
# find the indices starting a new distance value
# helpful to apply spline based predict function on a limited nuber of times
uniq.idx <- order(gene.dist)[!duplicated(gene.dist)]

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time to find unique indices in genomic distance: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}

#======================
# check for bias correction
#======================
if (opt$BiasCorr == 0) {
	#======================
	# bias correction is not enabled
	# same as FitHiC output
	#======================
	for (k in (1:length(uniq.idx))) {
		# all the locations from si to ei will have the same probability 
		# since the gene distance value is the same
		si <- uniq.idx[k]
		if (k < length(uniq.idx)) {
			ei <- uniq.idx[k+1] - 1
		} else {
			ei <- numPairs	# last read	
		}

		# compute the probability according to the genomic distance (of the start location - single element)
		# this probability value will be used for all the values within the interval (si to ei)
		# predict from the new spline fit - sourya
		pp <- predict(fit2_new, gene.dist[si])
		p <- pp$y

		# estimate the probability and the P value by applying binomial distribution
		curr_dbinom <- dbinom(interaction.data[si:ei, opt$cccol], size=TotContact, prob=p)
		curr_pbinom <- pbinom(interaction.data[si:ei, opt$cccol], size=TotContact, prob=p, lower.tail=FALSE)

		# append the probability and P values to the final vector
	 	Spline_Binom_Prob_CC <- c(Spline_Binom_Prob_CC, curr_dbinom)
	 	Spline_Binom_P_Val_CC <- c(Spline_Binom_P_Val_CC, (curr_dbinom + curr_pbinom))

	 	# store the prior probability as well
	 	Prob_Val[si:ei] <- p

	}	# end distance specific loop for significance computation

} else {
	#======================
	# bias correction is enabled
	#======================
	
	if (opt$MultBias == 1) {

		#===========================
		# case 1 - multiply the probabilities with the bias values, to assess the statistical significance
		# Note: we analyze the full set of input data even if peak to peak background is used
		# we start with spline based probability for individual distance values 
		#===========================
		for (k in (1:length(uniq.idx))) {
			# all the locations from si to ei will have the same probability 
			# since the gene distance value is the same
			si <- uniq.idx[k]
			if (k < length(uniq.idx)) {
				ei <- uniq.idx[k+1] - 1
			} else {
				ei <- numPairs	# last read	
			}

			# compute the probability according to the genomic distance 
			# this probability value will be used for all the values within the interval (si to ei)
			# predict from the new spline fit - sourya
			pp <- predict(fit2_new, gene.dist[si])
			p <- pp$y		

			result_binomdistr <- as.data.frame(parallel:::mclapply( si:ei , mc.cores = ncore , function(idx){
				cc <- interaction.data[idx, opt$cccol]
				# probability by multiplying the bias values
				pr <- p * interaction.data[idx, bias1.col] * interaction.data[idx, bias2.col]
				db <- dbinom(cc, size=TotContact, prob=pr)
				pb <- pbinom(cc, size=TotContact, prob=pr, lower.tail=FALSE)
				return(c(idx,pr,db,(db+pb)))
			} ))

			Prob_Val[si:ei] <- as.double(result_binomdistr[2,1:ncol(result_binomdistr)])
			Spline_Binom_Prob_CC[si:ei] <- as.double(result_binomdistr[3,1:ncol(result_binomdistr)])
			Spline_Binom_P_Val_CC[si:ei] <- as.double(result_binomdistr[4,1:ncol(result_binomdistr)])

		}	# end distance loop 

	} else {

		#===========================
		# Case 2: linear regression between the bias values and the observed contact count
		# Note: here we use the training data 
		# since the spline fit and the bin log was performed on the training data
		# for peak to peak background, the training data is a subset of the total input set of interactions
		#===========================

		# ******* condition for using either equal occupancy binning (recommended) 
		# or equal distance binning **********

		if (opt$EqOcc == 1) {
			cat(sprintf("\n ****** Start regression - Equal Occupancy binning ****\n")) 
			#====================
			# case 2.1: equal occupancy binning
			#====================

			time_start <- Sys.time()

			res_Row_Training <- as.data.frame(parallel:::mclapply( 1:length(avg_int_dist) , mc.cores = ncore , function(rownum){

				if (rownum == 1) {
					training_si <- 1
					training_ei <- no_distinct_loci[rownum]
				} else {
					training_si <- sum(no_distinct_loci[1:(rownum-1)]) + 1
					training_ei <- training_si + no_distinct_loci[rownum] - 1
				}
				training_prior_prob <- prior_contact_prob[rownum]
				training_avg_int_dist <- avg_int_dist[rownum]

				if (0) {
					cat(sprintf("\n\n === bias correction - equal occupancy binning -- analyzing the bins --- rownum: %s   training_si: %s  training_ei: %s  training_prior_prob: %s  ", rownum, training_si, training_ei, training_prior_prob))
				}

				# form the probability vector (having the length training_ei - training_si + 1)
				# obtained from the spline on training data
				# it contains gene distance specific probabilities 
				Prob_Vec_Curr_Row <- c()
				curr_training_data_gene_dist <- Trainingdata_GeneDist[training_si:training_ei]
				curr_si <- training_si
				while (curr_si <= training_ei) {
					curr_dist <- curr_training_data_gene_dist[curr_si - training_si + 1]
					idx_list <- which(curr_training_data_gene_dist == curr_dist)					
					if (0) {
						cat(sprintf("\n curr_si : %s curr_dist: %s  length idx_list : %s ", curr_si, curr_dist, length(idx_list)))
					}
					# estimate spline based probability for the current distance
					pp <- predict(fit2_new, curr_dist)
					prob_curr_dist <- pp$y
					if (0) {
						cat(sprintf("   prob_curr_dist: %s ", prob_curr_dist))
					}
					Prob_Vec_Curr_Row <- c(Prob_Vec_Curr_Row, rep(prob_curr_dist, length(idx_list)))
					curr_si <- curr_si + length(idx_list)
				}

				# filter some interactions (training data) according to the bias thresholds
				Curr_Int_Data <- TrainingData[training_si:training_ei, ]
				bias1_val_vec <- as.numeric(Curr_Int_Data[,bias1.col])
				bias2_val_vec <- as.numeric(Curr_Int_Data[,bias2.col])
				BiasThrIdx <- which((bias1_val_vec >= opt$BiasLowThr) & (bias1_val_vec <= opt$BiasHighThr) & (bias2_val_vec >= opt$BiasLowThr) & (bias2_val_vec <= opt$BiasHighThr))
				Curr_Int_Set_BiasThr <- Curr_Int_Data[BiasThrIdx, ]

				if (0) {
					cat(sprintf("\n === Regression model (training data) - rownum: %s   training_si: %s  training_ei: %s  training_prior_prob: %s  average distance: %s number of interactions: %s  number of bias threshold satisfying interactions: %s \n", rownum, training_si, training_ei, training_prior_prob, training_avg_int_dist, nrow(Curr_Int_Data), nrow(Curr_Int_Set_BiasThr)))
				}

				# subset of the probability vector, which pass the bias threshold 
				Prob_Vec_Curr_Row_BiasFilt <- Prob_Vec_Curr_Row[BiasThrIdx]

				# variance within the contact count
				TSS_CC <- sum((Curr_Int_Set_BiasThr[, opt$cccol] - mean(Curr_Int_Set_BiasThr[, opt$cccol]))^2)

				# expected contact count according to the spline fit probability
				# computed with respect to the training contacts satisfying bias thresholds
				distCCVec <- TotTrainingDataContact * Prob_Vec_Curr_Row_BiasFilt
				logdistCCVec <- log10(distCCVec)

				# debug - sourya
				if (0) {
					tempoutfile <- paste0(outdir, '/', rownum, '_', training_si, '_', training_ei, '_', training_avg_int_dist, '_', nrow(Curr_Int_Data), '_', nrow(Curr_Int_Set_BiasThr), '.txt')
					cat(sprintf("\n before writing regression data in file : %s ", tempoutfile))			
					bias_regression_data_temp <- data.frame(bias1=Curr_Int_Set_BiasThr[, bias1.col], bias2=Curr_Int_Set_BiasThr[, bias2.col], CC=Curr_Int_Set_BiasThr[, opt$cccol])
					write.table(bias_regression_data_temp, tempoutfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
					cat(sprintf("\n wrote regression data in file : %s  - waiting for linear regression ", tempoutfile))
				}
				# end debug - sourya

				if (opt$Resid == 0) {
					#=====================
					# case 2.1.1
					# model the regression between the observed contact count and the bias values
					#=====================
					bias_regression_data <- data.frame(logbias1=log10(as.numeric(Curr_Int_Set_BiasThr[, bias1.col])), logbias2=log10(as.numeric(Curr_Int_Set_BiasThr[, bias2.col])), logCC=log10(as.numeric(Curr_Int_Set_BiasThr[, opt$cccol])))
					# model the regression
					Linear_logbias_logCC <- lm(logCC ~ logbias1 + logbias2, data=bias_regression_data)
					if (0) {
						cat(sprintf("\n ~~~ modeled linear regression of bias"))
					}				
					# the model has three coefficients: 1) intercept, 2) for logbias1, 3) for logbias2
					coeff_intercept <- Linear_logbias_logCC$coefficients[1]
					coeff_logbias1 <- Linear_logbias_logCC$coefficients[2]
					coeff_logbias2 <- Linear_logbias_logCC$coefficients[3]

				} else {
					#=====================
					# case 2.1.2
					# model the regression with respect to the (observed contact count - spline fitted contact count) and the bias values
					#=====================
					bias_regression_data <- data.frame(logbias1=log10(as.numeric(Curr_Int_Set_BiasThr[, bias1.col])), logbias2=log10(as.numeric(Curr_Int_Set_BiasThr[, bias2.col])), logCC=(log10(as.numeric(Curr_Int_Set_BiasThr[, opt$cccol]) - logdistCCVec)))
					# model the regression
					Linear_logbias_logCC <- lm(logResidCC ~ logbias1 + logbias2, data=bias_regression_data)
					# the model has three coefficients: 1) intercept, 2) for logbias1, 3) for logbias2
					coeff_intercept <- Linear_logbias_logCC$coefficients[1]
					coeff_logbias1 <- Linear_logbias_logCC$coefficients[2]
					coeff_logbias2 <- Linear_logbias_logCC$coefficients[3]
				}
				
				# get the performance statistics of the regression model
				RSS_LM <- sum((10^residuals(Linear_logbias_logCC))^2)
				R.square_LM <- 1 - ((RSS_LM * 1.0) / TSS_CC)
				AIC_LM <- AIC(Linear_logbias_logCC)	

				# debug - sourya
				if (0) {
					cat(sprintf("\n ===>> rownum : %s coeff_intercept : %s coeff_logbias1 : %s coeff_logbias2 : %s RSS_LM : %s R.square_LM : %s AIC_LM : %s ", rownum, coeff_intercept, coeff_logbias1, coeff_logbias2, RSS_LM, R.square_LM, AIC_LM))
				}

				# form a vector of the current distance and regression model cofficients
				# currently the vector content is same for both residual and non-residual model
				currvec <- c(training_avg_int_dist, coeff_intercept, coeff_logbias1, coeff_logbias2, RSS_LM, R.square_LM, AIC_LM)

				# return the row number and the regression model vectors
				return(c(rownum, currvec))
			
			} ))

			# form a data frame with transpose of the contents of res_Row_Training
			# from row 2 to last row, and for all columns (transpose)
			# since the row 1 contains the row number values
			Regression_model_coeff.df <- t(res_Row_Training[2:nrow(res_Row_Training), 1:ncol(res_Row_Training)])
			colnames(Regression_model_coeff.df) <- c("Dist", "Intercept", "LogBias1", "LogBias2", "RSS_LM", "R.square_LM", "AIC_LM")

			# write the training data specific regression model 
			Estimated_Model_Coeff_File <- paste0(outdir, '/Regression_Model_Coeff.log')	
			write.table(Regression_model_coeff.df, Estimated_Model_Coeff_File, row.names=F, col.names=T, sep="\t", quote=F, append=F)

			time_end <- Sys.time()
			if (0) {
				cat(sprintf("\n\n\n ==>>> Processing bin rows - res_Row_Training - parallel optimization - new time: %s \n\n", (time_end - time_start)))
			}
			cat(sprintf("\n ****** Written regression model coefficients in the file : %s ****\n", Estimated_Model_Coeff_File)) 

			#============================
			# fit spline for individual distance values
			# vs three parameters:
			# 1) coefficient of intercept
			# 2) coefficient of logbias1
			# 3) coefficient of logbias2
			# first model the spline with respect to the available binned data
			# next extend the spline with all possible distance values (in the vector TestDistVal)
			#============================

			# add - sourya
			# check non-NA elements from average interaction distance and bias regression model coefficient intercept, and use them for spline fitting
			non_NA_idx <- which((!is.na(Regression_model_coeff.df[,1])) & (!is.na(Regression_model_coeff.df[,2])) & (is.finite(Regression_model_coeff.df[,1])) & (is.finite(Regression_model_coeff.df[,2])))
			if (length(non_NA_idx) <= 2) {
				cat(sprintf("\n\n ********* all NA entries in average interaction distance AND / OR bias regression model coefficient intercept -- spline fit is not possible - FitHiChIP quits !!!  **** \n\n "))
				return()
			}

			# fit_spline_coeff_Intercept <- smooth.spline(Regression_model_coeff.df[,1], Regression_model_coeff.df[,2], cv=TRUE)
			fit_spline_coeff_Intercept <- smooth.spline(Regression_model_coeff.df[non_NA_idx, 1], Regression_model_coeff.df[non_NA_idx, 2], cv=TRUE)
			pp.fit_spline_coeff_Intercept <- predict(fit_spline_coeff_Intercept, TestDistVal)
			fit_spline_coeff_Intercept_new <- smooth.spline(pp.fit_spline_coeff_Intercept$x, pp.fit_spline_coeff_Intercept$y, cv=TRUE)

			cat(sprintf("\n ===>>> using bias regression --- modeled fit_spline_coeff_Intercept"))

			if (1) {
				plotfile1 <- paste0(outdir, '/fit_spline_coeff_Intercept_new.pdf')	
				# a <- data.frame(group = paste("Original"), x = Regression_model_coeff.df[,1], y = Regression_model_coeff.df[,2])
				a <- data.frame(group = paste("Original"), x = Regression_model_coeff.df[non_NA_idx, 1], y = Regression_model_coeff.df[non_NA_idx, 2])
				b <- data.frame(group = paste("spline fit"), x = fit_spline_coeff_Intercept_new$x, y = fit_spline_coeff_Intercept_new$y)
				curr_plotA <- ggplot(rbind(a, b), aes(x=x, y=y, color=group, fill=group)) + geom_point(size=0.1) + geom_line(size=0.3) + xlab("Average interaction distance") + ylab("Regression coefficient - intercept") + scale_colour_manual(values = c("black", "red"))
				curr_plotA + ggtitle("fit_spline_coeff_Intercept")
				ggsave(plotfile1, plot = curr_plotA, width=5, height=5)				
			}

			# plot the fitted spline "fit_spline_coeff_Intercept_new" in log scale
			if (0) { #(opt$Draw) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Intercept_new_LOGScale.pdf')	
				# a <- data.frame(group = paste("Original"), x = log10(Regression_model_coeff.df[,1]), y = log10(Regression_model_coeff.df[,2]))
				a <- data.frame(group = paste("Original"), x = log10(Regression_model_coeff.df[non_NA_idx, 1]), y = log10(Regression_model_coeff.df[non_NA_idx, 2]))
				b <- data.frame(group = paste("spline fit"), x = log10(fit_spline_coeff_Intercept_new$x), y = log10(fit_spline_coeff_Intercept_new$y))
				curr_plotA <- ggplot(rbind(a, b), aes(x=x, y=y, color=group, fill=group)) + geom_point(size=0.1) + geom_line(size=0.3) + xlab("Average interaction distance (log)") + ylab("Regression coefficient - intercept (log)") + scale_colour_manual(values = c("black", "red"))
				curr_plotA + ggtitle("fit_spline_coeff_Intercept - log scale")
				ggsave(plotfile1, plot = curr_plotA, width=5, height=5)				
			}

			# add - sourya
			# check non-NA elements from average interaction distance and bias regression model logbias1, and use them for spline fitting
			non_NA_idx <- which((!is.na(Regression_model_coeff.df[,1])) & (!is.na(Regression_model_coeff.df[,3])) & (is.finite(Regression_model_coeff.df[,1])) & (is.finite(Regression_model_coeff.df[,3])))
			if (length(non_NA_idx) <= 2) {
				cat(sprintf("\n\n ********* all NA entries in average interaction distance AND / OR bias regression model logbias1 -- spline fit is not possible - FitHiChIP quits !!!  **** \n\n "))
				return()
			}

			# fit_spline_coeff_Logbias1 <- smooth.spline(Regression_model_coeff.df[,1], Regression_model_coeff.df[,3], cv=TRUE)
			fit_spline_coeff_Logbias1 <- smooth.spline(Regression_model_coeff.df[non_NA_idx, 1], Regression_model_coeff.df[non_NA_idx, 3], cv=TRUE)
			pp.fit_spline_coeff_Logbias1 <- predict(fit_spline_coeff_Logbias1, TestDistVal)
			fit_spline_coeff_Logbias1_new <- smooth.spline(pp.fit_spline_coeff_Logbias1$x, pp.fit_spline_coeff_Logbias1$y, cv=TRUE)

			cat(sprintf("\n ===>>> using bias regression --- modeled fit_spline_coeff_Logbias1"))

			# plot the fitted spline "fit_spline_coeff_Logbias1_new"
			if (1) { #(opt$Draw) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Logbias1_new.pdf')	
				# a <- data.frame(group = paste("Original"), x = Regression_model_coeff.df[,1], y = Regression_model_coeff.df[,3])
				a <- data.frame(group = paste("Original"), x = Regression_model_coeff.df[non_NA_idx, 1], y = Regression_model_coeff.df[non_NA_idx, 3])
				b <- data.frame(group = paste("spline fit"), x = fit_spline_coeff_Logbias1_new$x, y = fit_spline_coeff_Logbias1_new$y)
				curr_plotA <- ggplot(rbind(a, b), aes(x=x, y=y, color=group, fill=group)) + geom_point(size=0.1) + geom_line(size=0.3) + xlab("Average interaction distance") + ylab("Regression coefficient - logbias1") + scale_colour_manual(values = c("black", "red"))
				curr_plotA + ggtitle("fit_spline_coeff_logbias1")
				ggsave(plotfile1, plot = curr_plotA, width=5, height=5)				
			}

			# plot the fitted spline "fit_spline_coeff_Logbias1_new" in log scale
			if (0) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Logbias1_new_LOGScale.pdf')	
				# a <- data.frame(group = paste("Original"), x = log10(Regression_model_coeff.df[,1]), y = log10(Regression_model_coeff.df[,3]))
				a <- data.frame(group = paste("Original"), x = log10(Regression_model_coeff.df[non_NA_idx, 1]), y = log10(Regression_model_coeff.df[non_NA_idx, 3]))
				b <- data.frame(group = paste("spline fit"), x = log10(fit_spline_coeff_Logbias1_new$x), y = log10(fit_spline_coeff_Logbias1_new$y))
				curr_plotA <- ggplot(rbind(a, b), aes(x=x, y=y, color=group, fill=group)) + geom_point(size=0.1) + geom_line(size=0.3) + xlab("Average interaction distance (log)") + ylab("Regression coefficient - logbias1 (log)") + scale_colour_manual(values = c("black", "red"))
				curr_plotA + ggtitle("fit_spline_coeff_logbias1 - log scale")
				ggsave(plotfile1, plot = curr_plotA, width=5, height=5)
			}

			# add - sourya
			# check non-NA elements from average interaction distance and bias regression model logbias2, and use them for spline fitting
			non_NA_idx <- which((!is.na(Regression_model_coeff.df[,1])) & (!is.na(Regression_model_coeff.df[,4])) & (is.finite(Regression_model_coeff.df[,1])) & (is.finite(Regression_model_coeff.df[,4])))
			if (length(non_NA_idx) <= 2) {
				cat(sprintf("\n\n ********* all NA entries in average interaction distance AND / OR bias regression model logbias2 -- spline fit is not possible - FitHiChIP quits !!!  **** \n\n "))
				return()
			}

			# fit_spline_coeff_Logbias2 <- smooth.spline(Regression_model_coeff.df[, 1], Regression_model_coeff.df[, 4], cv=TRUE)
			fit_spline_coeff_Logbias2 <- smooth.spline(Regression_model_coeff.df[non_NA_idx, 1], Regression_model_coeff.df[non_NA_idx, 4], cv=TRUE)
			pp.fit_spline_coeff_Logbias2 <- predict(fit_spline_coeff_Logbias2, TestDistVal)
			fit_spline_coeff_Logbias2_new <- smooth.spline(pp.fit_spline_coeff_Logbias2$x, pp.fit_spline_coeff_Logbias2$y, cv=TRUE)

			cat(sprintf("\n ===>>> using bias regression --- modeled fit_spline_coeff_Logbias2"))

			# plot the fitted spline "fit_spline_coeff_Logbias2_new"
			if (1) { 
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Logbias2_new.pdf')	
				# a <- data.frame(group = paste("Original"), x = Regression_model_coeff.df[,1], y = Regression_model_coeff.df[,4])
				a <- data.frame(group = paste("Original"), x = Regression_model_coeff.df[non_NA_idx, 1], y = Regression_model_coeff.df[non_NA_idx, 4])
				b <- data.frame(group = paste("spline fit"), x = fit_spline_coeff_Logbias2_new$x, y = fit_spline_coeff_Logbias2_new$y)
				curr_plotA <- ggplot(rbind(a, b), aes(x=x, y=y, color=group, fill=group)) + geom_point(size=0.1) + geom_line(size=0.3) + xlab("Average interaction distance") + ylab("Regression coefficient - logbias2") + scale_colour_manual(values = c("black", "red"))
				curr_plotA + ggtitle("fit_spline_coeff_logbias2") 
				ggsave(plotfile1, plot = curr_plotA, width=5, height=5)				
			}

			# plot the fitted spline "fit_spline_coeff_Logbias2_new" in log scale
			if (0) { 
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Logbias2_new_LOGScale.pdf')	
				# a <- data.frame(group = paste("Original"), x = log10(Regression_model_coeff.df[,1]), y = log10(Regression_model_coeff.df[,4]))
				a <- data.frame(group = paste("Original"), x = log10(Regression_model_coeff.df[non_NA_idx, 1]), y = log10(Regression_model_coeff.df[non_NA_idx, 4]))
				b <- data.frame(group = paste("spline fit"), x = log10(fit_spline_coeff_Logbias2_new$x), y = log10(fit_spline_coeff_Logbias2_new$y))
				curr_plotA <- ggplot(rbind(a, b), aes(x=x, y=y, color=group, fill=group)) + geom_point(size=0.1) + geom_line(size=0.3) + xlab("Average interaction distance (log)") + ylab("Regression coefficient - logbias2 (log)") + scale_colour_manual(values = c("black", "red"))
				curr_plotA + ggtitle("fit_spline_coeff_logbias2 - log scale")
				ggsave(plotfile1, plot = curr_plotA, width=5, height=5)				
			}
			#============================

			# debug - sourya
			time_start <- Sys.time()

			# using the fitted splines for the distance / probability and the bias correction coefficients
			# we now use the full set of input interactions
			# to model the expected contact count and the corresponding bias corrected probability 
			for (k in (1:length(uniq.idx))) {
				# all the locations from si to ei will have the same probability 
				# with respect to fitted spline (FitHiC distance)
				si <- uniq.idx[k]
				if (k < length(uniq.idx)) {
					ei <- uniq.idx[k+1] - 1
				} else {
					ei <- numPairs	# last read	
				}

				# probability of spline fit for average distance vs avg prior probability
				pp <- predict(fit2_new, gene.dist[si])
				p <- pp$y

				# also predict the coefficients of intercept, logbias1, and logbias2
				# from different splines with respect to this distance
				pp_I <- predict(fit_spline_coeff_Intercept_new, gene.dist[si])
				coeff_Intcpt <- pp_I$y
				pp_L1 <- predict(fit_spline_coeff_Logbias1_new, gene.dist[si])
				coeff_LogBias1 <- pp_L1$y
				pp_L2 <- predict(fit_spline_coeff_Logbias2_new, gene.dist[si])
				coeff_LogBias2 <- pp_L2$y

				if (0) {
					cat(sprintf("\n *** Processing uniq distance idx: %s  si: %s  ei: %s num elem : %s  gene dist : %s  Spline prob (p): %s  coeff_Intcpt : %s coeff_LogBias1 : %s coeff_LogBias2 : %s ", k, si, ei, (ei-si+1), gene.dist[si], p, coeff_Intcpt, coeff_LogBias1, coeff_LogBias2))
				}

				# copy the spline based probability 
				Prob_Val[si:ei] <- p

				# #================
				# # add - sourya

				# # first compute the expected contact count
				# # using the regression coefficients
				# # and then filter out the NA entries
				# if (opt$Resid == 0) {
				# 	expCCVec <- 10^(coeff_Intcpt + coeff_LogBias1 * log10(interaction.data[si:ei, bias1.col]) + coeff_LogBias2 * log10(interaction.data[si:ei, bias2.col]))
				# } else {
				# 	expCCVec <- 10^(coeff_Intcpt + coeff_LogBias1 * log10(interaction.data[si:ei, bias1.col]) + coeff_LogBias2 * log10(interaction.data[si:ei, bias2.col]) + log10(TotContact * p))
				# }

				# # get indices having non-NA expected contact count
				# # Note: here we add the quantity (si-1) since the which operator returns indices from 1
				# # to get the indices mapped back on the range si:ei, we add this offset
				# non_NA_idx <- which(!is.na(expCCVec))
				# nonzero_biasidx_set <- non_NA_idx + (si-1)
				# # the zero bias values are those not in nonzero_biasidx_set
				# zero_biasidx_set <- setdiff(seq(si, ei), nonzero_biasidx_set)

				# # end add - sourya
				# #================

				#================
				# get indices having zero and non zero bias 
				# Note: here we add the quantity (si-1) since the which operator returns indices from 1
				# to get the indices mapped back on the range si:ei, we add this offset

				bias1_val_vec <- as.numeric(interaction.data[si:ei, bias1.col])
				bias2_val_vec <- as.numeric(interaction.data[si:ei, bias2.col])

				if (opt$BiasType == 1) {
					# coverage bias
					nonzero_biasidx_set <- which((bias1_val_vec > 0) & (bias2_val_vec > 0)) + (si-1)
				} else {
					# ICE bias
					# use bias values < 1000
					nonzero_biasidx_set <- which((bias1_val_vec > 0) & (bias1_val_vec < 1000) & (bias2_val_vec > 0) & (bias2_val_vec < 1000)) + (si-1)
				}

				# the zero bias values are those not in nonzero_biasidx_set
				zero_biasidx_set <- setdiff(seq(si, ei), nonzero_biasidx_set)

				if (0) {
					if (opt$BiasType == 1) {
						cat(sprintf("\n\n ****** Coverage bias regression - analyzing locus pairs between indices %s and %s - total elements : %s genomic distance : %s Spline fitted probability (p): %s number of locus pairs with non-zero bias values : %s number of locus pairs with zero bias in at least one end : %s ", si, ei, (ei-si+1), gene.dist[si], p, length(nonzero_biasidx_set), length(zero_biasidx_set)))
					} else {
						cat(sprintf("\n\n ****** ICE bias regression - analyzing locus pairs between indices %s and %s - total elements : %s genomic distance : %s Spline fitted probability (p): %s number of locus pairs with non-zero bias values (and within allowed thresholds) : %s number of locus pairs with zero bias in at least one end, or outside the allowed thresholds : %s ", si, ei, (ei-si+1), gene.dist[si], p, length(nonzero_biasidx_set), length(zero_biasidx_set)))
					}
				}

				# if one of the bias values are zero
				# then just use the spline predicted probability as the estimated distance
				if (length(zero_biasidx_set) > 0) {
					Exp_CC_BiasRegr[zero_biasidx_set] <- (TotContact * p)	
					if (0) {
						cat(sprintf("\n *** for locus pairs with zero bias in at least one end (or bias values outside very high thresholds) - assigning expected contact count value : %s ", (TotContact * p)))
					}
				}
				
				# for non-zero bias values, use the bias regression predicted contact count
				if (length(nonzero_biasidx_set) > 0) {
					if (opt$Resid == 0) {
						Exp_CC_BiasRegr[nonzero_biasidx_set] <- 10^(coeff_Intcpt + coeff_LogBias1 * log10(as.numeric(interaction.data[nonzero_biasidx_set, bias1.col])) + coeff_LogBias2 * log10(as.numeric(interaction.data[nonzero_biasidx_set, bias2.col])))
					} else {
						# here residual contact count is predicted
						# add the spline fit contact count to the predicted value
						Exp_CC_BiasRegr[nonzero_biasidx_set] <- 10^(coeff_Intcpt + coeff_LogBias1 * log10(as.numeric(interaction.data[nonzero_biasidx_set, bias1.col])) + coeff_LogBias2 * log10(as.numeric(interaction.data[nonzero_biasidx_set, bias2.col])) + log10(TotContact * p))
					}

					# condition add - sourya
					# for ICE bias, check if the expected contact count is going beyond a certain limit
					if (opt$BiasType == 2) {						
						limit_exceed_idx <- which(Exp_CC_BiasRegr > (2^15))
						if (length(limit_exceed_idx) > 0) {
							Exp_CC_BiasRegr[limit_exceed_idx] <- (TotContact * p)
							if (0) {
								cat(sprintf("\n *** ICE bias regression -- number of elements with very high expected CC : %s trimmed their expected contact count to the value : %s ", length(limit_exceed_idx), (TotContact * p)))
							}
						}
					}
					# end condition - sourya
					if (0) {
						cat(sprintf("\n *** processing nonzero_biasidx_set - maximum Exp_CC_BiasRegr in these elements : %s ", max(Exp_CC_BiasRegr[nonzero_biasidx_set])))
					}
				}

			}	# end distance loop 

			if (1) {
				cat(sprintf("\n\n *** modeled the bias regression based probability and expected contact count **** \n\n"))
			}

			time_end <- Sys.time()
			if (0) {
				cat(sprintf("\n\n\n ==>>> modeled the bias regression -- time: %s \n\n", (time_end - time_start)))
			}

			# now use the expected contact count from the bias correction 
			# to model the binomial distribution
			NumElem <- length(Exp_CC_BiasRegr)
			non_NA_Exp_CC_BiasRegr <- which(!is.na(Exp_CC_BiasRegr))
			ExpTotContact_1 <- sum(Exp_CC_BiasRegr)
			# check the maximum integer range and trim the sum of expected contact count
			ExpTotContact <- as.integer(sum(Exp_CC_BiasRegr))
			# ExpTotContact <- as.integer(min(sum(Exp_CC_BiasRegr), as.integer(.Machine$integer.max)))

			if (0) {
				cat(sprintf("\n *** NumElem: %s  length non_NA_Exp_CC_BiasRegr : %s  ExpTotContact_1 : %s ExpTotContact : %s ", NumElem, length(non_NA_Exp_CC_BiasRegr), ExpTotContact_1, ExpTotContact))
			}

			# debug - sourya
			time_start <- Sys.time()

			# parallel execution to model the binomial distribution
			# with respect to the expected contact count
			# and the expected sum of contact counts
			# such separate execution is required
			# to model the pdf integral sum as 1

			# process in subsets 
			for (k in (1:length(uniq.idx))) {
				# all the locations from si to ei will have the same probability 
				# with respect to fitted spline (FitHiC distance)
				si <- uniq.idx[k]
				if (k < length(uniq.idx)) {
					ei <- uniq.idx[k+1] - 1
				} else {
					ei <- numPairs	# last read	
				}

				if (0) {
					cat(sprintf("\n *** Modeling regression bias probability based P value --- uniq distance idx: %s  si: %s  ei: %s num elem : %s ", k, si, ei, (ei-si+1)))
				}

				# parallel processing
				resbias <- as.data.frame(parallel:::mclapply( si:ei , mc.cores = ncore , function(idx){
					curr_cc <- interaction.data[idx, opt$cccol]
					pr_bias <- (Exp_CC_BiasRegr[idx] * 1.0) / ExpTotContact
					db2 <- dbinom(curr_cc, size=ExpTotContact, prob=pr_bias)
					pb2 <- pbinom(curr_cc, size=ExpTotContact, prob=pr_bias, lower.tail=FALSE)
					if (0) {
						cat(sprintf("\n idx : %s curr_cc: %s  pr_bias : %s  db2 : %s  pb2: %s ", idx, curr_cc, pr_bias, db2, pb2))
					}
					return(c(idx,pr_bias,db2,(db2+pb2)))
				} ))

				# merge all the probability values
				# and the binomial distribution
				# in the final array
				Prob_BiasRegr[si:ei] <- as.double(resbias[2,1:ncol(resbias)])
				Spline_Binom_Prob_CC[si:ei] <- as.double(resbias[3,1:ncol(resbias)])
				Spline_Binom_P_Val_CC[si:ei] <- as.double(resbias[4,1:ncol(resbias)])
			}	# end distance loop

			if (1) {
				cat(sprintf("\n\n *** p-value estimation is complete for the bias regression **** \n\n"))
			}

			time_end <- Sys.time()
			if (0) {
				cat(sprintf("\n\n\n ==>>> modeled the binomial distribution with respect to bias regression -- time: %s \n\n", (time_end - time_start)))
			}
			
			# temporary plotting the probability and contact count plots
			if (0) {
				# drawing the comparison between original bias correction 
				# and the new bias correction (regression) method
				ProbComparePlotDir <- paste0(outdir, '/Prob_Compare_Plots')
				system(paste("mkdir -p", ProbComparePlotDir))
				CCComparePlotDir <- paste0(outdir, '/CC_Compare_Plots')
				system(paste("mkdir -p", CCComparePlotDir))

				# dump the distance specific probabilities and expected contact count values
				# of the spline fitted and bias regression based models
				for (k in (1:length(uniq.idx))) {
					si <- uniq.idx[k]
					if (k < length(uniq.idx)) {
						ei <- uniq.idx[k+1] - 1
					} else {
						ei <- numPairs	# last read	
					}
					# now plot the variation of original probability vs the new probability
					# and the observed contact count vs the estimated contact count
					# for each specific distance values
					ProbValuePlotFile <- paste0(ProbComparePlotDir,'/Prob_Compare_Bias_Dist_', gene.dist[si], '.pdf')
					ContactCountPlotFile <- paste0(CCComparePlotDir,'/CC_Compare_Bias_Dist_', gene.dist[si], '.pdf')

					a <- data.frame(group = "Orig_Prob_Mult", x = seq(1,(ei-si+1)), y = log10(Prob_Val[si:ei]))
					b <- data.frame(group = "Est_Prob_Regr", x = seq(1,(ei-si+1)), y = log10(Prob_BiasRegr[si:ei]))
					curr_plotA <- ggplot(rbind(a,b), aes(x=x, y=y, fill=group, colour=group)) + geom_point(size=0.01) + xlab('Locus pairs') + ylab('Probability value (log10)')
					curr_plotA + ggtitle("Change in probability estimation")
					ggsave(ProbValuePlotFile, plot = curr_plotA, width=8, height=6)		

					c <- data.frame(group = "Orig_CC", x = seq(1,(ei-si+1)), y = interaction.data[si:ei, opt$cccol])
					d <- data.frame(group = "Est_CC_Regr", x = seq(1,(ei-si+1)), y = Exp_CC_BiasRegr[si:ei])
					curr_plotB <- ggplot(rbind.data.frame(c,d), aes(x=x, y=y, fill=group, colour=group)) + geom_point(size=0.01) + xlab('Locus pairs') + ylab('Contact count')
					curr_plotB + ggtitle("Change in contact count estimation")
					ggsave(ContactCountPlotFile, plot = curr_plotB, width=8, height=6)
				} # end distance loop for plotting
			}	# end dummy if


		} else {
			#====================
			# case 2.2: equal distance binning
			#====================

			# indices having distinct genomic distance
			# from the sorted genomic distance vector
			# with respect to the training data
			uniq.Training.idx <- order(Trainingdata_GeneDist)[!duplicated(Trainingdata_GeneDist)]

			# comment - sourya
			# sequential processing - loop
			# for (k in (1:length(uniq.Training.idx))) {

			# add - sourya
			# parallel processing - sourya
			time_start <- Sys.time()	# debug - souya

			res_Row_Training <- as.data.frame(parallel:::mclapply( 1:length(uniq.Training.idx) , mc.cores = ncore , function(k){

				# all the locations from training_si to training_ei will have the same probability 
				# since the gene distance value is the same
				training_si <- uniq.Training.idx[k]
				if (k < length(uniq.Training.idx)) {
					training_ei <- uniq.Training.idx[k+1] - 1
				} else {
					training_ei <- numTrainingPairs	# last read	
				}

				# compute the probability according to the genomic distance (of the start location - single element)
				# this probability value will be used for all the values within the interval (training_si to training_ei)
				# predict from the new spline fit - sourya
				pp <- predict(fit2_new, Trainingdata_GeneDist[training_si])
				training_prior_prob <- pp$y

				# current probability value
				Prob_Vec_Curr_Row <- rep(training_prior_prob, (training_ei - training_si + 1))

				if (0) {
					cat(sprintf("\n\n === bias correction - equal distance binning -- analyzing the bins:   k: %s  distance: %s  training_si: %s  training_ei: %s  training_prior_prob: %s  ", k, Trainingdata_GeneDist[training_si], training_si, training_ei, training_prior_prob))
				}

				# filter some interactions (training data) according to the bias thresholds
				Curr_Int_Data <- TrainingData[training_si:training_ei, ]
				BiasThrIdx <- which((Curr_Int_Data[,bias1.col] >= opt$BiasLowThr) & (Curr_Int_Data[,bias1.col] <= opt$BiasHighThr) & (Curr_Int_Data[,bias2.col] >= opt$BiasLowThr) & (Curr_Int_Data[,bias2.col] <= opt$BiasHighThr))
				Curr_Int_Set_BiasThr <- Curr_Int_Data[BiasThrIdx, ]

				if (0) {
					cat(sprintf("\n === Regression model (training data) - number of interactions: %s  number of bias threshold satisfying interactions: %s \n", nrow(Curr_Int_Data), nrow(Curr_Int_Set_BiasThr)))
				}

				# subset of the probability vector
				# according to the bias threshold pass
				Prob_Vec_Curr_Row_BiasFilt <- Prob_Vec_Curr_Row[BiasThrIdx]

				# variance within the contact count
				TSS_CC <- sum((Curr_Int_Set_BiasThr[, opt$cccol] - mean(Curr_Int_Set_BiasThr[, opt$cccol]))^2)

				# expected contact count according to the spline fit probability
				# computed with respect to the training contacts satisfying bias thresholds
				distCCVec <- TotTrainingDataContact * Prob_Vec_Curr_Row_BiasFilt
				logdistCCVec <- log10(distCCVec)

				if (opt$Resid == 0) {
					# if residual option is 0, model the regression 
					# between the observed contact count  
					# and the bias values
					bias_regression_data <- cbind.data.frame(log10(Curr_Int_Set_BiasThr[, bias1.col]), log10(Curr_Int_Set_BiasThr[, bias2.col]), log10(Curr_Int_Set_BiasThr[, opt$cccol]))
					colnames(bias_regression_data) <- c('logbias1', 'logbias2', 'logCC')

					# model the regression
					Linear_logbias_logCC <- lm(logCC ~ logbias1 + logbias2, data=bias_regression_data)
				
					# the model has three coefficients: 1) intercept, 2) for logbias1, 3) for logbias2
					coeff_intercept <- Linear_logbias_logCC$coefficients[1]
					coeff_logbias1 <- Linear_logbias_logCC$coefficients[2]
					coeff_logbias2 <- Linear_logbias_logCC$coefficients[3]

				} else {
					# here model the regression with respect to the 
					# (observed contact count - spline fitted contact count)
					# and the bias values
					# modification - sourya
					# previously the expression was wrong as : log10((Curr_Int_Set_BiasThr[, opt$cccol] - distCCVec)))					
					bias_regression_data <- cbind.data.frame(log10(Curr_Int_Set_BiasThr[, bias1.col]), log10(Curr_Int_Set_BiasThr[, bias2.col]), (log10(Curr_Int_Set_BiasThr[, opt$cccol]) - logdistCCVec))	
					colnames(bias_regression_data) <- c('logbias1', 'logbias2', 'logResidCC')

					# model the regression
					Linear_logbias_logCC <- lm(logResidCC ~ logbias1 + logbias2, data=bias_regression_data)

					# the model has three coefficients: 1) intercept, 2) for logbias1, 3) for logbias2
					coeff_intercept <- Linear_logbias_logCC$coefficients[1]
					coeff_logbias1 <- Linear_logbias_logCC$coefficients[2]
					coeff_logbias2 <- Linear_logbias_logCC$coefficients[3]

				}
				
				# get the performance statistics of the regression model
				RSS_LM <- sum((10^residuals(Linear_logbias_logCC))^2)
				R.square_LM <- 1 - ((RSS_LM * 1.0) / TSS_CC)
				AIC_LM <- AIC(Linear_logbias_logCC)		

				# form a vector of the current distance and regression model cofficients
				# currently the vector content is same for both residual and non-residual model
				currvec <- c(Trainingdata_GeneDist[training_si], coeff_intercept, coeff_logbias1, coeff_logbias2, RSS_LM, R.square_LM, AIC_LM)

				# # comment - sourya
				# # was used in sequential processing - sourya
				# # append the vector to an existing (or new) data frame
				# if (rownum == 1) {
				# 	Regression_model_coeff.df <- currvec
				# } else {
				# 	Regression_model_coeff.df <- rbind(Regression_model_coeff.df, currvec)
				# }

				# add - Sourya
				# for parallel processing
				# we return the row number
				# and the regression model vectors
				return(c(k, currvec))
			
			# add - sourya
			# for parallel processing
			} ))

			# # comment - sourya
			# # was used in sequential processing - sourya
			# }	# end unique distance processing - regression model

			# add - sourya
			# for parallel processing
			# now form a data frame with transpose of the contents of res_Row_Training
			# from row 2 to last row, and for all columns (transpose)
			# since the row 1 contains the row number values
			# t() function returns the transpose
			Regression_model_coeff.df <- t(res_Row_Training[2:nrow(res_Row_Training), 1:ncol(res_Row_Training)])
			
			# write the training data specific regression model 
			# coefficients for individual distance values
			# currently the vector content is same for both residual and non-residual model
			Estimated_Model_Coeff_File <- paste0(outdir,'/Regression_Model_Coeff.log')	
			write.table(Regression_model_coeff.df, Estimated_Model_Coeff_File, row.names = FALSE, col.names = c("Dist", "Intercept", "LogBias1", "LogBias2", "RSS_LM", "R.square_LM", "AIC_LM"), sep = "\t", quote=FALSE, append=FALSE)

			# debug - sourya
			time_end <- Sys.time()

			if (0) {
				cat(sprintf("\n\n\n ==>>> Processing bin rows - res_Row_Training - parallel optimization - new time: %s \n\n", (time_end - time_start)))
			}

			#============================
			# fit spline for individual distance values
			# vs three parameters:
			# 1) coefficient of intercept
			# 2) coefficient of logbias1
			# 3) coefficient of logbias2
			# first model the spline with respect to the available binned data
			# next extend the spline with all possible distance values (in the vector TestDistVal)
			fit_spline_coeff_Intercept <- smooth.spline(Regression_model_coeff.df[,1], Regression_model_coeff.df[,2], cv=TRUE)
			pp.fit_spline_coeff_Intercept <- predict(fit_spline_coeff_Intercept, TestDistVal)
			fit_spline_coeff_Intercept_new <- smooth.spline(pp.fit_spline_coeff_Intercept$x, pp.fit_spline_coeff_Intercept$y, cv=TRUE)

			# plot the fitted spline "fit_spline_coeff_Intercept_new"
			if (0) { #(opt$Draw) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Intercept_new.pdf')	
				pdf(plotfile1, width=8, height=6)
				# par(mar=c(5,5,2,2)+0.1)
				plot(Regression_model_coeff.df[,1], Regression_model_coeff.df[,2], cex=0.5, col="black", xlab="Average interaction distance", ylab="Regression coefficient - intercept")
				lines(fit_spline_coeff_Intercept_new$x, fit_spline_coeff_Intercept_new$y, col="red", lwd=0.5)
				title("fit_spline_coeff_Intercept")
				dev.off()
			}

			# plot the fitted spline "fit_spline_coeff_Intercept_new" in log scale
			if (0) { #(opt$Draw) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Intercept_new_LOGScale.pdf')	
				pdf(plotfile1, width=8, height=6)
				# par(mar=c(5,5,2,2)+0.1)
				plot(log10(Regression_model_coeff.df[,1]), log10(Regression_model_coeff.df[,2]), cex=0.5, col="black", xlab="Average interaction distance (log)", ylab="Regression coefficient - intercept (log)")
				lines(log10(fit_spline_coeff_Intercept_new$x), log10(fit_spline_coeff_Intercept_new$y), col="red", lwd=0.5)
				title("fit_spline_coeff_Intercept - log scale")
				dev.off()
			}

			fit_spline_coeff_Logbias1 <- smooth.spline(Regression_model_coeff.df[,1], Regression_model_coeff.df[,3], cv=TRUE)
			pp.fit_spline_coeff_Logbias1 <- predict(fit_spline_coeff_Logbias1, TestDistVal)
			fit_spline_coeff_Logbias1_new <- smooth.spline(pp.fit_spline_coeff_Logbias1$x, pp.fit_spline_coeff_Logbias1$y, cv=TRUE)

			# plot the fitted spline "fit_spline_coeff_Logbias1_new"
			if (0) { #(opt$Draw) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Logbias1_new.pdf')	
				pdf(plotfile1, width=8, height=6)
				# par(mar=c(5,5,2,2)+0.1)
				plot(Regression_model_coeff.df[,1], Regression_model_coeff.df[,3], cex=0.5, col="black", xlab="Average interaction distance", ylab="Regression coefficient - logbias1")
				lines(fit_spline_coeff_Logbias1_new$x, fit_spline_coeff_Logbias1_new$y, col="red", lwd=0.5)
				title("fit_spline_coeff_logbias1")
				dev.off()
			}

			# plot the fitted spline "fit_spline_coeff_Logbias1_new" in log scale
			if (0) { #(opt$Draw) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Logbias1_new_LOGScale.pdf')	
				pdf(plotfile1, width=8, height=6)
				# par(mar=c(5,5,2,2)+0.1)
				plot(log10(Regression_model_coeff.df[,1]), log10(Regression_model_coeff.df[,3]), cex=0.5, col="black", xlab="Average interaction distance (log)", ylab="Regression coefficient - logbias1 (log)")
				lines(log10(fit_spline_coeff_Logbias1_new$x), log10(fit_spline_coeff_Logbias1_new$y), col="red", lwd=0.5)
				title("fit_spline_coeff_logbias1 - log scale")
				dev.off()
			}

			fit_spline_coeff_Logbias2 <- smooth.spline(Regression_model_coeff.df[,1], Regression_model_coeff.df[,4], cv=TRUE)
			pp.fit_spline_coeff_Logbias2 <- predict(fit_spline_coeff_Logbias2, TestDistVal)
			fit_spline_coeff_Logbias2_new <- smooth.spline(pp.fit_spline_coeff_Logbias2$x, pp.fit_spline_coeff_Logbias2$y, cv=TRUE)

			# plot the fitted spline "fit_spline_coeff_Logbias2_new"
			if (0) { #(opt$Draw) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Logbias2_new.pdf')	
				pdf(plotfile1, width=8, height=6)
				# par(mar=c(5,5,2,2)+0.1)
				plot(Regression_model_coeff.df[,1], Regression_model_coeff.df[,4], cex=0.5, col="black", xlab="Average interaction distance", ylab="Regression coefficient - logbias2")
				lines(fit_spline_coeff_Logbias2_new$x, fit_spline_coeff_Logbias2_new$y, col="red", lwd=0.5)
				title("fit_spline_coeff_logbias2")
				dev.off()
			}

			# plot the fitted spline "fit_spline_coeff_Logbias2_new" in log scale
			if (0) { #(opt$Draw) {
				plotfile1 <- paste0(outdir,'/fit_spline_coeff_Logbias2_new_LOGScale.pdf')	
				pdf(plotfile1, width=8, height=6)
				# par(mar=c(5,5,2,2)+0.1)
				plot(log10(Regression_model_coeff.df[,1]), log10(Regression_model_coeff.df[,4]), cex=0.5, col="black", xlab="Average interaction distance (log)", ylab="Regression coefficient - logbias2 (log)")
				lines(log10(fit_spline_coeff_Logbias2_new$x), log10(fit_spline_coeff_Logbias2_new$y), col="red", lwd=0.5)
				title("fit_spline_coeff_logbias2 - log scale")
				dev.off()
			}
			#============================

			# debug - sourya
			time_start <- Sys.time()

			# using the fitted splines for the distance / probability
			# and the bias correction coefficients
			# we now use the full set of input interactions
			# to model the expected contact count and the corresponding bias correcte probability 
			for (k in (1:length(uniq.idx))) {
				# all the locations from si to ei will have the same probability 
				# with respect to fitted spline (FitHiC distance)
				si <- uniq.idx[k]
				if (k < length(uniq.idx)) {
					ei <- uniq.idx[k+1] - 1
				} else {
					ei <- numPairs	# last read	
				}

				# probability of spline fit for average distance vs avg prior probability
				pp <- predict(fit2_new, gene.dist[si])
				p <- pp$y

				# also predict the coefficients of intercept, logbias1, and logbias2
				# from different splines with respect to this distance
				pp_I <- predict(fit_spline_coeff_Intercept_new, gene.dist[si])
				coeff_Intcpt <- pp_I$y
				pp_L1 <- predict(fit_spline_coeff_Logbias1_new, gene.dist[si])
				coeff_LogBias1 <- pp_L1$y
				pp_L2 <- predict(fit_spline_coeff_Logbias2_new, gene.dist[si])
				coeff_LogBias2 <- pp_L2$y

				if (0) {
					cat(sprintf("\n *** Processing uniq distance idx: %s  si: %s  ei: %s num elem : %s Spline prob: %s  coeff_Intcpt : %s  coeff_LogBias1 : %s   coeff_LogBias2 : %s ", k, si, ei, (ei-si+1), p, coeff_Intcpt, coeff_LogBias1, coeff_LogBias2))
				}

				# copy the spline based probability 
				Prob_Val[si:ei] <- p

				# get indices having zero and non zero bias 
				# Note: here we add the quantity (si-1) since the which operator returns indices from 1
				# to get the indices mapped back on the range si:ei, we add this offset
				nonzero_biasidx_set <- which((interaction.data[si:ei, bias1.col] > 0) & (interaction.data[si:ei, bias2.col] > 0)) + (si-1)
				# the zero bias values are those not in nonzero_biasidx_set
				zero_biasidx_set <- setdiff(seq(si, ei), nonzero_biasidx_set)
				cat(sprintf(" ---  length nonzero_biasidx_set: %s  length zero_biasidx_set : %s ", length(nonzero_biasidx_set), length(zero_biasidx_set)))

				# if one of the bias values are zero
				# then just use the spline predicted probability as the estimated distance
				if (length(zero_biasidx_set) > 0) {
					Exp_CC_BiasRegr[zero_biasidx_set] <- (TotContact * p)	
				}

				# for non-zero bias values, use the bias regression predicted contact count
				if (length(nonzero_biasidx_set) > 0) {
					if (opt$Resid == 0) {
						Exp_CC_BiasRegr[nonzero_biasidx_set] <- 10^(coeff_Intcpt + coeff_LogBias1 * log10(interaction.data[nonzero_biasidx_set, bias1.col]) + coeff_LogBias2 * log10(interaction.data[nonzero_biasidx_set, bias2.col]))
					} else {
						# here residual contact count is predicted
						# add the spline fit contact count to the predicted value
						Exp_CC_BiasRegr[nonzero_biasidx_set] <- 10^(coeff_Intcpt + coeff_LogBias1 * log10(interaction.data[nonzero_biasidx_set, bias1.col]) + coeff_LogBias2 * log10(interaction.data[nonzero_biasidx_set, bias2.col]) + log10(TotContact * p))
					}					
				}
				
			}	# end distance loop 

			if (0) {
				cat(sprintf("\n\n *** modeled the bias regression based probability and expected contact count **** \n\n"))
			}

			# debug - sourya
			time_end <- Sys.time()

			if (0) {
				cat(sprintf("\n\n\n ==>>> modeled the bias regression -- time: %s \n\n", (time_end - time_start)))
			}

			# now use the expected contact count from the bias correction 
			# to model the binomial distribution
			NumElem <- length(Exp_CC_BiasRegr)
			ExpTotContact <- as.integer(sum(Exp_CC_BiasRegr))

			# debug - sourya
			time_start <- Sys.time()

			# parallel execution to model the binomial distribution
			# with respect to the expected contact count
			# and the expected sum of contact counts
			# such separate execution is required
			# to model the pdf integral sum as 1

			# process in subsets 
			for (k in (1:length(uniq.idx))) {
				# all the locations from si to ei will have the same probability 
				# with respect to fitted spline (FitHiC distance)
				si <- uniq.idx[k]
				if (k < length(uniq.idx)) {
					ei <- uniq.idx[k+1] - 1
				} else {
					ei <- numPairs	# last read	
				}

				if (0) {
					cat(sprintf("\n *** Modeling regression bias probability based P value --- uniq distance idx: %s  si: %s  ei: %s num elem : %s ", k, si, ei, (ei-si+1)))
				}

				# parallel processing
				resbias <- as.data.frame(parallel:::mclapply( si:ei , mc.cores = ncore , function(idx){
					curr_cc <- interaction.data[idx, opt$cccol]
					pr_bias <- (Exp_CC_BiasRegr[idx] * 1.0) / ExpTotContact
					db2 <- dbinom(curr_cc, size=ExpTotContact, prob=pr_bias)
					pb2 <- pbinom(curr_cc, size=ExpTotContact, prob=pr_bias, lower.tail=FALSE)
					return(c(idx,pr_bias,db2,(db2+pb2)))
				} ))

				# merge all the probability values
				# and the binomial distribution
				# in the final array
				Prob_BiasRegr[si:ei] <- as.double(resbias[2,1:ncol(resbias)])
				Spline_Binom_Prob_CC[si:ei] <- as.double(resbias[3,1:ncol(resbias)])
				Spline_Binom_P_Val_CC[si:ei] <- as.double(resbias[4,1:ncol(resbias)])
			}	# end distance loop

			if (0) {
				cat(sprintf("\n\n *** modeled the binomial distribution with respect to bias regression probability **** \n\n"))
			}

			# debug - sourya
			time_end <- Sys.time()
			if (0) {
				cat(sprintf("\n\n\n ==>>> modeled the binomial distribution with respect to bias regression -- time: %s \n\n", (time_end - time_start)))
			}
			
			# temporary plotting the probability and contact count plots
			if (0) {
				# drawing the comparison between original bias correction 
				# and the new bias correction (regression) method
				ProbComparePlotDir <- paste0(outdir, '/Prob_Compare_Plots')
				system(paste("mkdir -p", ProbComparePlotDir))
				CCComparePlotDir <- paste0(outdir, '/CC_Compare_Plots')
				system(paste("mkdir -p", CCComparePlotDir))

				# dump the distance specific probabilities and expected contact count values
				# of the spline fitted and bias regression based models
				for (k in (1:length(uniq.idx))) {
					si <- uniq.idx[k]
					if (k < length(uniq.idx)) {
						ei <- uniq.idx[k+1] - 1
					} else {
						ei <- numPairs	# last read	
					}
					# now plot the variation of original probability vs the new probability
					# and the observed contact count vs the estimated contact count
					# for each specific distance values
					ProbValuePlotFile <- paste0(ProbComparePlotDir,'/Prob_Compare_Bias_Dist_', gene.dist[si], '.pdf')
					ContactCountPlotFile <- paste0(CCComparePlotDir,'/CC_Compare_Bias_Dist_', gene.dist[si], '.pdf')

					a <- data.frame(group = "Orig_Prob_Mult", x = seq(1,(ei-si+1)), y = log10(Prob_Val[si:ei]))
					b <- data.frame(group = "Est_Prob_Regr", x = seq(1,(ei-si+1)), y = log10(Prob_BiasRegr[si:ei]))
					curr_plotA <- ggplot(rbind(a,b), aes(x=x, y=y, fill=group, colour=group)) + geom_point(size=0.01) + xlab('Locus pairs') + ylab('Probability value (log10)')
					curr_plotA + ggtitle("Change in probability estimation")
					ggsave(ProbValuePlotFile, plot = curr_plotA, width=8, height=6)		

					c <- data.frame(group = "Orig_CC", x = seq(1,(ei-si+1)), y = interaction.data[si:ei, opt$cccol])
					d <- data.frame(group = "Est_CC_Regr", x = seq(1,(ei-si+1)), y = Exp_CC_BiasRegr[si:ei])
					curr_plotB <- ggplot(rbind.data.frame(c,d), aes(x=x, y=y, fill=group, colour=group)) + geom_point(size=0.01) + xlab('Locus pairs') + ylab('Contact count')
					curr_plotB + ggtitle("Change in contact count estimation")
					ggsave(ContactCountPlotFile, plot = curr_plotB, width=8, height=6)
				} # end distance loop for plotting
			}	# end dummy if

		}	# end condition of either distance or equal occupancy binning
	}	# end condition multiply probability by bias value
}	# end bias correction condition

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time for Spline based distribution compute (p value) (sorted interaction file): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}

# from the generated P values, obtain the Q value using BH correction
Spline_Binom_QVal <- p.adjust(Spline_Binom_P_Val_CC, method = "BH")	

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time to compute Q value (from P value) in spline: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}


# sourya - comment 
if (0) {

	# accumulate all results - also add header information
	if ((opt$BiasCorr == 0) | ((opt$BiasCorr == 1) & (opt$MultBias == 1))) {
		# either there is no bias correction
		# or bias correction is done by multiplying the probabilities
		FinalData <- cbind(interaction.data, Prob_Val, Spline_Binom_Prob_CC, Spline_Binom_P_Val_CC, Spline_Binom_QVal)
		colnames(FinalData) <- c(colnames(interaction.data), "p", "dbinom", "P-Value", "Q-Value")	
	} else {
		# here accumulate both the original spline fit probability
		# and also the probability obtained from the bias regression
		# and the expected contact count from the bias regression model
		# P and Q values are similar
		FinalData <- cbind(interaction.data, Prob_Val, Exp_CC_BiasRegr, Prob_BiasRegr, Spline_Binom_Prob_CC, Spline_Binom_P_Val_CC, Spline_Binom_QVal)	
		colnames(FinalData) <- c(colnames(interaction.data), "p", "exp_cc_Bias", "p_Bias", "dbinom_Bias", "P-Value_Bias", "Q-Value_Bias")	

		# append the Spline distribution probability and corresponding P value as separate columns
		# and write in a separate text file
		temp_outfile <- paste0(OutIntDir, '/', 'temp_out.bed')
		write.table(FinalData, temp_outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)	
	}

	# now sort the file contents and write that in the final specified output file
	system(paste('sort -k1,1 -k2,2n -k5,5n', paste0('-k',opt$cccol,',',opt$cccol,'nr'), temp_outfile, '>', opt$OutFile))

}	# end dummy if

# add - sourya
if (1) {

	if ((opt$BiasCorr == 0) | ((opt$BiasCorr == 1) & (opt$MultBias == 1))) {
		FinalData <- cbind.data.frame(Prob_Val, Spline_Binom_Prob_CC, Spline_Binom_P_Val_CC, Spline_Binom_QVal)
		colnames(FinalData) <- c("p", "dbinom", "P-Value", "Q-Value")	
	} else {
		FinalData <- cbind.data.frame(Prob_Val, Exp_CC_BiasRegr, Prob_BiasRegr, Spline_Binom_Prob_CC, Spline_Binom_P_Val_CC, Spline_Binom_QVal)
		colnames(FinalData) <- c("p", "exp_cc_Bias", "p_Bias", "dbinom_Bias", "P-Value_Bias", "Q-Value_Bias")
	}
	temp_outfile_2 <- paste0(OutIntDir, '/temp_out2.bed')
	write.table(FinalData, temp_outfile_2, row.names=F, col.names=T, sep="\t", quote=F, append=F)

	temp_outfile <- paste0(OutIntDir, '/temp_out.bed')
	system(paste("paste", opt$InpFile, temp_outfile_2, ">", temp_outfile))

	# delete the temporary output file
	system(paste("rm", temp_outfile_2))

	# now sort the file contents and write that in the final specified output file
	system(paste('sort -k1,1 -k2,2n -k5,5n', paste0('-k',opt$cccol,',',opt$cccol,'nr'), temp_outfile, '>', opt$OutFile))

}	# end dummy if

# end add - sourya

# delete the temporary output file
system(paste('rm', temp_outfile))

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time to write the interaction data (with Q value) for spline based distribution: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}

# plot the variation of q-value from 0 to 0.1, in the steps of 0.01
# and corresponding number of interaction count per range of q-values
qvalue_seq <- seq(0.01,0.1,by=0.01)
IntCount <- c()
LogIntCount <- c()
for (i in (1:length(qvalue_seq))) {
	cnt <- length(which(FinalData[, ncol(FinalData)] <= qvalue_seq[i]))
	IntCount <- c(IntCount, cnt)
	LogIntCount <- c(LogIntCount, log2(cnt+1))
}

plotfile1 <- paste0(outdir,'/Interaction_vs_qval.png')
a <- data.frame(group = paste("Interactions below FDR threshold"), x = qvalue_seq, y = IntCount)
curr_plotA <- ggplot(a, aes(x=x, y=y, fill=group, colour=group)) + geom_line(color="blue") + xlab('FDR threshold') + ylab('No of significant interactions') + xlim(min(qvalue_seq), max(qvalue_seq)) + ylim(0, (max(IntCount) + 10))
curr_plotA + ggtitle("Number of significant interactions below various FDR threshold values")
ggsave(plotfile1, plot = curr_plotA, width=8, height=6)		

plotfile1 <- paste0(outdir,'/Interaction_log2_vs_qval.png')
a <- data.frame(group = paste("Interactions below FDR threshold"), x = qvalue_seq, y = LogIntCount)
curr_plotA <- ggplot(a, aes(x=x, y=y, fill=group, colour=group)) + geom_line(color="blue") + xlab('FDR threshold') + ylab('No of significant interactions (log2)') + xlim(min(qvalue_seq), max(qvalue_seq)) + ylim(0, (max(LogIntCount) + 0.25))
curr_plotA + ggtitle("Number of significant interactions (log2) below various FDR threshold values")
ggsave(plotfile1, plot = curr_plotA, width=8, height=6)		

