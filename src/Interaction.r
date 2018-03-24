#!/usr/bin/env Rscript

#===========================================================
# R script for analyzing the statistical significance of interactions
# this is an implementation of the paper Ay et. al. 2014 (FitHiC)
# incorporates spline fit for raw contact count
# also implements the bias calculation and the bias based filtering and spline fit procedure

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================

library(splines)
library(fdrtool)
library(parallel)	# library for parallel processing

library(optparse)

# #===========================
# # parallel computation of binomial distribution on input vector
# parallel.BinomDistr <- function(inpvec) {
# 	idx <- inpvec[1]
# 	cc <- inpvec[2]
# 	tc <- inpvec[3]
# 	currprob <- inpvec[4]

# 	curr_dbinom <- dbinom(cc, size=tc, prob=currprob)
# 	curr_pbinom <- pbinom(cc, size=tc, prob=currprob, lower.tail=FALSE)
# 	return c(idx, curr_dbinom, (curr_dbinom + curr_pbinom))
# }

# #===========================

option_list = list(
  	make_option(c("--InpFile"), type="character", default=NULL, help="File with interactions + normalization features among individual genomic bins, preferably sorted with respect to interaction distance"),
  	make_option(c("--headerInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input interaction file has the header. Default FALSE."),
	make_option(c("--OutFile"), type="character", default=NULL, help="Output file name storing interactions + probability, P value and Q value"),
	make_option(c("--P2P"), type="integer", action="store", default=0, help="Can be 0 or 1. If 1, spline fit is modeled using only peak to peak interactions. Default 0."),
  	make_option(c("--EqLenBin"), type="logical", action="store_true", default=FALSE, help="If TRUE, Equal length bins are used. Default FALSE, means that equal occupancy bins are used. Users need not alter this parameter"),
  	make_option(c("--BiasCorr"), type="integer", action="store", default=0, help="If 0, raw contact count is used. Else, bias corrected contact count is used before applying FitHiC. Default 0."),
	make_option(c("--BiasLowThr"), type="numeric", default=0.2, help="Lower threshold of bias. Default 0.2"),
	make_option(c("--BiasHighThr"), type="numeric", default=5, help="Higher threshold of bias. Default 5"),  	
  	make_option(c("--nbins"), type="integer", action="store", default=200, help="Number of bins employed for FitHiC.", metavar="nbins"),
  	make_option(c("--Draw"), type="logical", action="store_true", default=FALSE, help="If TRUE, spline fit of FitHiC is plotted. Default FALSE."),
  	make_option(c("--cccol"), type="integer", action="store", default=7, help="Column number of the file storing the contact count"),
  	make_option(c("--TimeFile"), type="character", default=NULL, help="If specified, denotes the file which will contain time profiling information"),
  	make_option(c("--BiasFilt"), type="integer", action="store", default=0, help="If 1, interactions are filtered according to BiasLowThr and BiasHighThr, before processing. Default 0"),
	make_option(c("--ProbBias"), type="integer", action="store", default=0, help="If 1, probability values are multiplied with the bias values. Default 0.")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

nbins <- as.integer(opt$nbins) # number of bins used in FitHiC
if (is.null(opt$TimeFile)) {
	timeprof <- 0
} else {
	timeprof <- 1
}

# peak to peak background usage option
Peak2PeakBackg <- as.integer(opt$P2P)
cat(sprintf("\n Peak to peak background usage for spline fit: %s ", Peak2PeakBackg))

# number of processors within the system
ncore <- detectCores()
cat(sprintf("\n Number of cores in the system: %s ", ncore))

# directory containing the input interaction file
inpdir <- dirname(opt$InpFile)

# directory to contain the spline fitted output interaction file
OutIntDir <- dirname(opt$OutFile)

# this directory will store the spline graph according to the model
outdir <- paste0(OutIntDir,'/Results')
system(paste('mkdir -p', outdir))

# this file stores the spline fitted model
if (opt$EqLenBin) {
	plotfile <- paste0(outdir,'/','EqLenBin_SplinePass1.pdf')
} else {
	plotfile <- paste0(outdir,'/','EqOccBin_SplinePass1.pdf')	
}

if (timeprof == 1) {
	starttime <- Sys.time()
}

# load the interaction data matrix
# Note: the data has a header information
interaction.data <- read.table(opt$InpFile, header=opt$headerInp)
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

cat(sprintf("\n Number of unfiltered interactions: %s ", nrow(interaction.data)))

# return if the number of interaction is 0
if (nrow(interaction.data) == 0) {
	cat(sprintf("\n No of interactions : 0 - exit !!!"))
	return()
}

if (ncol(interaction.data) < opt$cccol) {
	cat(sprintf("\n There is no contact count column or the formatting of the input file has problems - exit !!!"))
	return()
}

# vector of initial (prior) probabilities for contact counts
Prob_Val <- c()

# probability of the observed contact count for a single locus pair, from the spline fit
Spline_Binom_Prob_CC <- c()

# P-value of the observed contact count for a single locus pair
# binomial distribution + spline based estimation
Spline_Binom_P_Val_CC <- c()

# for each bin, count the no of distinct locus pairs
no_distinct_loci <- c()

# for each bin, no of observed contacts
NumContact <- c()

# for each bin, stores the average contact count per locus pair
avg_contact <- c()

# prior contact probability for a specific locus pair falling within a bin
prior_contact_prob <- c()

# average interaction distance for all locus pairs falling within a bin
avg_int_dist <- c()

#======================================================
# if normalization using bias values is enabled 
# then discard the interactions where either loci has bias values outside specified thresholds
if ((opt$BiasCorr == 1) & (opt$BiasFilt == 1)) {
	
	# specifying the columns in the interaction data 
	# where the bias information is provided
	bias1.col <- 10
	bias2.col <- 16

	if (ncol(interaction.data) >= bias2.col) {
		interaction.data <- interaction.data[which((interaction.data[,bias1.col] >= opt$BiasLowThr) & (interaction.data[,bias1.col] <= opt$BiasHighThr) & (interaction.data[,bias2.col] >= opt$BiasLowThr) & (interaction.data[,bias2.col] <= opt$BiasHighThr)), ]
		cat(sprintf("\n Bias specific filtering is enabled - Number of interactions where both loci satisfy bias criterion: %s ", nrow(interaction.data)))
	}
}

#======================================================
# absolute genomic distance for an interaction instance
gene.dist <- abs(interaction.data[,2] - interaction.data[,5])

# no of interactions pairs
numPairs <- length(gene.dist)
cat(sprintf("\n *** Total number of interactions: %s ", numPairs))

# total number of contacts for all the interactions
TotContact <- sum(interaction.data[,opt$cccol])

if (timeprof == 1) {
	starttime <- Sys.time()
}

#======================================================
# add - sourya
# we form the training data out of the total number of interactions
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
	cat(sprintf("\n No training data for FitHiC spline fit - error - return !!! "))
	return()
}

cat(sprintf("\n *** Total number of training interactions: %s ", numTrainingPairs))

# interaction distance vector for the training data only
Trainingdata_GeneDist <- abs(TrainingData[,2] - TrainingData[,5])

# contact count sum for the training data
TotTrainingDataContact <- sum(TrainingData[,opt$cccol])

#=====================================================
# divide the genomic distance into b quantiles (number of bins)

# Note: for equal length bins (opt$EqLenBin is TRUE), each bin would have equal no of entries (pairs of loci)
# from the training data
# on the other hand, for equal occupancy bins (opt$EqLenBin is FALSE), each bin would have equal no of contact count
# from the training data
# In both cases, however, corresponding bin intervals will be variable

# error condition - sourya
# if the number of interactions is less than the argument 'nbins'
# then re-adjust the values of nbins and then assign the 
# no of entries in each of the bins

if (opt$EqLenBin) {
	# no of entries (occupancies) for each bin 
	# each bin have same no of locus pairs
	if (numTrainingPairs < nbins) {
		nbins <- numTrainingPairs
		nentry <- 1
	} else {
		nentry <- floor(numTrainingPairs / nbins)
	}
} else {
	# no of contacts for each bin
	# each bin have same no of contacts
	if (TotTrainingDataContact < nbins) {
		nbins <- numTrainingPairs
	}
	ncontactbin <- floor(TotTrainingDataContact / nbins)
}

#==============================================

if (opt$EqLenBin) {
	# equal length bin case (equal locus entries)
	for (bin_idx in (1:nbins)) {
		si <- (bin_idx - 1) * nentry + 1
		if (bin_idx < nbins) {
			ei <- si + nentry - 1		
		} else {
			ei <- numTrainingPairs	# last read
		}
		no_distinct_loci[bin_idx] <- (ei-si+1)
		NumContact[bin_idx] <- sum(TrainingData[si:ei,opt$cccol])
		avg_contact[bin_idx] <- NumContact[bin_idx] / no_distinct_loci[bin_idx]
		prior_contact_prob[bin_idx] <- (avg_contact[bin_idx] / TotTrainingDataContact)
		avg_int_dist[bin_idx] <- mean(Trainingdata_GeneDist[si:ei])
	}
} else {
	# equal occupancy (contact count) bin - preferred

	# keeps track of total no of contact counts for a particular bin
	cumContactCount <- 0
	nelem <- 0
	ei <- 0
	bin_idx <- 0

	# keeps track of the no of different genomic distance values encountered within this bin
	# and also the no of elements (locus pairs) having that particular genomic distance
	ValGeneDist <- c()
	nelemSameGeneDist <- c()

	while (ei < numTrainingPairs) {
		si <- ei + 1
		curr_gene_dist <- abs(TrainingData[si,2] - TrainingData[si,5])
		idx_list <- which(Trainingdata_GeneDist == curr_gene_dist)
		nelem <- nelem + length(idx_list)
		ei <- ei + length(idx_list)
		cumContactCount <- cumContactCount + sum(TrainingData[si:ei, opt$cccol])
		
		ValGeneDist <- c(ValGeneDist, curr_gene_dist)
		nelemSameGeneDist <- c(nelemSameGeneDist, length(idx_list))

		if (cumContactCount >= ncontactbin) {
			# current cumulative contact exceeds the expected average contact per bin
			# and also, all the contacts with this specified genomic distance is covered
			# so, we fix this bin 
			bin_idx <- bin_idx + 1
			no_distinct_loci[bin_idx] <- nelem
			NumContact[bin_idx] <- cumContactCount
			avg_contact[bin_idx] <- cumContactCount / nelem
			prior_contact_prob[bin_idx] <- (avg_contact[bin_idx] / TotTrainingDataContact)

			avg_int_dist[bin_idx] <- 0
			for (j in (1:length(ValGeneDist))) {
				avg_int_dist[bin_idx] <- avg_int_dist[bin_idx] + ((nelemSameGeneDist[j] * 1.0) / nelem) * ValGeneDist[j]
			}

			# reset the couners
			cumContactCount <- 0
		 	ValGeneDist <- c()
		 	nelemSameGeneDist <- c()
		 	nelem <- 0
		}
	}
	# now dump the data in a text file
	OutBinfile <- paste0(OutIntDir, '/', 'Bin_Info.log')
	write.table(cbind(avg_int_dist, no_distinct_loci, NumContact, avg_contact, prior_contact_prob), OutBinfile, row.names = FALSE, col.names = c("AvgDist","NumLoci","NumContact","AvgContact","PriorProb"), sep = "\t", quote=FALSE, append=FALSE)	
}

cat(sprintf("\n after printing the bin log file"))

#=====================================
# now prepare a spline fit where x axis is the 'avg_int_dist'
# and y axis is the 'prior_contact_prob'
#=====================================

# # the first spline works according to the specified degree of freedom (16 - for testing)
# fit <- smooth.spline(avg_int_dist, prior_contact_prob, df=16)

# if (timeprof == 1) {
# 	endtime <- Sys.time()
# 	fp <- file(opt$TimeFile, open="a")
# 	outstr <- paste('\n Time for spline (df=16): ', (endtime - starttime))
# 	write(outstr, file=fp, append=T)
# 	close(fp)
# 	starttime <- Sys.time()
# }

# the spline works according to the cross validation principle
fit2 <- smooth.spline(avg_int_dist, prior_contact_prob, cv=TRUE)

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time for spline (CV): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
	starttime <- Sys.time()
}

# # perform anti-tonic regression on the first spline
# fit.mr <- monoreg(fit$x, fit$y, type="antitonic")

# if (timeprof == 1) {
# 	endtime <- Sys.time()
# 	fp <- file(opt$TimeFile, open="a")
# 	outstr <- paste('\n Time for antitonic regression on spline (df = 16): ', (endtime - starttime))
# 	write(outstr, file=fp, append=T)
# 	close(fp)
# 	starttime <- Sys.time()
# }

# perform anti-tonic regression on the second spline
# antitonic = monotonic decreasing (as the probability decreases with increasing distance)
fit2.mr <- monoreg(fit2$x, fit2$y, type="antitonic")

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time for antitonic regression on spline (CV): ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}

# # plot the data, if plotting is enabled via the parameter opt$Draw
# if (opt$Draw) {
# 	pdf(plotfile, width=14, height=10)
# 	plot(avg_int_dist, prior_contact_prob, cex=0.5, col="darkgrey", xlab="Average interaction distance", ylab="Prior contact probability")
# 	title("Smooth spline - antitonic regression - pass 1")
# 	lines(fit.mr, col="red", lwd=2)
# 	lines(fit2.mr, col="blue", lwd=2)
# 	legend("topright",legend=c("16 DF", paste(as.integer(fit2$df), "df (cv)")), col=c("red", "blue"), lty=1, lwd=2, cex=0.8)
# 	dev.off()		
# }

#=========================
# now use this fitted spline
# to get predicted probability values for all distance values of the interactions
# in the unit of bin size employed in the current file)
binsize <- abs(interaction.data[1,2] - interaction.data[1,5])
TestDistVal <- seq(min(gene.dist), max(gene.dist), binsize)

# probability obtained from the spline
pp.fit2 <- predict(fit2, TestDistVal)

# use this computed probability and the input distance to fit a spline again
fit2_new <- smooth.spline(pp.fit2$x, pp.fit2$y, cv=TRUE)

# perform antitonic regression on this spline
fit2_new.mr <- monoreg(fit2_new$x, fit2_new$y, type="antitonic")

# plot the fitted spline and the smoothing regression
if (opt$Draw) {
	plotfile1 <- paste0(outdir,'/','EqOccBin_SplinePass1.pdf')	
	pdf(plotfile1, width=10, height=8)
	plot(avg_int_dist, prior_contact_prob, cex=0.5, col="black", xlab="Average interaction distance", ylab="Prior contact probability", xlim=c(0, max(gene.dist)))
	lines(fit2$x, fit2$y, col="yellow", lwd=0.5)
	lines(fit2.mr$x, fit2.mr$yf, col="blue", lwd=0.5)
	lines(fit2_new$x, fit2_new$y, col="green", lwd=0.5)
	lines(fit2_new.mr$x, fit2_new.mr$yf, col="red", lwd=0.5)
	legend("topright",legend=c(paste("spline (fit2) - ", as.integer(fit2$df), "df (cv)"), "MR on spline (fit2.mr)", paste("5 Kb uniform spline (fit2_new) - ", as.integer(fit2_new$df), "df (cv)"), "MR on 5 Kb uniform spline (fit2_new.mr)"), col=c("yellow", "blue", "green", "red"), lty=1, lwd=1, cex=0.8)
	title("Smooth spline - antitonic regression - pass 1")
	dev.off()
}

# plot the fitted spline and the smoothing regression
# in log scale
if (opt$Draw) {
	plotfile1 <- paste0(outdir,'/','EqOccBin_SplinePass1_LOGScale.pdf')	
	pdf(plotfile1, width=10, height=8)
	plot(avg_int_dist, log10(prior_contact_prob), cex=0.5, col="black", xlab="Average interaction distance", ylab="Prior contact probability (log)", xlim=c(0, max(gene.dist)))
	lines(fit2$x, log10(fit2$y), col="yellow", lwd=0.5)
	lines(fit2.mr$x, log10(fit2.mr$yf), col="blue", lwd=0.5)
	lines(fit2_new$x, log10(fit2_new$y), col="green", lwd=0.5)
	lines(fit2_new.mr$x, log10(fit2_new.mr$yf), col="red", lwd=0.5)
	legend("topright",legend=c(paste("spline (fit2) - ", as.integer(fit2$df), "df (cv)"), "MR on spline (fit2.mr)", paste("5 Kb uniform spline (fit2_new) - ", as.integer(fit2_new$df), "df (cv)"), "MR on 5 Kb uniform spline (fit2_new.mr)"), col=c("yellow", "blue", "green", "red"), lty=1, lwd=1, cex=0.8)
	title("Smooth spline - antitonic regression - pass 1 - log scale")
	dev.off()
}

# create a data frame to write the old spline and old MR regression results
old_sample.df <- cbind.data.frame(fit2$x, fit2$y, fit2.mr$x, fit2.mr$yf)
colnames(old_sample.df) <- c('fit2$x', 'fit2$y', 'fit2.mr$x', 'fit2.mr$yf')
write.table(old_sample.df, paste0(outdir,'/EqOccBin_OldSpline_Regression.txt'), row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

# create a data frame to write the old spline and old MR regression results
new_sample.df <- cbind.data.frame(fit2_new$x, fit2_new$y, fit2_new.mr$x, fit2_new.mr$yf)
colnames(new_sample.df) <- c('fit2_new$x', 'fit2_new$y', 'fit2_new.mr$x', 'fit2_new.mr$yf')
write.table(new_sample.df, paste0(outdir,'/EqOccBin_NewSpline_Regression.txt'), row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

#--------------------------------------------
# Paper: Ferhat Ay, 2014
# the probability is computed with respect to the spline plot
# which is then substituted to the binomial distribution, for P-value estimation

if (is.unsorted(gene.dist) == TRUE) {

	cat(sprintf("\n **** The interaction file is unsorted --- check ****"))

	# here the genomic distance is unsorted
	# so we employ the predict function and probability distribution computation
	# for the complete array

	# get the spline interpolation prediction from the already modeled smoothing spline object (CV spline)
	# with respect to the input genomic distance vector
	# the returned value is a component of 2 fields: (x,y)
	# where x is the input data and y is the predicted data
	
	# input will be the genomic distance
	# output is the corresponding contact probability 

	# note: the probability was computed with respect to 'numPairs' placed as the denominator
	# accordingly we modify the binomial distribution equation

	# predict from the new spline fit - sourya
	# pp <- predict(fit2, gene.dist)
	pp <- predict(fit2_new, gene.dist)

	for (k in (1:numPairs)) {
		# probability for this particular interaction is in the y field
		p <- pp$y[k]
		if (opt$BiasCorr == 0) {
			# model binomial distribution with this probability
			Spline_Binom_Prob_CC[k] <- dbinom(interaction.data[k, opt$cccol], size=TotContact, prob=p)
			# compute the p value with respect to the modified distribution
			Spline_Binom_P_Val_CC[k] <- pbinom(interaction.data[k, opt$cccol], size=TotContact, prob=p, lower.tail=FALSE) + Spline_Binom_Prob_CC[k]
			# store the prior probability
			Prob_Val[k] <- p
		} else {
			# model binomial distribution with the product of probability, and the bias values of two segments
			if (opt$ProbBias == 1) {
				curr_prob <- p * interaction.data[k, bias1.col] * interaction.data[k, bias2.col]	
			} else {
				curr_prob <- p
			}
			Spline_Binom_Prob_CC[k] <- dbinom(interaction.data[k, opt$cccol], size=TotContact, prob=curr_prob)
			# compute the p value with respect to the modified distribution
			Spline_Binom_P_Val_CC[k] <- pbinom(interaction.data[k, opt$cccol], size=TotContact, prob=curr_prob, lower.tail=FALSE) + Spline_Binom_Prob_CC[k]
			# store the prior probability
			Prob_Val[k] <- p
		}
	}
} else {

	if (timeprof == 1) {
		starttime <- Sys.time()
	}

	# here the genomic distance is sorted
	# so we selectively apply the predict function to only the distinct elements of the list
	uniq.idx <- order(gene.dist)[!duplicated(gene.dist)]

	if (timeprof == 1) {
		endtime <- Sys.time()
		fp <- file(opt$TimeFile, open="a")
		outstr <- paste('\n Time to find unique indices in genomic distance: ', (endtime - starttime))
		write(outstr, file=fp, append=T)
		close(fp)
	}

	if (timeprof == 1) {
		starttime <- Sys.time()
	}

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
		# pp <- predict(fit2, gene.dist[si])
		pp <- predict(fit2_new, gene.dist[si])
		p <- pp$y

		#==========================
		# compute the probability distribution for the current interval
		#==========================
		if (opt$BiasCorr == 0) {
			curr_dbinom <- dbinom(interaction.data[si:ei, opt$cccol], size=TotContact, prob=p)
			curr_pbinom <- pbinom(interaction.data[si:ei, opt$cccol], size=TotContact, prob=p, lower.tail=FALSE)
			# append the distribution to the final vector
		 	Spline_Binom_Prob_CC <- c(Spline_Binom_Prob_CC, curr_dbinom)
		 	Spline_Binom_P_Val_CC <- c(Spline_Binom_P_Val_CC, (curr_dbinom + curr_pbinom))	
		 	# store the prior probability as well
		 	Prob_Val[si:ei] <- p
		} else {
			# bias based probability correction

			# comment - sourya
			# for (idx in (si:ei)) {
			# 	# model binomial distribution with the product of probability, and the bias values of two segments
			# 	curr_prob <- p * interaction.data[idx, (opt$cccol + 1)] * interaction.data[idx, (opt$cccol + 3)]
			# 	curr_dbinom <- dbinom(interaction.data[idx, opt$cccol], size=TotContact, prob=curr_prob)
			# 	curr_pbinom <- pbinom(interaction.data[idx, opt$cccol], size=TotContact, prob=curr_prob, lower.tail=FALSE)
			# 	# append the distribution to the final vector
			#  	Spline_Binom_Prob_CC[idx] <- curr_dbinom
			#  	Spline_Binom_P_Val_CC[idx] <- (curr_dbinom + curr_pbinom)	
			#  	# store the prior probability as well
		 		# Prob_Val[idx] <- curr_prob
			# }

			# add - sourya

			# cat(sprintf("\n *** Before executing Parallel from si: %s to ei: %s ", si, ei))

			# now apply this matrix on the parallel version of binomial distribution computation
			# the second argument 1 means that individual rows of matrix is processed 
			# in the parallel computation
			result_binomdistr <- as.data.frame(parallel:::mclapply( si:ei , mc.cores = ncore , function(idx){
				cc <- interaction.data[idx, opt$cccol]
				if (opt$ProbBias == 1) {
					pr <- p * interaction.data[idx, bias1.col] * interaction.data[idx, bias2.col]
				} else {
					pr <- p
				}
				db <- dbinom(cc, size=TotContact, prob=pr)
				pb <- pbinom(cc, size=TotContact, prob=pr, lower.tail=FALSE)
				return(c(idx,pr,db,(db+pb)))
				} ))

			# cat(sprintf("\n *** rows: %s   columns: %s ", nrow(result_binomdistr), ncol(result_binomdistr)))

			# copy the results
			# columns : si to ei
			# rows: 4 
			# for (i in (1:ncol(result_binomdistr))) {
			# 	idx <- result_binomdistr[1,i]
			# 	Prob_Val[idx] <- result_binomdistr[2,i]
			# 	Spline_Binom_Prob_CC[idx] <- result_binomdistr[3,i]
			# 	Spline_Binom_P_Val_CC[idx] <- result_binomdistr[4,i]
			# }
			Prob_Val[si:ei] <- as.double(result_binomdistr[2,1:ncol(result_binomdistr)])
			Spline_Binom_Prob_CC[si:ei] <- as.double(result_binomdistr[3,1:ncol(result_binomdistr)])
			Spline_Binom_P_Val_CC[si:ei] <- as.double(result_binomdistr[4,1:ncol(result_binomdistr)])

			# end add - Sourya
		}
	}

	if (timeprof == 1) {
		endtime <- Sys.time()
		fp <- file(opt$TimeFile, open="a")
		outstr <- paste('\n Time for Spline based distribution compute (p value) (sorted interaction file): ', (endtime - starttime))
		write(outstr, file=fp, append=T)
		close(fp)
	}
}

if (timeprof == 1) {
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
}

if (timeprof == 1) {
	starttime <- Sys.time()
}

# accumulate all results - also add header information
FinalData <- cbind(interaction.data, Prob_Val, Spline_Binom_Prob_CC, Spline_Binom_P_Val_CC, Spline_Binom_QVal)
colnames(FinalData) <- c(colnames(interaction.data), "p", "dbinom", "P-Value", "Q-Value")

#===================================
# add - sourya - debug
# scan through individual genomic distance
# and plot the number of contacts / interactions for each distance
qval_thr <- 0.01
uniq_gene_dist <- sort(unique(gene.dist))

# distance specific contact count and no of interactions
unfilt_cc_dist <- c()
filt_cc_dist <- c()
unfilt_p2p_cc_dist <- c()
filt_p2p_cc_dist <- c()
unfilt_int_dist <- c()
filt_int_dist <- c()
unfilt_p2p_int_dist <- c()
filt_p2p_int_dist <- c()

for (i in 1:length(uniq_gene_dist)) {
	curr_dist <- uniq_gene_dist[i]
	curr_dist_unfilt_int_idx <- which(abs(FinalData[,2] - FinalData[,5]) == curr_dist)
	curr_dist_filt_int_idx <- which((abs(FinalData[,2] - FinalData[,5]) == curr_dist) & (FinalData[,ncol(FinalData)] < qval_thr))
	
	if (ncol(interaction.data) >= 15) {
		curr_dist_unfilt_p2p_idx <- which((abs(FinalData[,2] - FinalData[,5]) == curr_dist) & (FinalData[,9] == 1) & (FinalData[,15] == 1))
		curr_dist_filt_p2p_idx <- which((abs(FinalData[,2] - FinalData[,5]) == curr_dist) & (FinalData[,9] == 1) & (FinalData[,15] == 1) & (FinalData[,ncol(FinalData)] < qval_thr))		
	}

	unfilt_int_dist <- c(unfilt_int_dist, length(curr_dist_unfilt_int_idx))
	unfilt_cc_dist <- c(unfilt_cc_dist, sum(FinalData[curr_dist_unfilt_int_idx, 7]))
	filt_int_dist <- c(filt_int_dist, length(curr_dist_filt_int_idx))
	filt_cc_dist <- c(filt_cc_dist, sum(FinalData[curr_dist_filt_int_idx, 7]))

	if (ncol(interaction.data) >= 15) {
		unfilt_p2p_int_dist <- c(unfilt_p2p_int_dist, length(curr_dist_unfilt_p2p_idx))
		unfilt_p2p_cc_dist <- c(unfilt_p2p_cc_dist, sum(FinalData[curr_dist_unfilt_p2p_idx, 7]))
		filt_p2p_int_dist <- c(filt_p2p_int_dist, length(curr_dist_filt_p2p_idx))
		filt_p2p_cc_dist <- c(filt_p2p_cc_dist, sum(FinalData[curr_dist_filt_p2p_idx, 7]))	
	}
}

x <- sum(unfilt_cc_dist)
unfilt_cc_dist_frac <- unfilt_cc_dist / (x * 1.0)

x <- sum(filt_cc_dist)
filt_cc_dist_frac <- filt_cc_dist / (x * 1.0)

x <- sum(unfilt_int_dist)
unfilt_int_dist_frac <- unfilt_int_dist / (x * 1.0)

x <- sum(filt_int_dist)
filt_int_dist_frac <- filt_int_dist / (x * 1.0)

if (ncol(interaction.data) >= 15) {
	
	x <- sum(unfilt_p2p_cc_dist)
	unfilt_p2p_cc_dist_frac <- unfilt_p2p_cc_dist / (x * 1.0)

	x <- sum(filt_p2p_cc_dist)
	filt_p2p_cc_dist_frac <- filt_p2p_cc_dist / (x * 1.0)

	x <- sum(unfilt_p2p_int_dist)
	unfilt_p2p_int_dist_frac <- unfilt_p2p_int_dist / (x * 1.0)

	x <- sum(filt_p2p_int_dist)
	filt_p2p_int_dist_frac <- filt_p2p_int_dist / (x * 1.0)

	out.df <- cbind.data.frame(uniq_gene_dist, unfilt_int_dist, unfilt_int_dist_frac, unfilt_cc_dist, unfilt_cc_dist_frac, unfilt_p2p_int_dist, unfilt_p2p_int_dist_frac, unfilt_p2p_cc_dist, unfilt_p2p_cc_dist_frac, filt_int_dist, filt_int_dist_frac, filt_cc_dist, filt_cc_dist_frac, filt_p2p_int_dist, filt_p2p_int_dist_frac, filt_p2p_cc_dist, filt_p2p_cc_dist_frac)	

	colnamevec <- c("dist", "U-Num", "%U-Num", "U-CC", "%U-CC", "U-P2P-Num", "%U-P2P-Num", "U-P2P-CC", "%U-P2P-CC", "F-Num", "%F-Num", "F-CC", "%F-CC", "F-P2P-Num", "%F-P2P-Num", "F-P2P-CC", "%F-P2P-CC")	
} else {
	out.df <- cbind.data.frame(uniq_gene_dist, unfilt_int_dist, unfilt_int_dist_frac, unfilt_cc_dist, unfilt_cc_dist_frac, filt_int_dist, filt_int_dist_frac, filt_cc_dist, filt_cc_dist_frac)
	colnamevec <- c("dist", "U-Num", "%U-Num", "U-CC", "%U-CC", "F-Num", "%F-Num", "F-CC", "%F-CC")	
}

colnames(out.df) <- colnamevec
write.table(out.df, paste0(outdir, '/Count_Summary.txt'), row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

colorvec <- c("black", "yellow", "blue", "green", "red", "brown", "magenta", "cyan")

plotfile1 <- paste0(outdir,'/Count_Summary.pdf')	
pdf(plotfile1, width=10, height=8)
plot(out.df[,1], out.df[,3], cex=0.5, col=colorvec[1], type="l", xlab="Distance", ylab="Percentage")
lines(out.df[,1], out.df[,5], col=colorvec[2], lwd=0.5)
lines(out.df[,1], out.df[,7], col=colorvec[3], lwd=0.5)
lines(out.df[,1], out.df[,9], col=colorvec[4], lwd=0.5)

if (ncol(interaction.data) >= 15) {
	lines(out.df[,1], out.df[,11], col=colorvec[5], lwd=0.5)
	lines(out.df[,1], out.df[,13], col=colorvec[6], lwd=0.5)
	lines(out.df[,1], out.df[,15], col=colorvec[7], lwd=0.5)
	lines(out.df[,1], out.df[,17], col=colorvec[8], lwd=0.5)
}

if (ncol(interaction.data) >= 15) {
	legendvec <- c(colnamevec[3], colnamevec[5], colnamevec[7], colnamevec[9], colnamevec[11], colnamevec[13], colnamevec[15], colnamevec[17])
} else {
	legendvec <- c(colnamevec[3], colnamevec[5], colnamevec[7], colnamevec[9])
}

legend("topright",legend=legendvec, col=colorvec, lty=1, lwd=1, cex=0.8)
title("% of different quantities")
dev.off()
# end add - sourya - debug
#===================================


# append the Spline distribution probability and corresponding P value as separate columns
# and write in a separate text file
temp_outfile <- paste0(OutIntDir, '/', 'temp_out.bed')
write.table(FinalData, temp_outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE) 

# now sort the file contents and write that in the final specified output file
system(paste('sort -k1,1 -k2,2n -k5,5n', paste0('-k',opt$cccol,',',opt$cccol,'nr'), temp_outfile, '>', opt$OutFile))
# delete the temporary output file
system(paste('rm', temp_outfile))

# # sort the data according to the chromosome name, start positions of both the interacting segments
# # and finally with the decreasing contact count
# FinalData <- FinalData[order( FinalData[,1], FinalData[,2], FinalData[,5], -FinalData[,opt$cccol] ), ]

# # write the data with the header information
# write.table(FinalData, opt$OutFile, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)

if (timeprof == 1) {
	endtime <- Sys.time()
	fp <- file(opt$TimeFile, open="a")
	outstr <- paste('\n Time to write the interaction data (with Q value) for spline based distribution: ', (endtime - starttime))
	write(outstr, file=fp, append=T)
	close(fp)
}
