#!/usr/bin/env Rscript

#===========================================================
# R script for plotting variation of different normalization related features
# to check their relations and mutual dependency
# The contact count is plotted as a third feature, in a variable color scale
# to depict the variation of various normalization features with respect to the color map

# We also compute mean and SAD features for each plot

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================
library(optparse)
library(ggplot2)
# library(RColorBrewer)
library(data.table)

options(scipen = 999)
options(datatable.fread.datatable=FALSE)


# variables depicting the dimension of plots
PlotWidth <- 14
PlotHeight <- 10
cexthr <- 0.05

#===================================================
# this function is useful for plotting X-Y distribution
# subject to a binning employed
# the number of bins are also provided as an input argument
PlotBinnedDistr <- function(var1, var2, plotfile, xlab, ylab, titlestr, nbin=100) {
	colorvec <- c("blue", "cyan", "green", "yellow", "orange", "brown", "violet", "red")
	curr_plot <- ggplot(data.frame(var1,var2), aes(x=var1,y=var2)) + geom_bin2d(bins = nbin) + scale_fill_gradientn(colours=colorvec) + labs(x = xlab, y = ylab, title = titlestr)
	ggsave(plotfile, plot = curr_plot, width=PlotWidth, height=PlotHeight)
}

# #===================================================
# # function to compute the bin specific interactions
# # where the bins are formed according to the distribution of values in the input data
# # parameters:
# # var1 and var2 are the two variables (data of two columns to be used for X and Y axis)
# # colorvar: variable which is used for coloring the heatmap 
# # nbin: no of bins to be employed for each variable
# # xlab, ylab, colorlab: labels associated with individual variables
# # titlestr: title of the generated plot
# # UseLog: if 1, uses the log10 scale of the color quantity variable

# PlotBinnedCC <- function(var1, var2, colorvar, plotfile, xlab, ylab, colorlab, titlestr, UseLog=0, nbin=20) {
# 	quant_interval <- (1.0 / nbin)
# 	colorvec <- c("blue", "cyan", "green", "yellow", "orange", "brown", "violet", "red")
# 	# pointsize <- 0.01

# 	# quantile of two input variables
# 	# to place them in equal occupancy bins
# 	var1_Quant_values <- unique(quantile(var1, probs=seq(0, 1, quant_interval)))
# 	var2_Quant_values <- unique(quantile(var2, probs=seq(0, 1, quant_interval)))

# 	var1.quant <- as.integer(cut(var1, var1_Quant_values, include.lowest=TRUE))
# 	var2.quant <- as.integer(cut(var2, var2_Quant_values, include.lowest=TRUE))
	
# 	# create a data frame
# 	# Note: simple cbind operation does not produce the data frame object
# 	df <- data.frame(var1Q=var1.quant, var2Q=var2.quant, contact=as.numeric(colorvar))
	
# 	# the color variable (values) for the common bin pairs are averaged (using the FUN=mean argument) 
# 	# a new data frame needs to be assigned
# 	# its column  names also need to be adjusted
# 	df1 <- aggregate(df$contact, by=list(var1Q=df$var1Q,var2Q=df$var2Q), FUN=mean)
# 	colnames(df1) <- colnames(df)
	
# 	# plot the generated aggregated data frame
# 	if (UseLog == 1) {
# 		# curr_plot <- ggplot(df1, aes(x=var1Q, y=var2Q)) + geom_point(aes(colour=log10(contact)),size=pointsize) + scale_colour_gradientn(colours=colorvec) + labs(colour = colorlab, x = xlab, y = ylab, title = titlestr)
# 		curr_plot <- ggplot(df1, aes(x=var1Q, y=var2Q)) + geom_tile(aes(fill=log2(contact)),color="blue") + scale_fill_gradientn(colours=colorvec) + labs(x = xlab, y = ylab, title = titlestr) + scale_x_continuous(expand = c(0, 0), breaks=1:length(var1_Quant_values), labels=paste(var1_Quant_values)) + scale_y_continuous(expand = c(0, 0), breaks=1:length(var2_Quant_values), labels=paste(var2_Quant_values)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 	} else {
# 		# curr_plot <- ggplot(df1, aes(x=var1Q, y=var2Q)) + geom_point(aes(colour=contact),size=pointsize) + scale_colour_gradientn(colours=colorvec) + labs(colour = colorlab, x = xlab, y = ylab, title = titlestr)
# 		curr_plot <- ggplot(df1, aes(x=var1Q, y=var2Q)) + geom_tile(aes(fill=contact),color="blue") + scale_fill_gradientn(colours=colorvec) + labs(x = xlab, y = ylab, title = titlestr) + scale_x_continuous(expand = c(0, 0), breaks=1:length(var1_Quant_values), labels=paste(var1_Quant_values)) + scale_y_continuous(expand = c(0, 0), breaks=1:length(var2_Quant_values), labels=paste(var2_Quant_values)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 	}
# 	ggsave(plotfile, width=PlotWidth, height=PlotHeight)	
# }


#===================================================
# NEW function 
# to compute the bin specific interactions
# where the bins are formed according to the distribution of values in the input data
# parameters:
# var1 and var2 are the two variables (data of two columns to be used for X and Y axis)
# colorvar: variable which is used for coloring the heatmap (generally the contact count)
# nbin: no of bins to be employed for each variable
# xlab, ylab, colorlab: labels associated with individual variables
# titlestr: title of the generated plot
# UseLog: if 1, uses the log10 scale of the color quantity variable

PlotBinnedCC <- function(var1, var2, colorvar, plotfile, xlab, ylab, colorlab, titlestr, UseLog=0, nbin=20) {
	quant_interval <- (1.0 / nbin)
	colorvec <- c("blue", "cyan", "green", "yellow", "orange", "brown", "violet", "red")
	# pointsize <- 0.01

	# the range which will be used for displaying the heatmap
	# every map should have identical range for comparison of the color map
	minrange <- 0
	maxrange <- 10

	# quantile of two input variables
	# to place them in equal occupancy bins
	var1_Quant_values <- unique(quantile(var1, probs=seq(0, 1, quant_interval)))
	var2_Quant_values <- unique(quantile(var2, probs=seq(0, 1, quant_interval)))

	var1.quant <- as.integer(cut(var1, var1_Quant_values, include.lowest=TRUE))
	var2.quant <- as.integer(cut(var2, var2_Quant_values, include.lowest=TRUE))
	
	# now analyze the vector depicting the heatmap 
	TargetVec <- as.numeric(colorvar)

	# if log scale is used, first convert the input data into log2 scale
	if (UseLog == 1) {
		TargetVec <- log2(TargetVec)
	}

	# create a data frame
	# Note: simple cbind operation does not produce the data frame object
	df <- data.frame(var1Q=var1.quant, var2Q=var2.quant, contact=TargetVec)
	
	# the color variable (values) for the common bin pairs are averaged (using the FUN=mean argument) 
	# a new data frame needs to be assigned
	# its column  names also need to be adjusted
	df1 <- aggregate(df$contact, by=list(var1Q=df$var1Q,var2Q=df$var2Q), FUN=mean)
	colnames(df1) <- colnames(df)

	# we shift the range of contact count distribution to set a fixed color scale 
	range_mult <- (maxrange - minrange) / (max(df1$contact) - min(df1$contact))
	df1$contact <- df1$contact * range_mult
	offset <- minrange - min(df1$contact)
	df1$contact <- df1$contact + offset
	
	# compute the mean contact count (bin specific) and also the SAD of contact
	# with respect to all bins
	# after shifting the scale to the common color scale and range
	mean_contact <- mean(df1$contact)
	SAD_mean_contact <- sum(abs(df1$contact - mean_contact))
	titlestr <- paste0(titlestr, ' mean contact: ', mean_contact, ' SAD mean contact: ', SAD_mean_contact)

	# plot the generated aggregated data frame	
	# the TargetVec is used for filling the colors
	curr_plot <- ggplot(df1, aes(x=var1Q, y=var2Q)) + geom_tile(aes(fill=contact),color="blue") + scale_fill_gradientn(colours=colorvec, limits=c(minrange, maxrange)) + labs(x = xlab, y = ylab, title = titlestr, colour = colorlab) + scale_x_continuous(expand = c(0, 0), breaks=1:length(var1_Quant_values), labels=paste(var1_Quant_values)) + scale_y_continuous(expand = c(0, 0), breaks=1:length(var2_Quant_values), labels=paste(var2_Quant_values)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	
	ggsave(plotfile, plot = curr_plot, width=PlotWidth, height=PlotHeight)	
}

# #===================================================
# # generic function to plot 2D histograms
# Plot2DHist <- function(var1, var2, plotfile, xlab, ylab, colorlab, titlestr, nbin=200) {
# 	# quant_interval <- (1.0 / nbin)
# 	colorvec <- c("blue", "cyan", "green", "yellow", "orange", "brown", "violet", "red")
# 	# var1_Quant_values <- unique(quantile(var1, probs=seq(0, 1, quant_interval)))
# 	# var2_Quant_values <- unique(quantile(var2, probs=seq(0, 1, quant_interval)))	
# 	# var1.quant <- as.integer(cut(var1, var1_Quant_values, include.lowest=TRUE))
# 	# var2.quant <- as.integer(cut(var2, var2_Quant_values, include.lowest=TRUE))
# 	# df <- data.frame(var1Q=var1.quant, var2Q=var2.quant)
# 	df <- data.frame(var1Q=var1, var2Q=var2)
# 	curr_plot <- ggplot(df, aes(x=var1Q, y=var2Q)) + stat_bin2d(bins = nbin) + scale_fill_gradientn(colours=colorvec) + labs(x = xlab, y = ylab, title = titlestr, colour = colorlab)
# 	ggsave(plotfile, plot = curr_plot, width=PlotWidth, height=PlotHeight)
# }


#===================================================
option_list = list(
  	make_option(c("--IntFile"), type="character", default=NULL, help="Complete interaction file with normalization features between the genomic bins"),
  	make_option(c("--CommonDir"), type="character", default=NULL, help="Output directory for storing the non-bias specific plots"),
  	make_option(c("--BiasSpecificDir"), type="character", default=NULL, help="Output directory for storing the non-bias specific plots"),
	make_option(c("--MappThr"), type="numeric", default=0, help="The threshold of mappability score below which every segment (peak or non peak) is not considered for analysis. Default 0."),
	make_option(c("--GCThr"), type="numeric", default=0, help="The threshold of GC content below which every segment (peak or non peak) is not considered for analysis. Default 0."),
	make_option(c("--cccol"), type="integer", action="store", default=7, help="Column number of the file storing the contact count"),
	make_option(c("--OverWrite"), type="integer", action="store", default=0, help="if 1, overwrites the plots"),
	make_option(c("--NoMappGC"), type="logical", action="store_true", default=FALSE, help="If TRUE, there is no Mappability, GC content, or RE sites information in the feature file. Default FALSE.")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$IntFile)) {
	cat(sprintf("\n User did not provide any input interaction file - quit !!"))
	system('exit 1')
}

cat(sprintf("\n Input interaction file: %s Output directory: %s  ", opt$IntFile, opt$CommonDir))

if (opt$NoMappGC == FALSE) {
	cat(sprintf("\n Mappability threshold: %s  GC content threshold: %s ", opt$MappThr, opt$GCThr))	
}

# create the output directory
system(paste('mkdir -p', opt$CommonDir))
system(paste('mkdir -p', opt$BiasSpecificDir))

# the interaction matrix
# first 6 columns denote two chromosome intervals 
# then contact count, 
# then coverage1, isPeak1, Bias1, mappability1, GC content1, and number of cut sites(1) for the first interval
# then coverage2, isPeak2, Bias2, mappability2, GC content2, and number of cut sites(2) for the second interval

# Note that the file has one header line
# Interaction_Initial <- read.table(opt$IntFile, header=T, sep="\t", stringsAsFactors=F)
Interaction_Initial <- data.table::fread(opt$IntFile, header=T, sep="\t", stringsAsFactors=F)
colnames(Interaction_Initial) <- c("chr1", "s1", "e1", "chr2", "s2", "e2", "cc", "d1", "isPeak1", "Bias1", "mapp1", "gc1", "cut1", "d2", "isPeak2", "Bias2", "mapp2", "gc2", "cut2")
cat(sprintf('\n Total interactions: %s ', nrow(Interaction_Initial)))

# filter the interactions such that every pair of segments has a mappability score >= mappability threshold
# and GC content >= GC content threshold
# Note the comma after the condition - row index
Interaction_Filt <- Interaction_Initial[Interaction_Initial[,11] >= opt$MappThr & Interaction_Initial[,17] >= opt$MappThr & Interaction_Initial[,12] >= opt$GCThr & Interaction_Initial[,18] >= opt$GCThr, ] 
cat(sprintf('\n Interactions satisfying mappability and GC content thresholds: %s', nrow(Interaction_Filt)))

#====================================
# checking different interaction types
#====================================
IntList <- c("peaktopeak", "peaktononpeak", "peaktoall", "alltoall")

for (inttype in IntList) {

	# output directory storing the plots (for non bias related features)
	PlotDirCommon <- paste0(opt$CommonDir, '/', inttype)	
	system(paste('mkdir -p', PlotDirCommon))

	# output directory storing the plots (for bias related features)
	PlotDirBias <- paste0(opt$BiasSpecificDir, '/', inttype)	
	system(paste('mkdir -p', PlotDirBias))

	# In the following derivations
	# Interaction: filtered rows of interacion data according to the specified criterion
	# ContactCol: Data of absolute contact count
	# CoverageEnrichment: Ratio of contacts and the number of fragment pairs = contact count / (coverage1 * coverage2)
	if (inttype == "peaktopeak") {
		Interaction <- Interaction_Filt[which(Interaction_Filt[,9]==1 & Interaction_Filt[,15]==1), ]
		ContactCol <- Interaction[, opt$cccol]
		CoverageEnrichment <- (Interaction[, opt$cccol] / Interaction[, 8]) / Interaction[, 14]
	} else if (inttype == "peaktononpeak") {
		Interaction <- Interaction_Filt[which(Interaction_Filt[,9]==1 & Interaction_Filt[,15]==0), ]
		ContactCol <- Interaction[, opt$cccol]
		CoverageEnrichment <- (Interaction[, opt$cccol] / Interaction[, 8]) / Interaction[, 14]
	} else if (inttype == "peaktoall") {
		Interaction <- Interaction_Filt[which(Interaction_Filt[,9]==1), ]
		ContactCol <- Interaction[, opt$cccol]
		CoverageEnrichment <- (Interaction[, opt$cccol] / Interaction[, 8]) / Interaction[, 14]
	} else {
		Interaction <- Interaction_Filt
		ContactCol <- Interaction[, opt$cccol]
		CoverageEnrichment <- (Interaction[, opt$cccol] / Interaction[, 8]) / Interaction[, 14]
	}

	# find the indices (rows) where both intervals satisfy bias threshold 
	Bias_Low_Thr <- 0.2
	Bias_High_Thr <- 5
	biasidx <- which(Interaction[,10]>=Bias_Low_Thr & Interaction[,10]<=Bias_High_Thr & Interaction[,16]>=Bias_Low_Thr & Interaction[,16]<=Bias_High_Thr)

	#====================================
	# variation of product of a measure with the contact count
	#====================================
	if (1) {
		plotfile <- paste0(PlotDirCommon, '/ContactCount_vs_IntDist.pdf')
		# absolute genomic Distance for different interactions
		gene.dist <- abs(Interaction[,2] - Interaction[,5])	
		if ((file.exists(plotfile) == FALSE) | (opt$OverWrite == 1)) {
			PlotBinnedDistr(gene.dist, ContactCol, plotfile, "Genomic Distance", "Contact count", "Variation between genomic Distance and contact count")
		}

		plotfile <- paste0(PlotDirCommon, '/ProdCoverage_vs_ContactCount.pdf')
		# product of read depth values
		prod.coverage <- log10(Interaction[,8]) + log10(Interaction[,14])
		if ((file.exists(plotfile) == FALSE) | (opt$OverWrite == 1)) {
			PlotBinnedDistr(prod.coverage, log10(ContactCol), plotfile, "Prod Coverage (LOG)", "Contact count (LOG)", "Variation between product of read coverages and contact count (in log scale)", 50)
		}

		if (opt$NoMappGC == FALSE) {
			plotfile <- paste0(PlotDirCommon, '/ProdMappability_vs_ContactCount.pdf')
			# product of mappability values 
			prod.mappability <- Interaction[,11] * Interaction[,17]
			if ((file.exists(plotfile) == FALSE) | (opt$OverWrite == 1)) {
				PlotBinnedDistr(prod.mappability, ContactCol, plotfile, "Prod Mappability", "Contact count", "Variation between product of Mappabilities and contact count", 20)
			}

			plotfile <- paste0(PlotDirCommon, '/ProdGCContent_vs_ContactCount.pdf')
			# product of GC content values
			prod.GC <- Interaction[,12] * Interaction[,18]
			if ((file.exists(plotfile) == FALSE) | (opt$OverWrite == 1)) {
				PlotBinnedDistr(prod.GC, ContactCol, plotfile, "Prod GCContent", "Contact count", "Variation between product of GCContent and contact count", 20)
			}

			plotfile <- paste0(PlotDirCommon, '/ProdNumCutSites_vs_ContactCount.pdf')
			# product of the number of cut sites
			prod.cut_sites <- Interaction[,13] * Interaction[,19]
			if ((file.exists(plotfile) == FALSE) | (opt$OverWrite == 1)) {
				PlotBinnedDistr(prod.cut_sites, ContactCol, plotfile, "Prod NumCutSites", "Contact count", "Variation between product of NumCutSites and contact count", 20)
			}
		}

		plotfile <- paste0(PlotDirBias, '/ProdBias_vs_ContactCount.pdf')	
		# product of bias
		# first select the interactions whose bias values (both ends) are within a specified range
		prod.bias <- Interaction[biasidx,10] * Interaction[biasidx,16]
		titlestr <- paste0("Variation between product of Bias and contact count (subject to bias interval of ", Bias_Low_Thr, " - ", Bias_High_Thr, ")")
		if ((file.exists(plotfile) == FALSE) | (opt$OverWrite == 1)) {
			PlotBinnedDistr(prod.bias, Interaction[biasidx, opt$cccol], plotfile, "Prod Bias", "Contact count", titlestr, 10)
		}

		plotfile <- paste0(PlotDirBias, '/ProdBias_vs_ContactCount_Both_LOG.pdf')	
		# product of bias
		# first select the interactions whose bias values (both ends) are within a specified range
		log_prod.bias <- log10(Interaction[biasidx,10]) + log10(Interaction[biasidx,16])
		cc <- log10(Interaction[biasidx, opt$cccol])
		titlestr <- paste0("Variation between product of Bias (LOG) and contact count (LOG) (subject to bias interval of ", Bias_Low_Thr, " - ", Bias_High_Thr, ")")
		if ((file.exists(plotfile) == FALSE) | (opt$OverWrite == 1)) {
			PlotBinnedDistr(log_prod.bias, cc, plotfile, "Log(Prod Bias)", "Log(Contact count)", titlestr, 10)
		}
	}

	#====================================
	# then modeling variation of individual measures with the contact count
	#====================================

	if (0) {
		plotfile <- paste0(PlotDirCommon, '/ReadDepth_vs_AbsContactCount.pdf')
		PlotBinnedCC(Interaction[,8], Interaction[,14], ContactCol, plotfile, "Read depth 1", "Read depth 2", "Contact count", "Read depth vs contact count")
		plotfile <- paste0(PlotDirCommon, '/ReadDepth_vs_LOGContactCount.pdf')
		PlotBinnedCC(Interaction[,8], Interaction[,14], ContactCol, plotfile, "Read depth 1", "Read depth 2", "Log Contact count", "Read depth vs contact count (LOG)", 1)
		plotfile <- paste0(PlotDirCommon, '/ReadDepth_vs_CoverageEnrichment.pdf')
		PlotBinnedCC(Interaction[,8], Interaction[,14], CoverageEnrichment, plotfile, "Read depth 1", "Read depth 2", "Coverage Enrichment", "Read depth vs Coverage Enrichment", 1)


		if (opt$NoMappGC == FALSE) {
			plotfile <- paste0(PlotDirCommon, '/Mappability_vs_AbsContactCount.pdf')
			PlotBinnedCC(Interaction[,11], Interaction[,17], ContactCol, plotfile, "Mappability 1", "Mappability 2", "Contact count", "Mappability vs contact count")
			plotfile <- paste0(PlotDirCommon, '/Mappability_vs_LOGContactCount.pdf')
			PlotBinnedCC(Interaction[,11], Interaction[,17], ContactCol, plotfile, "Mappability 1", "Mappability 2", "Log Contact count", "Mappability vs contact count (LOG)", 1)
			plotfile <- paste0(PlotDirCommon, '/Mappability_vs_CoverageEnrichment.pdf')
			PlotBinnedCC(Interaction[,11], Interaction[,17], CoverageEnrichment, plotfile, "Mappability 1", "Mappability 2", "Coverage Enrichment", "Mappability vs Coverage Enrichment", 1)

			plotfile <- paste0(PlotDirCommon, '/GCContent_vs_AbsContactCount.pdf')
			PlotBinnedCC(Interaction[,12], Interaction[,18], ContactCol, plotfile, "GC content 1", "GC content 2", "Contact count", "GC content vs contact count")
			plotfile <- paste0(PlotDirCommon, '/GCContent_vs_LOGContactCount.pdf')
			PlotBinnedCC(Interaction[,12], Interaction[,18], ContactCol, plotfile, "GC content 1", "GC content 2", "Log Contact count", "GC content vs contact count (LOG)", 1)
			plotfile <- paste0(PlotDirCommon, '/GCContent_vs_CoverageEnrichment.pdf')
			PlotBinnedCC(Interaction[,12], Interaction[,18], CoverageEnrichment, plotfile, "GC content 1", "GC content 2", "Coverage Enrichment", "GC content vs Coverage Enrichment", 1)

			plotfile <- paste0(PlotDirCommon, '/NumCutSites_vs_AbsContactCount.pdf')
			PlotBinnedCC(Interaction[,13], Interaction[,19], ContactCol, plotfile, "Num cut site 1", "Num cut site 2", "Contact count", "Cut site vs contact count")
			plotfile <- paste0(PlotDirCommon, '/NumCutSites_vs_LOGContactCount.pdf')
			PlotBinnedCC(Interaction[,13], Interaction[,19], ContactCol, plotfile, "Num cut site 1", "Num cut site 2", "Log Contact count", "Cut site vs contact count (LOG)", 1)
			plotfile <- paste0(PlotDirCommon, '/NumCutSites_vs_CoverageEnrichment.pdf')
			PlotBinnedCC(Interaction[,13], Interaction[,19], CoverageEnrichment, plotfile, "Num cut site 1", "Num cut site 2", "Coverage Enrichment", "Cut site vs Coverage Enrichment", 1)
		}

		# incorporating the bias specific plots
		plotfile <- paste0(PlotDirBias, '/Bias_vs_AbsContactCount.pdf')
		titlestr <- paste0("Bias vs contact count - bias interval of ", Bias_Low_Thr, " - ", Bias_High_Thr)
		PlotBinnedCC(Interaction[biasidx,10], Interaction[biasidx,16], Interaction[biasidx, opt$cccol], plotfile, "Bias 1", "Bias 2", "Contact count", titlestr)
		plotfile <- paste0(PlotDirBias, '/Bias_vs_LOGContactCount.pdf')
		titlestr <- paste0("Bias vs contact count (LOG) - bias interval of ", Bias_Low_Thr, " - ", Bias_High_Thr)
		PlotBinnedCC(Interaction[biasidx,10], Interaction[biasidx,16], Interaction[biasidx, opt$cccol], plotfile, "Bias 1", "Bias 2", "Log Contact count", titlestr, 1)
		plotfile <- paste0(PlotDirBias, '/Bias_vs_CoverageEnrichment.pdf')
		titlestr <- paste0("Bias vs Coverage Enrichment - bias interval of ", Bias_Low_Thr, " - ", Bias_High_Thr)
		PlotBinnedCC(Interaction[biasidx,10], Interaction[biasidx,16], ((Interaction[biasidx, opt$cccol] / Interaction[biasidx, 8]) / Interaction[biasidx, 14]), plotfile, "Bias 1", "Bias 2", "Coverage Enrichment", titlestr, 1)		
	}

	if (inttype == "peaktoall") {
		if (0) {

			if (opt$NoMappGC == FALSE) {
				# peak part
				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_Mapp1_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,11], ContactCol, plotfile, "Read depth 1", "Mappability 1", "Contact count", "Read depth + Mappability (1) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_Mapp1_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,11], ContactCol, plotfile, "Read depth 1", "Mappability 1", "Log Contact count", "Read depth + Mappability (1) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_Mapp1_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,11], CoverageEnrichment, plotfile, "Read depth 1", "Mappability 1", "Coverage Enrichment", "Read depth + Mappability (1) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_GC1_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,12], ContactCol, plotfile, "Read depth 1", "GC content 1", "Contact count", "Read depth + GC content (1) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_GC1_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,12], ContactCol, plotfile, "Read depth 1", "GC content 1", "Log Contact count", "Read depth + GC content (1) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_GC1_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,12], CoverageEnrichment, plotfile, "Read depth 1", "GC content 1", "Coverage Enrichment", "Read depth + GC content (1) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_CutSite1_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,13], ContactCol, plotfile, "Read depth 1", "Cut Site 1", "Contact count", "Read depth + Cut Site (1) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_CutSite1_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,13], ContactCol, plotfile, "Read depth 1", "Cut Site 1", "Log Contact count", "Read depth + Cut Site (1) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/ReadDepth1_CutSite1_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,8], Interaction[,13], CoverageEnrichment, plotfile, "Read depth 1", "Cut Site 1", "Coverage Enrichment", "Read depth + Cut Site (1) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/Mapp1_GC1_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,11], Interaction[,12], ContactCol, plotfile, "Mappability 1", "GC content 1", "Contact count", "Mappability + GC content (1) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/Mapp1_GC1_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,11], Interaction[,12], ContactCol, plotfile, "Mappability 1", "GC content 1", "Log Contact count", "Mappability + GC content (1) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/Mapp1_GC1_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,11], Interaction[,12], CoverageEnrichment, plotfile, "Mappability 1", "GC content 1", "Coverage Enrichment", "Mappability + GC content (1) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/Mapp1_CutSite1_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,11], Interaction[,13], ContactCol, plotfile, "Mappability 1", "Cut site 1", "Contact count", "Mappability + Cut site (1) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/Mapp1_CutSite1_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,11], Interaction[,13], ContactCol, plotfile, "Mappability 1", "Cut site 1", "Log Contact count", "Mappability + Cut site (1) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/Mapp1_CutSite1_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,11], Interaction[,13], CoverageEnrichment, plotfile, "Mappability 1", "Cut site 1", "Coverage Enrichment", "Mappability + Cut site (1) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/GC1_CutSite1_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,12], Interaction[,13], ContactCol, plotfile, "GC content 1", "Cut site 1", "Contact count", "GC content + Cut site (1) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/GC1_CutSite1_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,12], Interaction[,13], ContactCol, plotfile, "GC content 1", "Cut site 1", "Log Contact count", "GC content + Cut site (1) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/GC1_CutSite1_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,12], Interaction[,13], CoverageEnrichment, plotfile, "GC content 1", "Cut site 1", "Coverage Enrichment", "GC content + Cut site (1) vs Coverage Enrichment", 1)

				# all part
				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_Mapp2_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,17], ContactCol, plotfile, "Read depth 2", "Mappability 2", "Contact count", "Read depth + Mappability (2) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_Mapp2_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,17], ContactCol, plotfile, "Read depth 2", "Mappability 2", "Log Contact count", "Read depth + Mappability (2) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_Mapp2_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,17], CoverageEnrichment, plotfile, "Read depth 2", "Mappability 2", "Coverage Enrichment", "Read depth + Mappability (2) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_GC2_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,18], ContactCol, plotfile, "Read depth 2", "GC content 2", "Contact count", "Read depth + GC content (2) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_GC2_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,18], ContactCol, plotfile, "Read depth 2", "GC content 2", "Log Contact count", "Read depth + GC content (2) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_GC2_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,18], CoverageEnrichment, plotfile, "Read depth 2", "GC content 2", "Coverage Enrichment", "Read depth + GC content (2) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_CutSite2_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,19], ContactCol, plotfile, "Read depth 2", "Cut Site 2", "Contact count", "Read depth + Cut Site (2) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_CutSite2_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,19], ContactCol, plotfile, "Read depth 2", "Cut Site 2", "Log Contact count", "Read depth + Cut Site (2) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/ReadDepth2_CutSite2_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,14], Interaction[,19], CoverageEnrichment, plotfile, "Read depth 2", "Cut Site 2", "Coverage Enrichment", "Read depth + Cut Site (2) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/Mapp2_GC2_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,17], Interaction[,18], ContactCol, plotfile, "Mappability 2", "GC content 2", "Contact count", "Mappability + GC content (2) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/Mapp2_GC2_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,17], Interaction[,18], ContactCol, plotfile, "Mappability 2", "GC content 2", "Log Contact count", "Mappability + GC content (2) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/Mapp2_GC2_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,17], Interaction[,18], CoverageEnrichment, plotfile, "Mappability 2", "GC content 2", "Coverage Enrichment", "Mappability + GC content (2) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/Mapp2_CutSite2_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,17], Interaction[,19], ContactCol, plotfile, "Mappability 2", "Cut site 2", "Contact count", "Mappability + Cut site (2) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/Mapp2_CutSite2_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,17], Interaction[,19], ContactCol, plotfile, "Mappability 2", "Cut site 2", "Log Contact count", "Mappability + Cut site (2) vs contact count (LOG)", 1)
				plotfile <- paste0(PlotDirCommon, '/Mapp2_CutSite2_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,17], Interaction[,19], CoverageEnrichment, plotfile, "Mappability 2", "Cut site 2", "Coverage Enrichment", "Mappability + Cut site (2) vs Coverage Enrichment", 1)

				plotfile <- paste0(PlotDirCommon, '/GC2_CutSite2_vs_AbsContactCount.pdf')
				PlotBinnedCC(Interaction[,18], Interaction[,19], ContactCol, plotfile, "GC content 2", "Cut site 2", "Contact count", "GC content + Cut site (2) vs contact count")
				plotfile <- paste0(PlotDirCommon, '/GC2_CutSite2_vs_LOGContactCount.pdf')
				PlotBinnedCC(Interaction[,18], Interaction[,19], ContactCol, plotfile, "GC content 2", "Cut site 2", "Log Contact count", "GC content + Cut site (2) vs contact count (LOG)", 1)		
				plotfile <- paste0(PlotDirCommon, '/GC2_CutSite2_vs_CoverageEnrichment.pdf')
				PlotBinnedCC(Interaction[,18], Interaction[,19], CoverageEnrichment, plotfile, "GC content 2", "Cut site 2", "Coverage Enrichment", "GC content + Cut site (2) vs Coverage Enrichment", 1)		
			}
		}
	}	# end peak to all interaction type check

}
