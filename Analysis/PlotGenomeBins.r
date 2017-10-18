#!/usr/bin/env Rscript

#===========================================================
# R script for plotting different features of individual genomic bins - 1D plots

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================
library(optparse)

library(ggplot2)
# library(RColorBrewer)

# variables depicting the dimension of plots
PlotWidth <- 14
PlotHeight <- 10
cexthr <- 0.05

#===================================================
# this function is useful for plotting X-Y distribution
# subject to a binning employed
# the number of bins are also provided as an input argument
PlotBinnedDistr <- function(var1, var2, plotfile, xlab, ylab, titlestr, nbin=200) {
	colorvec <- c("blue", "cyan", "green", "yellow", "orange", "brown", "violet", "red")
	curr_plot <- ggplot(data.frame(var1,var2), aes(x=var1,y=var2)) + geom_bin2d(bins = nbin) + scale_fill_gradientn(colours=colorvec) + labs(x = xlab, y = ylab, title = titlestr)
	ggsave(plotfile, plot = curr_plot, width=PlotWidth, height=PlotHeight)	
	# quant_interval <- (1.0 / nbin)
	# var1_Quant_values <- unique(quantile(var1, probs=seq(0, 1, quant_interval)))
	# var1.quant <- as.integer(cut(var1, var1_Quant_values, include.lowest=TRUE))
	# df <- data.frame(var1Q=var1.quant, var2Q=var2)
	# df1 <- aggregate(df$var2Q, by=list(var1Q=df$var1Q), FUN=mean)
	# colnames(df1) <- colnames(df)
	# pdf(plotfile, width=PlotWidth, height=PlotHeight)
	# plot(df1$var1Q, df1$var2Q, cex=cexthr, col="darkblue", xlab=xlab, ylab=ylab)
	# title(titlestr)
	# dev.off()
}

#===========================================================
# function for plotting distribution
PlotDistr <- function(inpvec, plotfile, titlestr, xlabelstr, ylabelstr, stepsize=0.1) {
	pdf(plotfile, width=14, height=10)
	# plot a histogram where X axis denotes the input vector and Y axis denotes the corresponding frequency
	# step size of the histogram is default 0.1 (unless specified as an input parameter)
	maxelem <- max(inpvec)
	no_of_breaks <- as.integer(maxelem / stepsize)
	hist(inpvec, main=titlestr, xlab=xlabelstr, ylab=ylabelstr, border="red", col="blue", breaks=no_of_breaks)
	dev.off()
}

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

#===================================================
option_list = list(
  	make_option(c("--GenomeBinFile"), type="character", default=NULL, help="File containing the genome bins and associated features", metavar="GenomeBinFile"),
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory", metavar="OutDir")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read the genome bins and associated normalization features
# Note: the file does not have header information
InpData <- read.table(opt$GenomeBinFile, header=F)	
# colnames(InpData) <- c("chr1","s1","e1","coverage","isPeak", "Bias", "Mapp", "GCContent", "numREsites")

system(paste('mkdir -p', opt$OutDir))

# find the indices corresponding to peak and non peaks
PeakIdx <- which(InpData[,5] == 1)
NonPeakIdx <- which(InpData[,5] == 0)

# note the indices having peak information and non zero coverage
PeakIDXNonZeroCoverage <- intersect(which(InpData[,4]> 0), which(InpData[,5]==1))

# note the indices of non peak and non zero coverage
nonPeakIDXNonZeroCoverage <- intersect(which(InpData[,4]> 0), which(InpData[,5]==0))

# Bias of the peak segments and non-zero coverage
BiasPeakVec <- InpData[PeakIDXNonZeroCoverage, 6]

# Bias of the non-peak segments and non-zero coverage
BiasNonPeakVec <- InpData[nonPeakIDXNonZeroCoverage, 6]

#=========================
# Plotting the distribution / histogram of bias with respect to the peak segments
#=========================
OutPlotFile <- paste0(opt$OutDir, '/PeakCoverageDistr.pdf') 
PlotDistr(BiasPeakVec, OutPlotFile, "Bias distribution for peaks", "Bias (ratio w.r.t mean)", "Peak frequency")

#=========================
# Plotting the distribution / histogram of bias with respect to the non-peak segments
#=========================
OutPlotFile <- paste0(opt$OutDir, '/NonPeakCoverageDistr.pdf') 
PlotDistr(BiasNonPeakVec, OutPlotFile, "Bias distribution for non-peaks", "Bias (ratio w.r.t mean)", "Non-Peak frequency")

#=========================
# Density Plot - GC content vs Coverage - for peaks, non peaks and for all segments
#=========================
plotfile <- paste0(opt$OutDir, '/GCContent_vs_Coverage_ALL.pdf')
# pdf(plotfile, width=PlotWidth, height=PlotHeight)
# plot(InpData[,8], InpData[,4], cex=cexthr, col="darkblue", xlab="GC content", ylab="Coverage")
# title("Variation between GC content and coverage for all segments")
# dev.off()	
PlotBinnedDistr(InpData[,8], InpData[,4], plotfile, "GC content", "Coverage", "Variation between GC content and coverage for all segments")

plotfile <- paste0(opt$OutDir, '/GCContent_vs_Coverage_Peaks.pdf')
# pdf(plotfile, width=PlotWidth, height=PlotHeight)
# plot(InpData[PeakIdx,8], InpData[PeakIdx,4], cex=cexthr, col="darkblue", xlab="GC content", ylab="Coverage")
# title("Variation between GC content and coverage for Peaks")
# dev.off()
PlotBinnedDistr(InpData[PeakIdx,8], InpData[PeakIdx,4], plotfile, "GC content", "Coverage", "Variation between GC content and coverage for Peaks")

plotfile <- paste0(opt$OutDir, '/GCContent_vs_Coverage_NonPeaks.pdf')
# pdf(plotfile, width=PlotWidth, height=PlotHeight)
# plot(InpData[NonPeakIdx,8], InpData[NonPeakIdx,4], cex=cexthr, col="darkblue", xlab="GC content", ylab="Coverage")
# title("Variation between GC content and coverage for Non peaks")
# dev.off()
PlotBinnedDistr(InpData[NonPeakIdx,8], InpData[NonPeakIdx,4], plotfile, "GC content", "Coverage", "Variation between GC content and coverage for Non peaks")

#=========================
# Density Plot - Mappability vs Coverage - for peaks, non peaks and for all segments
#=========================
plotfile <- paste0(opt$OutDir, '/Mappability_vs_Coverage_ALL.pdf')
# pdf(plotfile, width=PlotWidth, height=PlotHeight)
# plot(InpData[,7], InpData[,4], cex=cexthr, col="darkblue", xlab="Mappability", ylab="Coverage")
# title("Variation between Mappability and coverage for all segments")
# dev.off()	
PlotBinnedDistr(InpData[,7], InpData[,4], plotfile, "Mappability", "Coverage", "Variation between Mappability and coverage for all segments")

plotfile <- paste0(opt$OutDir, '/Mappability_vs_Coverage_Peaks.pdf')
# pdf(plotfile, width=PlotWidth, height=PlotHeight)
# plot(InpData[PeakIdx,7], InpData[PeakIdx,4], cex=cexthr, col="darkblue", xlab="Mappability", ylab="Coverage")
# title("Variation between Mappability and coverage for Peaks")
# dev.off()
PlotBinnedDistr(InpData[PeakIdx,7], InpData[PeakIdx,4], plotfile, "Mappability", "Coverage", "Variation between Mappability and coverage for Peaks")

plotfile <- paste0(opt$OutDir, '/Mappability_vs_Coverage_NonPeaks.pdf')
# pdf(plotfile, width=PlotWidth, height=PlotHeight)
# plot(InpData[NonPeakIdx,7], InpData[NonPeakIdx,4], cex=cexthr, col="darkblue", xlab="Mappability", ylab="Coverage")
# title("Variation between Mappability and coverage for Non peaks")
# dev.off()
PlotBinnedDistr(InpData[NonPeakIdx,7], InpData[NonPeakIdx,4], plotfile, "Mappability", "Coverage", "Variation between Mappability and coverage for Non peaks")

