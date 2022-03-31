#!/usr/bin/env Rscript

#===========================================================
# R script for plotting the distance vs contact count
# for significant interactions

# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
#===========================================================

library(optparse)
library(ggplot2)
library(data.table)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

#================================================
option_list = list(
  	make_option(c("--IntFile"), type="character", default=NULL, help="FitHiChIP significant interaction file."),
  	make_option(c("--OutFile"), type="character", default=NULL, help="Output file storing the distance vs contact count plot."),
  	make_option(c("--TitleStr"), type="character", default=NULL, help="Title string to be employed in the plot. Default NULL means that no title will be employed in the plot.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

OutDir <- dirname(opt$OutFile)
system(paste("mkdir -p", OutDir))

cat(sprintf("\n ---- Within function of plotting distance vs contact count ---- \n Input interaction file: %s ", opt$IntFile))
cat(sprintf("\n Output plot file: %s ", opt$OutFile))


# provide mean and median distance statistics 
# for these interactions
InpLoopData <- data.table::fread(opt$IntFile, header=T, sep="\t", stringsAsFactors=F)
CN <- colnames(InpLoopData)
if ("Dist" %in% CN) {
	## new loop file
	DistVec <- as.numeric(InpLoopData$Dist)
} else {
	DistVec <- abs(as.numeric(InpLoopData[,5] - InpLoopData[,2]))
}

median_dist <- median(DistVec)
mean_dist <- mean(DistVec)

textfile <-  paste0(OutDir, '/Statistics_Dist.txt')
fp_out <- file(textfile, "w")
outtext <- paste0("\n Total number of loops: ", nrow(InpLoopData))
writeLines(outtext, con=fp_out, sep="\n")
outtext <- paste0("\n Median distance of these loops: ", median_dist)
writeLines(outtext, con=fp_out, sep="\n")
outtext <- paste0("\n Mean distance of these loops: ", mean_dist)
writeLines(outtext, con=fp_out, sep="\n")
close(fp_out)

# plot histogram of interaction distance vs significant interaction count
# bin width is set as 10 Kb
min_dist <- min(20000, min(DistVec))
max_dist <- max(2000000, max(DistVec))
histdata <- data.frame(dist=DistVec)

max_hist_freq <- max(hist(histdata$dist, breaks=seq(min_dist, max_dist, by=10000), plot=FALSE)$counts)

curr_plotA <- ggplot(data=histdata, aes(x=dist)) + geom_histogram(col="skyblue", fill="green", alpha = .2, binwidth=10000) + labs(title="Loop Count vs Distance", x="Distance", y="Loop count") + xlim(c(min_dist,max_dist)) + geom_vline(aes(xintercept=median(dist)), color="red") + geom_text(aes(x=median(dist), y=(max_hist_freq / 2)),label=paste0("Median Dist: ", median_dist), hjust=-0.25, size=6)
ggsave(opt$OutFile, plot = curr_plotA, width=4, height=3)






