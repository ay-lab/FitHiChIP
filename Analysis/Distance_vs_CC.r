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

options(scipen = 999)
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

if (0) {

tempOutFile <- paste0(OutDir, '/temp_CC_Dist.bed')

# first print the distance vs interaction count in a text file
# where the first column is distance, and the second column is interaction count
system(paste("awk \'function abs(v) {return v < 0 ? -v : v} {print abs($5-$2)}\'", opt$IntFile, "| sort -k1,1n | uniq -c | awk \'{if (NR>1) {print $2\"\t\"$1}}\' - > ", tempOutFile))

# InpData <- read.table(tempOutFile, header=F, sep="\t", stringsAsFactors=F)
InpData <- data.table::fread(tempOutFile, header=F, sep="\t", stringsAsFactors=F)

a <- data.frame(group = paste("Distance vs Interaction"), x = InpData[,1], y = InpData[,2])
curr_plotA <- ggplot(a, aes(x=x, y=y, fill=group, colour=group)) + geom_line(color="blue") + xlab('Interaction Distance') + ylab('No of significant interactions')
curr_plotA + ggtitle("FitHiChIP - distance vs significant interactions")
ggsave(opt$OutFile, plot = curr_plotA, width=8, height=6)

# now remove the temporary file
system(paste("rm", tempOutFile))

}	# end dummy if


# provide mean and median distance statistics 
# for these interactions
# InpLoopData <- read.table(opt$IntFile, header=T)
InpLoopData <- data.table::fread(opt$IntFile, header=T)
DistVec <- abs(InpLoopData[,5] - InpLoopData[,2])
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

curr_plotA <- ggplot(data=histdata, aes(histdata$dist)) + geom_histogram(col="skyblue", fill="green", alpha = .2, binwidth=10000) + labs(title="Loop Count vs Distance", x="Distance", y="Loop count") + xlim(c(min_dist,max_dist)) + geom_vline(aes(xintercept=median(histdata$dist)), color="red") + geom_text(aes(x=median(histdata$dist), y=(max_hist_freq / 2)),label=paste0("Median Dist: ", median_dist), hjust=-0.25, size=6)
ggsave(opt$OutFile, plot = curr_plotA, width=4, height=3)






