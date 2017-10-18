#!/usr/bin/env Rscript

#===========================================================
# R script for summarizing the results of generalized FitHiC
# Here, two different plots are generated:
# 1) distribution of contact count vs interaction distance for two different sets of interactions (separated by Q value)
# 2) distribution of contact count (as boxplot) for two different sets of interactions (separated by Q value)

# Q threshold is set as 0.01

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript result_summary.r $inpfile
#===========================================================

library(ggplot2)

args <- commandArgs(TRUE)

# the 1st file is the file without any q-value based filtering
unfilt.file <- args[1] 
inpdir <- dirname(unfilt.file)

# column no having the absolute contact count
contactcol <- as.integer(args[2])

# threshold of Q value which is used for filtering the contacts
QTHR <- as.double(args[3])

outdir <- paste0(inpdir,'/Results')
system(paste('mkdir -p', outdir))

boxplotCCQvalfile <- paste0(outdir,'/','CC_Qval.png')
plotfile <- paste0(outdir,'/','CC_IntDist.png')

if ((file.exists(boxplotCCQvalfile) == FALSE) | (file.exists(plotfile) == FALSE)) {

	# load the interaction matrix of the unfiltered file
	# Note: this unfiltered interaction file has header information
	unfilt.data <- read.table(unfilt.file, header=T)

	# absolute genomic distance data
	gen.dist <- abs(unfilt.data[,5] - unfilt.data[,2])

	# indices corresponding to significant and insignificant interactions
	idx.qthr.pass <- which(unfilt.data[,ncol(unfilt.data)] < QTHR)
	idx.qthr.fail <- which(unfilt.data[,ncol(unfilt.data)] >= QTHR)

	# we only process and analyze if we have at least 2 success and 2 failure candidates
	if ((length(idx.qthr.pass) > 1) && (length(idx.qthr.fail) > 1)) {

		# get the genomic distance and contact count for the significant interactions (Q value < QTHR)
		contactcountcol <- unfilt.data[,contactcol]
		CC.qthr.pass <- contactcountcol[idx.qthr.pass]
		gen.dist.qthr.pass <- gen.dist[idx.qthr.pass]

		# get the genomic distance and contact count for the nonsignificant interactions (Q value >= QTHR)
		idx.qthr.fail <- which(unfilt.data[,ncol(unfilt.data)] >= QTHR)
		contactcountcol <- unfilt.data[,contactcol]
		CC.qthr.fail <- contactcountcol[idx.qthr.fail]
		gen.dist.qthr.fail <- gen.dist[idx.qthr.fail]

		#===========================================
		# plot the distributions of contact count for two different sets of Q values
		# The first set is the set of significant interactions (Q val < 0.01)
		# The second set is the set of significant interactions (Q val >= 0.01)
		if (file.exists(boxplotCCQvalfile) == FALSE) {
			a <- data.frame(group = paste0("CC_Qval<", QTHR), value = CC.qthr.pass)
			b <- data.frame(group = paste0("CC_Qval>=", QTHR), value = CC.qthr.fail)
			ccQVal <- rbind(a, b)
			#pdf(boxplotCCQvalfile, width=10, height=6)
			currplot <- ggplot(ccQVal, aes(x=group, y=value, fill=group)) + geom_boxplot()
			currplot + ggtitle("Box plot for contact counts - significant vs non-significant interactions")
			ggsave(boxplotCCQvalfile, width=10, height=6)
			#dev.off()		
		}
		#===========================================

		# first bind the genomic distance and the contact count information
		# for both Q value sets
		qthr.pass_GenDist_CC <- cbind(gen.dist.qthr.pass, CC.qthr.pass)
		qthr.fail_GenDist_CC <- cbind(gen.dist.qthr.fail, CC.qthr.fail)

		# now sort the data frames according to the genomic distance
		# (first column ascending order)
		# also remove the duplicate entries
		qthr.pass_GenDist_CC_uniq <- unique(qthr.pass_GenDist_CC[ order(qthr.pass_GenDist_CC[,1]), ])
		qthr.fail_GenDist_CC_uniq <- unique(qthr.fail_GenDist_CC[ order(qthr.fail_GenDist_CC[,1]), ])

		# now plot all the elements of the significant (pass) interactions
		# and selectively plot the elements of the non-significant (fail) interactions
		nelem <- length(qthr.pass_GenDist_CC_uniq[,1])

		qval.fail.meanCC <- c()
		AvgIntDist.fail <- c()
		num.qval.fail <- length(qthr.fail_GenDist_CC_uniq[,1])
		nentry <- floor(num.qval.fail / nelem)
		for (i in (1:nelem)) {
			si <- (i - 1) * nentry + 1
			if (i < nelem) {
				ei <- si + nentry - 1		
			} else {
				ei <- num.qval.fail
			}
			qval.fail.meanCC[i] <- mean(qthr.fail_GenDist_CC_uniq[si:ei,2])
			AvgIntDist.fail[i] <- mean(qthr.fail_GenDist_CC_uniq[si:ei,1])
		}

		a <- data.frame(group = "Qval<0.01", x = qthr.pass_GenDist_CC_uniq[,1], y = qthr.pass_GenDist_CC_uniq[,2])
		b <- data.frame(group = "Qval>=0.01", x = AvgIntDist.fail, y = qval.fail.meanCC)
		curr_plot <- ggplot(rbind(a,b), aes(x=x, y=y, fill=group, colour=group)) + geom_point(size=0.01) + xlab('Average interaction distance') + ylab('Contact count')
		curr_plot + ggtitle("Interaction distance vs Contact count - both Qvalue classes")
		ggsave(plotfile, width=12, height=8)

	#=========================================== 		
	} else {
		cat(sprintf("\n Either significant or insignificant interactions have count 0 - quit !!"))
	}
}

