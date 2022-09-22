#!/usr/bin/env Rscript

##======================
## filters input interactions according to the distance range
## allows circukler genome
##======================

args <- commandArgs(TRUE)
InpFile <- as.character(args[1])
OutFile <- as.character(args[2])
LowDistThres <- as.integer(args[3])
UppDistThres <- as.integer(args[4])
CircularGenome <- as.integer(args[5])
ChrSizeFile <- as.character(args[6])

if (CircularGenome == 0) {
	## default interaction filtering based on the absolute genomic distance
	system(paste0("awk -v l=", LowDistThres, " -v u=", UppDistThres, " -F[\'\t\'] \'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || (($1==$4) && ($7>0) && (abs($2-$5)>=l) && (abs($2-$5)<=u))) {print $0}}\' ", InpFile, " > ", OutFile))
} else {
	## use circular genome and chromosome size information
	system(paste0("awk \'(NR==1)\' ", InpFile, " > ", OutFile))

	ChrSizeData <- read.table(ChrSizeFile, header=F, sep="\t", stringsAsFactors=F)
	ChrList <- unique(ChrSizeData[,1])
	for (j in 1:length(ChrList)) {

		currchr <- ChrList[j]
		cat(sprintf("\n ==>> circular genome - genomic distance based filtering - processing chromosome : %s ", currchr))
		## chromosome size for this chromosome
		max_coord <- ChrSizeData[which(ChrSizeData[,1] == currchr), 2]
		cat(sprintf(" -->> max_coord (from ChrSizeFile) for this chromosome : %s ", max_coord))

		## analyze the circular genomic distance
		## include cis-pairs, with distance within the specified distance thresholds		
		system(paste0("awk -F\'[\t]\' -v R=", max_coord, " -v c=\"", currchr, "\" -v l=", LowDistThres, " -v u=", UppDistThres, " \'function abs(v) {return v < 0 ? -v : v}; function min(x,y) {return x < y ? x : y} {if ((NR>1) && ($1==c) && ($1==$4) && ($7>0)) {d=min(min(abs($5-$2),abs((R-$5)+$2)),abs((R-$2)+$5)); if ((d>=l) && (d<=u)) {print $0}}}\' ", InpFile, " >> ", OutFile))
	
	}	# end chromosome loop
}
