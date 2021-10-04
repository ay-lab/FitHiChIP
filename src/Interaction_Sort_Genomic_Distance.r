#!/usr/bin/env Rscript

#======================
# script to sort the input locus pairs according to their genomic distance (ascending order)
# required for spline fitting
## if the input reference genome is circular genome
## the distance calculation is also adjusted
#======================
args <- commandArgs(TRUE)

InpFile <- as.character(args[1])
OutFile <- as.character(args[2])
CircularGenome <- as.integer(args[3])
ChrSizeFile <- as.character(args[4])

## copy the first line (header) in the output file
## also append one field: genomic distance
system(paste("head ", InpFile, " | awk \'{if (NR==1) {print $0\"\tDist\"}}\' - > ", OutFile))

# ## also get the number of fields in the input file
# nfield <- as.integer(system(paste("cat", OutFile, " | awk \'{print NF}\' -"), intern = TRUE))
# cat(sprintf("\n Number of columns in the input file : %s ", nfield))

## file to store the genomic distances for individual locus pairs
temp_Dist_LocusPairFile <- paste0(dirname(OutFile), '/temp_Dist_LocusPair.bed')
if (file.exists(temp_Dist_LocusPairFile)) {
	system(paste("rm", temp_Dist_LocusPairFile))
}

## file to store the individual genomic distance values
tempDistFile <- paste0(dirname(OutFile), '/temp_Genomic_Distance.bed')

if (CircularGenome == 0) {
	## traditional reference genome
	system(paste0("awk -F\'[\t]\' \'function abs(v) {return v < 0 ? -v : v} {if (NR>1) {print $0\"\t\"abs($5-$2)}}\' ", InpFile, " > ", temp_Dist_LocusPairFile))
} else {
	## circular genome
	ChrSizeData <- read.table(ChrSizeFile, header=F, sep="\t", stringsAsFactors=F)
	ChrList <- unique(ChrSizeData[,1])
	for (j in 1:length(ChrList)) {
		currchr <- ChrList[j]
		cat(sprintf("\n ==>> circular genome - creating adjusted genomic distance - processing chromosome : %s ", currchr))
		## chromosome size for this chromosome
		max_coord <- ChrSizeData[which(ChrSizeData[,1] == currchr), 2]
		cat(sprintf(" -->> max_coord (from ChrSizeFile) for this chromosome : %s ", max_coord))
		system(paste0("awk -F\'[\t]\' -v R=", max_coord, " -v c=\"", currchr, "\" \'function abs(v) {return v < 0 ? -v : v}; function min(x,y) {return x < y ? x : y} {if ((NR>1) && ($1==c)) {print $0\"\t\" min(min(abs($5-$2),abs((R-$5)+$2)),abs((R-$2)+$5))}}\' ", InpFile, " >> ", temp_Dist_LocusPairFile))
	}
}

system(paste0("awk -F\'[\t]\' \'{print $NF}\' ", temp_Dist_LocusPairFile, " | sort -k1,1n | uniq > ", tempDistFile))
## read the unique distance values
DistData <- read.table(tempDistFile, header=F, sep="\t", stringsAsFactors=F)
DistValues <- sort(as.integer(DistData[,1]))
for (i in 1:length(DistValues)) {
	currDist <- DistValues[i]
	cat(sprintf("\n ---- Processing distance value for sorting locus pairs based on interaction distance: %s ", currDist))
	## old code - don't use the genomic distance (last field)
	# system(paste0("awk -F[\'\t\'] \'{if ($NF == ", currDist, ") {print $0}}\' ", temp_Dist_LocusPairFile, " | sort -k7,7nr | cut -f1-", nfield, " >> ", OutFile))
	## new code - use the genomic distance (last field)
	system(paste0("awk -F[\'\t\'] \'{if ($NF == ", currDist, ") {print $0}}\' ", temp_Dist_LocusPairFile, " | sort -k7,7nr >> ", OutFile))
}	# end distance value loop

if (file.exists(tempDistFile)) {
	system(paste("rm", tempDistFile))
}
if (file.exists(temp_Dist_LocusPairFile)) {
	system(paste("rm", temp_Dist_LocusPairFile))
}
