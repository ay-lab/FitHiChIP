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
cccol <- as.integer(args[5])

## temporary file to store the genomic distances for individual locus pairs
temp_Dist_LocusPairFile <- paste0(dirname(OutFile), '/temp_Dist_LocusPair.bed')
if (file.exists(temp_Dist_LocusPairFile)) {
	system(paste("rm", temp_Dist_LocusPairFile))
}

## copy the first line (header) in the output file
## also append one field: genomic distance
system(paste("head ", InpFile, " | awk \'{if (NR==1) {print $0\"\tDist\"}}\' - > ", OutFile))

## get the number of fields in the input file
nfield <- as.integer(system(paste("awk \'{if (NR==1) {print NF}}\' ", OutFile), intern = TRUE))
cat(sprintf("\n Number of columns in the output file : %s ", nfield))

if (CircularGenome == 0) {
	## traditional reference genome
	## sort by the genomic distance and then by the contact count
	system(paste0("awk -F\'[\t]\' \'function abs(v) {return v < 0 ? -v : v} {if (NR>1) {print $0\"\t\"abs($5-$2)}}\' ", InpFile, " | sort -k", nfield, ",", nfield, "n -k", cccol, ",", cccol, "nr - >> ", OutFile))
} else {
	## circular genome
	## construct locus pairs with genomic distance for individual chromosomes
	ChrSizeData <- read.table(ChrSizeFile, header=F, sep="\t", stringsAsFactors=F)
	ChrList <- unique(ChrSizeData[,1])
	for (j in 1:length(ChrList)) {
		currchr <- ChrList[j]
		cat(sprintf("\n ==>> circular genome - creating adjusted genomic distance - processing chromosome : %s ", currchr))
		## chromosome size for this chromosome
		max_coord <- ChrSizeData[which(ChrSizeData[,1] == currchr), 2]
		cat(sprintf(" -->> max_coord (from ChrSizeFile) for this chromosome : %s ", max_coord))
		system(paste0("awk -F\'[\t]\' -v R=", max_coord, " -v c=\"", currchr, "\" \'function abs(v) {return v < 0 ? -v : v}; function min(x,y) {return x < y ? x : y} {if ((NR>1) && ($1==c)) {print $0\"\t\"min(min(abs($5-$2),abs((R-$5)+$2)),abs((R-$2)+$5))}}\' ", InpFile, " >> ", temp_Dist_LocusPairFile))
	}
	## now sort by the genomic distance and then by the contact count
	system(paste0("sort -k", nfield, ",", nfield, "n -k", cccol, ",", cccol, "nr ", temp_Dist_LocusPairFile, " >> ", OutFile))	
}

if (file.exists(temp_Dist_LocusPairFile)) {
	system(paste("rm", temp_Dist_LocusPairFile))
}

