#!/usr/bin/env Rscript

library(optparse)
library(tools)
library(data.table)

options(scipen=999)
options(datatable.fread.datatable=FALSE)

#==============================
# program to estimate HiChIP contact map from HiC contact map
# using input ChIP coverage data 
# fixed sized bin is used
#==============================

#================================================
option_list = list(
	make_option(c("--chr"), type="character", default=NULL, help="Comma or colon seperated list of chromosome names for which the simulation would be performed. Default NULL."),
	make_option(c("--OutDir"), type="character", default=NULL, help="Base Output directory. Mandatory parameter."),
	make_option(c("--ChrSizeFile"), type="character", default=NULL, help="File containing size of chromosomes for reference genome. Mandatory parameter."),
	make_option(c("--ChIPAlignFile"), type="character", default=NULL, help="ChIP-seq file, either alignment (BAM) format, or in bedgraph (4 column) format. Default is NULL. Mandatory parameter."),
	make_option(c("--HiCMapFile"), type="character", default=NULL, help="File containing HiC locus pairs and contact counts (7 columns), with respect to the current chromosome."),
	make_option(c("--TotalRead"), type="integer", action="store", default=1000000, help="Total number of reads for the simulated HiChIP contact matrix. Default = 1000000 or 1M")	
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#==============
# check input files
#==============
if (is.null(opt$OutDir)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! Output directory is not provided - check the option --OutDir \n", call.=FALSE)
}
if (is.null(opt$ChrSizeFile)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! chromosome size file for the reference genome is not provided - check the option --ChrSizeFile \n", call.=FALSE)
}
if (is.null(opt$ChIPAlignFile)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! ChIP-seq data is not provided - check the option --ChIPAlignFile \n", call.=FALSE)
}
if (is.null(opt$HiCMapFile)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! HiC contact map data is not provided - check the option --HiCMapFile \n", call.=FALSE)
}


#=========
# process input parameters
ChrSizeData <- data.table::fread(opt$ChrSizeFile)

# target chromosome name 
if (!is.null(opt$chr)) {
	chrnames <- as.character(unlist(strsplit(opt$chr,"[,:]")))	
	cat(sprintf("\n ** specified list of target chromosomes for HiChIP simulation: %s ", chrnames))
} else {
	chrnames <- as.vector(ChrSizeData[,1])
	cat(sprintf("\n *** As no target chromosome name is provided explicitly, HiChIP map for the following chromosomes would be tried : %s  ****** \n", paste(chrnames, sep=" ")))
}

# create base output directory
system(paste("mkdir -p", opt$OutDir))

# bin size
if (tools::file_ext(opt$HiCMapFile) == "gz") {
	BinSize <- as.integer(system(paste("zcat", opt$HiCMapFile, " | awk \'{if (NR==2) {print ($3-$2)}}\' - "), intern = TRUE))
} else {
	BinSize <- as.integer(system(paste("awk \'{if (NR==2) {print ($3-$2)}}\'", opt$HiCMapFile), intern = TRUE))	
}

cat(sprintf("\n\n *** Parameter list ****** \n\n"))
cat(sprintf("\n *** Output directory : %s  ****** \n", opt$OutDir))
cat(sprintf("\n *** BinSize : %s  ****** \n", BinSize))
cat(sprintf("\n *** File containing chromosome size of reference genome : %s  ****** \n", opt$ChrSizeFile))
cat(sprintf("\n *** File containing ChIP-seq data : %s  ****** \n", opt$ChIPAlignFile))
cat(sprintf("\n *** HiC chromosome contact map file : %s  ****** \n", opt$HiCMapFile))
cat(sprintf("\n *** Target read count for the simulated HiChIP map: %s  ****** \n", opt$TotalRead))

# cat(sprintf("\n *** ChIP coverage file : %s  ****** \n", opt$ChIPCovFile))

# first get the binned distribution of specified chromosome size file
TargetBinnedChrFile <- paste0(opt$OutDir, '/Binned_Chromosome_Size_RefGenome.bed')
system(paste("bedtools makewindows -g", opt$ChrSizeFile, "-w", BinSize, ">", TargetBinnedChrFile))

# process input ChIP-seq alignment file and get the ChIP seq coverage 
# using the TargetBinnedChrFile
MergedChIPCovFile <- paste0(opt$OutDir, '/ChIP_Coverage_ALL.bed')
if (tools::file_ext(opt$ChIPAlignFile) == "bam") {
	cat(sprintf("\n ===>> processing ChIP-seq alignment file in BAM format : %s ====>> \n", opt$ChIPAlignFile))
} else if (tools::file_ext(opt$ChIPAlignFile) == "gz") {
	cat(sprintf("\n ===>> processing ChIP-seq alignment file in gzipped bedgraph (4 column) format : %s ====>> \n", opt$ChIPAlignFile))
} else {
	cat(sprintf("\n ===>> processing ChIP-seq alignment file in bedgraph (4 column) format : %s ====>> \n", opt$ChIPAlignFile))
}
if (tools::file_ext(opt$ChIPAlignFile) == "gz") {
	# gzipped bedgraph format
	system(paste("zcat", opt$ChIPAlignFile, "| bedtools coverage -a", TargetBinnedChrFile, "-b stdin -counts >", MergedChIPCovFile))
} else {
	# either BAM file or plain bedgraph format
	system(paste("bedtools coverage -a", TargetBinnedChrFile, "-b", opt$ChIPAlignFile, "-counts >", MergedChIPCovFile))
}

# read the ChIP-seq coverage for the input ChIP-seq data
# subject to the specified bin size
Merged_ChIPCovData <- data.table::fread(MergedChIPCovFile, header=F, sep="\t", stringsAsFactors=F)

# now process individual chromosomes and simulate HiChIP data
for (chrIdx in 1:length(chrnames)) {
	CurrChrName <- chrnames[chrIdx]
	cat(sprintf("\n\n ===>> Starting processing the chromosome : %s  ==>>> \n\n", CurrChrName))

	# extract ChIP-seq coverage for the current chromosome
	ChIPCoverageData <- Merged_ChIPCovData[which(Merged_ChIPCovData[,1] == CurrChrName), ]
	if (nrow(ChIPCoverageData) == 0) {
		next
	}

	# extract HiC contact map for the current chromosome	
	# generate a text file of three columns
	# 1. bin 1 index 2. bin 2 index 3. raw contact.count
	HiCInp_BinPair_File <- paste0(opt$OutDir, '/temp_HiC_LocusPair_Dumped_CC.bed')
	if (tools::file_ext(opt$HiCMapFile) == "gz") {
		system(paste0("zcat ", opt$HiCMapFile, " | awk \'{print ($3/", BinSize, ")\"\t\"($6/", BinSize, ")\"\t\"$7}\' - > ", HiCInp_BinPair_File))
	} else {
		system(paste0("awk \'{print ($3/", BinSize, ")\"\t\"($6/", BinSize, ")\"\t\"$7}\' ", opt$HiCMapFile, " > ", HiCInp_BinPair_File))
	}
	nline <- as.integer(system(paste("cat", HiCInp_BinPair_File, "| wc -l"), intern = TRUE))
	if (nline == 0) {
		next
	}

	TargetOutDir <- paste0(opt$OutDir, '/Simulate_HiChIP_chr_', CurrChrName, '_B_', BinSize)
	system(paste("mkdir -p", TargetOutDir))

	# file to dump the normalized regression factors
	NormRegrFactorFile <- paste0(TargetOutDir, '/Dumped_Balanced_CC_RegrFactor.bed')

	# statistics regarding the iterations
	StatFile <- paste0(TargetOutDir, '/Dumped_Balanced_CC_Statistics.bed')

	# number of bins in the current chromosome
	# used to determine the dimension of output matrix
	N <- nrow(ChIPCoverageData)

	# Store the ChIP coverage vector
	ChIPCoverageVec <- as.vector(ChIPCoverageData[, ncol(ChIPCoverageData)])

	# define a symmetric matrix of N X N dimension
	# to store the regression factors
	RegrFactorMat <- matrix(0, nrow=N, ncol=N)
	cat(sprintf("\n Initialized regression matrix of %s X %s dimension : ", nrow(RegrFactorMat), ncol(RegrFactorMat)))	

	HiCInp_BinPair_Data <- data.table::fread(HiCInp_BinPair_File, header=F, sep="\t", stringsAsFactors=F)
	cat(sprintf("\n Read HiC input Contact map (bin pair) information for the current chromosome"))

	# important - sourya
	# in the HiCInp_BinPair_Data, first column is matrix row indices
	# second column is matrix column indices
	# and the third column is contact count (raw / normalized)
	rowIdxVec <- HiCInp_BinPair_Data[, 1]
	colIdxVec <- HiCInp_BinPair_Data[, 2]	
	ValVec <- HiCInp_BinPair_Data[, 3]

	# use the following link
	# http://eamoncaddigan.net/r/programming/2015/10/22/indexing-matrices/
	# to assign values to matrix locations
	# wrong method: RegrFactorMat[rowIdxVec, colIdxVec] <- ValVec
	# right method:
	RegrFactorMat[rowIdxVec + nrow(RegrFactorMat) * (colIdxVec-1)] <- ValVec
	cat(sprintf("\n Loaded data in regression matrix of %s X %s dimension : ", nrow(RegrFactorMat), ncol(RegrFactorMat)))

	# define a symmetric matrix - Note: casting in form of a matrix is essential
	# comment - sourya
	# RegrFactorMat <- as.matrix(Matrix::forceSymmetric(RegrFactorMat, uplo="U"))
	# add - sourya
	# using rows as columns and columns as rows
	RegrFactorMat[colIdxVec + nrow(RegrFactorMat) * (rowIdxVec-1)] <- ValVec
	cat(sprintf("\n Made regression coefficient matrix - symmetric - %s X %s dimension : ", nrow(RegrFactorMat), ncol(RegrFactorMat)))

	# normalize the ChIP coverage vector
	# such that its mean is 1
	m <- mean(ChIPCoverageVec)
	if (m > 0) {
		ChIPCoverageVec_Norm <- (ChIPCoverageVec / (1.0 * m))
	}

	# normalize the RegrFactorMat row wise
	# such that mean of each row is 1
	r <- rowMeans(RegrFactorMat)
	r[r==0] <- 1	# to avoid division by 0 error
	RegrFactorMat <- sweep(RegrFactorMat, 1, r, "/")	# row wise division by mean values

	#===============
	# iteration - scale the matrix to have row and column sum = ChIP coverage
	# alternatively scale the row and column vectors
	#===============
	MAX_ITER <- 2	#10	#500	# sourya - modify later
	Thr <- 0.01		# 0.00001	# iteration stop threshold
	ThrDiff <- 0.02 	# threshold to track changes of SAD in successive iterations
	newMat <- matrix(0, nrow=N, ncol=N)	# matrix, to store the current output
	flag <- FALSE	# flag variable indicating convergence

	for (iter_cnt in c(1:MAX_ITER)) {
		cat(sprintf("\n\n iteration count : %s ", iter_cnt))

		if ((iter_cnt %% 2) == 1) {
			# row wise normalization
			s <- rowSums(RegrFactorMat)
			s[s==0] <- 1	# to avoid division by 0 error
			# first multiply with the target ChIP coverage (row wise)
			newMat <- sweep(RegrFactorMat, 1, ChIPCoverageVec_Norm, "*")
			# then normalize with the rowSums (row wise)
			newMat <- sweep(newMat, 1, s, "/")
		} else {
			# column wise normalization
			s <- colSums(RegrFactorMat)
			s[s==0] <- 1	# to avoid division by 0 error
			# first multiply with the target ChIP coverage (row wise)
			newMat <- sweep(RegrFactorMat, 2, ChIPCoverageVec_Norm, "*")
			# then normalize with the colSums (column wise)
			newMat <- sweep(newMat, 2, s, "/")		
		}
		cat(sprintf("\n Copied to newmat entries "))

		# newMat <- as.matrix(Matrix::forceSymmetric(newMat))
		# cat(sprintf("\n Made newmat symmetric "))

		# check the difference between earlier matrix and the new matrix
		# if it is less than the allowed threshold, quit
		# diff_val <- ((sum(abs(newMat - RegrFactorMat)) * 1.0) / (N * N))
		diff_val <- sum(abs(newMat - RegrFactorMat))
		cat(sprintf("\n difference between original and new matrix (SAD) : %s ", diff_val))
		if (diff_val < Thr) {
			cat(sprintf("\n\n\n *** Convergence of iterations is achieved !!!"))
			flag <- TRUE
		}

		# following check is performed after 5 iterations
		# check if the change of SAD in successive iterations is also low
		if (iter_cnt > 5) {
			if (abs(prev_diff_val - diff_val) < ThrDiff) {
				cat(sprintf("\n\n\n *** prev_diff_val : %s  current diff_val : %s  --- Convergence of iterations is achieved !!!", prev_diff_val, diff_val))
				flag <- TRUE
			}
		}

		# store the current SAD 
		prev_diff_val <- diff_val

		# assign it to the regression factor matrix
		RegrFactorMat <- newMat

		# clear newMat
		newMat[newMat>0] <- 0

		# find the correlation between the generated RegrFactorMat
		# in terms of its row wise coverage
		# and column wise coverage
		# with the reference ChIP seq coverage
		RS <- rowSums(RegrFactorMat)
		CS <- colSums(RegrFactorMat)
		idx1 <- which(!is.na(RS) & !is.na(ChIPCoverageVec))
		R_corr <- cor(ChIPCoverageVec[idx1], RS[idx1])
		idx2 <- which(!is.na(CS) & !is.na(ChIPCoverageVec))
		C_corr <- cor(ChIPCoverageVec[idx2], CS[idx2])

		# cat(sprintf("\n R_corr : %s  C_corr : %s ", R_corr, C_corr))

		# dump the iteration wise statistics
		currdf <- data.frame(Itn=iter_cnt, Diff=diff_val, R_corr=R_corr, C_corr=C_corr)
		if (iter_cnt == 1) {
			write.table(currdf, StatFile, sep="\t", quote=F, row.names=F, col.names=T, append=F)	
		} else {
			write.table(currdf, StatFile, sep="\t", quote=F, row.names=F, col.names=F, append=T)
		}
		if (flag == TRUE) {
			break
		}
	}	# end iterations 

	# # # make the matrix RegrFactorMat symmetric
	# # RegrFactorMat <- as.matrix(Matrix::forceSymmetric(RegrFactorMat))
	# ValVec <- RegrFactorMat[rowIdxVec + nrow(RegrFactorMat) * (colIdxVec-1)]
	# RegrFactorMat[colIdxVec + nrow(RegrFactorMat) * (rowIdxVec-1)] <- ValVec

	# dump the contact map entries after matrix balancing (basically the output RegrFactorMat)
	# subject to the locus pairs having nonzero HiC contact counts
	# (specified in rowIdxVec and colIdxVec)
	OutRegrFactor_Norm <- RegrFactorMat[rowIdxVec + nrow(RegrFactorMat) * (colIdxVec-1)]

	# now modify (scale) the entries in RegrFactorMat
	# according to the specified ChIP coverage in the original ChIPCoverageVec
	s <- rowSums(RegrFactorMat)
	s[s==0] <- 1	# to avoid division by 0 error
	RegrFactorMat <- sweep(RegrFactorMat, 1, s, "/")
	RegrFactorMat <- sweep(RegrFactorMat, 1, ChIPCoverageVec, "*")
	OutRegrFactor_Scaled <- RegrFactorMat[rowIdxVec + nrow(RegrFactorMat) * (colIdxVec-1)]

	# scale the contact count to the total sum of HiChIP reads
	s1 <- as.integer(opt$TotalRead)
	s2 <- sum(OutRegrFactor_Scaled)
	scaling_factor <- (s1 * 1.0) / s2
	# cat(sprintf("\n\n *** raw HiChIP count sum s1 : %s estimated count sum s2 : %s scaling_factor employed : %s ***** \n\n", s1, s2, scaling_factor))

	# contact count estimated after scaling
	NewCCVec <- as.integer(round(OutRegrFactor_Scaled * scaling_factor))

	# create the final data frame of scaled contact count for this chromosome
	OutDF <- cbind.data.frame(rep(CurrChrName, length(ValVec)), ((rowIdxVec * BinSize) - BinSize), (rowIdxVec * BinSize), rep(CurrChrName, length(ValVec)), ((colIdxVec * BinSize) - BinSize), (colIdxVec * BinSize), ValVec, NewCCVec)
	colnames(OutDF) <- c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'OrigCC', 'EstimatedCC')
	write.table(OutDF, NormRegrFactorFile, sep="\t", quote=F, append=F, col.names=T, row.names=F)

}	# end chromosome loop

if (file.exists(HiCInp_BinPair_File) == TRUE) {
	system(paste("rm", HiCInp_BinPair_File))
}
if (file.exists(TargetBinnedChrFile) == TRUE) {
	system(paste("rm", TargetBinnedChrFile))
}

