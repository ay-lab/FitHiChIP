
FitHiChIP_HiCPro
----------------

Software for analyzing HiChIP data and deriving the statistical significant interactions.
Uses outputs of HiC-Pro pipeline (the valid paired end reads), and applies the package FitHiC (custom implemented to 
process HiChIP data) to derive the significant interactions.
Also employs a basic normalization method, by using the bias correction technique mentioned in FitHiC.
		The bias is computed by taking the ratio of coverage (of the current interval) with the mean coverage of all the intervals.
		Bias values for both of the interacting regions are accounted to filter the interactions, and compute the relative probability of a particular interaction.

Prerequisites
-------------
FitHiChIP_HiCPro requires the following packages / libraries to be installed:

1) The package HiC-pro (obtained from "https://github.com/nservant/HiC-Pro")
	The package is used to align input fasta files to produce a valid pairs file (containing paired end reads). 
	This file is applied to another utility named "build_matrix" to produce the interaction matrix at a specified resolution.
	FitHiChIP_HiCPro requires the valid pairs file and (optionally) the interaction matrix files as the input.
2) Python (version 2.7 was used for development)
	Python libraries to be installed:
		OptionParser (from the library optparse), gzip
3) R (version 3.3.3 was used for development)
	R Libraries to be installed:
		optparse, ggplot2, splines, fdrtool, parallel

FitHiChIP_HiCPro is tested in linux environment, and requires bash for executing the main script.

Execution
-----------
Upon extracting the archieve, the file "FitHiChIP_HiCPro.sh" is the main executable.
A text file "configfile" lists the configuration parameters (described below in detail).
The folder "TestData" contains a sample valid pairs file, and associated files with respect to the reference genome hg19.

The tool is executed by typing the following command in a bash terminal:

FitHiChIP_HiCPro.sh -C configfile -b 0/1 -e 0/1 -l 0.2 -h 5

Where:
	configfile: configuration file listing input data and various input options, output directory, etc.
	-b 0/1: A binary option. If 1, indicates that the output interactions will be filtered according to the bias values.
	-e 0/1: A binary option. If 1, indicates that the probability of interactions (used in FitHiC) will be multiplied with the bias values.
	-l NUM1: A number (integer or fraction) depicting the lower cutoff threshold of bias. Default 0.2
	-h NUM2: A number (integer or fraction) depicting the higher cutoff threshold of bias. Default 5
		If -b option is 1, the interactions whose both ends (reads) have bias values within the thresholds specified in -l and -h options, are processed for finding the significant interactions using FitHiC.
		If -b option is 0, these thresholds do not have any impact.


Setting the configuration file
-------------------------------

The configuration file has the following format:
Param=ParamValue
where "Param" indicates one parameter and "ParamValue" is the corresponding value (numeric or string format).

Note: The file and directory paths are to be either mentioned in the absolute path format, or the relative path with respect to the configuration file itself.

Different parameters of th current pipeline, and their default values are summarized below:

A) ValidPairs: 
	Output of HiC-pro pipeline, as a .txt or .txt.gz file. It contains the paired end reads (one line for 
	each read) involving cis or trans interactions. Current pipeline extracts the CIS interactions from the 
	given validpairs file.
	Note: If not specified, the parameters (B) and (C) need to be specified.

B) Interval:
	The "scripts/Build_Matrix" utility in the HiC-pro package generates interaction matrix from the given valid pairs file. Such a matrix is represented by two files. The first file, also known as the interval file, by default ends with a suffix '_abs.bed'. It contains the binned intervals used for a particular execution. The binning is performed according to the bin size specified by the user. Individual bins are also assigned a distinct number.

C) Matrix:
	The interaction matrix between the bins (with respect to the bins mentioned in the interval file).

	The user may issue the following command to create interaction matrix from the given validpairs file:
	zcat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 

	Where,
		$InpValidPairsFile: input valid pairs file (in .txt.gz format)
		$MatrixBuildExec: "scripts/Build_Matrix" executable in the HiC-pro package.
		$BIN_SIZE: bin size to be specified (integer). For instance, 5000 means 5 Kb bin size.
		$ChrSizeFile: File depicting the size of individual chromosomes, corresponding to the reference chromosome employed to create the valid pairs file from HiC-pro package.
		$OutPrefix: output prefix (including the output directory) used in the names of interval and matrix files.

	Note: If the user does not pre-compute these matrices, he can leave the fields (B) and (C) blank. The current pipeline computes the matrix from the given valid pairs file (parameter A).

D) PeakFile:
	Peak detection output from the reference ChIP seq signal (generally, MACS2 is used to compute the reference ChIP seq peaks).
	Note: this is a mandatory parameter. The reference peaks are used to compute the interactions involving peak segments.

E) OutDir:
	 Output base directory under which all results will be stored.
	 Default: present working directory.

F) RefGenome:
	Reference genome for the current alignment (and valid pairs file).
	Default: hg19

G) ChrSizeFile:
	Path of the file containing the information of reference chromosome size
	Default: chrom_hg19.sizes (provided within the TestData folder)

H) RefFasta:
	Fasta formatted sequence of the reference chromosome.
	Default: hg19.fa
	Note: this is used to compute the mappability and GC content information

I) MappabilityFile:
	Absolute path of the Reference mappability file (according to the reference genome)
	(may be downloaded from the site  http://hgdownload.cse.ucsc.edu/goldenPath/)
	Should be provided in the bedgraph format
	If user has bigWig file, the following UCSC utility is to be used for conversion to bedGraph format
		BigWigToBedgraph inp.bw out.bedGraph

	A default mappability file "wgEncodeCrgMapabilityAlign50mer.bedGraph" is provided in the TestData folder, corresponding to the hg19 reference chromosome.

J) REFragFile:
	File containing the restriction fragment information, with respect to the reference chromosome, and a restriction site.
	The file is of the following format:
				chr     interval_start  	interval_end
	By default, MboI restriction fragment file "MboI_hg19_RE_Fragments.bed" is provided in the TestData folder (most common used restriction fragment in various HiChIP pipelines).
	Creation of a restriction fragment file can be found in the HiC-Pro faqs. 

K) GCSize:
	Size of the window upstream and downstream the restriction site. 
	Used to compute the GC content information. 
	Default value is 200 (according to the package HiCNorm)

L) MappSize:
	Size of the window upstream and downstream the restriction site. 
	Used to calculate the mappability information.
	Default value is 500 (according to the package HiCNorm)

M) BINSIZE:
	Size of the bins, depicting the resolution employed.
	Default= 5000, indicating 5 Kb resolution.	

N) LowDistThr:
	Lower distance threshold of interaction between two intervals (CIS)
	Default: 20000 (indicates 20 Kb)

O) UppDistThr:
	Upper distance threshold of interaction between two intervals (CIS)
	Default: 2000000 (indicates 2 Mb)

P) QVALUE:
	Minimum FDR (q-value) cutoff for detecting significant interactions.
	Default: 0.01

Q) NBins:
	Max number of equal occupancy bins employed in the FitHiC
	Default: 200

R) HiCProBasedir:
	The path of HiC-pro installation. The base path is required to find the utility 'build_matrix' for constructing the interaction 
	matrix from the input valid pairs file.
	Note: If the matrices are already computed and provided as the input (using the parameters (B) and (C), this option is not required)

S) PREFIX:
	Prefix string used before any output file name.

T) Draw:
	A binary variable (1/0). If 1, draws various statistical and analysis plots. Default 0.


Describing the Output
---------------------

Within the specified "OutDir", following files and folders exist:

1) Parameters.txt: a file listing the parameters used for current execution.

2) HiCPro_Matrix_BinSize"BINSIZE": If the parameters (B) and (C) are empty, this folder contains the files "MatrixHiCPro_abs.bed" and "MatrixHiCPro.matrix" which denote the interaction matrices generated from the HiC-Pro utility function "build_matrix".
	The folder also contains the CIS interactions within the file "*.interactions.initial.bed" and further filters these interactions (according to the distance thresholds specified in parameters (N) and (O) of the configuration file) in the file "*.cis.interactions.DistThr.bed".

3) The folder "NormFeatures" contains features computed for individual genomic intervals (with respect to the specified bin size). 
Specifically, the file "*.AllBin_CompleteFeat.bed" is of the following format:
	Columns 1 to 3 denote the chromosome interval.
	Column 4 denote its coverage
	Column 5 denotes whether this interval belongs to a peak segment (1) or non-peak segment (0) where the reference peak information is provided as the parameter (D) of the configuration file.
	Column 6 denotes the bias of this interval, obtained by coverage normalization technique.
	Columns 7 and 8 denote the mappability and GC contents, respectively, for the current interval.

4) The folders FitHiChIP_"INTTYPE"_b"BINSIZE"_L"LowDistThr"_U"UppDistThr" contain the results for different types of interactions, subject to the specified bin size, and the distance thresholds considered for an interaction.
	INTTYPE: interaction type: can be either all to all (interaction among all segments), peak to peak, peak to non peak, and peak to all (includes both peak and non peak segments).
	Within each directory, following files and folders exist:

	4.1) Interactions.bed: File depicting the interactions among different segments (according to the criteria of peak / non-peak / all employed). 
	
	4.2) Interactions.sortedGenDist.bed: The interactions are sorted in increasing order of distance between the interacting segments.
	
	4.3) FitHiC_EqOccBin: If -b and -e options are set 0 (please refer to the section "Execution"), this folder contains the significant interactions generated by FitHiC when no bias correction is used.

		4.3.1) Bin_Info.log: File listing the bins employed in FitHiC
		4.3.2) *.interactions_FitHiCPass1.bed: Lists all interactions along with their probabilities, P and Q values.
		4.3.3) *.interactions_FitHiCPass1_FILTER.bed: Lists the significant interactions by applying an FDR threshold (Q value threshold) of 0.01
		4.3.4) *.interactions_FitHiCPass1_FILTER_WashU.bed: The significant interactions are converted to this text file, which can be used for display on WashU epigenome browser.

	4.4) FitHiC_EqOccBin_BiasCorr_"NUM1"_"NUM2"_b"0/1"_e"0/1":
		If any of the options -b and -e is 1, this folder is created with the input values described in the section "Execution".
		Contents of this folder are similar to the folder "FitHiC_EqOccBin"


Contact
--------

For any query, please contact:

Sourya Bhattacharyya (sourya@lji.org)
Ferhat Ay (ferhatay@lji.org)




