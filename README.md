
FitHiChIP_HiCPro
----------------

FitHiChIP_HiCPro targets for analyzing HiChIP data and deriving the statistical significant interactions, subject to a 
distance threshold. 

Generally, mid to long range CIS interactions (20 Kb to 2 Mb) are sought, although this range can be supplied by user. 
We also experimented with the distance thresholds 10 Kb to 3 Mb.

The package uses output of HiC-Pro pipeline (paired end reads) and applies the package FitHiC (custom implemented to 
process HiChIP data) to derive the significant interactions.

User can opt for one of the following two outputs (statistically significant interactions):

1) Basic FitHiC with spline fit based p-value computation and subsequent FDR filtering (with a custom FDR threshold - default 0.01)

2) Or, FitHiC applied with a bias correction method. The bias for a genomic interval is defined as the ratio of its coverage with the mean coverage computed for all the intervals. Interactions are filtered according to the bias values for both of the interacting regions. Probability of a particular interaction also depends on their bias values. These probability values are used to compute the p-values and determine the statistical significance of a particular interaction.

Choice of a particular option can be provided through the configuration parameters (described below).

Prerequisites
-------------

FitHiChIP_HiCPro requires the following packages / libraries to be installed:

1) The package HiC-pro (obtained from "https://github.com/nservant/HiC-Pro")
	
	The package is used to align two input fasta files (read pairs). The output is a list of valid paired end reads. 
	Details of HiC-pro pipeline are described in the link: https://github.com/nservant/HiC-Pro

	If the complete HiC_Pro pipeline is executed, an interaction matrix is generated from the paired end reads. The matrix is computed according to the bin size (resolution) specified.

	Otherwise, user may use HiC-pro pipeline upto the alignment stage. In such a case, generated paired end reads need to be applied to the utility "build_matrix" (placed under the directory "scripts" of the HiC-Pro installation directory) to produce such interaction matrix.

	An interaction matrix is represented by two files: 1) The first file lists the intervals used for the matrix. Size of an interval equals the input bin size. 2) The second file lists the interactions (contacts) among different intervals. Both CIS and TRANS interactions are reported.

	**** Note: For convenience of the user, FitHiChIP_HiCPro can be executed with the following two input settings:

	1) Either user provides both the valid pairs and interaction matrix files.
	2) Or, user only provides the valid pairs file, but does not provide the interaction matrix files. In such a case, FitHiChIP_HiCPro itself computes the interaction matrix. User only needs to provide the installation directory of HiC-pro package in this case.

2) Python (version 2.7 was used for development)
	
	Python libraries to be installed:
		OptionParser (from the library optparse), gzip, networkx

3) R (version 3.3.3 was used for development)
	
	R Libraries to be installed:
		optparse, ggplot2, splines, fdrtool, parallel

FitHiChIP_HiCPro is tested in linux environment, and requires bash for executing the main script.

Extracting the archieve
-----------------------

A) Download the source code and test data from GitHub.

B) The folder "TestData" contains the following files:
		
		1) Sample_ValidPairs.txt.gz: Sample valid pairs file, an output from HiC-Pro pipeline.
		
		2) Sample.Peaks.gz: Reference peak information

			Unzip the peak file, by applying the following command:

				gunzip Sample.Peaks.gz
		
		3) MboI_hg19_RE_Fragments.bed.gz: Restriction fragments generated using MboI restriction enzyme on the reference genome 'hg19'

			Extract the contents of this peak file, by applying the command:

				gunzip MboI_hg19_RE_Fragments.bed.gz

Note: For test data involving different reference genome and restriction fragments, user needs to first download or arrange for such reference genome, restriction fragment, and reference peak information. Subsequently, location (path) of these files needs to be mentioned via the configuration file (described below in detail).

Download data (mandatory step)
------------------------------

	Given the above mentioned valid pairs and peak files, user needs to first download / copy a few files and place them within the "TestData" folder:
		
		1) chromosome size file (of the reference chromosome hg19): Download by using the following link:
			
			http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
		
		2) Reference genome (hg19) fasta sequence and its index:

			2.1) First, user may download the reference genome in 2 bit format, by using the following link: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

			2.2) Subsequently, an utility program "twoBitToFa" (UCSC basic tools) may be used to extract .fa file(s) from this archieve. The program can be downloaded from the link: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/. For a detailed description, please refer to http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/.
			
			2.3) Index the converted fasta file (to create .fai file)

			2.4) Place the generated .fa and .fai files within the folder "TestData" (or to any directory as per choice)

		3) Download the reference mappability information:

			3.1) In BigWig format: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig

	 		3.2) Convert it to bedgraph format, using the UCSC utilty "BigWigToBedgraph" by using the following command:
	 		
	 			BigWigToBedgraph wgEncodeCrgMapabilityAlign50mer.bigWig wgEncodeCrgMapabilityAlign50mer.bedGraph

	 			Place this mappability file within the "TestData" folder (or to any other folder as per choice).

	****** Note: if user wants to test another dataset (different paired end read file and different peak file), user needs to download 
	the appropriate genome reference and mappability related files.


Execution
---------

The shell script "FitHiChIP_HiCPro.sh" is the main executable.

A text file "configfile" lists the configuration parameters (described below in detail). It also lists the locations of input files and reference genome files (assuming all the files are located in the folder "TestData").

The tool is executed by typing the following command in a bash terminal:

FitHiChIP_HiCPro.sh -C configfile

Where, configfile denotes the configuration file of the current execution.
	

Setting the configuration file
-------------------------------

Each entry of the configuration file has the following format:

Param=ParamValue

where "Param" indicates one parameter and "ParamValue" is the corresponding value (numeric or string format).

**** Note: File and directory paths in the configuration file need to be mentioned either by absolute paths, or by relative paths (with respect to the directory having the configuration file).

FitHiChIP_HiCPro employs the following parameters (and their default values):

#===================================
A) Input File name related parameters:
#======================================

A.1) ValidPairs: Output of HiC-pro pipeline, as a .txt or .txt.gz file. 

	Contains the paired end reads (one line for each read) involving cis or trans interactions. 

	Current pipeline extracts the CIS interactions from the given validpairs file.
	
	**** Note: If not specified, the parameters (B) and (C) (mentioned below) must be provided.

A.2) Interval: File depicting the bins of the interaction matrix. Size of an interval depends on the bin size. Individual bins are also assigned a distinct number. By default, this file name ends with a suffix '_abs.bed'.

A.3) Matrix: File listing the number of interactions (contacts) among the bins listed in the "Interval" file.

	#=========================
	
	The user may issue the following command to create interaction matrix from the given validpairs file:

	zcat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 

	Where,

		$InpValidPairsFile: input valid pairs file (in .txt.gz format)
		
		$MatrixBuildExec: "scripts/Build_Matrix" executable in the HiC-pro package.
		
		$BIN_SIZE: bin size to be specified (integer). For instance, 5000 means 5 Kb bin size.
		
		$ChrSizeFile: File depicting the size of individual chromosomes, corresponding to the reference chromosome employed to create the valid pairs file from HiC-pro package.
		
		$OutPrefix: output prefix (including the output directory) used in the names of interval and matrix files.

	Note: If the user does not pre-compute these matrices, he can leave the fields (B) and (C) blank. In such a case, the valid pairs file (parameter A) must be provided. FitHiChIP_HiCPro computes the matrix from the given valid pairs file.

	#=========================

A.4) PeakFile: Peak detection output from the reference ChIP seq signal (generally, MACS2 is used to compute the reference ChIP seq peaks).
	
		**** Note: this is a mandatory parameter. The reference peaks are used to compute the interactions involving peak segments. 
		
		**** User may separately compute the peaks, or download the reference peaks from ENCODE.

A.5) OutDir: Output directory containing all the results. Default: present working directory.

#===================================
B) Reference genome related parameters:
#======================================

B.1) RefGenome: Reference genome name for the current alignment (valid pairs file). Default: hg19

B.2) ChrSizeFile: File containing the reference chromosome size information. 

	Default: chrom_hg19.sizes (already provided within the TestData folder).

B.3) RefFasta: Fasta formatted sequence of the reference chromosome. Default: hg19.fa

	User needs to download the reference chromosome fasta file, using the above mentioned description.

	Note: Reference fasta sequence is used to compute the mappability and GC content information

B.4) MappabilityFile: Reference mappability file (according to the reference genome).

	Downloaded from the site  http://hgdownload.cse.ucsc.edu/goldenPath/
	
	Should be provided in the bedgraph format
	
	If user has bigWig file, the following UCSC utility is to be used for conversion in bedGraph format
		
		BigWigToBedgraph inp.bw out.bedGraph

B.5) REFragFile: File containing the restriction fragment information, with respect to the reference genome and the restriction site.
	
		The file is of the following format:
				
				chr     interval_start  	interval_end
	
	The MboI restriction fragment file "MboI_hg19_RE_Fragments.bed" (most commonly used restriction fragment in various HiChIP pipelines) is provided in the TestData folder. For other restriction fragment files, please refer to the HiC_Pro manual to create them.
	
B.6) GCSize: Size of the window upstream and downstream the restriction site. 

	Used to compute the GC content information. Default value is 200 (as per the specification in the package HiCNorm [PMID: 23023982])

B.7) MappSize: Size of the window upstream and downstream the restriction site. 

	Used to calculate the mappability information. Default value is 500 (according to the package HiCNorm)

#===================================
C) FitHiC related parameters:
#======================================

C.1) BINSIZE: Size of the bins, depicting the resolution employed. Default= 5000 (means 5 Kb resolution)

C.2) LowDistThr: Lower distance threshold of interaction between two intervals (CIS). Default: 20000 (indicates 20 Kb). Interactions below this distance threshold will not be considered for statistical significance.

C.3) UppDistThr: Upper distance threshold of interaction between two intervals (CIS). Default: 2000000 (indicates 2 Mb). Interactions above this distance threshold will not be considered for statistical significance.

C.4) QVALUE: Minimum FDR (q-value) cutoff for detecting significant interactions. Default: 0.01

C.5) NBins: Max number of equal occupancy bins employed in the FitHiC. Default: 200

	*** Note: FitHiC with equal occupancy bin is used since it is recommended in the FitHiC settings.

C.6) BeginBiasFilter: A binary option. If 1, indicates that the interactions will be pre-filtered according to their bias values. Default 1.

C.7) EndBiasFilter: A binary option. If 1, indicates that the bias correction is enabled. Here, probability of interactions (used in FitHiC) will be multiplied with the bias values.

C.8) biaslowthr: A number (integer or fraction) depicting the lower cutoff of bias. Default 0.2. Applied when the parameter "BeginBiasFilter" is 1.

C.9) biashighthr: A number (integer or fraction) depicting the higher cutoff of bias. Default 5. Applied when the parameter "BeginBiasFilter" is 1.

**** Note: If "BeginBiasFilter" option is 1, interactions having bias values (both reads) within the thresholds "biaslowthr" and "biashighthr" are only considered for statistical significance test. If "BeginBiasFilter" option is 0, these thresholds do not have any impact.

#===================================
D) Miscellaneous parameters:
#======================================

D.1) HiCProBasedir: The installation directory of HiC-pro package. Required if the user does not supply the interaction matrix, and the package requires to compute them using the input valid pairs file.

	*** Note: If the matrices are already computed and provided as the input (using the parameters (B) and (C), this option is not required)

D.2) PREFIX: Prefix string used before any output file name.

D.3) Draw: A binary variable (1/0). If 1, draws various statistical and analysis plots. Default 0.

D.4) TimeProf: A binary variable (1/0). If 1, logs the time elapsed during individual steps. Default 0.


Describing the Output
---------------------

Within the specified output directory "OutDir", following files and folders exist:

1) Parameters.txt: a file listing the parameters used for current execution.

2) HiCPro_Matrix_BinSize"BINSIZE": If the parameters (B) and (C) are empty, this folder contains the interaction matrix files generated from the input valid pairs file.

	The interaction matrix files are "MatrixHiCPro_abs.bed" and "MatrixHiCPro.matrix".

	CIS interactions are also dumped separately in the file "*.interactions.initial.bed".

	According to the distance thresholds "l" and "h" specified in the parameters (C.2) and (C.3), a directory L_"l"_U"u" is created within this folder. A file "*.cis.interactions.DistThr.bed" within this directory contains the CIS interactions filtered within the specified distance thresholds.

3) The folder "NormFeatures" contains features computed for individual genomic intervals (with respect to the specified bin size). 

	Specifically, the file "*.AllBin_CompleteFeat.bed" is of the following format:
		Columns 1 to 3 denote the chromosome interval.
		
		Column 4 denote its coverage
		
		Column 5 denotes whether this interval belongs to a peak segment (1) or non-peak segment (0) where the reference peak information is provided as the parameter (D) of the configuration file.
		
		Column 6 denotes the bias of this interval, obtained by coverage normalization technique.
		
		Columns 7 and 8 denote the mappability and GC contents, respectively, for the current interval.

4) The folders FitHiChIP_"INTTYPE"_b"BINSIZE"_L"LowDistThr"_U"UppDistThr" contain the results for different types of interactions, subject to the specified bin size, and the distance thresholds specified.
	
		INTTYPE: interaction type: can be either ALL2ALL (interaction among all segments), Peak2Peak (interactions between two peaks), Peak2NonPeak (interactions between a peak and a non-peak), and Peak2ALL (interactions involving a peak segment - the other end can be a peak or a non-peak).
	
		Within each directory, following files and folders exist:

			4.1) Interactions.bed: File depicting the interactions among different segments (according to the criteria of peak / non-peak / all employed). 
	
			4.2) Interactions.sortedGenDist.bed: The interactions are sorted in increasing order of distance between the interacting segments.
	
			4.3) FitHiC_EqOccBin: If the options (C.6) and (C.7) are set 0 (please refer to the section "Execution"), this folder contains the significant interactions generated by FitHiC when no bias correction is used.

				4.3.1) Bin_Info.log: File listing the bins employed in FitHiC
				4.3.2) *.interactions_FitHiCPass1.bed: Lists all interactions along with their probabilities, P and Q values.
				4.3.3) *.interactions_FitHiCPass1_FILTER.bed: Lists the significant interactions by applying an FDR threshold (Q value threshold) of 0.01
				4.3.4) *.interactions_FitHiCPass1_FILTER_WashU.bed: The significant interactions are converted to this text file, which can be used for display on WashU epigenome browser.

			4.4) FitHiC_EqOccBin_BiasCorr_"biaslowthr"_"biashighthr"_b"0/1"_e"0/1": This directory is generated when any of the options (C.6) and (C.7) is 1.

					If option (C.6) = 1, the folder name has '-b1' as a substring.
					If option (C.7) = 1, the folder name has '-e1' as a substring.

					Contents of this folder are similar to the folder "FitHiC_EqOccBin"


Contact
--------

For any query, please contact:

Sourya Bhattacharyya (sourya@lji.org)
Ferhat Ay (ferhatay@lji.org)




