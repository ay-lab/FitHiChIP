
FitHiChIP_HiCPro
----------------

FitHiChIP_HiCPro analyzes HiChIP/PLAC-seq data and derives the statistical significant interactions, given a specific distance range (default = 20 Kb to 2 Mb for 5kb-binned data, resolution and distance range can be user defined as well).

The package starts from valid read pairs from the HiC-Pro pipeline (Servant et. al. 2015), a commonly used HiC data processing pipeline. The valid read pairs are then processed through different steps of FitHiChIP for finding statistically significant interactions. FitHiChIP implements a HiChIP-specific normalization and peak calling routine followed by a p-value calculation similar to FitHiC (Ay and Noble, 2014) but augments with a bias correction technique to find statistically significant interactions.

User can check the configuration parameters (described below) to run FitHiChIP.

Prerequisites
-------------

FitHiChIP_HiCPro requires the following packages / libraries to be installed:

1) Python (version 2.7 was used for development)
	
	Python libraries to be installed:
		OptionParser (from the library optparse), gzip, networkx

2) R (version 3.3.3 / 3.4.0 is recommended)
	
	R Libraries to be installed:
		optparse, ggplot2, splines, fdrtool, parallel

3) The package HiC-pro (obtained from "https://github.com/nservant/HiC-Pro").
	**** We recommend downloading and installing the latest version (2.9.0 or above) of HiC-pro pipeline.

	The package is used to align paired end reads (fasta files). 

	*** After installing HiC-pro, user needs to provide the base installation directory of HiC-pro as 
	an input to the FitHiChIP pipeline (user may check the configuration file of FitHiChIP for details).
 
	Outputs of HiC-pro (to be provided as an input to the FitHiChIP pipeline):

	A) A collection of valid paired end reads (filename: *.validPairs). MUST BE PROVIDED AS AN INPUT TO FitHiChIP.
	
	B) Interaction matrix generated from the input paired-end reads, according to the specified bin size. The matrix is 
	represented by two files: a) Interval file (filename: *_abs.bed) listing the binning intervals, and b) interactions (contacts) 
	among different intervals (filename: *.matrix). Both CIS and TRANS interactions are reported. 

	--- Files related to "Interaction matrix" can be optionally provided by the user to the FitHiChIP pipeline. Otherwise,
	FitHiChIP itself computes the interaction matrix by using the installation directory of the HiC-pro package.

FitHiChIP_HiCPro is tested in linux environment, and requires bash for executing the main script.

Extracting the archieve
-----------------------

A) Download the source code and test data from GitHub.

B) The folder "TestData" contains the following files:
		
		1) Sample_ValidPairs.txt.gz: Sample valid pairs file, an output from HiC-Pro pipeline.
		
		2) Sample.Peaks.gz: Reference peak information

			Unzip the peak file by applying the following command:

				gunzip Sample.Peaks.gz
		
		3) MboI_hg19_RE_Fragments.bed.gz: Restriction fragments generated using MboI enzyme on the reference genome 'hg19'

			Extract the contents of this file by applying the command:

				gunzip MboI_hg19_RE_Fragments.bed.gz

**** Note: For test data involving different reference genome and restriction fragments, 
the user needs to first download or generate files for such reference genome, restriction fragment, 
and reference peak information (example commands are provided below). 

**** Subsequently, location (absolute path) of these files needs to 
be mentioned in the configuration file (as described below).

Download data (mandatory step)
------------------------------

	Given the above mentioned files, user needs to first download / copy a few files and place them within the "TestData" folder:
		
		1) chromosome size file (of the reference chromosome hg19): Download by using the following link:
			
			http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
		
		2) Reference genome (hg19) fasta sequence and its index:

			2.1) First, user may download the reference genome in 2 bit format, by using the following link: 

			http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

			2.2) Subsequently, an utility program "twoBitToFa" (UCSC basic tools) may be used to extract .fa file(s) 
			from this archieve. The program can be downloaded from the link: 

			http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/ 

			For a detailed description, please refer to http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/
			
			2.3) Index the converted fasta file (to create .fai file)

			2.4) Place the generated .fa and .fai files within the folder "TestData" (or to any directory as per choice)

		3) Download the reference mappability information:

			3.1) In BigWig format: 
		http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig

	 		3.2) Convert it to bedgraph format, using the UCSC utilty "BigWigToBedgraph" by using the following command:
	 		
	 			BigWigToBedgraph wgEncodeCrgMapabilityAlign50mer.bigWig wgEncodeCrgMapabilityAlign50mer.bedGraph

	 		3.3) Place this mappability file within the "TestData" folder (or to any other folder as per choice).

	****** Note: For testing a different dataset (different paired end read file and different peak file), 
	user needs to download the appropriate genome reference and mappability related files.

	****** There is no restriction of placing all the files in the TestData folder. The only requirement is to 
	provide the absolute (and correct) paths of different files in the respective configuration options.


Execution
---------

The shell script "FitHiChIP_HiCPro.sh" is the main executable.

Four configuration files are also provided, indicating a choice of coverage / ICE bias correction, 
and use of peak to peak background or not.

The tool is executed by typing the following command in a bash terminal:

FitHiChIP_HiCPro.sh -C configuration_file_name
	

Setting the configuration file
-------------------------------

Each entry of the configuration file has the following format:

Param=ParamValue

where "Param" indicates one parameter and "ParamValue" is the corresponding value (numeric or string format).

**** Note: User should mention the file and directory paths in the configuration file in absolute paths.

#===================================
A) Input File name related parameters:
#======================================

A.1) ValidPairs: Valid pairs generated by HiC-pro pipeline (either simple text format, or can be in gzipped format). 

	**** Note: Mandatory parameter.

A.2) Interval: File depicting the bins of the interaction matrix. Size of an interval depends on the bin size. 
Individual bins are also assigned a distinct number. By default, this file name ends with a suffix '_abs.bed'.

A.3) Matrix: File listing the number of interactions (contacts) among the bins listed in the "Interval" file.

	#=========================
	User may leave the fields A.2) and A.3) blank, if wishes to not pre-compute these matrices. In such a case, 
	FitHiChIP computes the matrix from the parameter A.1.

	Otherwise, user can run the following command to create interaction matrix from the given validpairs file (assumed to be 
		in gzipped format):

	zcat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 

	Where,

		$InpValidPairsFile: input valid pairs file (in .txt.gz format)
		
		$MatrixBuildExec: "scripts/Build_Matrix" executable in the HiC-pro package.
		
		$BIN_SIZE: bin size to be specified (integer). For instance, 5000 means 5 Kb bin size.
		
		$ChrSizeFile: File depicting the size of individual chromosomes, corresponding to the reference chromosome employed to create the valid pairs file from HiC-pro package.
		
		$OutPrefix: output prefix (including the output directory) used in the names of interval and matrix files.

	#=========================

A.4) PeakFile: Reference ChIP-seq peak file (recommended to use standard ENCODE peaks. Otherwise, User may 
separately compute the peaks, using macs2, or even employ the package hichipper to compute the peaks).
	
		**** Mandatory parameter.
		
A.5) OutDir: Output directory which would contain all the results. Default: present working directory.

A.6) HiCProBasedir: Installation directory of HiC-pro package. (example: /custom_path/HiCPro/HiC-Pro_2.9.0/) 

	*** Mandatory parameter.

#===================================
B) Reference genome related parameters:
#======================================

B.1) RefGenome: Reference genome name for the current valid pairs. Can be hg19, mm9, etc. Default: hg19

B.2) ChrSizeFile: File having the reference chromosome size. 

	Default: chrom_hg19.sizes (already provided within the TestData folder).

B.3) RefFasta: Fasta formatted sequence of the reference chromosome (such as hg19, mm9). 
	 Used to compute the mappability and GC content information (downloading instructions are mentioned above).

		Default: hg19.fa

B.4) MappabilityFile: Reference mappability file (according to the reference genome). 
		Download from the site: http://hgdownload.cse.ucsc.edu/goldenPath/

	*** Should be provided in the bedgraph format. For bigWig file, user should use the following UCSC utility for converting to bedGraph format: BigWigToBedgraph inp.bw out.bedGraph

B.5) REFragFile: File containing the restriction fragment information, with respect to the reference 
	genome and the restriction site. The file format (tab delimited):		chr     interval_start  	interval_end
	
	The MboI restriction fragment file "MboI_hg19_RE_Fragments.bed" 
	(most commonly used restriction fragment in various HiChIP pipelines) 
	is provided in the TestData folder. For other restriction fragment files, 
	please refer to the HiC Pro manual for how to create them.
	
B.6) GCSize: Size of the window upstream and downstream the restriction site, for computing the GC content. 
	Default = 200 (as per the specification in the package HiCNorm [PMID: 23023982])

B.7) MappSize: Size of the window upstream and downstream the restriction site, to calculate the mappability information. 
	Default = 500 (as per the specification of HiCNorm [PMID: 23023982])

#===================================
C) Generic parameters for profiling interactions:
#======================================

C.1) IntType: Type of interaction to be computed. Options are:

		(1) peak to peak: contacts between all pairs of peak segments (subject to fixed size binning)

		(2) peak to non peak: contacts between pairs of segments such that one is a peak and the other is a non-peak segment

		(3) peak to all (default): here one interacting segment is a peak, while the other can be a peak or a non-peak. Encapsulates 
		the options 1 and 2. This is the default setting, as HiChIP protocol is targeted to find protein centric interactions.								
		(4) all to all: Interactions between every possible pairs of segments are computed. This resembles Hi-C mode.

		(5) Everything from 1 to 4. That is, all of the above mentioned interactions are computed.

C.2) BINSIZE: Size of the bins, depicting the resolution employed. Default= 5000 (means 5 Kb resolution)

C.3) LowDistThr: Lower distance threshold of interaction between two intervals (CIS). Default: 20000 (indicates 20 Kb). Interactions below this distance threshold will not be considered for statistical significance.

C.4) UppDistThr: Upper distance threshold of interaction between two intervals (CIS). Default: 2000000 (indicates 2 Mb). Interactions above this distance threshold will not be considered for statistical significance.

C.5) QVALUE: Minimum FDR (q-value) cutoff for detecting significant interactions. Default: 0.01

#===================================
D) FitHiC related parameters:
#======================================

D.1) NBins: Max number of equal occupancy bins employed in the FitHiC. Default: 200

	*** Note: Equal occupancy binning (recommended in FitHiC) is used for contact probability estimation.

D.2) UseP2PBackgrnd: Can be 0 or 1. Applicable specifically for peak to all interactions.

	A value of 1 means that peak to peak background will be used for contact probability estimation. That is, 
	contacts between peak segments are only used for background probability estimation. This 
	ensures highly stringent significance estimation.

	Value of 0 means that all of the interactions (of peak to all category) are used for background estimation. 
	This indicates lower stringency of significance estimation.

	User may opt for the value 1 (recommended) if very stringent output (only very highly significant interactions would be reported) 
	is targeted. Otherwise, value of 0 means less stringency (higher number of interactions would be reported).

#===================================
E) Bias correction related parameters:
#======================================

E.1) BiasCorrection: Indicates if the bias correction is enabled (1) or not (0). Default 1 (recommended).

	User should enable the bias correction, since without bias correction, many false positive interactions 
	are reported.

E.2) BiasType: Can be 1 or 2. This parameter signifies the type of bias correction used. 

	A value of 1 means that coverage specific bias is used. Value of 0 means that ICE specific bias is computed. 
	Experiments show that 1 is recommended (higher performance) (default).


#===================================
F) Merging nearby interactions:
#======================================

F.1) MergeInt: Has the value 0 or 1. If 1, means that interactions close to each other (both of their source and 
target bin pairs are very close) are merged and reported as a single interaction. If 0, no merging step is performed. 

 	**** Note: Currently this feature is under development - so recommended 0 for the moment.

#===================================
G) Miscellaneous parameters:
#======================================

G.1) PREFIX: Prefix string used before any output file name.

G.2) Draw: A binary variable (1/0). If 1, draws various statistical and analysis plots. Default 0.

G.3) TimeProf: A binary variable (1/0). If 1, logs the time elapsed during individual steps. Default 0.

G.4) OverWrite: A binary variable (1/0). If 1, overwrites previous FitHiChIP output file in the current execution. Otherwise, if 0, 
skips computation of already existing files.


Describing the Output
---------------------

Within the specified output directory "OutDir", following files and folders exist:

1) Parameters.txt: a file listing the parameters used for current execution.

2) TimingProfile.txt: If the parameter (G.3) is 1, this file lists the time taken for individual steps.

3) HiCPro_Matrix_BinSize"BINSIZE": If the parameters (A.2) and (A.3) are empty, this folder 
contains the interaction matrix files generated from the input valid pairs file. The interaction matrix 
files are "MatrixHiCPro_abs.bed" and "MatrixHiCPro.matrix".

	CIS interactions are also dumped separately in the file "*.interactions.initial.bed".

	According to the distance thresholds "l" and "h" specified in the parameters (C.3) and (C.4), 
	a directory L_"l"_U"u" is created within this folder. 
	A file "*.cis.interactions.DistThr.bed" within this directory contains the 
	CIS interactions filtered within the specified distance thresholds.

3) "NormFeatures": This folder contains features computed for individual genomic intervals / bins  
(with respect to the specified bin size).

	User may look at the file "*.AllBin_CompleteFeat.bed" specifically. The file is of the following format:
		
		Columns 1 to 3: chromosome interval. 		

		Column 4: Coverage of this interval (number of reads mapped onto it)

		Column 5: denotes whether this interval belongs to a peak segment (1) or non-peak segment (0).
		Here, the peak information is with respect to the parameter (A.4) of the configuration file.

		Column 6: Bias value of this interval. Computed either by coverage specific bias correction, or 
		the ICE specific bias normalization, according to the value of the parameter (E.2).
		
		Columns 7 and 8 denote the mappability and GC contents, respectively, for the current interval.

4) The folders FitHiChIP_"INTTYPE"_b"BINSIZE"_L"LowDistThr"_U"UppDistThr" contain the results for different types of interactions, subject to the specified bin size, and the distance thresholds specified.
	
		INTTYPE: interaction type: Have one of the following values:

				ALL2ALL (when the parameter (C.1) = 4)

				Peak2Peak (when the parameter (C.1) = 1)

				Peak2NonPeak (when the parameter (C.1) = 2)

				Peak2ALL (when the parameter (C.1) = 3)

		Within this folder, a directory structure of the following naming pattern is created:

			P2PBckgr_[UseP2PBackgrnd]/[BiasType]

			Where: 

				[UseP2PBackgrnd]: Value of the parameter (D.2)

				[BiasType]: String of either "Coverage_Bias" or "ICE_Bias", according to the value of parameter (E.2).


			4.1) Within this generic directory structure, following files and folders are created:

				4.1.1) Interactions.bed: All contacts among the selected pairs of segments (depending upon the value of INTTYPE). 
	
				4.1.2) Interactions.sortedGenDist.bed: Above contacts are sorted by increasing genomic distance. Used as an input to 
			the custom FitHiC implementation.
	
				4.1.3) FitHiC/BinomDistr: If the option (E.1) = 0, this folder contains the significant interactions generated by Naive FitHiC method (without any bias correction technique).

				4.1.4) FitHiC_BiasCorr_Resid_0_EqOcc_1/BinomDistr: this folder is created when (E.1) = 1.


				Directory structure within the folders (4.1.3) and (4.1.4) are similar. Following files and folders reside within this directory:

						A) Bin_Info.log: File listing the bins employed in FitHiChIP (equal occupancy bins)
						
						B) configuration.txt: list of configuration parameters for FitHiChIP
						
						C) PREFIX.interactions_FitHiC.bed: Lists all interactions along with their probabilities, p and q values, computed 
						using the FitHiChIP method.

						D) PREFIX.interactions_FitHiC_Q*.bed: Lists the significant interactions according to the specified q-value significance threshold.

						E) PREFIX.interactions_FitHiC_Q*_WashU.bed: Significant interactions are converted to this text file, for 
						display on WashU epigenome browser.

						F) Merge_Nearby_Interactions: Directory containing the output of merging nearby interactions (if the parameter 
							(F.1) = 1).
 


Contact
--------

For any query, please contact:

Sourya Bhattacharyya (sourya@lji.org)
Ferhat Ay (ferhatay@lji.org)




