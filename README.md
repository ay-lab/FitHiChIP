
FitHiChIP_HiCPro
----------------

FitHiChIP_HiCPro package ia targeted for analyzing HiChIP data and deriving the statistical significant interactions, subject to a 
distance threshold. 

Generally, mid to long range CIS interactions (10 Kb to 3 Mb) are sought, although this range can be supplied by user.

Uses output of HiC-Pro pipeline (paired end reads) and applies the package FitHiC (custom implemented to 
process HiChIP data) to derive the significant interactions.

User can opt for one of the following two outputs (statistically significant interactions):

1) Basic FitHiC with spline fit based p-value computation and subsequent FDR filtering (with a custom FDR threshold - default 0.01)

2) Or, FitHiC applied with a bias correction method. The bias for a genomic interval is defined as the ratio of its coverage with the mean coverage computed for all the intervals. Interactions are filtered according to the bias values for both of the interacting regions. Probability of a particular interaction also depends on their bias values. These probability values are used to compute the p-value and determine the statistical significance of a particular interaction.

Prerequisites
-------------

FitHiChIP_HiCPro requires the following packages / libraries to be installed:

1) The package HiC-pro (obtained from "https://github.com/nservant/HiC-Pro")
	
	The package is used to align two input fasta files. The output is a list of valid paired end reads. If HiC_Pro is executed upto all its steps (please see the reference manual of HiC-Pro), the paired end reads are used to generate an interaction matrix, according to the bin size (resolution) specified.

	Alternatively, user may apply the paired end reads to the utility "build_matrix" (placed under the directory "scripts" of the HiC-Pro installed package) to produce such interaction matrix (two output files).

	FitHiChIP_HiCPro requires either the valid pairs file or the interaction matrix files as the input. If no interaction matrix is computed, user must provide the paired end read file to FitHiChIP_HiCPro, which computes the interaction matrix itself.
	
2) Python (version 2.7 was used for development)
	
	Python libraries to be installed:
		OptionParser (from the library optparse), gzip

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

Note: For test data involving different reference genome and restriction fragment, user needs to first download or arrange for such reference genome, restriction fragment, and reference peak information. Subsequently, location of those data needs to be mentioned via the configuration file (see below for a detailed description of configuration file).
	

Download data (mandatory step)
------------------------------

	With respect to the above mentioned valid pairs and peak files, user needs to download / copy a few files and place them within the "TestData" folder:
		
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


Execution
---------

The shell script "FitHiChIP_HiCPro.sh" is the main executable.

A text file "configfile" lists the configuration parameters (described below in detail). It also lists the locations of input files and reference genome related files (assuming all the files are located in the folder "TestData").

The tool is executed by typing the following command in a bash terminal:

FitHiChIP_HiCPro.sh -C configfile -b 0/1 -e 0/1 -l 0.2 -h 5

Where:
	
	configfile: configuration file listing input data and various input options, output directory, etc.
	
	-b 0/1: A binary option. If 1, indicates that the bias correction is enabled. Here, interactions will be filtered according to the bias values.
	
	-e 0/1: A binary option. If 1, indicates that the bias correction is enabled. Here, probability of interactions (used in FitHiC) will be multiplied with the bias values.
	
	-l NUM1: A number (integer or fraction) depicting the lower cutoff of bias. Default 0.2
	
	-h NUM2: A number (integer or fraction) depicting the higher cutoff of bias. Default 5
	
		If -b option is 1, interactions having bias values (both reads) within the thresholds of -l and -h options, are only considered for analyzing using FitHiC. If -b option is 0, these thresholds do not have any impact.


Setting the configuration file
-------------------------------

The configuration file has the following format:

Param=ParamValue

where "Param" indicates one parameter and "ParamValue" is the corresponding value (numeric or string format).

Note: File and directory paths in the configuration file are to be mentioned either by their absolute path, or a relative path with respect to the directory having the configuration file.

FitHiChIP_HiCPro employs the following parameters (and their default values):

A) ValidPairs: Output of HiC-pro pipeline, as a .txt or .txt.gz file. It contains the paired end reads (one line for 
	each read) involving cis or trans interactions. 

	Current pipeline extracts the CIS interactions from the given validpairs file.
	
	Note: If not specified, the parameters (B) and (C) (mentioned below) must be provided.

B) Interval: The utility "scripts/Build_Matrix" of the HiC-pro package generates interaction matrix from the given valid pairs file. Such a matrix is represented by two files. The first file is called the "Interval file" which contains the binned intervals (according to the specified bin size). Individual bins are also assigned a distinct number. By defaultm, this file name ends with a suffix '_abs.bed'.

C) Matrix: This is the second file representing the interaction matrix. It lists the number of interactions between the bins (represented by the bin numbers provided in the Interval file).

	The user may issue the following command to create interaction matrix from the given validpairs file:
	zcat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 

	Where,
		$InpValidPairsFile: input valid pairs file (in .txt.gz format)
		$MatrixBuildExec: "scripts/Build_Matrix" executable in the HiC-pro package.
		$BIN_SIZE: bin size to be specified (integer). For instance, 5000 means 5 Kb bin size.
		$ChrSizeFile: File depicting the size of individual chromosomes, corresponding to the reference chromosome employed to create the valid pairs file from HiC-pro package.
		$OutPrefix: output prefix (including the output directory) used in the names of interval and matrix files.

	Note: If the user does not pre-compute these matrices, he can leave the fields (B) and (C) blank. In such a case, the valid pairs file (parameter A) must be provided. FitHiChIP_HiCPro computes the matrix from the given valid pairs file.

D) PeakFile: Peak detection output from the reference ChIP seq signal (generally, MACS2 is used to compute the reference ChIP seq peaks).
	Note: this is a mandatory parameter. The reference peaks are used to compute the interactions involving peak segments. 

		User may separately compute the peaks, or download the reference peaks from ENCODE.

E) OutDir: Output directory containing all the results. Default: present working directory.

F) RefGenome: Reference genome name for the current alignment (valid pairs file). Default: hg19

G) ChrSizeFile: File containing the reference chromosome size information. Default: chrom_hg19.sizes (already provided within the TestData folder).

H) RefFasta: Fasta formatted sequence of the reference chromosome. Default: hg19.fa

	User needs to download the reference chromosome fasta file, using the procedure mentioned above.

	Note: Reference fasta sequence is used to compute the mappability and GC content information

I) MappabilityFile: Reference mappability file (according to the reference genome).

	Downloaded from the site  http://hgdownload.cse.ucsc.edu/goldenPath/
	
	Should be provided in the bedgraph format
	
	If user has bigWig file, the following UCSC utility is to be used for conversion in bedGraph format
		
		BigWigToBedgraph inp.bw out.bedGraph

J) REFragFile: File containing the restriction fragment information, with respect to the reference genome and the restriction site.
	
		The file is of the following format:
				
				chr     interval_start  	interval_end
	
	The MboI restriction fragment file "MboI_hg19_RE_Fragments.bed" (most commonly used restriction fragment in various HiChIP pipelines) is provided in the TestData folder. For other restriction fragment files, please refer to the HiC_Pro manual to create such a reference restriction fragment file.
	
K) GCSize: Size of the window upstream and downstream the restriction site. Used to compute the GC content information. Default value is 200 (according to the package HiCNorm)

L) MappSize: Size of the window upstream and downstream the restriction site. Used to calculate the mappability information. Default value is 500 (according to the package HiCNorm)

M) BINSIZE: Size of the bins, depicting the resolution employed. Default= 5000, indicating 5 Kb resolution.	

N) LowDistThr: Lower distance threshold of interaction between two intervals (CIS). Default: 20000 (indicates 20 Kb)

O) UppDistThr: Upper distance threshold of interaction between two intervals (CIS). Default: 2000000 (indicates 2 Mb)

P) QVALUE: Minimum FDR (q-value) cutoff for detecting significant interactions. Default: 0.01

Q) NBins: Max number of equal occupancy bins employed in the FitHiC. Default: 200

R) HiCProBasedir: The path of HiC-pro installation. The base path is required to find the utility 'build_matrix' for constructing the interaction matrix from the input valid pairs file.
	
		Note: If the matrices are already computed and provided as the input (using the parameters (B) and (C), this option is not required)

S) PREFIX: Prefix string used before any output file name.

T) Draw: A binary variable (1/0). If 1, draws various statistical and analysis plots. Default 0.


Describing the Output
---------------------

Within the specified output directory "OutDir", following files and folders exist:

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

			4.4) FitHiC_EqOccBin_BiasCorr_"NUM1"_"NUM2"_b"0/1"_e"0/1": If any of the options -b and -e is 1, this folder is created with the input values described in the section "Execution". Contents of this folder are similar to the folder "FitHiC_EqOccBin"


Contact
--------

For any query, please contact:

Sourya Bhattacharyya (sourya@lji.org)
Ferhat Ay (ferhatay@lji.org)




