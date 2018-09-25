
FitHiChIP
----------------

Developers: Sourya Bhattacharyya, Ferhat Ay

La Jolla Institute for Allergy and Immunology

La Jolla, CA 92037, USA

**************************

FitHiChIP analyzes HiChIP / PLAC-seq data and derives the statistical significant CIS interactions, given a specific distance range (default = 20 Kb to 2 Mb: can be user defined).

Input of FitHiChIP is the valid CIS read pairs generated from the HiC-Pro pipeline (Servant et. al. 2015), the most commonly 
used Hi-C data processing pipeline. These read pairs are analyzed to find the statistically significant loops within the given distance range. FitHiChIP uses the distance decay model of the interaction count, proposed in the method FitHiC (Ay and Noble, 2014) and also applies a bias specific regression model to find the statistical significant interactions. For details, please check the manuscript (link provided below).

Prerequisites
-------------

FitHiChIP requires the following packages / libraries to be installed:

1) HiC-pro (obtained from "https://github.com/nservant/HiC-Pro")
	
		**** User should install HiC-pro version 2.9.0 or above.

		**** Base installation directory of HiC-pro, along with the following files, are to be provided within the configuration parameters of FitHiChIP (details of the configuration parameters are provided in the next section):

		A) valid paired end reads (filename: *.validPairs). MUST BE PROVIDED AS AN INPUT TO FitHiChIP.
	
		B) Interaction matrix generated from the input paired end reads, according to the specified bin size. 
	
		The matrix is represented by two files: 
		a) Interval file (filename: *_abs.bed) listing the binning intervals, and 
		b) interactions (contacts) among different intervals (filename: *.matrix). 
		Both CIS and TRANS interactions are reported. 

		--- Files related to "Interaction matrix" are optional. If not provided, FitHiChIP computes those matrices using the .validPairs file, and the installed HiC-pro package. 
		
		

2) Python 	

		**** version 2.7 was used for development
	
		Python libraries to be installed:
			1) OptionParser (from the library optparse), 
			2) gzip, 
			3) networkx (https://networkx.github.io/)

3) R 		

		**** we have used version 3.4.3 
	
		R Libraries to be installed:
		1) optparse, 
		2) ggplot2, 
		3) splines, 
		4) fdrtool, 
		5) parallel, 
		6) tools, 
		7) GenomicRanges, 
		8) dplyr
		
		In addition, the R package EdgeR (https://bioconductor.org/packages/release/bioc/html/edgeR.html) is to be installed if user performs differential analysis of FitHiChIP loops.
		
		
4) bedtools (http://bedtools.readthedocs.io/en/latest/)

		**** we have used bedtools version 2.26.0 
	
5) samtools (http://www.htslib.org/)

		**** At least version 1.6 is required. We have tested with version 1.9

6) macs2 for peak calling (https://github.com/taoliu/MACS)


**** FitHiChIP is tested in linux environment, and requires bash for executing the main script.


Extracting the source code archieve
-----------------------------------

A) The source codes are placed within the directory code/

		A.1) The file "sample_script.sh" within this folder contains the basic execution script of this pipeline. The script invokes a configuration file.

		A.2) Four different configuration files are also placed within this folder. Those files vary according to the FitHiChIP bias correction (coverage or ICE bias) method and whether peak to all (loose) background or peak to peak (stringent) background are employed. 

		** All of these configuration files are initialized with the location of the testing data employed (mentioned below). The user only needs to check the parameter "HiCProBasedir" and mention the directory containing HiC-pro installed package.


B) The folder data/ contains the following files:
		
		1) Sample_ValidPairs.txt.gz: Sample valid pairs file, an output from HiC-Pro pipeline.
		
		2) Sample.Peaks.gz: Reference peak information

			Unzip the peak file, by applying the following command:

			gunzip Sample.Peaks.gz
				
			**** Note: The peak file is applicable for the given test data. For any other test data, user needs to download the reference peak files from a reference site such as ENCODE. 

			**** User may also infer peaks from HiChIP data itself (details are described afterwards)
		
		3) MboI_hg19_RE_Fragments.bed.gz: Restriction fragments generated using MboI restriction enzyme on the reference genome 'hg19'. 
		
			**** Note: applying the restriction fragment file as a parameter is now optional !!!!
			**** User can keep the option blank (please refer to the configuration parameters below)
			**** however, if user requires computing the number of restriction sites involved in bins, he can provide the restriction fragment file as an input. 
		
			For the current restriction fragment file, user can extract its contents by applying the command:

			gunzip MboI_hg19_RE_Fragments.bed.gz

			**** Note: For test data involving different reference genome and restriction fragments, user needs to first download or arrange for such reference genome, restriction fragment, and reference peak information (example commands are mentioned below). 

			**** Subsequently, location (preferably absolute path) of these files needs to be mentioned in the configuration file (as described below).

		4) chrom_hg19.sizes: length of individual chromosomes corresponding to the reference genome hg19


Downloading data / parameter files (mandatory step)
----------------------------------------------------

Given the above mentioned valid pairs and peak files, user needs to first download / copy a few files and place them within the "TestData" folder. 

**** Note: some of the following parameters are mentioned as optional. 

**** User can skip downloading the optional parameters and keep the corresponding parameter (in the configuration file) as blank.
		
	1) chromosome size file (such as for the reference chromosome hg19): 

		**** This parameter is MANDATORY.
			
		Download this file by using the following link:

		http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

		For any other reference chromosome, download the corresponding chromosome size information from the UCSC site.
		
	2) Reference genome (hg19) fasta sequence and its index:
		
		***** This parameter is OPTIONAL.
		***** If provided, used to compute the GC content information of respective bins. 

		2.1) First, user may download the reference genome in 2 bit format, by using the following link: 

		http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

		2.2) Subsequently, an utility program "twoBitToFa" (UCSC basic tools) 
		may be used to extract .fa file(s) 
		from this archieve. The program can be downloaded from the link: 

		http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/ 

		For a detailed description, please refer to http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/

		2.3) Index the converted fasta file (to create .fai file)

		2.4) Place the generated .fa and .fai files within the folder "TestData" 
		(or to any directory as per choice)

	3) Reference mappability file (corresponding to the reference genome):
		
		***** This parameter is OPTIONAL.
		***** If provided, used to compute the mappability for interacting bins.
		***** Downloading this file involves the following steps:

		3.1) In BigWig format: 
			http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig

		3.2) Convert it to bedgraph format, using the UCSC utilty "BigWigToBedgraph" 
		by using the following command:
	 		
	 	BigWigToBedgraph wgEncodeCrgMapabilityAlign50mer.bigWig wgEncodeCrgMapabilityAlign50mer.bedGraph

	 	3.3) Place this mappability file within the "TestData" folder (or to any other folder as per choice).

		****** Note: For testing a different dataset (different paired end read file and different peak file), 
		user needs to download the appropriate genome reference (fasta file) and mappability related files.

		****** There is no restriction of placing all the files in the TestData folder. The only requirement is to 
		provide the absolute (and correct) paths of different files in the respective configuration options.

		****** The optional parameters and corresponding files can be skipped by the user, and corresponding configuration options in the text file (described below) can be left as blank !!!


Execution
---------

The shell script "FitHiChIP_HiCPro.sh" is the main executable.

Four configuration files are also provided along with this repository, indicating a choice of coverage / ICE bias correction, 
and use of peak to peak background or not. 

	configfile_BiasCorrection_CoverageBias: when coverage bias correction and loose (peak to all) background is used.
	configfile_BiasCorrection_ICEBias: when ICE bias correction and loose (peak to all) background is used.
	configfile_P2P_BiasCorrection_CoverageBias: when coverage bias correction and stringent (peak to peak) background is used.
	configfile_P2P_BiasCorrection_ICEBias: when ICE bias correction and stringent (peak to peak) background is used.

FitHiChIP is executed by typing the following command in a bash terminal (assuming the executable is in current directory):

sh FitHiChIP_HiCPro.sh -C configuration_file_name
	

Setting the configuration file
-------------------------------

Each entry of the configuration file has the following format:

	Param=ParamValue

where "Param" indicates one parameter and "ParamValue" is the corresponding value (numeric or string format).

**** Note: We recommend users to mention the file and directory paths in the configuration file in absolute paths.

#===================================

A) Input File name related parameters:

#======================================

	A.1) ValidPairs: 
		Valid pairs generated by HiC-pro pipeline (either simple text format, or can be in gzipped format). 
		**** Note: Mandatory parameter.

	A.2) Interval: 
		File depicting the bins of the interaction matrix. 
		Size of an interval depends on the bin size. 
		Individual bins are also assigned a distinct number. 
		By default, this file name ends with a suffix '_abs.bed'.

	A.3) Matrix: File listing the number of interactions (contacts) among the bins listed in the "Interval" file.

	#=========================	
	User may leave the fields A.2) and A.3) blank, if wishes to not pre-compute these matrices. In such a case, FitHiChIP computes the matrix from the parameter A.1.
	Otherwise, user can run the following command to create an interaction matrix from the given validpairs file (assumed to be in .txt.gz format):
		zcat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 
		Where,
			$InpValidPairsFile: input valid pairs file (in .txt.gz format)
			$MatrixBuildExec: executable of the name "Build_Matrix" placed under the folder "scripts" within the HiC-pro installation directory.	
			$BIN_SIZE: bin size (integer) in base pairs. For instance, 5000 means 5 Kb bin size.
			$ChrSizeFile: File containing the size of individual chromosomes, with respect to the reference genome (such as hg19). This file should be pre-downloaded from UCSC genome browser. Chromosome size file for the reference chromosome hg19 is already placed within the test data folder.
			$OutPrefix: output prefix string, including the output directory path, used in the names of output interval and contact matrix files.
	#=========================

	A.4) PeakFile: Reference ChIP-seq peak file (recommended to use standard ENCODE peaks). User may use pre-computed ChIP-seq peaks, or pre-computed HiChIP peaks (procedure to generate HiChIP peaks is described in the next section).
	
		**** Mandatory parameter.
		**** Note: user may employ MACS2, or a latest package hichipper (Lareau et al 2018) to compute the ChIP-seq / HiChIP peaks, and use them (bed formatted) as an input to FitHiChIP.
	
	
	A.5) OutDir: Output directory which would contain all the results. Default: present working directory.

	A.6) HiCProBasedir: Installation directory of HiC-pro package. (example: /custom_path/HiCPro/HiC-Pro_2.9.0/) 

		*** Mandatory parameter.
	
		*** Using this base path, following executables are tracked:
		$HiCProBasedir/scripts/build_matrix
			(to compute the inetraction matrix from the valid pairs file)
		$HiCProBasedir/scripts/ice
			(to compute the ICE bias, if ICE bias correction is enabled)
	

#===================================

B) Reference genome related parameters:

#======================================

	B.1) ChrSizeFile: File having the reference chromosome size. 

		Default: chrom_hg19.sizes (already provided within the TestData folder).

	B.2) RefFasta: Fasta formatted sequence of the reference chromosome (such as hg19, mm9). 

		***** Optional parameter. Default value is blank (empty string).
		***** If provided by the user (like the file hg19.fa), used to compute the mappability and GC content information (instructions for downloading this file are already mentioned above).

	B.3) MappabilityFile: Reference mappability file (according to the reference genome). 
		
		***** Optional parameter. Default value is blank (empty string).
		***** If provided by the user, used to compute the mappability information. In such a case, user needs to download this file from the site: http://hgdownload.cse.ucsc.edu/goldenPath/
		*** Should be provided in the bedgraph format. 
		*** For bigWig file, user should use the following UCSC utility for converting to bedGraph format: BigWigToBedgraph inp.bw out.bedGraph

	B.4) REFragFile: File containing the restriction fragment information, with respect to the reference genome and the restriction site. 
	
		**** Optional parameter. Default value is blank (empty string).
		**** If provided, the file format (tab delimited) is:
			chr     interval_start  	interval_end
	
		The MboI restriction fragment file "MboI_hg19_RE_Fragments.bed" (most commonly used restriction fragment in various HiChIP pipelines) is provided as an example file, in the TestData folder. 
		**** For other restriction fragment files, please refer to the HiC-Pro manual (https://github.com/nservant/HiC-Pro) to know about their generation.
	
	B.5) GCSize: Size of the window upstream and downstream the restriction site, for computing the GC content. 
		Default = 200 (as per the specification in the package HiCNorm [PMID: 23023982])
		**** This value is not used unless reference fasta sequence, mappability file, and the restriction fragment file is not provided.

	B.6) MappSize: Size of the window upstream and downstream the restriction site, to calculate the mappability information. 
		Default = 500 (as per the specification of HiCNorm [PMID: 23023982])
		**** This value is not used unless reference fasta sequence, mappability file, and the restriction fragment file is not provided.

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

	D.2) UseP2PBackgrnd: Can be 0 or 1. Applicable only for peak to all interactions (IntType = 3).

		A value of 1 means that peak to peak background will be used for contact probability estimation. Also termed as the stringent background of FitHiChIP. Here, contacts between the peak segments are only used for background probability estimation. This model is thus highly stringent, and reports very highly significant loops.
		
		Value of 0 means that all the interactions (both peak to peak and peak to non-peak) are used for background estimation. This indicates lower stringency of significance estimation.
		
		**** We recommend using the value 1. 

#===================================

E) Bias correction related parameters:

#======================================

	E.1) BiasCorrection: Indicates if the bias correction is enabled (1) or not (0). Default 1 (we recommend as MUST).

		User should enable the bias correction, since without bias correction, many false positive interactions are reported.

	E.2) BiasType: Can be 1 or 2. This parameter signifies the type of bias correction used. 

		A value of 1 means that coverage specific bias is used. 
		Value of 2 means that ICE specific bias is computed. 
		Experiments show that 1 produces higher performance (default).


#===================================

F) Merging adjacent interactions:

#======================================

	F.1) MergeInt: Has the value 0 or 1. 
		If 1, means that adjacent interactions are merged via connected component modeling (refer to Manuscript). 
		If 0, no merging step is performed. 
		*** We recommend setting this value as 1.

#===================================

G) Miscellaneous parameters:

#======================================

	G.1) PREFIX: Prefix string used before any output file name. Default = "FitHiChIP"

	G.2) Draw: A binary variable (1/0). If 1, draws various statistical and analysis plots. Default 0.

	G.3) TimeProf: A binary variable (1/0). If 1, logs the time elapsed during individual steps. Default 0.

	G.4) OverWrite: A binary variable (1/0). If 1, overwrites existing FitHiChIP output files. Otherwise, if 0, skips re-computation of existing outputs.


Describing the Output
---------------------

***** Execution of FitHiChIP pipeline generates an HTML file "Summary_results_FitHiChIP.html" within the output directory "OutDir". 
This file lists all output files generated from the FitHiChIP pipeline, according to the given input parameters. 


***** For convenience of the users, here we describe all the files and folders existing within the specified output directory "OutDir":

1) Parameters.txt: a file listing the parameters used for current execution.

2) TimingProfile.txt: If the parameter (G.3) is 1, this file lists the time taken for individual steps.

3) HiCPro_Matrix_BinSize"BINSIZE": If the parameters (A.2) and (A.3) are empty, this folder contains the interaction matrix files generated from the input valid pairs file. The interaction matrix files are "MatrixHiCPro_abs.bed" and "MatrixHiCPro.matrix".

	"*.interactions.initial.bed": CIS interactions with respect to this contact matrix file (interacting regions, and their contact counts) are also dumped separately in this file.
	
	L_"l"_U"u": folder name according to the distance thresholds "l" and "h" specified in the parameters (C.3) and (C.4).

		"*.cis.interactions.DistThr.bed": File within the above mentioned directory contains the CIS interactions filtered within the specified distance thresholds l and u.

3) "NormFeatures": Folder containing various statistical features computed for individual genomic intervals / bins (with respect to the specified bin size).

	User may look at the file "*.AllBin_CompleteFeat.bed" specifically. The file is of the following format:
		
		Columns 1 to 3: chromosome interval (according to the specified bin size).
		
		Column 4: Coverage of this interval (number of reads mapped onto it)
		
		Column 5: denotes whether this interval belongs to a peak segment (1) or non-peak segment (0).
		Here, the peak information is with respect to the parameter (A.4) of the configuration file.

		Column 6: Bias value of this interval. Computed either by coverage specific bias correction, or 
		the ICE specific bias normalization, according to the value of the parameter (E.2).
		
		Columns 7 and 8 denote the mappability and GC contents, respectively, for the current interval. 
		
		***** These values are zero if no mappability, restriction fragment and reference genome fasta sequence are provided.

4) The folders FitHiChIP_"INTTYPE"_b"BINSIZE"_L"l"_U"u" contain the results for different types of interactions, subject to the specified bin size, and the distance thresholds specified.
	
		INTTYPE: interaction type: Have one of the following values:

			ALL2ALL (when the parameter (C.1) = 4)

			Peak2Peak (when the parameter (C.1) = 1)

			Peak2NonPeak (when the parameter (C.1) = 2)

			Peak2ALL (when the parameter (C.1) = 3)

		Within this folder, a directory structure of the following naming pattern is created:

			P2PBckgr_[UseP2PBackgrnd]/[BiasType]

			Where: 

				[UseP2PBackgrnd]: Value of the parameter (D.2)

				[BiasType]: String of either "Coverage_Bias" or "ICE_Bias", 
				according to the value of parameter (E.2).

			4.1) Within this generic directory structure, following files and folders are created:

				4.1.1) Interactions.bed: All contacts among the selected pairs of segments 
				(depending upon the value of INTTYPE). 
	
				4.1.2) Interactions.sortedGenDist.bed: Above contacts are sorted by 
				increasing genomic distance. Used as an input to the custom FitHiC implementation.
	
				4.1.3) FitHiC: If the option (E.1) = 0, this folder contains the significant interactions 
				generated by Naive FitHiC method (without any bias correction technique).

				4.1.4) FitHiC_BiasCorr: this folder is created when (E.1) = 1.

					Directory structure within the folders (4.1.3) and (4.1.4) are similar. 
					Following files and folders reside within this directory:

					A) Bin_Info.log: File listing the bins employed in FitHiChIP (equal occupancy bins)
						
					B) configuration.txt: list of configuration parameters for FitHiChIP
						
					C) PREFIX.interactions_FitHiC.bed: Lists all interactions along with their 
					probabilities, p and q values, computed using the FitHiChIP method.

					D) PREFIX.interactions_FitHiC_Q"$qval".bed: Lists the significant interactions 
					according to the q-value ($qval) significance threshold specified as an 
					input argument. For example, if q value threshold is 0.01, the file name 
					becomes PREFIX.interactions_FitHiC_Q0.01.bed.
					
					Note: for a custom FDR threshold (say 0.05), user may use the following awk script 
					to filter the file (C) based on its last column.
					
					awk '$NF<0.05' PREFIX.interactions_FitHiC.bed > FitHiC_out_Q0.05.bed

					E) PREFIX.interactions_FitHiC_Q"$qval"_WashU.bed: Significant interactions 
					(as provided in file (D)) are converted to WashU epigenome browser 
					compatible format, and saved in this text file. The file has the following 
					format:
					
						chr1,start1,end1	chr2,start2,end2	-log10(q-value)
						
					F) PREFIX.interactions_FitHiC_Q"$qval"_Dist_CC.pdf: Plots the contact count 
					vs genomic distance, corresponding to the significant loops.
					
					F) Merge_Nearby_Interactions: Directory containing the output of 
					merging adjacent interactions (if the parameter (F.1) = 1).
					
						F.1) PREFIX.interactions_FitHiC_Q"$qval"_MergeNearContacts.bed: 
						file containing the loops after merging the adjacent loops in (D)
						
						F.2) PREFIX.interactions_FitHiC_Q"$qval"_MergeNearContacts_WashU.bed: 
						WashU epigenome browser compatible loops.
						
						F.3) PREFIX.interactions_FitHiC_Q"$qval"_MergeNearContacts_Dist_CC.pdf:
						Plot between contact count and interaction distance, corresponding 
						to the significant and merged loops.




Various Utility functions 
-------------------------

In addition to the main HiChIP data processing pipeline discussed above, FitHiChIP also includes various utility functions and related scripts, which are described in this section.

Utility 1 - generating contact matrices of different resolution (bin size)
---------------------------------------------------------------------------

The script "ContactMatDiffRes.sh" within the folder "Imp_Scripts" generates contact matrix of desired bin size given a valid pairs file generated from the HiC-pro pipeline. Parameters associated with this script is:

	-V ValidPairsFile 	Name of the valid pairs file generated from the HiC-pro pipeline
	-B BinSize         	Size of the bin (in bytes: target resolution): default 5000 (5 Kb)
        -C ChrSizeFile          File containing the size of the chromosomes corresponding to the 
				reference genome (such as hg19.chrom.sizes)
	-D OutDir               Output directory which will contain the contact matrix of the target resolution
	-H HiCProDir            Directory containing the HiC-pro installed package.
	
	
	A folder "HiCPro_Matrix_BinSize${BinSize}" will be created under the specified output directory. The bin specific 
	interval file name is MatrixHiCPro_abs.bed and the corresponding contact matrix file name is MatrixHiCPro.matrix

Utility 2 - inferring peaks from HiChIP data (for use in the HiChIP pipeline)
-----------------------------------------------------------------------------

The script "PeakInferHiChIP.sh" within the folder "Imp_Scripts" is used to infer peaks from HiChIP data. The script can be 
used if HiC-pro pipeline is already executed on a given pair of reads (such as .fastq.gz read pairs). The script uses MACS2 
for inferring the peaks. Parameters associated with this script are as follows:

	-H  HiCProDir		Directory containing the reads generated by HiC-pro pipeline. Within this directory, 
				files of the formats .ValidPairs, .DEPairs, .REPairs, and .SCPairs are present, 
				which corresponds to different categories of reads generated by the HiC-pro pipeline.
	-D  OutDir		Directory to contain the output set of peaks. Default: current directory
	-R  refGenomeStr        Reference genome specific string used for MACS2. Default is 'hs' for human 
				chromosome. For mouse, specify 'mm'.
	-M  MACS2ParamStr	String depicting the parameters for MACS2. Default: "--nomodel --extsize 147 -q 0.01"
	-L  ReadLength		Length of reads for the HiC-pro generated reads. Default 75
	
	The script uses all of the DE, SC, RE and validpairs reads generated from the HiC-pro pipeline to infer peaks. 
	The folder "MACS2_ExtSize" within the specified output directory contains the MACS2 generated peaks.
	
	
Utility 3 - merging a set of ChIP-seq alignment files
-----------------------------------------------------

The script "MergeBAMInferPeak.sh" within the folder "Imp_Scripts" processes a set of ChIP-seq alignment files, to generate a merged alignment, infer peak (using MACS2) from the generated alignment, and also determines the coverage of individual bins (according to a specified fixed bin size) with respect to individual input ChIP-seq alignment files. Output of this script is used as an input to the differential analysis module (mentioned next). Parameters associated with this script are as follows:

	-I  InpFile                     A list of ChIP-seq BAM files. Multiple BAM files are to be mentioned in 
					the format: -I bamfile1.bam -I bamfile2.bam -I bamfile3.bam and so on
	-D  OutBaseDir                  Directory containing the output set of peaks. Default: current directory
	-R  refGenome                   Reference genome string used for MACS2. Default is 'hs' for human 
					chromosome. For mouse, specify 'mm'
	-M  MACS2ParamStr           	String depicting the parameters for MACS2. Default: "--nomodel --extsize 147 -q 0.01"
	-C  ChrSizeFile                 Filename containing the chromosome size information for the reference genome
	-b  BinSize                     BinSize in base pair. Used to compute the ChIP-seq coverage for individual 
					input BAM files. Default = 5000 (means 5 Kb)

	Output files generated from this script:
		1) ${OutBaseDir}/merged_input.bam: merged ChIP-seq alignment file
		2) ${OutBaseDir}/MACS2_Out/out.macs2_peaks.narrowPeak: Peak (derived by MACS2) 
		file generated by merged ChIP-seq alignment.
		3) ${OutBaseDir}/ChIPCoverage1.bed, ${OutBaseDir}/ChIPCoverage2.bed, ..., ${OutBaseDir}/ChIPCoverageN.bed, 
		where N is the number of input ChIP-seq alignment files. Each of these output files contain the 
		coverage of corresponding input alignment file (with respect to the specified bin size).
	
Utility 4 - Differential analysis of HiChIP loops - two categories, each with multiple replicates
-------------------------------------------------------------------------------------------------

The R script "DiffAnalysisHiChIP.r" within the folder "Imp_Scripts" is the code for differential analysis of two categories of 
FitHiChIP loops, having M and N replicates, respectively. These two categories may correspond to two different cell types 
or cell lines. However, bin size employed for all of these interactions should be identical.

Parameters associated with this script are as follows:

	--AllLoopList	FILELIST	Comma or colon separated list of FitHiChIP loop files from all the replicates 
					of both input categories, where individual interaction files are of the 
					format: PREFIX.interactions_FitHiC.bed. That is, all significant and 
					non-significant interactions of individual categories and replicates are 
					used as the input (without applying any FDR threshold). MANDATORY PARAMETER.
					
	--FDRThrLoop	Threshold	FDR threshold used for FitHiChIP loops. Default = 0.01 (same used in 
					FitHiChIP implementation)

	--OutDir	DirName		Base Output directory under which all results will be stored. 
					MANDATORY PARAMETER.
					
	--UseRawCC	0/1		If 1, uses the raw contact count for differential analysis. Else, uses both 
					raw and expected contact count (obtained from the bias regression model) values 
					for differential analysis. Default = 0 (Recommended).
					
	--PeakFileCat1	filename	ChIP-seq peak file obtained by merging ChIP-seq replicates of the first category. 
					User may use the script "MergeBAMInferPeak.sh" for producing such a file. 
					MANDATORY PARAMETER.
					
	--PeakFileCat2	filename	ChIP-seq peak file obtained by merging ChIP-seq replicates of the second category. 
					User may use the script "MergeBAMInferPeak.sh" for producing such a file. 
					MANDATORY PARAMETER.
					
	--CategoryList	Names		Comma or colon separated list of strings depicting the names of two categories. 
					User may provide the names of two cell lines or cell types. 
					Default: Category1:Category2
					
	--ReplicaCount	Counts		Comma or colon separated list of two integer values - 
					the number of replicates belonging to individual 
					input categories. Default: 1:1 meaning that both categories have single 
					replicate.
					
	--ReplicaLabels1 Names		Comma or colon separated list of the label of replicates for the first 
					category. Default: R1:R2:R3 etc (as per the replicate counts)
					
	--ReplicaLabels2 Names		Comma or colon separated list of the label of replicates for the second 
					category. Default: R1:R2:R3 etc (as per the replicate counts)
					
	--BinCoverageList FILELIST	List of files storing the ChIP-seq coverage of individual alignment files (of 
					different replicates of both categories) according to the specified bin size. 
					User may use the script "MergeBAMInferPeak.sh" for producing these files. 
					MANDATORY PARAMETER.
					
	--InpTSSFile	FileName	Name of file containing TSS information of the reference genome. 
					Please check the description mentioned below to know how to 
					generate this file. MANDATORY PARAMETER.
					
	--GeneExprFileList TwoFileNames	Comma or colon separated list of two files, storing the gene expression 
					values for individual categories. THIS PARAMETER is OPTIONAL; if provided, 
					gene expression of the differential loops are also analyzed.
					
	--GeneNameColList Counts	Comma or colon separated list of two integer values - the column numbers of 
					individual gene expression files (mentioned in the parameter 
					--GeneExprFileList) storing the gene expression values. 
					THIS PARAMETER is OPTIONAL; required only if the parameter 
					--GeneExprFileList is provided.
					
	--ExprValColList  Counts	Comma or colon separated list of two integer values - the column numbers of 
					individual gene expression files (mentioned in the parameter 
					--GeneExprFileList) storing the name of corresponding genes. 
					THIS PARAMETER is OPTIONAL; required only if the parameter 
					--GeneExprFileList is provided.
					
	--FoldChangeThr	  integer	Fold change threshold employed in EdgeR (log2 scale). Default = 3, meaning 
					log2(3) is used as the fold change threshold.
					
	--FDRThr 	threshold	FDR threshold for determining the significance of EdgeR. Default is 0.05, 
					means that loops with FDR < 0.05, and fold change >= log2(FoldChangeThr) 
					would be considered as differential.
	
	--bcv		threshold	If EdgeR is used with single samples (replica count = 1 for any of the 
					categories), this value is the square-root-dispersion.  
					For datasets arising from well-controlled experiments are 0.4 for human data, 
					0.1 for data on genetically identical model organisms or 
					0.01 for technical replicates. For details, see the edgeR manual. 
					By default, the value is set as 0.4. Used only if a category contains 
					single replicate.


**** Example of differential loop calling using the above mentioned parameters, and with respect to the test data provided, 
is described in the file "DiffAnalysisHiChIP_script.sh" placed within the folder "Imp_Scripts".

**** This file also contains sample scripts to generate the TSS containing file (used in the parameter --InpTSSFile) 
with respect to different reference genomes (such as hg19, mm9, mm10, etc.) Users are requested to check the script 
in details for understanding the parameters.


Sample logs from the console, corresponding to the TestData
--------------------------------------------------------------

Upon executing the given sample script, test data, and the configuration files, texts logged in the console are provided in the file "SampleConsoleOutput.txt". 

User can compare this text file with his / her console output, to check whether the installation has been successful.

Note: Only the file names and relative paths should vary from the given sample output and the user's, according to the target environment.


Citation
---------

If you are using FitHiChIP, please cite:

Sourya Bhattacharyya, Vivek Chandra, Pandurangan Vijayanand, and Ferhat Ay, 

"FitHiChIP: Identification of significant chromatin contacts from HiChIP data", 

preprint at https://www.biorxiv.org/content/early/2018/09/10/412833, 

DOI: https://doi.org/10.1101/412833


Contact
--------

For any query, please contact:

Sourya Bhattacharyya (sourya@lji.org)

Ferhat Ay (ferhatay@lji.org)