FitHiChIP
----------------

Developers: Sourya Bhattacharyya, Ferhat Ay

La Jolla Institute for Allergy and Immunology

La Jolla, CA 92037, USA

**************************

FitHiChIP analyzes HiChIP / PLAC-seq data and derives the statistical significant CIS interactions, given a specific distance range (default = 20 Kb to 2 Mb: can be user defined).

Prerequisites
-------------

FitHiChIP requires the following packages / libraries to be installed:

1) HiC-pro (obtained from "https://github.com/nservant/HiC-Pro") (>= version 2.9.0)
	
2) Python2 (>= version 2.7) along with libraries: 1) OptionParser (from the library optparse), 2) gzip, 3) networkx (https://networkx.github.io/)

3) R (>= 3.4.3) with the libraries: 1) optparse, 2) ggplot2, 3) splines, 4) fdrtool, 5) parallel, 6) tools, 7) GenomicRanges, 8) dplyr, 9) EdgeR (https://bioconductor.org/packages/release/bioc/html/edgeR.html) (if differential analysis of FitHiChIP loops is to be performed).	
		
4) bedtools (http://bedtools.readthedocs.io/en/latest/) (>= version 2.26.0)
	
5) samtools (http://www.htslib.org/) (>= version 1.6), 

	==== (New from version 6.0) Also install htslib (version >= 1.6) specifically for the utilities "bgzip" and "tabix"

6) macs2 for peak calling (https://github.com/taoliu/MACS)


**** User can check if the above packages are installed by typing each package name in the command line terminal and checking if the package exists in the system.


**** FitHiChIP is tested in linux environment, and requires bash for executing the main script.


Execution
---------

"sample_script.sh": sample script to execute FitHiChIP, having the following command.

	bash FitHiChIP_HiCPro.sh -C configuration_file_name

where, 
		
	"FitHiChIP_HiCPro.sh": main executable code of FitHiChIP pipeline.
		
Four configuration files are provided: 

	configfile_BiasCorrection_CoverageBias: FitHiChIP(L) & FitHiChIP(L+M) with coverage bias correction.
	configfile_BiasCorrection_ICEBias: FitHiChIP(L) & FitHiChIP(L+M) with ICE bias correction.
	configfile_P2P_BiasCorrection_CoverageBias: FitHiChIP(S) & FitHiChIP(S+M) with coverage bias correction.
	configfile_P2P_BiasCorrection_ICEBias: FitHiChIP(S) & FitHiChIP(S+M) with ICE bias correction.
	
User may look into any one of the configuration files or the detailed README file for understanding the configuration parameters. Some of the parameters are kept blank in the configuration files, meaning those parameters are optional.

User should provide the following mandatory parameters:

	1) ValidPairs: HiC-pro pipeline generated valid pairs file (gzipped or uncompressed tab delimited text format)
	2) PeakFile: Either reference ChIP-seq peak file, or peaks inferred from HiChIP data (for details, please see the detailed README file, utilities section)
	3) OutDir: Output directory containing all the results
	4) HiCProBasedir: HiC-pro installation directory
	5) ChrSizeFile: size of individual chromosomes corresponding to the reference genome (one example file "chrom_hg19.sizes" is provided within the TestData folder. Can be downloaded from UCSC genome browser).
	6) IntType: we recommend the value 3 (peak to all foreground / output loops, meaning both peak to peak and peak to non-peak loops will be reported)
	7) BINSIZE: Size of the bins (in bp). We have used 5000 (means 5 Kb) as a bin size.
	8)  LowDistThr: Lower distance threshold of interaction between two intervals (CIS). Default: 20000 (indicates 20 Kb). Interactions below this distance threshold will not be considered for statistical significance.
	9) UppDistThr: Upper distance threshold of interaction between two intervals (CIS). Default: 2000000 (indicates 2 Mb). Interactions above this distance threshold will not be considered for statistical significance.
	10) QVALUE: Minimum FDR (q-value) cutoff for detecting significant interactions. Default: 0.01
	11) UseP2PBackgrnd: if 0, FitHiChIP(L) or loose background is employed; else if 1, FitHiChIP(S) or stringent background is employed. This parameter is however applicable only for peak to all interactions (IntType = 3). In such a case, we recommend using 1.
	12) BiasCorrection: This value MUST be set as 1
	13) BiasType: Can be 1 (coverage bias) or 2 (ICE bias). We recommend 1.
	14) MergeInt: Boolean variable. If 1, loops with adjacent bins in both sides are merged together via connected component modeling. We recommed 1.  
	15) OverWrite: Overwrites the existing input files. User may set this value 1 if previous execution was interrupted or a new version of code is run. 

Output
-------

***** Execution of FitHiChIP pipeline generates an HTML file "Summary_results_FitHiChIP.html" within the output directory "OutDir". This file lists all important output files generated from the FitHiChIP pipeline, according to the given input parameters. 

Detailed description of FitHiChIP output files are provided in the detailed README file.

Utility functions 
------------------

FitHiCHIP offers various utility programs such as generating contact matrices of different bin size from the given valid pairs file, inferring peaks directly from the HiChIP data (for using them in the FitHiChIP pipeline instead of using the reference ChIP-seq peaks), and most importantly, differential analysis of FitHiChIP loops between two categories (such as two different cell types) with one or more replicates, such that the differential loops are solely caused by 3D (loop) changes and not by 1D (ChIP) changes. Detailed description of these utility programs are provided in the detailed version of README file.

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
