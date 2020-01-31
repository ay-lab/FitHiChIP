FitHiChIP
----------------

Developers: Sourya Bhattacharyya, Ferhat Ay

La Jolla Institute for Allergy and Immunology

La Jolla, CA 92037, USA

**************************

FitHiChIP analyzes HiChIP / PLAC-seq data and derives the statistical significant CIS interactions, given a specific distance range (default = 20 Kb to 2 Mb: can be user defined).


Documentation of FitHiChIP is provided in 

https://ay-lab.github.io/FitHiChIP/


Citation
-----------
*FitHiChIP is now published at Nature Communications (<https://www.nature.com/articles/s41467-019-11950-y>)*

If you are using FitHiChIP, please cite:

Sourya Bhattacharyya, Vivek Chandra, Pandurangan Vijayanand, and Ferhat Ay, *Identification of significant chromatin contacts from HiChIP data by FitHiChIP*, Nature Communications, Vol 10, No 4221, 2019, DOI: <https://doi.org/10.1038/s41467-019-11950-y>


Release notes corresponding to version 8.0 (January 30, 2020)
--------------------------------------------------------------

	1. Minor bug fix in differential analysis module
	2. Using ggplot2 for plotting spline fit and regression coefficients, in the significant interaction module. Useful for running in computational cluster environment.

Release notes corresponding to version 7.1 (24th October 2019) 
------------------------------------------------------------------

	Major changes: 

	1. Updated differential analysis module - support processing ChIP-seq coverage files in BEDGraph format (in addition to process ChIP-seq alignment .bam files). Improved code with much lower running time. Improved documentation.
		** Requires additional R libraries to be installed : 1) data.table, 2) plyr
		** For details, please check the installation procedure in the main documentation
	2. Thoroughly tested support of HiC data. User now can provide HiC contact matrices in either HiC-pro based validpairs / matrix format, or simple bed format (peak file is not required) and specify ALL-to-ALL (4) interaction type to obtain HiC data specific significant interactions.
	3. Included code for simulating HiChIP data from input HiC and ChIP-seq (as published in our Nature Communication Paper). Note that this implementation is far from optimal, and was mainly to show the robustness of FitHiChIP. The simulation needs to be much improved.
	3. Lower running time in finding statistically significant interactions, by faster data reading and processing.
	4. Updated nextflow and docker installation steps.

	Minor changes:

	1. Added: differential analysis using gzipped input files.
	2. Removed: dependency of specifying HiC pro installation directory in the configuration file.
	2. Bug fix: differential loops using input files with a subset of chromosomes (even one chromosome).
	3. Bug fix: parsing input parameters - error in checking q-value range (invalid arithmetic operator)
	
	*** for older releases and corresponding release notes, please check the file "Release_Notes.txt"


Contact
--------

Please use the GitHub issues page for reporting any issues / suggestions (recommended). Alternatively, you can mail us:

- Sourya Bhattacharyya <sourya@lji.org>
- Ferhat Ay <ferhatay@lji.org>

