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


Data repository
-----------------

All the results in FitHiChIP, like the significant loops, HiChIP peak calling, performance analysis is now available in Zenodo

https://doi.org/10.5281/zenodo.3255048


Major release - version 9.1 (October, 2021)
---------------------------------------------------------

	1. Upgraded merge filtering routine to support Python3, instead of deprecated python2
	2. Added support for circular genome. In such a case, the genomic distance between interacting fragments are computed using the minimum of their linear distance and their circular genome specific distance. The configuration file now has a parameter "CircularGenome" which, if 1, denotes that the reference genome is circular. **Note** User must provide the appropriate chromosome size file in the configuration options.
	3. Updated the genomic distance based filtering of interactions according to the circular genome.
	4. The FitHiChIP output file (*fithic*.bed) has now one extra field, namely "Dist" (20th column) which explicitly mentions the genomic distance between the interacting fragments. Useful for the circular genome.
	5. Updated the HiChIP peak calling routine as well. Previously, all the reads from DE, SC, RE and valid pairs (from HiC-pro output) were required as the input set of reads. Now we've relaxed such constraints. If the user does not provide any DE, RE, or SC reads, but only provides the valid pairs, HiChIP peaks will be estimated from the valid pairs itself.	
	6. Updated README and manual

Release notes corresponding to version 9.0 (January 10, 2021)
----------------------------------------------------------------

	1. Updated singularity, Docker installation, incorporating the latest code and dependencies.
	2. Incorporated simultaneous generation of WashU browser compatible tracks for differential analysis.
	3. Minor warning fix regarding 0 size of input peaks for a given chromosome.
	4. Updated documentation.

Release notes corresponding to version 8.1 (May 10, 2020)
--------------------------------------------------------------

	1. Mandatory to use HiCPro version 2.11.4. This version automatically installs ICE via python package iced.
	2. Minor bug fix in testing the installed packages and versions before running FitHiChIP.


Contact
--------

Please use the GitHub issues page for reporting any issues / suggestions (recommended). Alternatively, you can mail us:

- Sourya Bhattacharyya <sourya@lji.org>
- Ferhat Ay <ferhatay@lji.org>

