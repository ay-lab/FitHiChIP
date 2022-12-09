FitHiChIP
----------------

Developers: Sourya Bhattacharyya, Ferhat Ay

La Jolla Institute for Immunology

La Jolla, CA 92037, USA

**************************

FitHiChIP analyzes HiChIP / PLAC-seq data and derives the statistical significant CIS interactions.


A comprehensive documentation of FitHiChIP is provided in 

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


Release notes - Version 11.0 (December 2022)
-----------------------------------------

	1. FitHiChIP now support HiChIP interactions in .hic and .cool / .mcool formats, in addition to the earlier formats.
	2. Updated configuration files corresponding to these new input options.
	3. Updated Docker and Singularity packages.
	4. Differential HiChIP loop calling does not require ChIP-seq alignment files as a mandatory option. If users do not have any ChIP-seq alignment file, they can just proceed with the differential analysis without considering the difference in 1D.
	5. FitHiChIP output loops are now converted to files compatible with WashU, UCSC and IGV epigenome browsers.


For the earlier release notes, please check the file *Release_Notes.txt*


Utility scripts for the manuscript
======================================

	Check the folder *UtilScript* and corresponding README file for the links to various utility scripts used to generate the figures in this manuscript.


Contact
--------

Please use the GitHub issues page for reporting any issues / suggestions (recommended). 

Alternatively, you can e-mail us:

- Sourya Bhattacharyya <sourya@lji.org>
- Ferhat Ay <ferhatay@lji.org>

