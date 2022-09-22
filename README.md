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


Release notes - Version 10.0 (April 2022)
-----------------------------------------

	1. HiC-pro installation directory is now checked only if user provides HiC-pro valid pairs output. If user provides matrix and bin interval files (obtained from another source), he/she does not require to install HiCPro.
	2. HiC-pro version checking is stopped. Although users are requested to use the latest version of HiCPro.
	3. Parallel processing is updated. Instead of using mclapply, we use lapply. Some users experienced halt / crash of parallel processing routine.


For the earlier release notes, please check the file *Release_Notes.txt*


Utility scripts for the manuscript
======================================

	Check the folder *UtilScript* and corresponding README file for the links to various utility scripts used to generate the figures in this manuscript.



Contact
--------

Please use the GitHub issues page for reporting any issues / suggestions (recommended). Alternatively, you can mail us:

- Sourya Bhattacharyya <sourya@lji.org>
- Ferhat Ay <ferhatay@lji.org>

