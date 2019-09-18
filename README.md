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


Release notes corresponding to version 7.0 (29th April 2019) (Major release)
---------------------------------------------------------------------------------------------------

	1. Included support for processing pre-computed set of locus pairs along with their observed contact count. The file should have the following seven fields in tab seperated format: 

		chr1	start1	end1	chr2	start2	end2	contactcount

	Such a file can be computed using any HiC or HiChIP data processing pipeline.

	In the sample configuration files, parameter (A.4) corresponds to this new parameter.

	2. Completely updated methodology of differential HiChIP loop finding using FitHiChIP loops of two categories, each with one or more replicates.

		The script "Differetial_Analysis_Script.sh" contains the updated commands.

		We have also included test data of differential loop finding under the folder "TestData/DiffLoopData"

	3. Modified HiChIP coverage computation for different types of inputs (i.e. either valid pairs file, or matrix + bin files, or pre-computed locus pairs file)

	4. Removed redundant and optional parameters from the command line options, so that users should enter only the essential parameters for executing FitHiChIP. Sample configuration files are updated.


	*** for older releases and corresponding release notes, please check the file "Release_Notes.txt"


Contact
-----------

- Sourya Bhattacharyya <sourya@lji.org>
- Ferhat Ay <ferhatay@lji.org>

