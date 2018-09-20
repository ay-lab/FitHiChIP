#!/bin/bash

# Rscript version
RscriptExec=`which Rscript`
echo 'Rscript version installed in the system : '$RscriptExec

# current working directory
currworkdir=`pwd`

# directory of this script
currscriptdir=`dirname $0`
echo 'current directory containing this script: '$currscriptdir

# go to the current script containing directory
cd $currscriptdir

# source code (R) of the differential analysis
DiffAnalysisCodeExec=$currscriptdir'/DiffAnalysisHiChIP.r'
echo 'R Code of differential analysis: '$DiffAnalysisCodeExec

# directory containing reference dataset
DataDir='../../data/'

# directory containing data specific to the differential analysis module
DiffAnalysisMainDir=$DataDir'DiffAnalysis/'

# directory containing all the results
ResultMainDir='../../results/'

#======================
# this script is required if user has to create a .gtf file containing promoters (TSS sites) 
# corresponding to the reference chromosome
# as we are providing the TSS sites (.gtf file) corresponding to the reference genome hg19
# we are disabling the script
# for other reference genomes, user should follow the subsequent instructions 
# and enable this script
#======================

# first download the reference genome annotation in .gtf format
# from GENCODE site:
# hg19: https://www.gencodegenes.org/releases/19.html
# hg38: https://www.gencodegenes.org/releases/28.html
# mm9: https://www.gencodegenes.org/mouse_releases/1.html
# mm10: https://www.gencodegenes.org/mouse_releases/18.html

# in each case, download the Comprehensive gene annotation (CHR) in .gtf format
# suppose the corresponding file name is "inp.gtf"

# now process this downloaded file to extract the TSS and 
# gene name information
# using the following scripts (awk + sed)
# suppose the output filename is "out.gtf"


# start of script

if [ 1 == 0 ]; then

awk -F'[;\t ]' '{if (substr($1,1,1)!="#")  {print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$16"\t"$22"\t"$13}}' inp.gtf | awk -F['\t'] '($6=="\"protein_coding\"")' - | awk -F['\t'] '{if ($4=="+") {print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7} else {print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}}' - | awk -F['\t'] '!seen[$1, $4, $5, $6, $7]++' - | awk -F['\t'] '{print $1"\t"$2"\t"$3"\t"NR"\t"$4"\t"$5"\t"$6"\t"$7}' - > out.gtf

sed -i 1i"Chr\\tStart\\tEnd\\tTSSIdx\\tStrand\\tGeneID\\tGeneType\\tGeneName" out.gtf

fi

# end of script






#======================
# this script is required if user has to merge the ChIP-seq alignments from individual replicates
# of a particular category
# as we are providing a sample merged alignment with the peaks as well
# we are disabling this script
# for any other dataset, user should enable this script
# and fill with the custom input and output file names
#=======================

# with respect to individual categories, first generate the peak file from the merged ChIP-seq alignments
# (of individual replicates)

# call the script "MergeBAMInferPeak.sh"
# with the ChIP seq alignment files of the first category (here Category1_R1.bam to Category1_R3.bam)
# output directory name: $DiffAnalysisMainDir'Category1_Peaks'
# MACS2 reference genome: 'hs'
# MACS2 parameters: '--nomodel --extsize 147 -q 0.01'
# output peak file name will be: $DiffAnalysisMainDir'Category1_Peaks/MACS2_Out/out.macs2_peaks.narrowPeak'
# bin size is mentioned as 5000 (5 Kb)
# ChIP-seq coverage corresponding to individual BAM files (and the specified bin size)
# will be placed in the folder $DiffAnalysisMainDir'Category1_Peaks'
# in the files ChIPCoverage1.bed to ChIPCoverage3.bed

if [ 1 == 0 ]; then

sh $currscriptdir'/MergeBAMInferPeak.sh' -I $DiffAnalysisMainDir'Category1_R1.bam' -I $DiffAnalysisMainDir'Category1_R2.bam' -I $DiffAnalysisMainDir'Category1_R3.bam' -D $DiffAnalysisMainDir'Category1_Peaks' -R 'hs' -M '--nomodel --extsize 147 -q 0.01' -C $DataDir'chrom_hg19.sizes' -b 5000

fi

# similarly call the script "MergeBAMInferPeak.sh"
# with the ChIP seq alignment files of the second category
# output directory name: $DiffAnalysisMainDir'Category2_Peaks'
# MACS2 reference genome: 'hs'
# MACS2 parameters: '--nomodel --extsize 147 -q 0.01'
# output peak file name will be: $DiffAnalysisMainDir'Category2_Peaks/MACS2_Out/out.macs2_peaks.narrowPeak'
# bin size is mentioned as 5000 (5 Kb)
# ChIP-seq coverage corresponding to individual BAM files (and the specified bin size)
# will be placed in the folder $DiffAnalysisMainDir'Category2_Peaks'
# in the files ChIPCoverage1.bed to ChIPCoverage3.bed

if [ 1 == 0 ]; then

sh $currscriptdir'/MergeBAMInferPeak.sh' -I $DiffAnalysisMainDir'Category2_R1.bam' -I $DiffAnalysisMainDir'Category2_R2.bam' -I $DiffAnalysisMainDir'Category2_R3.bam' -D $DiffAnalysisMainDir'Category2_Peaks' -R 'hs' -M '--nomodel --extsize 147 -q 0.01' -C $DataDir'chrom_hg19.sizes' -b 5000

fi





#==============================
# this is the main invoking of the differential analysis code
# all the parameters are with respect to the sample test data provided
# if user has to run this module with any other dataset
# then modify the corresponding file names accordingly

# Note: if multiple arguments are provided for a single parameter (such as --AllLoopList), those arguments are 
# separated by a colon (:) operator
#==============================

$RscriptExec $DiffAnalysisCodeExec --AllLoopList $DiffAnalysisMainDir'Category1_R1_FitHiC_All.bed':$DiffAnalysisMainDir'Category1_R2_FitHiC_All.bed':$DiffAnalysisMainDir'Category1_R3_FitHiC_All.bed':$DiffAnalysisMainDir'Category2_R1_FitHiC_All.bed':$DiffAnalysisMainDir'Category2_R2_FitHiC_All.bed':$DiffAnalysisMainDir'Category2_R3_FitHiC_All.bed' --FDRThrLoop 0.01 --UseRawCC 0 --OutDir $ResultMainDir'/DiffAnalysis_HiChIP_Output' --CategoryList 'CellType1':'CellType2' --ReplicaCount 3:3 --ReplicaLabels1 'Cat1_R1':'Cat1_R2':'Cat1_R3' --ReplicaLabels2 'Cat2_R1':'Cat2_R2':'Cat2_R3' --PeakFileCat1 $DiffAnalysisMainDir'Category1_Peaks/MACS2_Out/out.macs2_peaks.narrowPeak' --PeakFileCat2 $DiffAnalysisMainDir'Category2_Peaks/MACS2_Out/out.macs2_peaks.narrowPeak' --BinCoverageList $DiffAnalysisMainDir'Category1_Peaks/ChIPCoverage1.bed':$DiffAnalysisMainDir'Category1_Peaks/ChIPCoverage2.bed':$DiffAnalysisMainDir'Category1_Peaks/ChIPCoverage3.bed':$DiffAnalysisMainDir'Category2_Peaks/ChIPCoverage1.bed':$DiffAnalysisMainDir'Category2_Peaks/ChIPCoverage2.bed':$DiffAnalysisMainDir'Category2_Peaks/ChIPCoverage3.bed' --InpTSSFile $DiffAnalysisMainDir'hg19_Genes_Only_Protein_Coding_TSS.gtf' 




# now go back to the original working directory
cd $currworkdir















