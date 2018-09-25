#!/bin/bash 


currworkdir=`pwd`
currscriptdir=`dirname $0`
cd $currscriptdir
echo 'current directory containing this script: '$currscriptdir


#===============================================
# directory containing all the input data of differential analysis
DiffAnalysisMainDir=$currscriptdir'/data/DiffAnalysis/'

# directory to contain the output results of differential analysis
DiffAnalysisOutDir=$currscriptdir'/results/DiffAnalysis_HiChIP_Output/'

#===============================================

# Rscript version
RscriptExec=`which Rscript`
echo 'Rscript version installed in the system : '$RscriptExec

# source code (R) of the differential analysis
DiffAnalysisCodeExec=$currscriptdir'/Imp_Scripts/DiffAnalysisHiChIP.r'
echo 'R Code of differential analysis: '$DiffAnalysisCodeExec

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

# this out.gtf file would be used in the differential analysis module (as the .gtf file containing protein coding genes)

#===========
# start of script

if [ 1 == 0 ]; then

awk -F'[;\t ]' '{if (substr($1,1,1)!="#")  {print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$16"\t"$22"\t"$13}}' inp.gtf | awk -F['\t'] '($6=="\"protein_coding\"")' - | awk -F['\t'] '{if ($4=="+") {print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7} else {print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}}' - | awk -F['\t'] '!seen[$1, $4, $5, $6, $7]++' - | awk -F['\t'] '{print $1"\t"$2"\t"$3"\t"NR"\t"$4"\t"$5"\t"$6"\t"$7}' - > out.gtf

sed -i 1i"Chr\\tStart\\tEnd\\tTSSIdx\\tStrand\\tGeneID\\tGeneType\\tGeneName" out.gtf

fi

# end of script
#===========

#======================
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

sh $currscriptdir'/MergeBAMInferPeak.sh' -I $DiffAnalysisMainDir'Category1_R1.bam' -I $DiffAnalysisMainDir'Category1_R2.bam' -I $DiffAnalysisMainDir'Category1_R3.bam' -D $DiffAnalysisMainDir'Category1_Peaks' -R 'hs' -M '--nomodel --extsize 147 -q 0.01' -C $DiffAnalysisMainDir'chrom_hg19.sizes' -b 5000

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

sh $currscriptdir'/MergeBAMInferPeak.sh' -I $DiffAnalysisMainDir'Category2_R1.bam' -I $DiffAnalysisMainDir'Category2_R2.bam' -I $DiffAnalysisMainDir'Category2_R3.bam' -D $DiffAnalysisMainDir'Category2_Peaks' -R 'hs' -M '--nomodel --extsize 147 -q 0.01' -C $DiffAnalysisMainDir'chrom_hg19.sizes' -b 5000

fi


#=============================================
# now invoke the differential analysis code
# if multiple arguments are provided for a single parameter (such as --AllLoopList), those arguments are 
# separated by a colon (:) operator

# --AllLoopList: All the loops along with their FitHiC generated significance values (without any FDR thresholding) for different categories and replicates, separated by colon
# --FDRThrLoop: FDR Threshold used to determine loop significance
# --UseRawCC: 0 or 1: If 0, raw contact count is used for EdgeR. Otherwise, raw and estimated (by bias regression model) contact counts are used. Default and recommeded: 0
# --OutDir : output directory containing the differential analysis output
# --CategoryList: comma or colon separated list of names of two categories. Here we have specified the names as CellType1 and CellType2
# --ReplicaCount: comma or colon separated list of integers containing the number of replicates of individual categories
# --ReplicaLabels1: comma or colon separated list of names of replicates of individual categories. Here we have mentioned those names 
# 					as Cat1_R1, Cat1_R2, ...
# --PeakFileCat1: peak file corresponding to the first category (obtained by first merging the CHIP-seq alignments of individual 
# 				  replicates, and then inferring peaks using MACS2)
# --PeakFileCat2: peak file corresponding to the second category 
# --BinCoverageList: comma or colon separated list of filenames, where each file contains the ChIP-seq coverage of 
# 					 individual bins (according to the bin size of FitHiChIP loops), and for individual replicates. 
# --InpTSSFile: List of protein coding genes for the reference genome

#=============================================

$RscriptExec $DiffAnalysisCodeExec --AllLoopList $DiffAnalysisMainDir'Category1_R1_FitHiC_All.bed.gz':$DiffAnalysisMainDir'Category1_R2_FitHiC_All.bed.gz':$DiffAnalysisMainDir'Category1_R3_FitHiC_All.bed.gz':$DiffAnalysisMainDir'Category2_R1_FitHiC_All.bed.gz':$DiffAnalysisMainDir'Category2_R2_FitHiC_All.bed.gz':$DiffAnalysisMainDir'Category2_R3_FitHiC_All.bed.gz' --FDRThrLoop 0.01 --UseRawCC 0 --OutDir $DiffAnalysisOutDir --CategoryList 'CellType1':'CellType2' --ReplicaCount 3:3 --ReplicaLabels1 'Cat1_R1':'Cat1_R2':'Cat1_R3' --ReplicaLabels2 'Cat2_R1':'Cat2_R2':'Cat2_R3' --PeakFileCat1 $DiffAnalysisMainDir'Category1Peaks.narrowPeak.gz' --PeakFileCat2 $DiffAnalysisMainDir'Category2Peaks.narrowPeak.gz' --BinCoverageList $DiffAnalysisMainDir'Category1_R1_ChIPCoverage_5Kb.bed.gz':$DiffAnalysisMainDir'Category1_R2_ChIPCoverage_5Kb.bed.gz':$DiffAnalysisMainDir'Category1_R3_ChIPCoverage_5Kb.bed.gz':$DiffAnalysisMainDir'Category2_R1_ChIPCoverage_5Kb.bed.gz':$DiffAnalysisMainDir'Category2_R2_ChIPCoverage_5Kb.bed.gz':$DiffAnalysisMainDir'Category2_R3_ChIPCoverage_5Kb.bed.gz' --InpTSSFile $DiffAnalysisMainDir'hg19_Genes_Only_Protein_Coding_TSS.gtf.gz' 

#======================
# revert to the old directory
cd $currworkdir

















