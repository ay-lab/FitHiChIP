#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem-per-cpu=40G
# source ~/.bashrc
source ~/.bash_profile

cd /home/sourya/proj/GitHub_Codes/FitHiChIP_HiCPro

# code executable for differential loop analysis 
currworkdir=`pwd`
currscriptdir=$currworkdir	#`dirname $0`
cd $currscriptdir
echo 'current directory containing this script: '$currscriptdir

# Rscript version
RscriptExec=`which Rscript`
echo 'Rscript version installed in the system : '$RscriptExec

# source code (R) of the differential analysis
DiffAnalysisCodeExec=$currscriptdir'/Imp_Scripts/DiffAnalysisHiChIP.r'
echo 'R Code of differential analysis: '$DiffAnalysisCodeExec

BaseOutDir=$currscriptdir'/TestData/DiffLoopData/Results_DiffLoops_NEW_Feb15_2024'
mkdir -p $BaseOutDir

## input interaction files
cat1_file1=$currscriptdir'/TestData/DiffLoopData/Category1_Replicate1_FitHiChIP.interactions_FitHiC.bed.gz'
cat1_file2=$currscriptdir'/TestData/DiffLoopData/Category1_Replicate2_FitHiChIP.interactions_FitHiC.bed.gz'
cat2_file1=$currscriptdir'/TestData/DiffLoopData/Category2_Replicate1_FitHiChIP.interactions_FitHiC.bed.gz'
cat2_file2=$currscriptdir'/TestData/DiffLoopData/Category2_Replicate2_FitHiChIP.interactions_FitHiC.bed.gz'

## input ChIP files
chip1=$currscriptdir'/TestData/DiffLoopData/Category1_ChIPSeq_Align.bam'
chip2=$currscriptdir'/TestData/DiffLoopData/Category2_ChIPSeq_Align.bam'

$RscriptExec ${DiffAnalysisCodeExec} \
    --AllLoopList ${cat1_file1}:${cat1_file2}:${cat2_file1}:${cat2_file2} \
    --ChrSizeFile $currscriptdir'/TestData/chrom_hg19.sizes' \
    --FDRThr 0.01 \
    --CovThr 25 \
    --ChIPAlignFileList ${chip1}:${chip2} \
    --OutDir $BaseOutDir \
    --CategoryList 'C1':'C2' \
    --ReplicaCount 2:2 \
    --ReplicaLabels1 "R1":"R2" \
    --ReplicaLabels2 "R1":"R2" \
    --FoldChangeThr 2 \
    --DiffFDRThr 0.05 \
    --BackgroundFDRThr 0.1

# revert to the old directory
cd $currworkdir

