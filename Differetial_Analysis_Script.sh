#!/bin/bash 

# code executable for differential loop analysis 
currworkdir=`pwd`
currscriptdir=`dirname $0`
cd $currscriptdir
echo 'current directory containing this script: '$currscriptdir

# Rscript version
RscriptExec=`which Rscript`
echo 'Rscript version installed in the system : '$RscriptExec

# source code (R) of the differential analysis
DiffAnalysisCodeExec=$currscriptdir'/Imp_Scripts/DiffAnalysisHiChIP.r'
echo 'R Code of differential analysis: '$DiffAnalysisCodeExec

BaseOutDir=$currscriptdir'/TestData/DiffLoopData/Results_DiffLoops'
mkdir -p $BaseOutDir

$RscriptExec ${DiffAnalysisCodeExec} --AllLoopList $currscriptdir'/TestData/DiffLoopData/Category1_Replicate1_FitHiChIP.interactions_FitHiC.bed.gz':$currscriptdir'/TestData/DiffLoopData/Category1_Replicate2_FitHiChIP.interactions_FitHiC.bed.gz':$currscriptdir'/TestData/DiffLoopData/Category2_Replicate1_FitHiChIP.interactions_FitHiC.bed.gz':$currscriptdir'/TestData/DiffLoopData/Category2_Replicate2_FitHiChIP.interactions_FitHiC.bed.gz' --ChrSizeFile $currscriptdir'/TestData/chrom_hg19.sizes' --FDRThr 0.01 --CovThr 25 --ChIPAlignFileList $currscriptdir'/TestData/DiffLoopData/Category1_ChIPSeq_Align.bam':$currscriptdir'/TestData/DiffLoopData/Category2_ChIPSeq_Align.bam' --OutDir $BaseOutDir --CategoryList 'Category1':'Category2' --ReplicaCount 2:2 --ReplicaLabels1 "R1":"R2" --ReplicaLabels2 "R1":"R2" --FoldChangeThr 2 --DiffFDRThr 0.05 

# revert to the old directory
cd $currworkdir

