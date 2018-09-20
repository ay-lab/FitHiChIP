#!/bin/bash

#===============
# A stand alone executable which takes input of a contact matrix
# (either raw contact matrix or normalized contact matrix)
# generated from the HiC-pro pipeline
# along with a desired (bin) resolution
# and creates a contact matrix of the specified bin (resolution)

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============


#===============
# sample execution command for this script:
# ./ContactMatDiffRes.sh -V rawdata.validpairs -B 10000 -C hg19.chrom.sizes -D /home/SampleDir/ -H /home/HiC_Pro_2.9.0/
#===============

usage(){
cat << EOF

Options:
   	-V  ValidPairsFile		Name of the valid pairs file generated from the HiC-pro pipeline
   	-B 	BinSize 			Size of the bin (in bytes: target resolution): default 5000 (5 Kb)
   	-C  ChrSizeFile			File containing the size of the chromosomes (such as hg19.chrom.sizes)
   	-D 	OutDir 				Output directory which will contain the contact matrix of the target resolution
   	-H 	HiCProDir 			Directory containing the executable of HiC-pro pipeline
EOF
}

#==============
# default parameters
InpValidPairsFile=""
BIN_SIZE=5000
ChrSizeFile=""
OutDir=`pwd`'/'
HiCProBasedir=""
#==============

while getopts "V:B:C:D:H:" opt;
do
	case "$opt" in
		V) InpValidPairsFile=$OPTARG;;
		B) BIN_SIZE=$OPTARG;;
		C) ChrSizeFile=$OPTARG;;
		D) OutDir=$OPTARG;;
		H) HiCProBasedir=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

#===================
# verify the input parameters
#===================

if [[ -z $InpValidPairsFile ]]; then
	echo 'User did not provide any valid pairs file. - exit !!'
	exit 1
fi

if [[ -z $ChrSizeFile ]]; then
	echo 'User did not provide any file containing chromosome size information. - exit !!'
	exit 1
fi

if [[ -z $HiCProBasedir ]]; then
	echo 'User should provide the Base directory of HiC-pro installation path - exit !!'
	exit 1
fi

#===================
# create the output directory
#===================
mkdir -p $OutDir

# generate the matrix of Hi-C interactions (ALL)
# using HiC-pro pipeline
HiCProMatrixDir=$OutDir'/HiCPro_Matrix_BinSize'$BIN_SIZE
mkdir -p $HiCProMatrixDir

# executable of matrix building is to be obtained from the HiC-pro base directory
# provided as the input
MatrixBuildExecSet=( $(find $HiCProBasedir -type f -name 'build_matrix') )
len=${#MatrixBuildExecSet[@]}
# echo 'len: '$len
if [[ $len == 0 ]]; then
	echo 'Did not find HiC-pro package installation and the utility for matrix generation - quit !!'
	exit 1
fi
idx=`expr $len - 1`
# echo 'idx: '$idx
MatrixBuildExec=${MatrixBuildExecSet[$idx]}
echo -e '\n *** MatrixBuildExec: '$MatrixBuildExec

echo '*** Computing HiC-pro matrices from the input valid pairs file'

# This directory and prefix is used to denote the generated matrices
# using the HiC pro routine
OutPrefix=$HiCProMatrixDir'/MatrixHiCPro'

if [[ ! -f $OutPrefix'_abs.bed' || ! -f $OutPrefix'.matrix' ]]; then
	# check the extension of input valid pairs file
	# and extract accordingly
	if [[ $InpValidPairsFile == *.gz ]]; then
		echo "***** HiC-pro input valid pairs file in gzipped format"
		zcat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper'  
	else
		echo "***** HiC-pro input valid pairs file in simple text format"
		cat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 
	fi
fi

# now assign the matrix names to the designated variables
InpBinIntervalFile=$OutPrefix'_abs.bed'
InpMatrixFile=$OutPrefix'.matrix'

#============================
# this directory stores the features and associated data
#============================
FeatureDir=$OutDir'/NormFeatures_BinSize'$BIN_SIZE
mkdir -p $FeatureDir

#=======================================
# ICE specific bias is also computed
# using the ICE routine from the HiC-pro pipeline
# =======================================

AllFeatureDir=$FeatureDir'/ICE_Bias'
mkdir -p $AllFeatureDir

# here ICE specific bias is used
# compute the bias vector from the HiC-pro contact matrix
ICEExec=$HiCProBasedir'/scripts/ice'
echo -e '\n *** ICE computation Executable: '$ICEExec
echo '*** Computing ICE based bias vector from the HiC-pro contact matrix'

# run the ICE executable and store the normalized contact matrix 
# and the bias vector
NormContactMatrixFile=$AllFeatureDir'/ICE.norm.Contact.Matrix'
BiasVecFile=$NormContactMatrixFile'.biases'

if [[ ! -f $NormContactMatrixFile || ! -f $BiasVecFile ]]; then
	$ICEExec $InpMatrixFile --results_filename $NormContactMatrixFile --output-bias $BiasVecFile
	# replace the NANs of the derived bias vector with zero
	sed -i 's/nan/0/g' $BiasVecFile
fi


