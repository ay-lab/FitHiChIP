#!/bin/bash

#===============
# A stand alone executable of FitHiChIP module
# used to process the HiC-pro pipeline output (allvalidpairs.txt file)
# to generate the contact matrices and 
# the statistically significant interactions

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============

#===============
# Sourya - parse a configuration file
# with the appropriate parameters
#===============

# first read the configuration file name

usage(){
cat << EOF

Options:
   	-C  ConfigFile		Name of the configuration file storing the parameters of FitHiChIP.
EOF
}

# default values
biaslowthr=0.2
biashighthr=5
BeginBiasFilter=0
EndBiasFilter=0

while getopts "C:" opt;
do
	case "$opt" in
		C) ConfigFile=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

#======================
# default values of the parameters
#======================

# input files
InpValidPairsFile=""
InpBinIntervalFile=""
InpMatrixFile=""
PeakFILE=""

PREFIX='FitHiChIP'

# reference genome
RefGENOME='hg19'

# size of the chromosome that is to be provided
ChrSizeFile=""

# reference chromosome fasta file
RefFastaFile=""

# reference chromosome based mappability file
# (may be downloaded from the site  http://hgdownload.cse.ucsc.edu/goldenPath/)
MappabilityFile=""

# restriction fragment file compatible with the reference chromosome
REFragFile=""

# window size used to compute GC content
GCContentWindowSize=200

# window size used to compute the mappability
MappabilityWindowSize=500

# 5 Kb resolution
BIN_SIZE=5000

QVALUE=1e-2	# 0.01

# default value of output directory is the present working directory
OutDir=`pwd`'/'

# lower distance threshold for two cis interactions
LowDistThres=20000	# 20 Kb

# upper distance threshold for two cis interactions
UppDistThres=2000000 # 2 Mb

# number of bins employed for FitHiC
NBins=200

# default value of plotting analysis figures
DrawFig=0

# binning method for FitHiC technique
# 1 for equal occupancy bin (default)
FitHiCBinMethod=1

# option to note down the timing information
TimeProf=0

#=========================
# bias correction related parameters
#=========================

# default lower cutoff of bias value
biaslowthr=0.2

# default higher cutoff of bias value
biashighthr=5

# boolean variable for pre-filtering the interactions according to the bias value
BeginBiasFilter=0

# boolean variable for probability adjustment of the interactions according to the bias value
EndBiasFilter=0

# Merging interactions which are near to each other
MergeInteraction=1

#==============================
# read the configuration file and store various parameters
#==============================

# separator used in the config file
IFS="="
while read -r name value
do
	param=$name
	paramval=${value//\"/}
	if [[ -n $param ]]; then
		if [[ $param != \#* ]]; then
			#echo "Content of $param is $paramval"
			if [ $param == "ValidPairs" ]; then
				InpValidPairsFile=$paramval
			fi
			if [ $param == "Interval" ]; then
				InpBinIntervalFile=$paramval
			fi
			if [ $param == "Matrix" ]; then
				InpMatrixFile=$paramval
			fi
			if [ $param == "PeakFile" ]; then
				PeakFILE=$paramval
			fi
			if [ $param == "OutDir" ]; then
				if [[ ! -z $paramval ]]; then
					OutDir=$paramval
				fi
			fi
			if [ $param == "RefGenome" ]; then
				if [[ ! -z $paramval ]]; then
					RefGENOME=$paramval
				fi
			fi
			if [ $param == "ChrSizeFile" ]; then
				ChrSizeFile=$paramval
			fi
			if [ $param == "MappabilityFile" ]; then
				MappabilityFile=$paramval
			fi
			if [ $param == "RefFasta" ]; then
				RefFastaFile=$paramval
			fi
			if [ $param == "REFragFile" ]; then
				REFragFile=$paramval
			fi
			if [ $param == "GCSize" ]; then
				if [[ ! -z $paramval ]]; then
					GCContentWindowSize=$paramval
				fi
			fi
			if [ $param == "MappSize" ]; then
				if [[ ! -z $paramval ]]; then
					MappabilityWindowSize=$paramval
				fi
			fi
			if [ $param == "BINSIZE" ]; then
				if [[ ! -z $paramval ]]; then
					BIN_SIZE=$paramval
				fi
			fi
			if [ $param == "LowDistThr" ]; then
				if [[ ! -z $paramval ]]; then
					LowDistThres=$paramval
				fi
			fi
			if [ $param == "UppDistThr" ]; then
				if [[ ! -z $paramval ]]; then
					UppDistThres=$paramval
				fi
			fi
			if [ $param == "QVALUE" ]; then
				if [[ ! -z $paramval ]]; then
					QVALUE=$paramval
				fi
			fi
			if [ $param == "NBins" ]; then
				if [[ ! -z $paramval ]]; then
					NBins=$paramval
				fi
			fi
			if [ $param == "HiCProBasedir" ]; then
				HiCProBasedir=$paramval
			fi
			if [ $param == "PREFIX" ]; then
				if [[ ! -z $paramval ]]; then
					PREFIX=$paramval
				fi
			fi
			if [ $param == "Draw" ]; then
				if [[ ! -z $paramval ]]; then
					DrawFig=$paramval
				fi
			fi
			if [ $param == "TimeProf" ]; then
				if [[ ! -z $paramval ]]; then
					TimeProf=$paramval
				fi
			fi
			if [ $param == "BeginBiasFilter" ]; then
				if [[ ! -z $paramval ]]; then
					BeginBiasFilter=$paramval
				fi
			fi
			if [ $param == "EndBiasFilter" ]; then
				if [[ ! -z $paramval ]]; then
					EndBiasFilter=$paramval
				fi
			fi
			if [ $param == "biaslowthr" ]; then
				if [[ ! -z $paramval ]]; then
					biaslowthr=$paramval
				fi
			fi			
			if [ $param == "biashighthr" ]; then
				if [[ ! -z $paramval ]]; then
					biashighthr=$paramval
				fi
			fi
			if [ $param == "MergeInt" ]; then
				if [[ ! -z $paramval ]]; then
					MergeInteraction=$paramval
				fi
			fi
		fi
	fi
done < $ConfigFile

#===================
# verify the input parameters
#===================

if [[ -z $InpValidPairsFile ]]; then
	if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then
		echo 'User did not provide any valid pairs file. So, user needs to provide both of the interval and matrix files. But at least one of them is missing - exit !!'
		exit 1
	fi
fi

if [[ -z $PeakFILE ]]; then
	echo 'User should provide a reference peak detection file to compute the interactions involving peak segments - exit !!'
	exit 1
fi

if [[ -z $RefFastaFile ]]; then
	echo 'User should provide reference chromosome fasta file - quit !!'
	exit 1
fi

if [[ -z $REFragFile ]]; then
	echo 'User should provide a restriction fragment file - quit !!'
	exit 1
fi

if [[ -z $HiCProBasedir ]]; then
	if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then
		echo 'Input matrices are not provided and the Base directory of HiC-pro installation path is also not provided - exit !!'
		exit 1
	fi
fi

if [[ -z $ChrSizeFile ]]; then
	echo 'Chromosome size file is not specified - exit !!'
	exit 1
fi

if [[ -z $MappabilityFile ]]; then
	echo 'Reference mappability file is not specified - exit !!'
	exit 1
fi

#==============================
# here check if the configuration file has relative path names as the input
# in such a case, convert the relative path names (with respect to the location of the configuration file itself)
# in the absolute file

# directory of the configuration file
ConfigFileDir=$(dirname "${ConfigFile}")

# first go to the configuration file directory
cd $ConfigFileDir

if [[ ! -z $InpValidPairsFile ]]; then
	if [[ "${InpValidPairsFile:0:1}" != / && "${InpValidPairsFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		InpValidPairsFile="$(cd $(dirname $InpValidPairsFile); pwd)/$(basename $InpValidPairsFile)"
		echo 'Absolute converted path: InpValidPairsFile: '$InpValidPairsFile
	fi
fi

if [[ ! -z $InpBinIntervalFile ]]; then
	if [[ "${InpBinIntervalFile:0:1}" != / && "${InpBinIntervalFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		InpBinIntervalFile="$(cd $(dirname $InpBinIntervalFile); pwd)/$(basename $InpBinIntervalFile)"
		echo 'Absolute converted path: InpBinIntervalFile: '$InpBinIntervalFile
	fi
fi

if [[ ! -z $InpMatrixFile ]]; then
	if [[ "${InpMatrixFile:0:1}" != / && "${InpMatrixFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		InpMatrixFile="$(cd $(dirname $InpMatrixFile); pwd)/$(basename $InpMatrixFile)"
		echo 'Absolute converted path: InpMatrixFile: '$InpMatrixFile
	fi
fi

if [[ ! -z $PeakFILE ]]; then
	if [[ "${PeakFILE:0:1}" != / && "${PeakFILE:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		PeakFILE="$(cd $(dirname $PeakFILE); pwd)/$(basename $PeakFILE)"
		echo 'Absolute converted path: PeakFILE: '$PeakFILE
	fi
fi

if [[ ! -z $RefFastaFile ]]; then
	if [[ "${RefFastaFile:0:1}" != / && "${RefFastaFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		RefFastaFile="$(cd $(dirname $RefFastaFile); pwd)/$(basename $RefFastaFile)"
		echo 'Absolute converted path: RefFastaFile: '$RefFastaFile
	fi
fi

if [[ ! -z $REFragFile ]]; then
	if [[ "${REFragFile:0:1}" != / && "${REFragFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		REFragFile="$(cd $(dirname $REFragFile); pwd)/$(basename $REFragFile)"
		echo 'Absolute converted path: REFragFile: '$REFragFile
	fi
fi

if [[ ! -z $ChrSizeFile ]]; then
	if [[ "${ChrSizeFile:0:1}" != / && "${ChrSizeFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		ChrSizeFile="$(cd $(dirname $ChrSizeFile); pwd)/$(basename $ChrSizeFile)"
		echo 'Absolute converted path: ChrSizeFile: '$ChrSizeFile
	fi
fi

if [[ ! -z $MappabilityFile ]]; then
	if [[ "${MappabilityFile:0:1}" != / && "${MappabilityFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		MappabilityFile="$(cd $(dirname $MappabilityFile); pwd)/$(basename $MappabilityFile)"
		echo 'Absolute converted path: MappabilityFile: '$MappabilityFile
	fi
fi

if [[ ! -z $OutDir ]]; then
	if [[ "${OutDir:0:1}" != / && "${OutDir:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		OutDir="$(cd $(dirname $OutDir); pwd)/$(basename $OutDir)"
		echo 'Absolute converted path: OutDir: '$OutDir
	fi
fi

if [[ ! -z $HiCProBasedir ]]; then
	if [[ "${HiCProBasedir:0:1}" != / && "${HiCProBasedir:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		HiCProBasedir="$(cd $(dirname $HiCProBasedir); pwd)/$(basename $HiCProBasedir)"
		echo 'Absolute converted path: HiCProBasedir: '$HiCProBasedir
	fi
fi

# revert to the old directory
cd -

#===================
# create the output directory
#===================
mkdir -p $OutDir

#============================
# print the parameters and values
#============================
ConfFile=$OutDir/Parameters.txt

echo "InpValidPairsFile: $InpValidPairsFile " > $ConfFile
echo "InpBinIntervalFile: $InpBinIntervalFile " >> $ConfFile
echo "InpMatrixFile: $InpMatrixFile " >> $ConfFile
echo "PeakFILE: $PeakFILE " >> $ConfFile
echo "OutDir: $OutDir " >> $ConfFile
echo "RefGENOME: $RefGENOME " >> $ConfFile
echo "ChrSizeFile: $ChrSizeFile " >> $ConfFile
echo "MappabilityFile: $MappabilityFile " >> $ConfFile
echo "RefFastaFile: $RefFastaFile " >> $ConfFile
echo "REFragFile: $REFragFile " >> $ConfFile
echo "GCContentWindowSize: $GCContentWindowSize " >> $ConfFile
echo "MappabilityWindowSize: $MappabilityWindowSize " >> $ConfFile
echo "BIN_SIZE: $BIN_SIZE " >> $ConfFile
echo "LowDistThr: $LowDistThres " >> $ConfFile
echo "UppDistThr: $UppDistThres " >> $ConfFile
echo "QVALUE: $QVALUE " >> $ConfFile
echo "NBins: $NBins " >> $ConfFile
echo "HiCProBasedir: $HiCProBasedir " >> $ConfFile
echo "PREFIX: $PREFIX " >> $ConfFile
echo "DrawFig: $DrawFig " >> $ConfFile
echo "Timeprof: $TimeProf " >> $ConfFile
echo "Bias pre-filtering: $BeginBiasFilter " >> $ConfFile
echo "Prob Adjust due to bias: $EndBiasFilter " >> $ConfFile
echo "Bias lower cutoff: $biaslowthr " >> $ConfFile
echo "Bias higher cutoff: $biashighthr " >> $ConfFile
echo "Merging nearby interactions: $MergeInteraction " >> $ConfFile

#=======================================
# generate a file which will contain the timing profile
if [ $TimeProf == 1 ]; then
	OutTimeFile=$OutDir'/TimingProfile.txt'
	echo " ================ Time profiling =========== " > $OutTimeFile
	start=$(date +%s.%N)
fi

#==============================
# important - sourya
# first change the current working directory to the directory containing this script
# it is useful when the script is invoked from a separate directory
#==============================
currworkdir=`pwd`
currscriptdir=`dirname $0`
cd $currscriptdir

# generate the matrix of Hi-C interactions (ALL)
# using HiC-pro pipeline
HiCProMatrixDir=$OutDir'/HiCPro_Matrix_BinSize'$BIN_SIZE
mkdir -p $HiCProMatrixDir

#=================
# if the matrices are not provided and the validpairs text file is provided
# then compute the interaction matrices using the HiC-pro utility
#=================
if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then

	# this is an executable which builds matrix from the input valid pairs file
	# that's why we require the HiC pro executable directory as a command line option

	# executable of matrix building is to be obtained from the HiC-pro base directory
	# provided as the input
	# MatrixBuildExec=$HiCProBasedir'/scripts/build_matrix'
	MatrixBuildExecSet=( $(find $HiCProBasedir -type f -name 'build_matrix') )
	if [[ ${#MatrixBuildExecSet[@]} == 0 ]]; then
		echo 'Did not find HiC-pro package installation and the utility for matrix generation - quit !!'
		exit 1
	fi
	MatrixBuildExec=${MatrixBuildExecSet[0]}
	echo -e '\n *** MatrixBuildExec: '$MatrixBuildExec

	echo '*** Computing HiC-pro matrices from the input valid pairs file'

	# This directory and prefix is used to denote the generated matrices
	# using the HiC pro routine
	OutPrefix=$HiCProMatrixDir'/MatrixHiCPro'

	if [ ! -f $OutPrefix'_abs.bed' ]; then
		# check the extension of input valid pairs file
		# and extract accordingly
		if [[ $InpValidPairsFile == *.gz ]]; then
			zcat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 
		else
			cat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 
		fi
	fi

	# now assign the matrix names to the designated variables
	InpBinIntervalFile=$OutPrefix'_abs.bed'
	InpMatrixFile=$OutPrefix'.matrix'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for computing the interaction matrix using HiC-Pro build_matrix utility: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi

fi

#=======================
# Now generate the list of interactions from the HiC-pro matrix data
# ALL to ALL interactions
# Both cis and trans interactions are considered (with respect to the given bin size)
# No distance threshold based filtering is used
# Interaction format:
# chr1	start1	end1	chr2	start2	end2	cc
#=======================
Interaction_Initial_File=$HiCProMatrixDir/$PREFIX.interactions.initial.bed
if [ ! -f $Interaction_Initial_File ]; then
	Rscript ./src/InteractionHicPro.r $InpBinIntervalFile $InpMatrixFile $Interaction_Initial_File
	
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting CIS interactions: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#=======================
# generate filtered cis interactions 
# with respect to distance thresholds
# ALL to ALL interactions
#=======================
# create a directory for individual distance thresholds
InteractionThrDir=$HiCProMatrixDir'/L_'$LowDistThres'_U'$UppDistThres
mkdir -p $InteractionThrDir
Interaction_File=$InteractionThrDir/$PREFIX.cis.interactions.DistThr.bed

if [ ! -f $Interaction_File ]; then
	awk -v l="$LowDistThres" -v u="$UppDistThres" 'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || ($1==$4 && abs($2-$5)>=l && abs($2-$5)<=u)) {print $0}}' $Interaction_Initial_File > $Interaction_File

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting CIS interactions within distance thresholds: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#============================
# this directory stores the features and associated data
#============================
FeatureDir=$OutDir'/NormFeatures'
mkdir -p $FeatureDir

# file storing the RE fragments, mappability and GC content together
REFragMappGCFile=$FeatureDir'/REFrag_Mapp_GC_Merged.bed'

if [ ! -f $REFragMappGCFile ]; then

	#============================
	# generating the mappability information
	# first divide each RE fragment interval
	# from two ends
	# the offset size = $MappabilityWindowSize
	# upstream and downstream 
	# we do not cross the RE fragment boundary - safe for length overflow
	#============================
	echo 'Creating the fragment end (w.r.t window size) file -- to compute the mappability information!!'
	MappOffsetCutBedFile=$FeatureDir'/Temp_Fragment_Mapp_'$MappabilityWindowSize'bp.bed'
	if [ ! -f $MappOffsetCutBedFile ]; then
		awk -v s=$MappabilityWindowSize 'function max(x,y) {return x>y?x:y}; function min(x,y) {return x<y?x:y}; {printf "%s\t%d\t%d\n%s\t%d\t%d\n", $1, $2, min($2+s,$3), $1, max($3-s, $2), $3}' $REFragFile > $MappOffsetCutBedFile

		if [ $TimeProf == 1 ]; then
			duration=$(echo "$(date +%s.%N) - $start" | bc)
			echo " ++++ Time (in seconds) for computing the fragment file to compute the mappability: $duration" >> $OutTimeFile
			start=$(date +%s.%N)
		fi		
	fi

	#============================
	# computation of the mappability scores
	# the map utility of bedtools function checks the 
	# overlap from the 2nd file to the 1st file
	# the 4th column of the second file is used as the score
	# mean indicates the average score to be used in the final output
	# Note: Also replace the non-number entries with 0
	#============================
	echo 'Creating the mappability file !!'
	MappabilityOutFile=$FeatureDir'/Mappability_RE_Fragments.bed'
	if [ ! -f $MappabilityOutFile ]; then
		bedtools map -a $MappOffsetCutBedFile -b $MappabilityFile -c 4 -o mean | awk '{if ($4=="." || $4=="NA" || $4=="NaN") {$4=0}; print $0}' - > $MappabilityOutFile

		if [ $TimeProf == 1 ]; then
			duration=$(echo "$(date +%s.%N) - $start" | bc)
			echo " ++++ Time (in seconds) for computing the mappability: $duration" >> $OutTimeFile
			start=$(date +%s.%N)
		fi		
	fi

	#============================
	# generating the GC content information
	# first divide each RE fragment interval
	# from two ends
	# the offset size = $GCContentWindowSize
	# upstream and downstream 
	# we do not cross the RE fragment boundary - safe for length overflow
	#============================
	echo 'Creating the fragment end (w.r.t window size) file -- to compute the GC content information!!'
	GCOffsetCutBedFile=$FeatureDir'/Temp_Fragment_GC_'$GCContentWindowSize'bp.bed'
	if [ ! -f $GCOffsetCutBedFile ]; then
		awk -v s=$GCContentWindowSize 'function max(x,y) {return x>y?x:y}; function min(x,y) {return x<y?x:y}; {printf "%s\t%d\t%d\n%s\t%d\t%d\n", $1, $2, min($2+s,$3), $1, max($3-s, $2), $3}' $REFragFile > $GCOffsetCutBedFile

		if [ $TimeProf == 1 ]; then
			duration=$(echo "$(date +%s.%N) - $start" | bc)
			echo " ++++ Time (in seconds) for computing the fragment file to compute the GC content: $duration" >> $OutTimeFile
			start=$(date +%s.%N)
		fi		
	fi

	#============================
	# generation of %GC from the reference fasta file
	# bedtools suite is used
	# Note: Also replace the non-number entries with 0
	#============================
	echo 'Creating the GC content file !!'
	GCOutFile=$FeatureDir'/GC_Content_RE_Fragments.bed'
	if [ ! -f $GCOutFile ]; then
		nucBed -fi $RefFastaFile -bed $GCOffsetCutBedFile | awk '{if ($4=="." || $4=="NA" || $4=="NaN") {$4=0}; print $0}' - > $GCOutFile

		if [ $TimeProf == 1 ]; then
			duration=$(echo "$(date +%s.%N) - $start" | bc)
			echo " ++++ Time (in seconds) for computing the GC content: $duration" >> $OutTimeFile
			start=$(date +%s.%N)
		fi		
	fi

	#============================
	# examine the 4th field of generated mappability file
	# and 5th field of the generated GC content (%) file
	# contents of consecutive lines need to be averaged and dumped 
	#============================
	# the mappability fragment file does not have any header information
	# average of the 4th field
	Temp_Mapp_File=$FeatureDir'/Mappability_Dump.bed'
	awk '{sum+=$4} NR%2==0 {print sum/2; sum=0}' $MappabilityOutFile > $Temp_Mapp_File

	# the GC content fragment file has header
	# average of the 5th field
	Temp_GC_File=$FeatureDir'/GC_Dump.bed'
	awk '{if (NR>1) {sum+=$5}}; {if (NR%2!=0 && NR>1) {print sum/2; sum=0}}' $GCOutFile > $Temp_GC_File

	#============================
	# Now combine the average mappability and GC content values with the 
	# restriction fragments 
	#============================
	Rscript ./src/CombineREFragMappGC.r $REFragFile $Temp_Mapp_File $Temp_GC_File $REFragMappGCFile

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for generating the final mappability and GC content of intervals: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi

	# remove the temporary files
	rm $MappOffsetCutBedFile
	rm $MappabilityOutFile
	rm $GCOffsetCutBedFile
	rm $GCOutFile
	rm $Temp_Mapp_File
	rm $Temp_GC_File
fi

#=================
# From the input valid paired end read file, and the given input bin size parameter
# compute the coverage of individual genome segments (bins)
# the output is a list of chromosome and bins, their coverage, and a boolean indicator 
# whether the segment overlaps with a peak
#=================
CoverageFile=$FeatureDir'/'$PREFIX'.coverage.bed'
if [ ! -f $CoverageFile ]; then
	python ./src/CoverageBin.py -i $InpValidPairsFile -p $PeakFILE -b $BIN_SIZE -o $CoverageFile -c $ChrSizeFile
	echo 'Computed initial coverage information for individual genomic bins'
	
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting coverage of individual bins: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

# =======================================
# Compute the bias for peaks and non-peaks separately 
# Curerntly, the bias values are written in the coverage file itself (overwriting)
# and it uses only the coverage values information
# later, more complex models can be incorporated
# =======================================
CoverageBiasFile=$FeatureDir'/'$PREFIX'.coverage_Bias.bed'
if [ ! -f $CoverageBiasFile ]; then
	Rscript ./src/BiasCalc.r --CoverageFile $CoverageFile --OutFile $CoverageBiasFile
	echo 'Appended bias information for individual genomic bins'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for computing bias of individual bins: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi	
fi

#=================
# templates of filenames used throughout the execution
#=================
InteractionFileName='Interactions.bed'
InteractionSortedDistFileName='Interactions.sortedGenDist.bed'

#=================
# The coverage file (along with the fixed size chromosome bin) and the boolean peak information
# needs to be extended to contain the mappability and GC content features for individual bins
# here the RE fragment file (file b) has mappability information in 4th column 
# and GC content information in 5th column
# mean operation for multiple overlaps is used
# the option "-header" prints the header information of the file (a)
# missing values (no overlap) are indicated by 0

# the coverage file has header information, which needs to be discarded before this function
# after the operation, a file with 8 columns is produced
# chr start end coverage ispeak mappability GCcontent NoCutSites
#=================
AllFeatFile=$FeatureDir'/'$PREFIX'.AllBin_CompleteFeat.bed'
if [ ! -f $AllFeatFile ]; then
	AllFeatFile_temp1=$FeatureDir'/'$PREFIX'.AllBin_CompleteFeat_temp1.bed'
	AllFeatFile_temp2=$FeatureDir'/'$PREFIX'.AllBin_CompleteFeat_temp2.bed'
	awk 'NR>1' $CoverageBiasFile | bedtools map -a /dev/stdin -b $REFragMappGCFile -c 4 -o mean -null '0' > $AllFeatFile_temp1
	bedtools map -a $AllFeatFile_temp1 -b $REFragMappGCFile -c 5 -o mean -null '0' > $AllFeatFile_temp2
	bedtools map -a $AllFeatFile_temp2 -b $REFragMappGCFile -c 4 -o count -null '0' > $AllFeatFile
	rm $AllFeatFile_temp1
	rm $AllFeatFile_temp2

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for computing bin specific features: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi		
fi

#=================
# now we plot various features for individual genomic bins
#=================
if [ $DrawFig == 1 ]; then
	Rscript ./Analysis/PlotGenomeBins.r --GenomeBinFile $AllFeatFile --OutDir $FeatureDir'/Plots'
	echo 'Plotted the distribution of coverage values among peak and non peak segments'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for plotting bin specific features: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi			
fi

#============================
# create all to all interaction file
# with the features like read depth, mappability, GC content, 
# and number of RE sites
#============================
DirALLtoALL=$OutDir'/FitHiChIP_ALL2ALL_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
mkdir -p $DirALLtoALL
IntFileALLtoALL=$DirALLtoALL'/'$InteractionFileName

if [ ! -f $IntFileALLtoALL ]; then
	Rscript ./src/Significance_Features.r -I $Interaction_File -E $AllFeatFile -O $IntFileALLtoALL
fi

# derive the contact count column
cccol=`cat $Interaction_File | tail -n 1 | awk '{print NF}' -`
echo 'Contact count col: '$cccol

# derive the number of columns in the interaction file with normalization 
# related features
totcol=`cat $IntFileALLtoALL | tail -n 1 | awk '{print NF}' -`
echo 'Total number of columns for the complete feature interactions: '$totcol

if [ $TimeProf == 1 ]; then
	duration=$(echo "$(date +%s.%N) - $start" | bc)
	echo " ++++ Time (in seconds) for computing pairwise interactions among all segments: $duration" >> $OutTimeFile
	start=$(date +%s.%N)
fi		

#===================
# using the interaction file among all binned intervals
# associated with the normalization features
# plot the variation among different features
#===================
if [ $DrawFig == 1 ]; then
	Rscript ./Analysis/InteractionPlots.r --IntFile $IntFileALLtoALL --OutDir $OutDir'/Plots_Norm' --MappThr 0.5 --GCThr 0.2 --cccol $cccol
	echo 'Plotted the distribution of normalization features among peak and non peak segments'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for plotting normalization related features for different types of interactions: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#===================
# Create the interaction files for other types of interactions
# peak to peak, peak to non peak, and peak to all
#===================
DirPeaktoPeak=$OutDir'/FitHiChIP_Peak2Peak_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
DirPeaktoNonPeak=$OutDir'/FitHiChIP_Peak2NonPeak_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
DirPeaktoALL=$OutDir'/FitHiChIP_Peak2ALL_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
mkdir -p $DirPeaktoPeak
mkdir -p $DirPeaktoNonPeak
mkdir -p $DirPeaktoALL

IntFilePeaktoPeak=$DirPeaktoPeak'/'$InteractionFileName
IntFilePeaktoNonPeak=$DirPeaktoNonPeak'/'$InteractionFileName
IntFilePeaktoALL=$DirPeaktoALL'/'$InteractionFileName

if [ ! -f $IntFilePeaktoPeak ]; then
	# peak to peak interactions 
	# 9th and 15th fields are 1
	awk '(NR==1) || ($9==1 && $15==1)' $IntFileALLtoALL > $IntFilePeaktoPeak
fi

if [ ! -f $IntFilePeaktoNonPeak ]; then
	# peak to non peak interactions
	# 9th field is 1, but 15th field is 0
	awk '(NR==1) || ($9==1 && $15==0)' $IntFileALLtoALL > $IntFilePeaktoNonPeak
fi

if [ ! -f $IntFilePeaktoALL ]; then
	# peak to all interactions
	# just check if 9th field is 1
	awk '(NR==1) || ($9==1)' $IntFileALLtoALL > $IntFilePeaktoALL
fi

if [ $TimeProf == 1 ]; then
	duration=$(echo "$(date +%s.%N) - $start" | bc)
	echo " ++++ Time (in seconds) for assigning different types of interactions: $duration" >> $OutTimeFile
	start=$(date +%s.%N)
fi

#==============================
# navigate through individual types of interactions (corresponding folders)
# and apply FitHiC for individual interaction types
#==============================
# $DirALLtoALL is commented for the moment - sourya
for dirname in $DirPeaktoPeak $DirPeaktoNonPeak $DirPeaktoALL $DirALLtoALL; do
	
	echo 'Processing the directory: '$dirname

	#==============
	# first create interaction files with sorted genomic distance
	#==============
	CurrIntFile=$dirname'/'$InteractionFileName
	CurrIntFileSortDist=$dirname'/'$InteractionSortedDistFileName

	coln=`expr $totcol + 1`
	if [ ! -s $CurrIntFileSortDist ]; then
		awk -v OFS='\t' 'function abs(v) {return v < 0 ? -v : v} {print $0"\t"abs($5-$2)}' $CurrIntFile | sort -k${coln},${coln}n -k${cccol},${cccol}nr | cut -f1-${totcol} - > $CurrIntFileSortDist	
	fi

	echo 'Created sorted genomic distance based interaction file'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for sorting the interactions (according to genomic distance): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi

	#==============
	# now apply FitHiC on the sorted gene distance based interaction matrix
	#==============
	GenFitHiCDir=$dirname'/FitHiC_EqOccBin'
	if [[ $BeginBiasFilter == 0 && $EndBiasFilter == 0 ]]; then
		BiasCorr=0
	else		
		GenFitHiCDir=$GenFitHiCDir'_BiasCorr_'$biaslowthr'_'$biashighthr'_b'$BeginBiasFilter'_e'$EndBiasFilter
		BiasCorr=1
	fi
	mkdir -p $GenFitHiCDir

	FitHiC_Pass1_outfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1.bed
	FitHiC_Pass1_Filtfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER.bed
	FitHiC_Pass1_Filt_PeakCountfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER.Peakcount.bed
	FitHiC_Pass1_LogQ_file=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_WashU.bed
	FitHiC_Pass1_PeakCCDistr_Text=$GenFitHiCDir/$PREFIX.anchorPeakCCDistr.bed
	FitHiC_Pass1_Filt_MergedIntfile=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_MERGED.bed
	FitHiC_Pass1_Filt_MergedInt_LogQ_file=$GenFitHiCDir/$PREFIX.interactions_FitHiCPass1_FILTER_MERGED_WashU.bed

	# Modeling the statistical significance by FitHiC - main function
	if [ ! -f $FitHiC_Pass1_outfile ]; then
		Rscript ./src/Interaction.r --InpFile $CurrIntFileSortDist --OutFile $FitHiC_Pass1_outfile --Norm $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter
	fi

	echo 'Applied FitHiC'

	# Filter the interaction file with respect to significance (Q value < $QVALUE)
	# also print the header line
	if [ ! -f $FitHiC_Pass1_Filtfile ]; then
		awk -v q="$QVALUE" '{if ((NR==1) || ($NF < q && $NF > 0)) {print $0}}' $FitHiC_Pass1_outfile > $FitHiC_Pass1_Filtfile
	fi

	echo 'Extracted significant interactions from FitHiC'

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for applying FitHiC (significant interactions): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi		

	# no of significant interactions (FitHiC)
	nsigFitHiC=`cat $FitHiC_Pass1_Filtfile | wc -l`

	# Check the no of significant contacts associated with each short peak segment
	if [[ $dirname != $DirALLtoALL ]]; then
		if [ ! -f $FitHiC_Pass1_Filt_PeakCountfile ]; then
			# At least 10 significant interactions are required (empirical threshold)
			# considering the 1st line as header
			if [[ $nsigFitHiC -gt 11 ]]; then
				# skip the 1st header line
				awk 'NR>1' $FitHiC_Pass1_Filtfile | cut -f1-3,${cccol} | sort -k1,1 -k2,2n -k3,3n | awk -v OFS='\t' '{a[$1" "$2" "$3]+=$4}END{for (i in a){print i,a[i]}}' - > $FitHiC_Pass1_Filt_PeakCountfile
			else
				echo 'number of significant interactions for spline distribution < 10 - skip the peak count distribution function'
			fi
		fi
	fi

	if [ $DrawFig == 1 ]; then
		# the R file takes the spline fitted interaction file (without q-value based filtering)
		# and returns the contact count distribution for two different sets of interactions
		# separated by the Q value threshold of 0.01
		# check for non empty interactions file
		if [[ $nsigFitHiC -gt 1 ]]; then
			Rscript ./Analysis/result_summary.r $FitHiC_Pass1_outfile $cccol $QVALUE
		else
			echo 'Number of significant spline interaction <= 1 - no result summary'
		fi
	fi

	# the filtered interaction (with respect to the spline) file is used to create a session file
	# for applying in WashU epigenome browser
	# for that, a special file containing only the interacting chromosome intervals 
	# and the log of Q value is created
	if [ ! -s $FitHiC_Pass1_LogQ_file ]; then
		# check for non empty interactions file
		if [[ $nsigFitHiC -gt 2 ]]; then
			awk '{if (NR > 1) {print $1","$2","$3"\t"$4","$5","$6"\t"(-log($NF)/log(10))}}' $FitHiC_Pass1_Filtfile > $FitHiC_Pass1_LogQ_file
		else
			echo 'There is no significant interaction - so no WashU specific session file is created !!'
		fi
	fi

	# if [[ $DrawFig == 1 && $dirname != $DirALLtoALL ]]; then
	# 	res_outdir=$GenFitHiCDir'/Results'
	# 	mkdir -p $res_outdir
	# 	if [ -f $FitHiC_Pass1_Filt_PeakCountfile ]; then
	# 		# Distribution of significant contact counts (FitHiC) for individual peaks
	# 		if [ ! -f $res_outdir/$PREFIX.PeakCCDistr_FiltSpline.pdf ] || [ ! -f $FitHiC_Pass1_PeakCCDistr_Text ]; then
	# 			Rscript ./Analysis/ContactCountDistr.r $PeakFILE $FitHiC_Pass1_Filt_PeakCountfile $res_outdir/$PREFIX.PeakCCDistr_FiltSpline.pdf $FitHiC_Pass1_PeakCCDistr_Text
	# 		fi
	# 	fi
	# fi

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for post processing FitHiC results: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi

	# if merging nearby interactions are enabled
	# then we merge the nearby interactions from the earlier generated significant interactions
	# and also create a washu browser generated compatible file
	if [ $MergeInteraction == 1 ]; then
		if [ ! -f $FitHiC_Pass1_Filt_MergedIntfile ]; then
			python ./src/CombineNearbyInteraction.py --InpFile $FitHiC_Pass1_Filtfile --OutFile $FitHiC_Pass1_Filt_MergedIntfile --headerInp 1 --binsize $BIN_SIZE
		fi
		if [ ! -s $FitHiC_Pass1_Filt_MergedInt_LogQ_file ]; then
			nint=`cat $FitHiC_Pass1_Filt_MergedIntfile | wc -l`
			if [[ $nint -gt 2 ]]; then
				# 9th field stores the Q value
				awk '{if (NR > 1) {print $1","$2","$3"\t"$4","$5","$6"\t"(-log($9)/log(10))}}' $FitHiC_Pass1_Filt_MergedIntfile > $FitHiC_Pass1_Filt_MergedInt_LogQ_file
			else
				echo 'There is no significant interaction - so no WashU specific session file is created !!'
			fi
		fi
		echo 'Merged nearby significant interactions - created washu browser compatible file for these merged interactions!!!'
	fi

done 	# end of directory traversal loop

#============================
# sourya - now go back to the original working directory
#============================
cd $currworkdir
