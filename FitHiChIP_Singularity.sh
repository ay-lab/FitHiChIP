#!/bin/bash

usage(){
cat << EOF
Options:
   	-C  ConfigFile		Name of the configuration file storing the parameters of FitHiChIP.
EOF
}

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
# 3 different ways to provide input to FitHiChIP
# option 1: provide a valid pairs file generated from HiC pro pipeline
InpValidPairsFile=""
# option 2: provide the bin intervakl file and the matrix file generated from HiC pro pipeline
InpBinIntervalFile=""
InpMatrixFile=""
# option 3: provide the set of locus pairs (6 columns storing the interacting bins and the 7th column storing the contact count)
InpInitialInteractionBedFile=""

# reference ChIP-seq / HiChIP peak file
PeakFILE=""

# prefix string used for every output file
PREFIX='FitHiChIP'

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

QVALUE=0.01

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

# option to note down the timing information
TimeProf=0

# boolean variable indicating that previous existing output files
# will be overwritten (1) or not (0 - default)
OverWrite=0

#=========================
# Significant interaction call 
# and the bias correction related parameters
#=========================
# denotes the type of Interaction
# 1: peak to peak, 2: peak to non peak
# 3: peak to all (default) 4: all to all
# 5: accounting for all of 1 to 4
IntType=3

# variable indicating bias correction (1: On, 0: off)
# recommended = 1
BiasCorr=1

# type of bias vector (if bias correction is employed)
# 1: coverage specific bias
# 2: ICE specific bias (default)
BiasType=1

# temporary variable (binary)
# used to model FitHiChIP peak to all interactions
# using peak to peak background only
# applicable for only peak to all interactions
UseP2PBackgrnd=1

# Merging interactions which are near to each other
MergeInteraction=1

#========================
# following parameters are now obsolete
# so keeping them in a default state
# user cannot alter these values
#========================

# # reference genome
# # sourya - this parameter can be removed
# RefGENOME='hg19'

# # binning method for FitHiC technique
# # 1 for equal occupancy bin (default)
# FitHiCBinMethod=1

# # type of distribution for modeling the P value of FitHiC
# # 1: binomial distribution (employed in FitHiC - default)
# # 2: negative binomial distribution
# DistrType=1

# default lower cutoff of bias value
biaslowthr=0.2

# default higher cutoff of bias value
biashighthr=5

# boolean variable for pre-filtering the interactions according to the bias value
BeginBiasFilter=0

# boolean variable for probability adjustment of the interactions according to the bias value
EndBiasFilter=0

# temporary variable (binary)
# if 1, includes only nonzero contact based locus pairs for
# FitHiC spline fit implementation
# default : 0
UseNonzeroContacts=0

# two variables used for bias correction
# modeling the regression between observed contact count
# and the bias variables
resid_biascorr=0
eqocc_biascorr=1

# boolean variable indicating whether for bias correction
# multiplicative bias value would be used
# defult 0
MultBias=0

#========================

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
			# if there are multiple parameter values (separated by # - old values are kept)
			# then the following operation selects the current one
			paramval=$(echo "$paramval" | awk -F['#\t'] '{print $1}' | tr -d '[:space:]');
			if [ $param == "ValidPairs" ]; then
				InpValidPairsFile=$paramval
			fi
			if [ $param == "Interval" ]; then
				InpBinIntervalFile=$paramval
			fi
			if [ $param == "Matrix" ]; then
				InpMatrixFile=$paramval
			fi
			if [ $param == "Bed" ]; then
				InpInitialInteractionBedFile=$paramval
			fi
			if [ $param == "PeakFile" ]; then
				PeakFILE=$paramval
			fi
			if [ $param == "OutDir" ]; then
				if [[ ! -z $paramval ]]; then
					OutDir=$paramval
				fi
			fi
			if [ $param == "ChrSizeFile" ]; then
				ChrSizeFile=$paramval
			fi

			# these parameters are optional
			# if [ $param == "RefGenome" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		RefGENOME=$paramval
			# 	fi
			# fi			
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
			# end optional parameters

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
			if [ $param == "MergeInt" ]; then
				if [[ ! -z $paramval ]]; then
					MergeInteraction=$paramval
				fi
			fi

			# these four variables are added - sourya
			if [ $param == "BiasType" ]; then
				if [[ ! -z $paramval ]]; then
					BiasType=$paramval
				fi
			fi
			if [ $param == "UseP2PBackgrnd" ]; then
				if [[ ! -z $paramval ]]; then
					UseP2PBackgrnd=$paramval
				fi
			fi			
			if [ $param == "OverWrite" ]; then
				if [[ ! -z $paramval ]]; then
					OverWrite=$paramval
				fi
			fi
			if [ $param == "IntType" ]; then
				if [[ ! -z $paramval ]]; then
					IntType=$paramval
				fi
			fi
			if [ $param == "BiasCorrection" ]; then
				if [[ ! -z $paramval ]]; then
					BiasCorr=$paramval
				fi
			fi
			# end add - sourya

		fi
	fi
done < $ConfigFile

#===================
# verify the input parameters
#===================


if [[ -z $PeakFILE ]]; then
	echo 'User should provide a reference peak file (either ChIP or HiChIP) to compute the interactions involving peak segments - exit !!'
	exit 1
fi

#======================
# these parameters are not mandatory 
# so warning messages are removed

# if [[ -z $RefFastaFile ]]; then
# 	echo 'User should provide reference chromosome fasta file - quit !!'
# 	exit 1
# fi

# if [[ -z $REFragFile ]]; then
# 	echo 'User should provide a restriction fragment file - quit !!'
# 	exit 1
# fi

# if [[ -z $MappabilityFile ]]; then
# 	echo 'Reference mappability file is not specified - exit !!'
# 	exit 1
# fi
#======================


if [[ -z $ChrSizeFile ]]; then
	echo 'ERROR ====>>> Chromosome size file is not specified - exit !!'
	exit 1
fi

#*****************************
# error checking - 
# check if the required libraries are all installed or not
#*****************************
if [ 1 == 1 ]; then

# boolean variable indicating error condition
errcond=0




# directory of the configuration file
ConfigFileDir=$(dirname "${ConfigFile}")
BINDARG="$ConfigFileDir,$TMPDIR"

# first go to the configuration file directory
cd $ConfigFileDir

if [[ ! -z $InpValidPairsFile ]]; then
	if [ ! -f $InpValidPairsFile ]; then
		echo 'ERROR ===>>>> Valid pairs file provided as input : '$InpValidPairsFile
		echo 'However, no such valid pairs file exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi
	if [[ "${InpValidPairsFile:0:1}" != / && "${InpValidPairsFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		InpValidPairsFile="$(cd $(dirname $InpValidPairsFile); pwd)/$(basename $InpValidPairsFile)"
	fi
		BINDARG+=",$(dirname $InpValidPairsFile)"
fi

if [[ ! -z $InpBinIntervalFile ]]; then
	if [ ! -f $InpBinIntervalFile ]; then
		echo 'ERROR ===>>>> Bin interval file provided as input : '$InpBinIntervalFile
		echo 'However, no such bin interval file exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi
	if [[ "${InpBinIntervalFile:0:1}" != / && "${InpBinIntervalFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		InpBinIntervalFile="$(cd $(dirname $InpBinIntervalFile); pwd)/$(basename $InpBinIntervalFile)"
	fi
		BINDARG+=",$(dirname $InpBinIntervalFile)"
fi

if [[ ! -z $InpMatrixFile ]]; then
	if [ ! -f $InpMatrixFile ]; then
		echo 'ERROR ===>>>> HiChIP matrix file provided as input : '$InpMatrixFile
		echo 'However, no such matrix file exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi	
	if [[ "${InpMatrixFile:0:1}" != / && "${InpMatrixFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		InpMatrixFile="$(cd $(dirname $InpMatrixFile); pwd)/$(basename $InpMatrixFile)"
	fi
		BINDARG+=",$(dirname $InpMatrixFile)"
fi

if [[ ! -z $InpInitialInteractionBedFile ]]; then
	if [ ! -f $InpInitialInteractionBedFile ]; then
		echo 'ERROR ===>>>> Pre-computed locus pairs (along with their contact counts) are provided as input : '$InpInitialInteractionBedFile
		echo 'However, no such file containing locus pairs exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi	
	if [[ "${InpInitialInteractionBedFile:0:1}" != / && "${InpInitialInteractionBedFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		InpInitialInteractionBedFile="$(cd $(dirname $InpInitialInteractionBedFile); pwd)/$(basename $InpInitialInteractionBedFile)"
	fi
		BINDARG+=",$(dirname $InpInitialInteractionBedFile)"
fi

if [[ ! -z $PeakFILE ]]; then
	if [ ! -f $PeakFILE ]; then
		echo 'ERROR ===>>>> Peak file provided as input : '$PeakFILE
		echo 'However, no such peak file exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi
	if [[ "${PeakFILE:0:1}" != / && "${PeakFILE:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		PeakFILE="$(cd $(dirname $PeakFILE); pwd)/$(basename $PeakFILE)"
	fi
		BINDARG+=",$(dirname $PeakFILE)"
fi

if [[ ! -z $RefFastaFile ]]; then
	if [ ! -f $RefFastaFile ]; then
		echo 'ERROR ===>>>> Fasta file of reference genome provided as input : '$RefFastaFile
		echo 'However, no such fasta file exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi
	if [[ "${RefFastaFile:0:1}" != / && "${RefFastaFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		RefFastaFile="$(cd $(dirname $RefFastaFile); pwd)/$(basename $RefFastaFile)"
	fi
		BINDARG+=",$(dirname $RefFastaFile)"
fi

if [[ ! -z $REFragFile ]]; then
	if [ ! -f $REFragFile ]; then
		echo 'ERROR ===>>>> Restriction fragment file provided as input : '$REFragFile
		echo 'However, no such restriction fragment file exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi
	if [[ "${REFragFile:0:1}" != / && "${REFragFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		REFragFile="$(cd $(dirname $REFragFile); pwd)/$(basename $REFragFile)"
	fi
		BINDARG+=",$(dirname $REFragFile)"
fi

if [[ ! -z $ChrSizeFile ]]; then
	if [ ! -f $ChrSizeFile ]; then
		echo 'ERROR ===>>>> File containing size of chromosomes for the reference genome, as provided : '$ChrSizeFile
		echo 'However, no such chromosome size containing file exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi
	if [[ "${ChrSizeFile:0:1}" != / && "${ChrSizeFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		ChrSizeFile="$(cd $(dirname $ChrSizeFile); pwd)/$(basename $ChrSizeFile)"
	fi
		BINDARG+=",$(dirname $ChrSizeFile)"
fi

if [[ ! -z $MappabilityFile ]]; then
	if [ ! -f $MappabilityFile ]; then
		echo 'ERROR ===>>>> Mappability information containing file provided as input : '$MappabilityFile
		echo 'However, no such Mappability related file exists - check input file and path - FitHiChIP quits !!'
		exit 1
	fi
	if [[ "${MappabilityFile:0:1}" != / && "${MappabilityFile:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		MappabilityFile="$(cd $(dirname $MappabilityFile); pwd)/$(basename $MappabilityFile)"
	fi
		BINDARG+=",$(dirname $MappabilityFile)"
fi

if [[ ! -z $OutDir ]]; then
	if [[ "${OutDir:0:1}" != / && "${OutDir:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		OutDir="$(cd $(dirname $OutDir); pwd)/$(basename $OutDir)"
	fi
		BINDARG+=",$OutDir"
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

singularity exec -B "$BINDARG" library://tuvan/default/fithichip bash "/FitHiChIP/FitHiChIP_HiCPro.sh" -C "$ConfigFile" -s
fi

