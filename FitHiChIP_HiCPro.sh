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

while getopts "sC:" opt;
do
	case "$opt" in
		C) ConfigFile=$OPTARG;;
		s) Singularity=true;;
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

echo -e "\n ================ Parsing input configuration file ================="

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
			echo -e "Content of $param is $paramval"
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

			# these parameters are not required for the moment
			# if [ $param == "BeginBiasFilter" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		BeginBiasFilter=$paramval
			# 	fi
			# fi
			# if [ $param == "EndBiasFilter" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		EndBiasFilter=$paramval
			# 	fi
			# fi
			# if [ $param == "biaslowthr" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		biaslowthr=$paramval
			# 	fi
			# fi			
			# if [ $param == "biashighthr" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		biashighthr=$paramval
			# 	fi
			# fi
			# if [ $param == "DistrType" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		DistrType=$paramval
			# 	fi
			# fi			
			# if [ $param == "MultBias" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		MultBias=$paramval
			# 	fi
			# fi	
			# if [ $param == "BiasCorrResid" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		resid_biascorr=$paramval
			# 	fi
			# fi	
			# if [ $param == "BiasCorrEqOcc" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		eqocc_biascorr=$paramval
			# 	fi
			# fi	
			# if [ $param == "UseNonzeroContacts" ]; then
			# 	if [[ ! -z $paramval ]]; then
			# 		UseNonzeroContacts=$paramval
			# 	fi
			# fi			
			# end non-required parameters

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

echo -e "\n ================ Verifying input configuration parameters ================="

if [[ -z $InpValidPairsFile && -z $InpInitialInteractionBedFile ]]; then
	if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then
		echo -e 'There are three ways to provide FitHiChIP input: \n 1) provide valid pairs file (from HiC pro pipeline), \n 2) Provide both bin interval and matrix files (from HiC pro pipeline), or \n 3) provide a set of locus pairs along with their contact counts (7 columns). \n But user did not provide any of these input configurations. - exit !!'
		exit 1
	fi
fi

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
if [[ ! -z "$Singularity" ]]; then
	HiCProBasedir="/HiC-Pro-2.11.1/"
fi

if [[ -z $HiCProBasedir && -z "$Singularity" && -z $InpInitialInteractionBedFile ]]; then
	if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then
		echo 'ERROR ====>>> As user did not provide any pre-computed locus pairs (along with their contact count) in BED input, neither provided the HiC-pro base installation directory, FitHiChIP quits - exit !!'
		exit 1
	fi
fi

if [[ -z $ChrSizeFile ]]; then
	echo 'ERROR ====>>> Chromosome size file is not specified - exit !!'
	exit 1
fi

echo '***** Specified output directory : '$OutDir

#*****************************
# error checking - 
# check if the required libraries are all installed or not
#*****************************
if [ 1 == 1 ]; then

# boolean variable indicating error condition
errcond=0

# first check the HiC-pro installation
HiCProVersion=$(HiC-Pro -v | head -n 1 | awk -F[" "] '{print $3}' -)
if [[ -z "$HiCProVersion" ]]; then
    echo "ERROR ====>>> HiC-pro is not installed in the system !!! FitHiChIP quits !!!" 
    errcond=1
    # exit 1
else
	echo "HiC-pro is installed in the system"
fi
parsedVersion=$(echo "${HiCProVersion//./}")
if [[ "$parsedVersion" -lt "290" ]]; then 
    echo "ERROR ====>>> HiC-pro should have version >=  2.9.0 !!! FitHiChIP quits !!!"
    errcond=1
    # exit 1
else
	echo "Installed HiC-pro version: "${HiCProVersion}
fi

# first check the python version
# check if python is installed and its version is > 2.7.0
# pythonversion=$(python -V 2>&1 | grep -Po '(?<=Python )(.+)')
pythonversion=$(python --version 2>&1 | head -n 1 | awk '{print $2}' -)
if [[ -z "$pythonversion" ]]; then
    echo "ERROR ====>>> No Python installation is detected in the system !!! FitHiChIP quits !!!" 
    errcond=1
    # exit 1
fi
parsedVersion=$(echo "${pythonversion//./}")
# echo "parsedVersion : "$parsedVersion
if [[ "$parsedVersion" -lt "300" && "$parsedVersion" -gt "270" ]]; then 
    echo "*** Valid python version is detected - installed version: "$pythonversion
elif [[ "$parsedVersion" -lt "3000" && "$parsedVersion" -gt "2700" ]]; then 
	echo "*** Valid python version is detected - installed version: "$pythonversion
else
    echo "ERROR ====>>> Invalid python version : "$pythonversion
    echo " --- should be python 2 with version > 2.7.0 !!! FitHiChIP quits !!!"
    errcond=1
    # exit 1
fi

# check if python libraries are also installed
if python -c "import gzip"; then
    echo '*** Python library gzip is installed'
else
    echo 'ERROR ====>>> Python library gzip is not installed !!! FitHiChIP quits !!!'
    errcond=1
    # exit 1
fi
if python -c "from optparse import OptionParser"; then
    echo '*** Python module OptionParser (from the package optparse) is installed'
else
    echo 'ERROR ====>>> Python module OptionParser (from the package optparse) is not installed !!! FitHiChIP quits !!!'
    errcond=1
    # exit 1
fi
if python -c "import networkx"; then
    echo '*** Python package networkx is installed'
else
    echo 'ERROR ====>>> Python package networkx is not installed !!! FitHiChIP quits !!!'
    errcond=1
    # exit 1
fi

# check if MACS2 package is installed
macs2version=$(macs2 --version 2>&1 |  head -n 1 | awk -F[" "] '{print $2 }' -)
if [[ -z "$macs2version" ]]; then
    echo "ERROR ====>>> MACS2 peak calling package is not detected in the system !!! FitHiChIP quits !!!" 
    errcond=1
    # exit 1
else
	echo "*** Found MACS2 package (for peak calling) installed in the system -  "$macs2version
fi

# check the R version
Rversion=$(R --version | head -n 1 | awk -F[" "] '{print $3}' -)
if [[ -z "$Rversion" ]]; then
    echo "ERROR ====>>> No R installation is detected in the system !!! FitHiChIP quits !!!" 
    errcond=1
    # exit 1
fi
parsedVersion=$(echo "${Rversion//./}")
if [[ "$parsedVersion" -ge "330" && "$parsedVersion" -lt "1000" ]]; then 
    echo "*** Valid R version is detected - installed R version: "$Rversion
elif [[ "$parsedVersion" -ge "3300" && "$parsedVersion" -lt "10000" ]]; then 
    echo "*** Valid R version is detected - installed R version: "$Rversion
elif [[ "$parsedVersion" -ge "31000" ]]; then 
    echo "*** Valid R version is detected - installed R version: "$Rversion      
else
    echo "ERROR ====>>> Invalid R version: "$Rversion
    echo " -- should be at least R 3.3 !!! FitHiChIP quits !!!"
    errcond=1
    # exit 1
fi

# check the samtools version
samtoolsversion=$(samtools 2>&1 | grep "Version" | awk -F[" "] '{print $2}' -)
if [[ -z "$samtoolsversion" ]]; then
    echo "ERROR ====>>> No samtools installation is detected in the system !!! FitHiChIP quits !!!" 
    errcond=1
    # exit 1
fi
parsedsamtoolsVersion=$(echo "${samtoolsversion//./}")
if [[ ${samtoolsversion:0:1} == "0" ]]; then
	# samtools version is 0.*
	echo "ERROR ====>>> Invalid samtools version : "$samtoolsversion 
	echo " - should be at least 1.6 !!! FitHiChIP quits !!!"
	errcond=1
elif [[ "$parsedsamtoolsVersion" -ge "16" && "$parsedsamtoolsVersion" -lt "99" ]]; then
    echo "*** Valid samtools version is detected - installed version: "$samtoolsversion
elif [[ "$parsedsamtoolsVersion" -ge "160" && "$parsedsamtoolsVersion" -lt "999" ]]; then
	echo "*** Valid samtools version is detected - installed version: "$samtoolsversion
else
    echo "ERROR ====>>> Invalid samtools version : "$samtoolsversion 
    echo " - should be at least 1.6 !!! FitHiChIP quits !!!"
    errcond=1
    # exit 1
fi

# check if bgzip is installed in the system
bgzipnum=`bgzip -h 2>&1 | wc -l`
if [[ $bgzipnum -gt 5 ]]; then
	echo "*** bgzip utility is installed in the system"
else
	echo "ERROR =====>> bgzip utility is NOT installed in the system -- please install it from htslib (associated with samtools) version >= 1.6"
	errcond=1
	# exit 1
fi

# check if tabix is installed in the system
tabixnum=`tabix -h 2>&1 | wc -l`
if [[ $tabixnum -gt 5 ]]; then
	echo "*** tabix utility is installed in the system"
else
	echo "ERROR =====>> tabix utility is NOT installed in the system -- please install it from htslib (associated with samtools) version >= 1.6"
	errcond=1
	# exit 1
fi

# check the bedtools version
# command 1
bedtoolsversion1=$(bedtools --version | head -n 1 | awk -F[" "] '{print substr($2,2); }' -)
# echo "bedtoolsversion1: "$bedtoolsversion1
# command 2
bedtoolsversion2=$(bedtools 2>&1 | grep "Version" | awk -F[" "] '{print substr($NF,2); }' -)
# echo "bedtoolsversion2: "$bedtoolsversion2
# if both commands fail then surely bedtools is not installed in the system
if [[ -z "$bedtoolsversion1" && -z "$bedtoolsversion2" ]]; then
    echo "ERROR ====>>> No bedtools installation is detected in the system !!! FitHiChIP quits !!!" 
    errcond=1
    # exit 1
fi
if [[ ! -z "$bedtoolsversion1" ]]; then
    parsedbedtoolsVersion=$(echo "${bedtoolsversion1//./}")
    if [[ "$parsedbedtoolsVersion" -ge "2260" && "$parsedbedtoolsVersion" -lt "9999" ]]; then 
        echo "*** Valid bedtools version is detected - installed version: "$bedtoolsversion1
    elif [[ "$parsedbedtoolsVersion" -ge "300" && "$parsedbedtoolsVersion" -lt "999" ]]; then
        echo "*** Valid bedtools version is detected - installed version: "$bedtoolsversion1
    else
        echo "ERROR ====>>> Invalid bedtools version : "$bedtoolsversion1
        echo " - should be at least 2.26.0 !!! FitHiChIP quits !!!"
        errcond=1
        # exit 1
    fi
fi
if [[ -z "$bedtoolsversion1" && ! -z "$bedtoolsversion2" ]]; then
    parsedbedtoolsVersion=$(echo "${bedtoolsversion2//./}")
    if [[ "$parsedbedtoolsVersion" -ge "2260" && "$parsedbedtoolsVersion" -lt "9999" ]]; then 
        echo "*** Valid bedtools version is detected - installed version: "$bedtoolsversion2
    elif [[ "$parsedbedtoolsVersion" -ge "300" && "$parsedbedtoolsVersion" -lt "999" ]]; then
        echo "*** Valid bedtools version is detected - installed version: "$bedtoolsversion2
    else
        echo "ERROR ====>>> Invalid bedtools version : "$bedtoolsversion2
        echo " - should be at least 2.26.0 !!! FitHiChIP quits !!!"
        errcond=1
        # exit 1
    fi
fi

# final evaluation
if [[ $errcond == 1 ]]; then
	echo " --------->>>>>> ERROR - current system has one or more packages missing" 
	echo " ------- please check the above logs and install the missing packages before executing FitHiChIP"
	echo " ------- exiting for the moment !!"
	exit 1
fi

#*****************************
fi 	# end dummy if - testing installed packages
#*****************************


#*****************************
# here check if the configuration file has relative path names as the input
# in such a case, convert the relative path names (with respect to the location of the configuration file itself)
# in the absolute file

# also perform Error checking - check if input files provided are correct
# i.e. these files exist

echo -e "\n ================ Changing relative pathnames of the input files to their absolute path names ================="


# directory of the configuration file
ConfigFileDir=$(dirname "${ConfigFile}")

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
		echo 'Absolute converted path: InpValidPairsFile: '$InpValidPairsFile
	fi
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
		echo 'Absolute converted path: InpBinIntervalFile: '$InpBinIntervalFile
	fi
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
		echo 'Absolute converted path: InpMatrixFile: '$InpMatrixFile
	fi
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
		echo 'Absolute converted path: InpInitialInteractionBedFile: '$InpInitialInteractionBedFile
	fi
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
		echo 'Absolute converted path: PeakFILE: '$PeakFILE
	fi
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
		echo 'Absolute converted path: RefFastaFile: '$RefFastaFile
	fi
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
		echo 'Absolute converted path: REFragFile: '$REFragFile
	fi
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
		echo 'Absolute converted path: ChrSizeFile: '$ChrSizeFile
	fi
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

if [[ -z "$Singularity" && ! -z $HiCProBasedir ]]; then
	if [ ! -d $HiCProBasedir ]; then
		echo 'ERROR ===>>>> Base directory of HiC-Pro package as provided : '$HiCProBasedir
		echo 'However, no such directory exists in the system - check input file and path - FitHiChIP quits !!'
		exit 1
	fi	
	if [[ "${HiCProBasedir:0:1}" != / && "${HiCProBasedir:0:2}" != ~[/a-z] ]]; then
		# relative path - convert to absolute path
		HiCProBasedir="$(cd $(dirname $HiCProBasedir); pwd)/$(basename $HiCProBasedir)"
		echo 'Absolute converted path: HiCProBasedir: '$HiCProBasedir
	fi
fi

# revert to the old directory
cd -

#===================
# further error checking - check for false values
#===================

if [[ $IntType -lt 1 || $IntType -gt 5 ]]; then
	echo 'ERROR ===>>>> Parameter IntType should have value between 1 to 5: -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $BIN_SIZE -lt 0 ]]; then
	echo 'ERROR ===>>>> Parameter Bin size is specified negative: -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $LowDistThres -lt 0 ]]; then
	echo 'ERROR ===>>>> Parameter low distance threshold is specified negative: -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $UppDistThres -lt 0 ]]; then
	echo 'ERROR ===>>>> Parameter upper distance threshold is specified negative: -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $UppDistThres -lt $LowDistThres ]]; then
	echo 'ERROR ===>>>> Parameter upper distance threshold is specified lower than the lower distance threshold -- check both values -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $UseP2PBackgrnd -lt 0 || $UseP2PBackgrnd -gt 1 ]]; then
	echo 'ERROR ===>>>> Parameter UseP2PBackgrnd should have values either 0 or 1: -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $BiasType -lt 1 || $BiasType -gt 2 ]]; then
	echo 'ERROR ===>>>> Parameter BiasType should have values either 1 (indicating coverage bias regression) or 2 (ICE bias regression): -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $MergeInteraction -lt 0 || $MergeInteraction -gt 1 ]]; then
	echo 'ERROR ===>>>> Parameter MergeInt should have values either 0 or 1: -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $QVALUE -lt 0 || $QVALUE -ge 1 ]]; then
	echo 'ERROR ===>>>> Parameter q-value (FDR) should be positive but close to 0 (default 0.01) - but here invalid values are provided - check the parameters -- FitHiChIP quits !!!'
	exit 1
fi

if [[ $OverWrite -lt 0 || $OverWrite -gt 1 ]]; then
	echo 'ERROR ===>>>> Parameter OverWrite should have values either 0 or 1: -- FitHiChIP quits !!!'
	exit 1
fi


#===================
# create the output directory
#===================
mkdir -p $OutDir

#============================
# print the parameters and values
#============================
ConfFile=$OutDir/Parameters.txt

if [[ ! -f $ConfFile || $OverWrite == 1 ]]; then
	echo "Summarizing the parameters employed in this execution" > $ConfFile

	echo "Listing the input files: " >> $ConfFile
	if [[ ! -z $InpValidPairsFile ]]; then
		echo "Input file containing the valid pairs (generated from HIC-PRO pipeline): $InpValidPairsFile " >> $ConfFile
	fi
	if [[ ! -z $InpBinIntervalFile ]]; then
		echo "Input file containing the fixed size bin intervals (generated from HIC-PRO pipeline): $InpBinIntervalFile " >> $ConfFile
	fi
	if [[ ! -z $InpMatrixFile ]]; then
		echo "Input file storing the contact matrix between the fixed size bins (generated from HIC-PRO pipeline): $InpMatrixFile " >> $ConfFile
	fi
	if [[ ! -z $InpInitialInteractionBedFile ]]; then
		echo "Input file containing pre-computed locus pairs and contact count: $InpInitialInteractionBedFile " >> $ConfFile
	fi
	echo "Input ChIP-seq peak file: $PeakFILE " >> $ConfFile
	echo "File containing chromosome size information corresponding to the reference chromosome: $ChrSizeFile " >> $ConfFile

	if [[ ! -z $MappabilityFile ]]; then
		echo "Mappability File: $MappabilityFile " >> $ConfFile
	fi
	if [[ ! -z $RefFastaFile ]]; then
		echo "Reference Fasta sequence File: $RefFastaFile " >> $ConfFile
	fi
	if [[ ! -z $REFragFile ]]; then
		echo "Restriction Fragment File: $REFragFile " >> $ConfFile
	fi

	echo "Output directory which will store all the results: $OutDir " >> $ConfFile
	if [[ ! -z $HiCProBasedir ]]; then
		if [[ -z "$Singularity" ]]; then
			echo "Base directory containing the HIC-PRO executable: $HiCProBasedir " >> $ConfFile
		fi
		if [[ ! -z "$Singularity" ]]; then
			echo "Base directory containing the HIC-PRO executable: /HiC-Pro-2.11.1/ " >> $ConfFile
		fi
	fi

	echo "Low distance threshold: $LowDistThres " >> $ConfFile
	echo "Higher distance threshold: $UppDistThres " >> $ConfFile

	# echo "Genome specific parameters: " >> $ConfFile
	# echo "RefGENOME: $RefGENOME " >> $ConfFile
	
	if [[ ! -z $MappabilityFile && ! -z $RefFastaFile && ! -z $REFragFile ]]; then
		echo "GCContentWindowSize: $GCContentWindowSize " >> $ConfFile
		echo "MappabilityWindowSize: $MappabilityWindowSize " >> $ConfFile
	fi

	echo "BIN_SIZE: $BIN_SIZE " >> $ConfFile
	echo "PREFIX (string of output files): $PREFIX " >> $ConfFile
	# echo "Timeprof (1 means timing will be noted): $TimeProf " >> $ConfFile
	echo "OverWrite (1 means existing files will be overwritten): $OverWrite " >> $ConfFile
	# echo "DrawFig (1 means various summary statistics will be plotted, depending on the availability of reference fasta sequence file and mappability file): $DrawFig " >> $ConfFile
fi

#=======================================
# generate a file which will contain the timing profile
if [ $TimeProf == 1 ]; then
	OutTimeFile=$OutDir'/TimingProfile.txt'
	if [[ ! -f $OutTimeFile || $OverWrite == 1 ]]; then
		echo " ================ Time profiling =========== " > $OutTimeFile
	else 
		echo " ================ Time profiling (append) =========== " >> $OutTimeFile
	fi
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

#==============================
# note down the executables of python and Rscript
PythonExec=`which python`
echo 'Executable of python (2): '$PythonExec

RScriptExec=`which Rscript`
echo 'Executable of R : '$RScriptExec

# #************************
# # error checking
# # check if R libraries are also installed
# # otherwise, install those packages
# $RScriptExec ./Analysis/PackageTest.r

#================================

# generate the matrix of Hi-C interactions (ALL)
# using HiC-pro pipeline
HiCProMatrixDir=$OutDir'/HiCPro_Matrix_BinSize'$BIN_SIZE
mkdir -p $HiCProMatrixDir

# if pre-computed locus pairs are provided as input
# then use this interaction file
if [[ ! -z $InpInitialInteractionBedFile ]]; then

	echo -e "\n ================ Processing input pre-computed locus pairs along with the contact count ================="

	# check whether the file is gzippe or not
	if [[ $InpInitialInteractionBedFile == *.gz ]]; then
		gzipInt=1
		echo -e "\n The input locus pair file is in gzipped format --- "
	else
		gzipInt=0
	fi

	Interaction_Initial_File=$HiCProMatrixDir/$PREFIX.interactions.initial.bed
	
	# dump the lines having 2nd, 3rd, 5th, 6th and 7th fields as integer
	if [ $gzipInt == 1 ]; then
		zcat $InpInitialInteractionBedFile | awk '{if (($2 ~ /^[0-9]+$/) && ($3 ~ /^[0-9]+$/) && ($5 ~ /^[0-9]+$/) && ($6 ~ /^[0-9]+$/) && ($7 ~ /^[0-9]+$/)) {print $0}}' - | cut -f1-7 | awk '($7>0)' - > $Interaction_Initial_File
	else
		awk '{if (($2 ~ /^[0-9]+$/) && ($3 ~ /^[0-9]+$/) && ($5 ~ /^[0-9]+$/) && ($6 ~ /^[0-9]+$/) && ($7 ~ /^[0-9]+$/)) {print $0}}' $InpInitialInteractionBedFile | cut -f1-7 | awk '($7>0)' - > $Interaction_Initial_File
	fi
	
	# insert a header line
	sed -i '1ichr\ts1\te1\tchr2\ts2\te2\tcc' $Interaction_Initial_File

	echo -e "\n Dumped the set of input interactions along with the header --- "

else 

	#=================
	# if the matrices are not provided and the validpairs text file is provided
	# then compute the interaction matrices using the HiC-pro utility
	#=================
	echo -e "\n ================ Processing HiC-pro contact matrices ================="

	if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then

		# this is an executable which builds matrix from the input valid pairs file
		# thats why we require the HiC pro executable directory as a command line option

		# executable of matrix building is to be obtained from the HiC-pro base directory
		# provided as the input
		# MatrixBuildExec=$HiCProBasedir'/scripts/build_matrix'
		MatrixBuildExecSet=( $(find $HiCProBasedir -type f -name 'build_matrix') )
		len=${#MatrixBuildExecSet[@]}
		# echo 'len: '$len
		
		#==========================
		# error in the reference HiC-pro package
		# if such script is not found
		if [[ $len == 0 ]]; then		
			echo 'ERROR ===>>>> could not find the executable build_matrix within HiC-pro installed package ----'
			echo 'Check the directory of HiC-pro Package or check its version (recommended >= 2.9.0)'
			echo 'Matrix from the input valid pairs file could not be generated - FitHiChIP is quiting !!!'			
			exit 1
		fi
		#==========================
	
		idx=`expr $len - 1`
		# echo 'idx: '$idx
		MatrixBuildExec=${MatrixBuildExecSet[$idx]}
		echo -e '\n *** MatrixBuildExec: '$MatrixBuildExec

		echo '*** Computing HiC-pro matrices from the input valid pairs file'

		# This directory and prefix is used to denote the generated matrices
		# using the HiC pro routine
		OutPrefix=$HiCProMatrixDir'/MatrixHiCPro'

		if [[ ! -f $OutPrefix'_abs.bed' || $OverWrite == 1 ]]; then
			# check the extension of input valid pairs file
			# and extract accordingly
			if [[ $InpValidPairsFile == *.gz ]]; then
				echo "***** HiC-pro input valid pairs file in gzipped format"
				zcat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper'  
			else
				echo "***** HiC-pro input valid pairs file in simple text format"
				cat $InpValidPairsFile | $MatrixBuildExec --binsize $BIN_SIZE --chrsizes $ChrSizeFile --ifile /dev/stdin --oprefix $OutPrefix --matrix-format 'upper' 
			fi
			if [ $TimeProf == 1 ]; then
				duration=$(echo "$(date +%s.%N) - $start" | bc)
				echo " ++++ Time (in seconds) for computing the interaction matrix using HiC-Pro build_matrix utility: $duration" >> $OutTimeFile
				start=$(date +%s.%N)
			fi		
		fi

		# now assign the matrix names to the designated variables
		InpBinIntervalFile=$OutPrefix'_abs.bed'
		InpMatrixFile=$OutPrefix'.matrix'

	fi

	#=======================
	# Now generate the list of interactions from the HiC-pro matrix data
	# ALL to ALL interactions
	# Both cis and trans interactions are considered (with respect to the given bin size)
	# No distance threshold based filtering is used
	# Interaction format:
	# chr1	start1	end1	chr2	start2	end2	cc
	#=======================
	echo -e "\n ================ Creating input interactions ================="

	Interaction_Initial_File=$HiCProMatrixDir/$PREFIX.interactions.initial.bed

	if [[ ! -f $Interaction_Initial_File || $OverWrite == 1 ]]; then
		$RScriptExec ./src/InteractionHicPro.r $InpBinIntervalFile $InpMatrixFile $Interaction_Initial_File
		
		if [ $TimeProf == 1 ]; then
			duration=$(echo "$(date +%s.%N) - $start" | bc)
			echo " ++++ Time (in seconds) for getting CIS interactions: $duration" >> $OutTimeFile
			start=$(date +%s.%N)
		fi
	fi
fi

# check the number of interactions
numLoop=`cat $Interaction_Initial_File | wc -l`
echo 'Number of locus pairs with nonzero contact count (without any distance thresholding): '$numLoop
if [ $numLoop == 0 ]; then
	echo 'Number of locus pairs with nonzero contact count is zero - FitHiChIP is quiting !!!'			
	exit 1
fi


#=======================
# generate filtered cis interactions 
# with respect to distance thresholds
# ALL to ALL interactions
#=======================

echo -e "\n ================ Limiting input interactions to the specified distance ranges ${LowDistThres} to ${UppDistThres} ================="

# create a directory for individual distance thresholds
InteractionThrDir=$HiCProMatrixDir'/L_'$LowDistThres'_U'$UppDistThres
mkdir -p $InteractionThrDir
Interaction_File=$InteractionThrDir/$PREFIX.cis.interactions.DistThr.bed

if [[ ! -f $Interaction_File || $OverWrite == 1 ]]; then
	awk -v l="$LowDistThres" -v u="$UppDistThres" -F['\t'] 'function abs(v) {return v < 0 ? -v : v} {if ((NR==1) || ($1==$4 && ($7>0) && abs($2-$5)>=l && abs($2-$5)<=u)) {print $0}}' $Interaction_Initial_File > $Interaction_File

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting CIS interactions within distance thresholds: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

# check the number of interactions
numLoop=`cat $Interaction_File | wc -l`
echo 'Number of locus pairs with nonzero contact count (after distance thresholding): '$numLoop
if [ $numLoop == 0 ]; then
	echo 'Number of locus pairs with nonzero contact count (after distance thresholding) is zero - FitHiChIP is quiting !!!'			
	exit 1
fi

#============================
# this directory stores the features and associated data
#============================
FeatureDir=$OutDir'/NormFeatures'
mkdir -p $FeatureDir

# if all of these three files exist
# then we compute the RE fragment, GC content, and mappability information
# for individual bins
if [[ ! -z $MappabilityFile && ! -z $RefFastaFile && ! -z $REFragFile ]]; then

	echo -e "\n ================ Processing input mappability, GC content files ================="

	# file storing the RE fragments, mappability and GC content together
	REFragMappGCFile=$FeatureDir'/REFrag_Mapp_GC_Merged.bed'

	if [[ ! -f $REFragMappGCFile || $OverWrite == 1 ]]; then

		#============================
		# generating the mappability information
		# first divide each RE fragment interval
		# from two ends
		# the offset size = $MappabilityWindowSize
		# upstream and downstream 
		# we do not cross the RE fragment boundary - safe for length overflow
		# Important - sourya
		# the final generated file is applied on bedtools map function
		# do, it should be sorted by genome coordinate, using the function sort -k1,1 -k2,2n
		#============================
		echo '======== Creating the fragment end (w.r.t window size) file -- to compute the mappability information!!'
		MappOffsetCutBedFile=$FeatureDir'/Temp_Fragment_Mapp_'$MappabilityWindowSize'bp.bed'
		
		if [[ ! -f $MappOffsetCutBedFile || $OverWrite == 1 ]]; then
			awk -v s="$MappabilityWindowSize" -F['\t'] 'function max(x,y) {return x>y?x:y}; function min(x,y) {return x<y?x:y}; {printf "%s\t%d\t%d\n%s\t%d\t%d\n", $1, $2, min($2+s,$3), $1, max($3-s, $2), $3}' $REFragFile | sort -k1,1 -k2,2n > $MappOffsetCutBedFile

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
		echo '======== Creating the mappability file !!'
		MappabilityOutFile=$FeatureDir'/Mappability_RE_Fragments.bed'
		
		if [[ ! -f $MappabilityOutFile || $OverWrite == 1 ]]; then
			bedtools map -a $MappOffsetCutBedFile -b $MappabilityFile -c 4 -o mean | awk -F['\t'] '{if ($4=="." || $4=="NA" || $4=="NaN") {$4=0}; print $0}' - > $MappabilityOutFile

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
		# Important - sourya
		# the final generated file is applied on bedtools map function
		# do, it should be sorted by genome coordinate, using the function sort -k1,1 -k2,2n
		#============================
		echo '======== Creating the fragment end (w.r.t window size) file -- to compute the GC content information!!'
		GCOffsetCutBedFile=$FeatureDir'/Temp_Fragment_GC_'$GCContentWindowSize'bp.bed'
		
		if [[ ! -f $GCOffsetCutBedFile || $OverWrite == 1 ]]; then
			awk -v s="$GCContentWindowSize" -F['\t'] 'function max(x,y) {return x>y?x:y}; function min(x,y) {return x<y?x:y}; {printf "%s\t%d\t%d\n%s\t%d\t%d\n", $1, $2, min($2+s,$3), $1, max($3-s, $2), $3}' $REFragFile | sort -k1,1 -k2,2n > $GCOffsetCutBedFile

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
		echo '======== Creating the GC content file !!'
		GCOutFile=$FeatureDir'/GC_Content_RE_Fragments.bed'
		
		if [[ ! -f $GCOutFile || $OverWrite == 1 ]]; then
			nucBed -fi $RefFastaFile -bed $GCOffsetCutBedFile | awk -F['\t'] '{if ($4=="." || $4=="NA" || $4=="NaN") {$4=0}; print $0}' - > $GCOutFile

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
		awk -F['\t'] '{sum+=$4} NR%2==0 {print sum/2; sum=0}' $MappabilityOutFile > $Temp_Mapp_File

		# the GC content fragment file has header
		# average of the 5th field
		Temp_GC_File=$FeatureDir'/GC_Dump.bed'
		awk -F['\t'] '{if (NR>1) {sum+=$5}}; {if (NR%2!=0 && NR>1) {print sum/2; sum=0}}' $GCOutFile > $Temp_GC_File

		#============================
		# Now combine the average mappability and GC content values with the 
		# restriction fragments 
		#============================
		$RScriptExec ./src/CombineREFragMappGC.r $REFragFile $Temp_Mapp_File $Temp_GC_File $REFragMappGCFile

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

fi 	# end if mappability, GC content and RE frag files exist

#=================
# modifications - sourya
# from the generated initial CIS interaction file (without distance thresholding)
# get the HiChIP 1D bin coverage values, and also whether those bins overlap with
# input peak segments (either ChIP peaks or HiChIP peaks)

# previously a python script was used
# now a better optimized R script is used
#=================
echo -e "\n ================ Generating coverage statistics for individual bins ================="

CoverageFile=$FeatureDir'/'$PREFIX'.coverage.bed'

if [[ ! -f $CoverageFile || $OverWrite == 1 ]]; then

	# comment - sourya
	# $PythonExec ./src/CoverageBin.py -i $InpValidPairsFile -p $PeakFILE -b $BIN_SIZE -o $CoverageFile -c $ChrSizeFile

	# add - sourya
	$RScriptExec ./src/CoverageBin.r --InpFile $Interaction_Initial_File --PeakFile $PeakFILE --BinSize $BIN_SIZE --ChrSizeFile $ChrSizeFile --OutFile $CoverageFile

	echo '======== Computed initial coverage information for individual genomic bins'
	
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for getting coverage of individual bins: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi



#=======================================
# condition: value of BiasType - 1 means coverage specific bias
# and 2 means ICE specific bias
# if BiasType = 1, bias is computed solely from the coverage. 
# peaks and non-peaks are analyzed separately to compute the bias
# if BiasType = 2, ICE routine from the HiC-pro pipeline is used to compute the ICE bias
# bias information is appended to the coverage and peak information for individual bins
# =======================================

echo -e "\n ================ Merging coverage with bias statistics ================="

if [ $BiasType == 2 ]; then
	AllFeatureDir=$FeatureDir'/ICE_Bias'
	mkdir -p $AllFeatureDir
	CoverageBiasFile=$AllFeatureDir'/'$PREFIX'.coverage_ICE_Bias.bed'
else
	AllFeatureDir=$FeatureDir'/Coverage_Bias'
	mkdir -p $AllFeatureDir
	CoverageBiasFile=$AllFeatureDir'/'$PREFIX'.coverage_Bias.bed'
fi

if [[ ! -f $CoverageBiasFile || $OverWrite == 1 ]]; then
	
	if [ $BiasType == 1 ]; then
		# compute the bias from the coverage of individual bins
		# use peaks and non-peaks separately for bias computation
		$RScriptExec ./src/BiasCalc.r --CoverageFile $CoverageFile --OutFile $CoverageBiasFile
		echo '======== Appended coverage bias information for individual genomic bins'
	else
		# here ICE specific bias is used
		# compute the bias vector from the HiC-pro contact matrix
		ICEExec=$HiCProBasedir'/scripts/ice'
		echo -e '\n *** ICE computation Executable: '$ICEExec
		echo '*** Computing ICE based bias vector from the HiC-pro contact matrix'

		# run the ICE executable and store the normalized contact matrix 
		# and the bias vector
		NormContactMatrixFile=$AllFeatureDir'/'$PREFIX'.norm.Contact.Matrix'
		BiasVecFile=$NormContactMatrixFile'.biases'
		if [[ ! -f $BiasVecFile || $OverWrite == 1 ]]; then
			$ICEExec $InpMatrixFile --results_filename $NormContactMatrixFile --output-bias $BiasVecFile
			# replace the NANs of the derived bias vector with zero
			sed -i 's/nan/0/g' $BiasVecFile
		fi

		# the $BiasVecFile is basically a column vector (without any header information)
		# containing bias values for individual bins
		# these bins correspond to the bins specified in the $InpBinIntervalFile
		temp_ICEBias_binfile=$AllFeatureDir'/'$PREFIX'.temp_bin_bias.txt'
		if [[ ! -f $temp_ICEBias_binfile || $OverWrite == 1 ]]; then
			paste $InpBinIntervalFile $BiasVecFile | cut -f1,2,3,5 > $temp_ICEBias_binfile
		fi

		# merge the files containing coverage + peak information of individual bins
		# and the generated bias containing files
		# Note: these two files may have different ordering of chromosomes
		# so we do not use the "paste" function
		$RScriptExec ./src/BiasCalc.r --CoverageFile $CoverageFile --BiasFile $temp_ICEBias_binfile --OutFile $CoverageBiasFile

		echo '======== Appended ICE bias information for individual genomic bins'
	fi

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for computing bias of individual bins: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi	
fi


#==================================
# merge the bin specific coverage, peak, bias values 
# with the mappability and GC content values
# Note: provided if these mappability, GC content, and RE sites
# are computed (and the corresponding input files are provided as the inputs)

# here the RE fragment file (file b) has mappability information in 4th column 
# and GC content information in 5th column
# mean operation for multiple overlaps is used
# the option "-header" prints the header information of the file (a)
# missing values (no overlap) are indicated by 0

# after the operation, a file with 8 columns is produced
# chr start end coverage ispeak mappability GCcontent NoCutSites

# important - sourya
# before applying bedtools map, check whether the input is sorted by position
#==================================

echo -e "\n ================ Merging coverage + bias with mappability, GC content, and number of cut sites - creating all feature file ================="

# Use of ICE / coverage bias results in different feature files
if [ $BiasType == 2 ]; then
	AllFeatFile=$AllFeatureDir'/'$PREFIX'.AllBin_CompleteFeat_ICE.bed'
else
	AllFeatFile=$AllFeatureDir'/'$PREFIX'.AllBin_CompleteFeat.bed'
fi

if [[ ! -f $AllFeatFile || $OverWrite == 1 ]]; then
	
	#======================
	# first ensure that inputs to bed operation are sorted by chromosome name and coordinate
	#======================
	# the coverage / bias file has header information - discard the header before processing
	temp_CoverageBiasFile=$AllFeatureDir'/'$PREFIX'.coverage_Bias1.bed'
	awk -F['\t'] 'NR>1' $CoverageBiasFile | sort -k1,1 -k2,2n > $temp_CoverageBiasFile
	
	if [[ ! -z $MappabilityFile && ! -z $RefFastaFile && ! -z $REFragFile ]]; then

		# here mappability, GC content, and RE sites information are available
		# merge with the coverage and bias statistics

		temp_REFragMappGCFile=$AllFeatureDir'/REFrag_Mapp_GC_Merged1.bed'
		sort -k1,1 -k2,2n $REFragMappGCFile > $temp_REFragMappGCFile
		
		#======================
		# then apply the map function
		#======================
		# first merge the mappability and GC content information
		AllFeatFile_temp1=$AllFeatureDir'/'$PREFIX'.AllBin_CompleteFeat_temp1.bed'
		bedtools map -c 4,5 -o mean -null '0' -a $temp_CoverageBiasFile -b $temp_REFragMappGCFile > $AllFeatFile_temp1
		
		# then merge the number of RE sites
		bedtools map -a $AllFeatFile_temp1 -b $temp_REFragMappGCFile -c 4 -o count -null '0' > $AllFeatFile
	
		#======================
		# finally remove the temporary files
		#======================
		rm $temp_CoverageBiasFile
		rm $temp_REFragMappGCFile
		rm $AllFeatFile_temp1

	else
		# here we do not have any mappability, RE sites or GC content information
		# so use all 0s in their slots (3 columns to be appended)
		awk -F['\t'] '{if (NR>1) {print $0"\t0\t0\t0"}}' $temp_CoverageBiasFile > $AllFeatFile
	fi

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

	echo -e "\n ================ Plotting bias and non-bias related feature distributions for peaks and non-peaks ================="

	# only bias specific plots would be placed in different directories, based on 
	# whether coverage specific bias or ICE specific bias is used
	
	# condition based on having mappability, GC content, or RE sites 
	# information or not
	if [[ ! -z $MappabilityFile && ! -z $RefFastaFile && ! -z $REFragFile ]]; then
		$RScriptExec ./Analysis/PlotGenomeBins.r --GenomeBinFile $AllFeatFile --CommonDir $FeatureDir'/Plots_Common' --BiasSpecificDir $AllFeatureDir'/Plots_Bias' --OverWrite $OverWrite
	else
		$RScriptExec ./Analysis/PlotGenomeBins.r --GenomeBinFile $AllFeatFile --CommonDir $FeatureDir'/Plots_Common' --BiasSpecificDir $AllFeatureDir'/Plots_Bias' --OverWrite $OverWrite --NoMappGC
	fi

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for plotting bin specific features: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#=================
# templates of filenames used throughout the execution
#=================
InteractionFileName='Interactions.bed'
InteractionSortedDistFileName='Interactions.sortedGenDist.bed'

#============================
# create all to all interaction file
# with the features like read depth, mappability, GC content, bias (coverage / ICE bias)
# and number of RE sites
# depending on the bias type, two different directories are created 
# for each category of interactions
#============================
echo -e "\n ================ Generating interactions + features ================="

DirALLtoALLBase=$OutDir'/FitHiChIP_ALL2ALL_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
if [ $BiasType == 2 ]; then
	DirALLtoALL=$DirALLtoALLBase'/ICE_Bias'
else
	DirALLtoALL=$DirALLtoALLBase'/Coverage_Bias'
fi
mkdir -p $DirALLtoALL
IntFileALLtoALL=$DirALLtoALL'/'$InteractionFileName

if [[ ! -f $IntFileALLtoALL || $OverWrite == 1 ]]; then
	# merge the input interactions (chromosome interval + contact count)
	# with the feature file	
	$RScriptExec ./src/Significance_Features.r -I $Interaction_File -E $AllFeatFile -O $IntFileALLtoALL
	
	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for computing all to all interactions: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

# derive the contact count column
cccol=`cat $Interaction_File | tail -n 1 | awk -F['\t'] '{print NF}' -`
echo 'Contact count col: '$cccol

# derive the number of columns in the interaction file with normalization 
# related features
totcol=`cat $IntFileALLtoALL | tail -n 1 | awk -F['\t'] '{print NF}' -`
echo 'Total number of columns for the complete feature interactions: '$totcol

#===================
# using the merged interaction + feature file for all bins
# plot the variation among different features
#===================
if [ $DrawFig == 1 ]; then

	echo -e "\n ================ Plotting the distribution of normalization features among peak and non peak segments =============="

	# two different directories are employed for plotting
	# 1) common dir for plotting non-bias related features
	# 2) bias specific directory for plotting bias related features
	if [[ ! -z $MappabilityFile && ! -z $RefFastaFile && ! -z $REFragFile ]]; then
		$RScriptExec ./Analysis/InteractionPlots.r --IntFile $IntFileALLtoALL --CommonDir $DirALLtoALLBase'/Plots_Norm' --BiasSpecificDir $DirALLtoALL'/Plots_Norm' --MappThr 0.5 --GCThr 0.2 --cccol $cccol --OverWrite $OverWrite
	else
		$RScriptExec ./Analysis/InteractionPlots.r --IntFile $IntFileALLtoALL --CommonDir $DirALLtoALLBase'/Plots_Norm' --BiasSpecificDir $DirALLtoALL'/Plots_Norm' --MappThr 0.5 --GCThr 0.2 --cccol $cccol --OverWrite $OverWrite --NoMappGC
	fi

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for plotting normalization related features for different types of interactions: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi
fi

#===================
# depending on the input parameter "IntType"
# Create the interaction files for other types of interactions
# peak to peak, peak to non peak, and peak to all
#===================

# boolean variable indicating error condition
errcond=0

if [[ $IntType -ge 1 && $IntType -le 4 ]]; then
	IntLow=$IntType
	IntHigh=$IntType
else
	IntLow=1
	IntHigh=4
fi

echo "Specified IntType: "$IntType
echo "Derived IntLow: "$IntLow
echo "Derived IntHigh: "$IntHigh

#===============
# loop through different types of interactions specified in the input parameters
#===============
# for CurrIntType in $(seq $IntLow $IntHigh); do

CurrIntType=$IntLow
while [[ $CurrIntType -le $IntHigh ]]; do

	echo -e "\n\n **** Start of while Loop ----- current interaction type: $CurrIntType  ****** \n\n"
	if [[ $CurrIntType == 1 ]]; then
		DirPeaktoPeakBase=$OutDir'/FitHiChIP_Peak2Peak_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
		if [ $BiasType == 2 ]; then
			DirPeaktoPeak=$DirPeaktoPeakBase'/ICE_Bias'
		else
			DirPeaktoPeak=$DirPeaktoPeakBase'/Coverage_Bias'
		fi
		mkdir -p $DirPeaktoPeak
		IntFilePeaktoPeak=$DirPeaktoPeak'/'$InteractionFileName
		if [[ ! -f $IntFilePeaktoPeak || $OverWrite == 1 ]]; then
			# peak to peak interactions 
			# 9th and 15th fields are 1
			awk -F['\t'] '((NR==1) || ($9==1 && $15==1))' $IntFileALLtoALL > $IntFilePeaktoPeak
		fi
		currdirname=$DirPeaktoPeak

	elif [[ $CurrIntType == 2 ]]; then
		DirPeaktoNonPeakBase=$OutDir'/FitHiChIP_Peak2NonPeak_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
		if [[ $BiasType == 2 ]]; then
			DirPeaktoNonPeak=$DirPeaktoNonPeakBase'/ICE_Bias'
		else
			DirPeaktoNonPeak=$DirPeaktoNonPeakBase'/Coverage_Bias'
		fi
		mkdir -p $DirPeaktoNonPeak
		IntFilePeaktoNonPeak=$DirPeaktoNonPeak'/'$InteractionFileName
		if [[ ! -f $IntFilePeaktoNonPeak || $OverWrite == 1 ]]; then
			# peak to non peak interactions
			# 9th field is 1, but 15th field is 0
			# or 15th field is 1 and 9th field is 0
			awk -F['\t'] '((NR==1) || ($9==1 && $15==0) || ($9==0 && $15==1))' $IntFileALLtoALL > $IntFilePeaktoNonPeak
		fi
		currdirname=$DirPeaktoNonPeak

	elif [[ $CurrIntType == 3 ]]; then
		# for peak to all interactions, two subdirectories are created
		# depending on the usage of peak to peak background
		DirPeaktoALLBase=$OutDir'/FitHiChIP_Peak2ALL_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres'/P2PBckgr_'$UseP2PBackgrnd
		if [ $BiasType == 2 ]; then
			DirPeaktoALL=$DirPeaktoALLBase'/ICE_Bias'
		else
			DirPeaktoALL=$DirPeaktoALLBase'/Coverage_Bias'
		fi
		mkdir -p $DirPeaktoALL
		IntFilePeaktoALL=$DirPeaktoALL'/'$InteractionFileName
		if [[ ! -f $IntFilePeaktoALL || $OverWrite == 1 ]]; then
			# peak to all interactions
			# just check if 9th field is 1 or 15th field is 1
			awk -F['\t'] '((NR==1) || ($9==1) || ($15==1))' $IntFileALLtoALL > $IntFilePeaktoALL
		fi
		currdirname=$DirPeaktoALL
	else
		currdirname=$DirALLtoALL
	fi

	if [[ $TimeProf == 1 ]]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for assigning the interactions: $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi

	echo -e "\n ============ Calling significant interactions =============== \n"

	# first create interaction files with sorted genomic distance
	CurrIntFile=$currdirname'/'$InteractionFileName
	CurrIntFileSortDist=$currdirname'/'$InteractionSortedDistFileName

	coln=`expr $totcol + 1`
	if [[ ! -f $CurrIntFileSortDist || $OverWrite == 1 ]]; then
		awk -F['\t'] -v OFS='\t' 'function abs(v) {return v < 0 ? -v : v} {print $0"\t"abs($5-$2)}' $CurrIntFile | sort -k${coln},${coln}n -k${cccol},${cccol}nr | cut -f1-${totcol} - > $CurrIntFileSortDist	

		if [[ $TimeProf == 1 ]]; then
			duration=$(echo "$(date +%s.%N) - $start" | bc)
			echo " ++++ Time (in seconds) for sorting the interactions (according to genomic distance): $duration" >> $OutTimeFile
			start=$(date +%s.%N)
		fi
		echo 'Created sorted genomic distance based interaction file'
	fi

	#==============
	# now apply FitHiC on the sorted gene distance based interaction matrix
	# create folders based on the use of bias correction method (or not)
	# and also the use of binomial / negative binomial distribution
	# there are two different modeling used - FitHiC
	# 1) when every candidate interaction is used for spline fitting
	# 2) when only peak to peak interactions are used for spline fitting
	# different output directories are created to contain the results
	#==============
	GenFitHiCDir=$currdirname'/FitHiC'
	# if [[ $UseNonzeroContacts == 1 ]]; then
	# 	GenFitHiCDir=$GenFitHiCDir'_NonZeroCnt'
	# fi

	if [ $BiasCorr == 1 ]; then
		GenFitHiCDir=$GenFitHiCDir'_BiasCorr'
		# if [[ $BeginBiasFilter == 1 ]]; then
		# 	GenFitHiCDir=$GenFitHiCDir'_'$biaslowthr'_'$biashighthr'_b'$BeginBiasFilter
		# fi
		# if [[ $MultBias == 1 ]]; then
		# 	GenFitHiCDir=$GenFitHiCDir'_Mult_'$MultBias
		# else 
		# 	# name is appended with the bias correction regression parameters
		# 	GenFitHiCDir=$GenFitHiCDir'_Resid_'$resid_biascorr'_EqOcc_'$eqocc_biascorr
		# fi
	fi

	# if [[ $DistrType == 2 ]]; then
	# 	GenFitHiCDir=$GenFitHiCDir'/NegBinomDistr'
	# else
	# 	GenFitHiCDir=$GenFitHiCDir'/BinomDistr'
	# fi
	echo '============== ************* Current directory for FitHiChIP significant interaction calling: '$GenFitHiCDir
	mkdir -p $GenFitHiCDir

	#====================================
	# write the configuration in a text file
	outtext=$GenFitHiCDir'/configuration.txt'

	echo "Configurations used for FitHiC execution: " > $outtext	
	echo "Interaction type (1: peak to peak, 2: peak to non peak, 3: peak to all, 4: all to all, 5: all of 1 to 4): $IntType " >> $outtext
	echo "Current interaction type: $CurrIntType " >> $outtext
	# if [[ $DistrType == 1 ]]; then
	# 	echo "Binomial distribution is employed" >> $outtext
	# else
	# 	echo "Negative Binomial distribution is employed" >> $outtext
	# fi
	echo "FitHiC with equal occupancy binning" >> $outtext
	# echo "Using non zero contacts only: $UseNonzeroContacts " >> $outtext
	echo "Lower distance threshold (for significant interaction calling): $LowDistThres " >> $outtext
	echo "Upper distance threshold (for significant interaction calling): $UppDistThres " >> $outtext
	echo "FDR threshold: $QVALUE " >> $outtext
	echo "Number of bins for FitHiC spline modeling: $NBins " >> $outtext
	echo "Bias correction ? (1: yes, 0: no) : $BiasCorr " >> $outtext
	echo "Bias type employed (1: coverage specific, 2: ICE bias): $BiasType " >> $outtext
	
	if [ $BiasCorr == 1 ]; then
		# if [ $BeginBiasFilter == 1 ]; then
		# 	echo "Filtering interactions at beginning based on bias values: bias low threshold: "$biaslowthr"  bias high threshold: "$biashighthr >> $outtext
		# fi
		# if [ $MultBias == 1 ]; then
		# 	echo "Bias correction - Multiplying the probabilities with the bias values " >> $outtext
		# else
		echo "Bias correction - Regression model using the observed contact count and the bias values " >> $outtext
		# echo "Modeling residual contacts for regression: $resid_biascorr" >> $outtext
		# echo "Modeling equal occupancy bins for regression: $eqocc_biascorr" >> $outtext
		# fi
	fi
	if [ $UseP2PBackgrnd == 1 ]; then
		if [[ $currdirname == $DirPeaktoALL ]]; then
			echo "Peak to all interactions - and the background is peak to peak for spline fitting and finding contact significance" >> $outtext
		fi
	fi
	echo "Merging nearby interactions (1: yes, 0: no): $MergeInteraction " >> $outtext

	#====================================

	# files storing FitHiC interactions (significant + all)
	# along with the WashU compatible interactions
	FitHiC_Pass1_outfile=$GenFitHiCDir'/'$PREFIX'.interactions_FitHiC.bed'
	FitHiC_Pass1_Filtfile=$GenFitHiCDir'/'$PREFIX'.interactions_FitHiC_Q'${QVALUE}'.bed'	
	FitHiC_Pass1_LogQ_file=$GenFitHiCDir'/'$PREFIX'.interactions_FitHiC_Q'${QVALUE}'_WashU.bed'
	# FitHiC_Pass1_Filt_PeakCountfile=$GenFitHiCDir'/'$PREFIX'.interactions_FitHiC_Q'${QVALUE}'.PeakSpecificContact.bed'
	# FitHiC_Pass1_PeakCCDistr_Text=$GenFitHiCDir/$PREFIX.anchorPeakCCDistr.bed

	# directory containing the interactions created by merging close contacts
	MergeIntDir=$GenFitHiCDir'/Merge_Nearby_Interactions'
	mkdir -p $MergeIntDir

	FitHiC_Pass1_Filt_MergedIntfile=$MergeIntDir'/'$PREFIX'.interactions_FitHiC_Q'${QVALUE}'_MergeNearContacts.bed'
	FitHiC_Pass1_Filt_MergedInt_LogQ_file=$MergeIntDir'/'$PREFIX'.interactions_FitHiC_Q'${QVALUE}'_MergeNearContacts_WashU.bed'

	# Modeling the statistical significance by FitHiC - main function
	# depending on whether using zero contacts and corresponding pairs of loci
	# two different FitHiC versions are used
	if [[ ! -f $FitHiC_Pass1_outfile || $OverWrite == 1 ]]; then
		
		#===========================
		# comment - sourya
		#===========================
		# if [ $UseZeroCount == 0 ]; then
		# 	$RScriptExec ./src/Interaction.r --InpFile $CurrIntFileSortDist --headerInp --OutFile $FitHiC_Pass1_outfile --BiasCorr $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter --P2P $UseP2PBackgrnd
		# else
		# 	# only for peak to all interactions
		# 	# consider peak to peak background option
		# 	# for all other interaction type, use 0 for this option
		# 	if [ $currdirname == $DirPeaktoALL ]; then
		# 		# comment - sourya
		# 		# $RScriptExec ./src/FitHiC_new.r --InpFile $CurrIntFileSortDist --headerInp --OutFile $FitHiC_Pass1_outfile --CoverageFile $CoverageBiasFile --BinSize $BIN_SIZE --P2P $UseP2PBackgrnd --BiasCorr $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter --IntType $CurrIntType

		# 		# add - sourya
		# 		$RScriptExec ./src/FitHiC_new2.r --InpFile $CurrIntFileSortDist --headerInp --OutFile $FitHiC_Pass1_outfile --CoverageFile $CoverageBiasFile --BinSize $BIN_SIZE --P2P $UseP2PBackgrnd --BiasCorr $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter --IntType $CurrIntType --Resid $resid_biascorr --EqOcc $eqocc_biascorr
		# 	else
		# 		# comment - sourya
		# 		# $RScriptExec ./src/FitHiC_new.r --InpFile $CurrIntFileSortDist --headerInp --OutFile $FitHiC_Pass1_outfile --CoverageFile $CoverageBiasFile --BinSize $BIN_SIZE --P2P 0 --BiasCorr $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter --IntType $CurrIntType

		# 		# add - sourya
		# 		$RScriptExec ./src/FitHiC_new2.r --InpFile $CurrIntFileSortDist --headerInp --OutFile $FitHiC_Pass1_outfile --CoverageFile $CoverageBiasFile --BinSize $BIN_SIZE --P2P 0 --BiasCorr $BiasCorr --Draw --cccol $cccol --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --BiasFilt $BeginBiasFilter --ProbBias $EndBiasFilter --IntType $CurrIntType --Resid $resid_biascorr --EqOcc $eqocc_biascorr
		# 	fi
		# fi
		#===========================
		# end comment - sourya
		#===========================

		#===========================
		# add - sourya
		if [[ $currdirname == $DirPeaktoALL ]]; then			
			$RScriptExec ./src/FitHiC_SigInt.r --InpFile $CurrIntFileSortDist --headerInp --OutFile $FitHiC_Pass1_outfile --CoverageFile $CoverageBiasFile --BinSize $BIN_SIZE --P2P $UseP2PBackgrnd --IntType $CurrIntType --UseNonzeroContacts $UseNonzeroContacts --BiasCorr $BiasCorr --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --Draw --cccol $cccol --BiasFilt $BeginBiasFilter --MultBias $MultBias --Resid $resid_biascorr --EqOcc $eqocc_biascorr
		else
			$RScriptExec ./src/FitHiC_SigInt.r --InpFile $CurrIntFileSortDist --headerInp --OutFile $FitHiC_Pass1_outfile --CoverageFile $CoverageBiasFile --BinSize $BIN_SIZE --P2P 0 --IntType $CurrIntType --UseNonzeroContacts $UseNonzeroContacts --BiasCorr $BiasCorr --BiasLowThr $biaslowthr --BiasHighThr $biashighthr --Draw --cccol $cccol --BiasFilt $BeginBiasFilter --MultBias $MultBias --Resid $resid_biascorr --EqOcc $eqocc_biascorr
		fi
		# end add - sourya
		#===========================
		echo '******** FINISHED calling significant interactions'
	fi

	# error checking condition - 
	# check if FitHiChIP significance is computed at all !!
	if [ ! -f $FitHiC_Pass1_outfile ]; then
		echo "ERROR !!!!!!!! FitHiChIP could not compute the statistical significance of input interactions"
		echo "Check the input parameters, or check if the number of input nonzero contact locus pairs are too few !!!"
		echo "If you are using peak to peak background only (UseP2PBackgrnd=1 or stringent background), check the number of nonzero peak-to-peak locus pairs (peak-to-peak) !!! In such a case, you should go for loose background (UseP2PBackgrnd=0) for modeling interaction significance!!! Specially, if the sequencing depth of your data is very low...."

		errcond=1	# indicating errors
		# exit 1
	fi

	# Filter the interaction file with respect to significance (Q value < $QVALUE)
	# also print the header line
	# if [[ ! -f $FitHiC_Pass1_Filtfile || $OverWrite == 1 ]]; then
		# comment - sourya
		# awk -v q="$QVALUE" '{if ((NR==1) || ($NF < q)) {print $0}}' $FitHiC_Pass1_outfile > $FitHiC_Pass1_Filtfile
		# add - sourya
		# due to strange awk error - possibly due to format conversion error between awk and R
		# also we check whether the field is not NA
		awk -F['\t'] -v q="$QVALUE" '{if ((NR==1) || (($NF != "NA") && (sprintf("%0.400f",$NF) < q))) {print $0}}' $FitHiC_Pass1_outfile > $FitHiC_Pass1_Filtfile
		echo '----- Extracted significant interactions ---- FDR threshold lower than: '$QVALUE
	# fi

	if [ $TimeProf == 1 ]; then
		duration=$(echo "$(date +%s.%N) - $start" | bc)
		echo " ++++ Time (in seconds) for applying FitHiC (significant interactions): $duration" >> $OutTimeFile
		start=$(date +%s.%N)
	fi		

	# no of significant interactions (FitHiC)
	nsigFitHiC=`cat $FitHiC_Pass1_Filtfile | wc -l`

	# error checking condition - 
	# check if there is no significant loops at all
	if [[ $nsigFitHiC -le 1 ]]; then
		echo "SORRY !!!!!!!! FitHiChIP could not find any statistically significant interactions"
		echo "Check the input parameters, or check if the number of input nonzero contact locus pairs are too few !!!"
		echo "If you are using peak to peak background only (UseP2PBackgrnd=1 or stringent background), check the number of nonzero peak-to-peak locus pairs (peak-to-peak) !!! In such a case, you should go for loose background (UseP2PBackgrnd=0) for modeling interaction significance!!! Specially, if the sequencing depth of your data is very low...."

		errcond=1	# indicating errors
		# exit 1
	fi

	# generate distance vs CC plots
	if [[ $nsigFitHiC -gt 1 ]]; then
		DistPlotfile=$GenFitHiCDir'/'$PREFIX'.interactions_FitHiC_Q'${QVALUE}'_Dist_CC.png'	
		$RScriptExec ./Analysis/Distance_vs_CC.r --IntFile $FitHiC_Pass1_Filtfile --OutFile $DistPlotfile	
	fi

	# # Check the no of significant contacts associated with each peak segment
	# # except all to all interactions
	# if [ $currdirname != $DirALLtoALL ]; then		
	# 	if [[ ! -f $FitHiC_Pass1_Filt_PeakCountfile || $OverWrite == 1 ]]; then
	# 		# At least 10 significant interactions are required (empirical threshold)
	# 		# considering the 1st line as header
	# 		if [[ $nsigFitHiC -gt 11 ]]; then
	# 			# skip the 1st header line
	# 			awk -F['\t'] 'NR>1' $FitHiC_Pass1_Filtfile | cut -f1-3,${cccol} | sort -k1,1 -k2,2n -k3,3n | awk -F['\t'] -v OFS='\t' '{a[$1" "$2" "$3]+=$4}END{for (i in a){print i,a[i]}}' - > $FitHiC_Pass1_Filt_PeakCountfile
	# 			echo 'Extracted significant interactions associated with each peak segment'
	# 		else
	# 			echo 'number of significant interactions for spline distribution < 10 - skip the peak count distribution function'
	# 		fi
	# 	fi
	# fi

	if [ $DrawFig == 1 ]; then
		# the R file takes the spline fitted interaction file (without q-value based filtering)
		# and returns the contact count distribution for two different sets of interactions
		# separated by the Q value threshold of 0.01
		# check for non empty interactions file
		if [[ $nsigFitHiC -gt 1 ]]; then
			$RScriptExec ./Analysis/result_summary.r $FitHiC_Pass1_outfile $cccol $QVALUE
			echo 'derived contact count distribution for significant and non-significant interactions'
		else
			echo 'Number of significant spline interaction <= 1 - no result summary'
		fi
	fi

	# the filtered interaction (with respect to the spline) file is used to create a session file
	# for applying in WashU epigenome browser
	# for that, a special file containing only the interacting chromosome intervals 
	# and the log of Q value is created

	# comment - sourya
	# if [[ ! -f $FitHiC_Pass1_LogQ_file || $OverWrite == 1 ]]; then
		
		# check for non empty interactions file
		if [[ $nsigFitHiC -gt 1 ]]; then
			# print the midpoints of the interactions
		    # for better display in the WashU genome browser
		    # check: if Q score is 0, we replace the log of Q score by a default value of 100

		    # old code - sourya
		    # compatible with the old WashU browser 
		    # awk -F['\t'] '{if (NR > 1) {if ($NF > 0) {print $1","(($2+$3)/2-1)","(($2+$3)/2+1)"\t"$4","(($5+$6)/2-1)","(($5+$6)/2+1)"\t"(-log($NF)/log(10))} else {print $1","(($2+$3)/2-1)","(($2+$3)/2+1)"\t"$4","(($5+$6)/2-1)","(($5+$6)/2+1)"\t500"}}}' $FitHiC_Pass1_Filtfile > $FitHiC_Pass1_LogQ_file

		    # new code - sourya
		    # compatible with the new WashU browser (or generally applicable)
		    # converts into a tabix formatted gzipped file
		    awk -F['\t'] '{if (NR > 1) {if ($NF > 0) {print $1"\t"(($2+$3)/2-1)"\t"(($2+$3)/2+1)"\t"$4":"(($5+$6)/2-1)"-"(($5+$6)/2+1)","(-log($NF)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"(($2+$3)/2-1)"\t"(($2+$3)/2+1)"\t"$4":"(($5+$6)/2-1)"-"(($5+$6)/2+1)",500\t"(NR-1)"\t."}}}' $FitHiC_Pass1_Filtfile | sort -k1,1 -k2,2n > $FitHiC_Pass1_LogQ_file
		    if [ -f $FitHiC_Pass1_LogQ_file'.gz' ]; then
		    	rm $FitHiC_Pass1_LogQ_file'.gz'
	    	fi
		    bgzip $FitHiC_Pass1_LogQ_file
		    tabix -f -p bed $FitHiC_Pass1_LogQ_file'.gz'
			echo 'generated WashU epigenome browser compatible significant interactions'
		else
			echo 'There is no significant interaction - so no WashU specific session file is created !!'
		fi

	# comment - sourya
	# fi

	# if [[ $DrawFig == 1 && $currdirname != $DirALLtoALL ]]; then
	# 	res_outdir=$GenFitHiCDir'/Results'
	# 	mkdir -p $res_outdir
	# 	if [ -f $FitHiC_Pass1_Filt_PeakCountfile ]; then
	# 		# Distribution of significant contact counts (FitHiC) for individual peaks
	# 		if [ ! -f $res_outdir/$PREFIX.PeakCCDistr_FiltSpline.pdf ] || [ ! -f $FitHiC_Pass1_PeakCCDistr_Text ]; then
	# 			$RScriptExec ./Analysis/ContactCountDistr.r $PeakFILE $FitHiC_Pass1_Filt_PeakCountfile $res_outdir/$PREFIX.PeakCCDistr_FiltSpline.pdf $FitHiC_Pass1_PeakCCDistr_Text
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

	if [[ $MergeInteraction == 1 && $nsigFitHiC -gt 1 ]]; then
		if [[ ! -f $FitHiC_Pass1_Filt_MergedIntfile || $OverWrite == 1 ]]; then
			# create merged interactions - connected component based analysis
			$PythonExec ./src/CombineNearbyInteraction.py --InpFile $FitHiC_Pass1_Filtfile --OutFile $FitHiC_Pass1_Filt_MergedIntfile --headerInp 1 --binsize $BIN_SIZE --percent 100 --Neigh 2
			echo '-------- *** Merged filtering option is true'
			echo '----- Applied merged filtering (connected component model) on the adjacent loops of FitHiChIP'

			nint=`cat $FitHiC_Pass1_Filt_MergedIntfile | wc -l`

			# error checking condition - 
			# check if there is no significant loops at all
			if [[ $nint -le 1 ]]; then
				echo "SORRY !!!!!!!! FitHiChIP could not find any statistically significant interactions after applying merge filtering on the generated set of loops !!"
				echo "Option 1: use significant loops without merge filtering"
				echo "Option 2: If the number of significant loops (without merge filtering) is also very low, Check the input parameters, or check if the number of input nonzero contact locus pairs are too few !!!"
				echo "If you are using peak to peak background only (UseP2PBackgrnd=1 or stringent background), check the number of nonzero peak-to-peak locus pairs (peak-to-peak) !!! In such a case, you should go for loose background (UseP2PBackgrnd=0) for modeling interaction significance!!! Specially, if the sequencing depth of your data is very low...."

				errcond=1	# indicating errors
				# exit 1
			fi

			if [[ $nint -gt 1 ]]; then
				# 9th field stores the Q value
				# print the midpoints of the interactions
			    # for better display in the WashU genome browser
			    # check: if Q score is 0, we replace the log of Q score by a default value of 100

			    # old code - sourya
		    	# compatible with the old WashU browser
			    # awk -F['\t'] '{if (NR > 1) {if ($9 > 0) {print $1","(($2+$3)/2-1)","(($2+$3)/2+1)"\t"$4","(($5+$6)/2-1)","(($5+$6)/2+1)"\t"(-log($9)/log(10))} else {print $1","(($2+$3)/2-1)","(($2+$3)/2+1)"\t"$4","(($5+$6)/2-1)","(($5+$6)/2+1)"\t500"}}}' $FitHiC_Pass1_Filt_MergedIntfile > $FitHiC_Pass1_Filt_MergedInt_LogQ_file

			    # new code - sourya
			    # compatible with the new WashU browser (or generally applicable)
			    # converts into a tabix formatted gzipped file
			    awk -F['\t'] '{if (NR > 1) {if ($9 > 0) {print $1"\t"(($2+$3)/2-1)"\t"(($2+$3)/2+1)"\t"$4":"(($5+$6)/2-1)"-"(($5+$6)/2+1)","(-log($9)/log(10))"\t"(NR-1)"\t."} else {print $1"\t"(($2+$3)/2-1)"\t"(($2+$3)/2+1)"\t"$4":"(($5+$6)/2-1)"-"(($5+$6)/2+1)",500\t"(NR-1)"\t."}}}' $FitHiC_Pass1_Filt_MergedIntfile | sort -k1,1 -k2,2n > $FitHiC_Pass1_Filt_MergedInt_LogQ_file
			    if [ -f $FitHiC_Pass1_Filt_MergedInt_LogQ_file'.gz' ]; then
			    	rm $FitHiC_Pass1_Filt_MergedInt_LogQ_file'.gz'
		    	fi		    
			    bgzip $FitHiC_Pass1_Filt_MergedInt_LogQ_file
			    tabix -f -p bed $FitHiC_Pass1_Filt_MergedInt_LogQ_file'.gz'

				# generate distance vs CC plots
				DistPlotfile=$MergeIntDir'/'$PREFIX'.interactions_FitHiC_Q'${QVALUE}'_MergeNearContacts_Dist_CC.png'
				$RScriptExec ./Analysis/Distance_vs_CC.r --IntFile $FitHiC_Pass1_Filt_MergedIntfile --OutFile $DistPlotfile	

			else
				echo 'There is no significant interaction - so no WashU specific session file is created !!'
			fi
			echo 'Merged filtering significant interactions - created washu browser compatible file for these interactions!!!'
		fi
	fi

	# increment the counter
	CurrIntType=`expr $CurrIntType + 1`
	echo "Updated CurrIntType: "$CurrIntType

done 	# end of different types of interaction traversal loop

#=================================
# call a separate script which summarizes all the output files and their location
echo -e "\n **** Now summarizing FitHiChIP results in the HTML file *** \n"
bash ./Analysis/SummaryFitHiChIP.sh $OutDir $BIN_SIZE $LowDistThres $UppDistThres $IntType $BiasCorr $BiasType $UseP2PBackgrnd $MergeInteraction $PREFIX ${QVALUE}

#=================================

if [ $errcond == 1 ]; then
	echo "***** There are one or more errors / notifications (such as zero interaction count) associated with this execution !!!"
	echo "***** Please check the console to know about them !!"
else
	echo "***** FitHiChIP pipeline is completely executed - congratulations !!! *****"
fi

#============================
# sourya - now go back to the original working directory
#============================
cd $currworkdir
