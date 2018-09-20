#!/bin/bash

#===============================
# used to merge different ChIP seq samples (bam files)
# and then call MACS2 specific peaks
# also it generates the ChIP-seq coverage of individual BAM files
# according to the specified bin size
# and the reference genome specific chromosome size information

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#===============================

usage(){
cat << EOF

Options:
	-I  InpFile        	 	A list of ChIP-seq BAM files. Multiple BAM files are to be mentioned in the format: 
							-I bamfile1.bam -I bamfile2.bam -I bamfile3.bam and so on
   	-D  OutBaseDir			Directory containing the output set of peaks. Default: current directory
   	-R 	refGenome			Reference genome string used for MACS2. Default is 'hs' for human chromosome. For mouse, specify 'mm'
   	-M 	MACS2ParamStr		String depicting the parameters for MACS2. Default: "--nomodel --extsize 147 -q 0.01"
   	-C  ChrSizeFile 		Filename containing the chromosome size information for the reference genome
   	-b 	BinSize				BinSize in base pair. Used to compute the ChIP-seq coverage for this bin size. Default 5000 (means 5 Kb)
   	
EOF
}

#==============
# default parameters
OutBaseDir=`pwd`'/'
refGenome='hs'
MACS2ParamStr='--nomodel --extsize 147 -q 0.01'
BinSize=5000
#==============

while getopts "I:D:R:M:C:b:" opt;
do
	case "$opt" in
		I) InpFile+=($OPTARG);;
		D) OutBaseDir=$OPTARG;;
		R) refGenome=$OPTARG;;
		M) MACS2ParamStr=$OPTARG;;
		C) ChrSizeFile=$OPTARG;;
		b) BinSize=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $InpFile ]]; then
	echo 'User did not provide any input BAM file - exit for the moment !!'
	exit 1
fi

if [[ -z $ChrSizeFile ]]; then
	echo 'User did not provide any file containing size of different chromosomes for the reference genome - exit for the moment !!'
	exit 1
fi

# number of input files provided
nsample=${#InpFile[@]}
echo 'Number of input ChIP-seq BAM files : '$nsample

# generate the output directory
mkdir -p $OutBaseDir

echo 'Sample BAM files which are merged and used for peak calling' > $OutBaseDir'/SampleNames.txt'

# notes down the argument bam files
argstr1=''
for (( i=0; i<${nsample}; i++ ))
do
	echo "Sample: "${InpFile[$i]} >> $OutBaseDir'/SampleNames.txt'
	argstr1=$argstr1' '${InpFile[$i]}
done

# merging input bam files
samtools merge $OutBaseDir'/merged_input.bam' $argstr1 

# first ensure that the bam file contains only the valid chromosomes
inpfile_reduced=$OutBaseDir'/merged_input_validChr.bam'
samtools view -h $OutBaseDir'/merged_input.bam' | awk '(substr($1, 1, 1)=="@") || (( $3 ~ /^chr([1-9]|2[0-2]|1[0-9])$/ ) && ( ( $7 ~ /^chr([1-9]|2[0-2]|1[0-9])$/ ) || ($7=="=") || ($7=="*") ))' - | samtools view -bhS - > $inpfile_reduced

# now call MACS2 peaks
MACS2_OutDir=$OutBaseDir'/MACS2_Out'
MACS2PeakOutFile=$MACS2_OutDir'/out.macs2_peaks.narrowPeak'

# employ the macs2 parameters mentioned in input settings
# Reference genome
# and MACS2 string
# this peak file (with respect to the merged alignment)
macs2 callpeak -t $inpfile_reduced -f BAM -g $refGenome -n 'out.macs2' --outdir $MACS2_OutDir $MACS2ParamStr

# also infer the ChIP-seq coverage for the input alignment (BAM) files 
# (different replicates)
# according to the reference genome chromosomes size 
# and the specified bin size

# first create the bin specific distribution of the reference chromosome
RefBinFile=$OutBaseDir'/RefBin.bed'
awk -v b="$BinSize" '{for (i=0;i<$2;i=i+b) {print $1"\t"i"\t"(i+b)}}' $ChrSizeFile > $RefBinFile

# then for individual input BAM files,
# compute the ChIP-seq coverage
for (( i=0; i<${nsample}; i++ ))
do
	echo "Computing ChIP-seq coverage for the file: "${InpFile[$i]}
	fileno=`expr $i + 1`
	bedtools coverage -a $RefBinFile -b ${InpFile[$i]} -counts > $OutBaseDir'/ChIPCoverage'$fileno'.bed'
done

# # also generate the ChIP-seq coverage of the 
# # merged alignment (bam) file
# bedtools coverage -a $RefBinFile -b $OutBaseDir'/merged_input.bam' -counts > $OutBaseDir'/ChIPCoverage_MERGE.bed'




