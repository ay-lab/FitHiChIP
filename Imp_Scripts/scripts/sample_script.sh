#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=40GB
#PBS -l walltime=48:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR


#============================================

if [ 1 == 0 ]; then

basedir='/mnt/BioAdHoc/Groups/vd-vijay/sourya/DATA/HiChIP/Merged_Vivek_Feb_March_2018/'

CodeExec='/home/sourya/proj/2017_FitHiChIP_HiCPro/Imp_Scripts/ContactMatDiffRes.sh'

ChrSizeFile='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/chrom_hg19.sizes'

HiCProBasedir='/home/sourya/packages/HiCPro/HiC-Pro_2.9.0/'

# variation for different cell types
#for celltype in 'CD4_Naive' 'CD8_Naive' 'Monocyte' 'B_Cell' 'NK'; do
#for celltype in 'CD8_Naive' 'Monocyte' 'B_Cell' 'NK'; do
for celltype in 'NB' 'NK'; do

	inpvalidpairsfile=$basedir'Merged_'$celltype'/rawdata_allValidPairs.txt.gz'
	BaseOutDir=$basedir'Merged_'$celltype'/Contact_Matrix_Varying_Resolutions'

	# variation for different resolutions
	for targetres in 500000 50000 25000 5000 1000; do
		$CodeExec -V $inpvalidpairsfile -C $ChrSizeFile -B $targetres -H $HiCProBasedir -D $BaseOutDir
	done
done


fi 
#============================================

sh /home/sourya/proj/GitHub_Codes/FitHiChIP_HiCPro/Imp_Scripts/DiffAnalysisHiChIP_Script.sh

