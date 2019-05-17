#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=20GB
#PBS -l walltime=48:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

basedir='/mnt/BioAdHoc/Groups/vd-vijay/sourya/DATA/HiChIP/Merged_Vivek_Feb_March_2018/'

CodeExec='/home/sourya/proj/2017_FitHiChIP_HiCPro/Imp_Scripts/PeakInferHiChIP.sh'

# variation for different cell types
#for celltype in 'CD4_Naive' 'CD8_Naive' 'Monocyte' 'NB' 'NK'; do
for celltype in 'NB' 'NK'; do
	hicprodir=$basedir'Merged_'$celltype'/HiCPro_Reads_MergedAlignment'
	BaseOutDir=$basedir'Merged_'$celltype'/HiChIP_Inferred_Peaks_KeepDupALL_B'
	# MACS2 parameters are default mentioned in the code
	$CodeExec -H $hicprodir -D $BaseOutDir
done


