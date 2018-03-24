#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -l walltime=48:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

# executing FitHiChIP according to the settings of condifuration file
# any modification of input parameters should be done via the configuration file 
#./FitHiChIP_HiCPro.sh -C configfile_BiasCorrection_CoverageBias
#./FitHiChIP_HiCPro.sh -C configfile_BiasCorrection_ICEBias
./FitHiChIP_HiCPro.sh -C configfile_P2P_BiasCorrection_CoverageBias
#./FitHiChIP_HiCPro.sh -C configfile_P2P_BiasCorrection_ICEBias


