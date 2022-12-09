#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=20GB
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

# executing FitHiChIP according to the settings of configuration file
# any modification of input parameters should be done via the configuration file 

# executing FitHiChIP(L) with coverage bias regression
bash FitHiChIP_HiCPro.sh -C configfile_BiasCorrection_CoverageBias

# executing FitHiChIP(L) with ICE bias regression
#bash FitHiChIP_HiCPro.sh -C configfile_BiasCorrection_ICEBias

# executing FitHiChIP(S) with coverage bias regression
#bash FitHiChIP_HiCPro.sh -C configfile_P2P_BiasCorrection_CoverageBias

# executing FitHiChIP(S) with ICE bias regression
#bash FitHiChIP_HiCPro.sh -C configfile_P2P_BiasCorrection_ICEBias


# this command is commented
# user may uncomment this command for differential analysis of FitHiChIP loops
# for two categories each with multiple replicates 
# a sample test data is provided in this repository
# bash Differetial_Analysis_Script.sh


