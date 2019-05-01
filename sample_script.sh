#!/bin/bash

# executing FitHiChIP according to the settings of configuration file
# any modification of input parameters should be done via the configuration file 

# executing FitHiChIP(L) with coverage bias regression
bash FitHiChIP_HiCPro.sh -C configfile_BiasCorrection_CoverageBias

# executing FitHiChIP(L) with ICE bias regression
bash FitHiChIP_HiCPro.sh -C configfile_BiasCorrection_ICEBias

# executing FitHiChIP(S) with coverage bias regression
bash FitHiChIP_HiCPro.sh -C configfile_P2P_BiasCorrection_CoverageBias

# executing FitHiChIP(S) with ICE bias regression
bash FitHiChIP_HiCPro.sh -C configfile_P2P_BiasCorrection_ICEBias


# this command is commented
# user may uncomment this command for differential analysis of FitHiChIP loops
# for two categories each with multiple replicates 
# a sample test data is provided in this repository
# bash Differetial_Analysis_Script.sh


