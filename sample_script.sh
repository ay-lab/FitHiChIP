#!/bin/bash

# executing FitHiChIP according to the settings of condifuration file
# any modification of input parameters should be done via the configuration file 
./FitHiChIP_HiCPro.sh -C configfile_BiasCorrection_CoverageBias
./FitHiChIP_HiCPro.sh -C configfile_BiasCorrection_ICEBias
./FitHiChIP_HiCPro.sh -C configfile_P2P_BiasCorrection_CoverageBias
./FitHiChIP_HiCPro.sh -C configfile_P2P_BiasCorrection_ICEBias


