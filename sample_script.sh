#!/bin/bash

# executing FitHiChIP without any bias correction
./FitHiChIP_HiCPro.sh -C configfile

# executing FitHiChIP with bias correction
# -b 1 means that interactions having bias values within the default thresholds (0.2 and 5) will be only considered as an input to FitHiC
# -e 1 means that probability of an interaction (computed using FitHiC) will be multiplied with their bias values
# accordingly, computation of the P values will be adjusted
./FitHiChIP_HiCPro.sh -C configfile -b 1 -e 1

