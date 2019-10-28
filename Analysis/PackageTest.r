#!/usr/bin/env Rscript

#===========================================================
# R script for checking whether the required packages are installed or not
#===========================================================

#specify the packages of interest
packages = c("optparse","data.table","ggplot2","splines","fdrtool","parallel","tools","GenomicRanges","dplyr","EdgeR")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
    	# install the package
        install.packages(x, dependencies = TRUE)
        # loading the library
        # library(x, character.only = TRUE)
    }
})