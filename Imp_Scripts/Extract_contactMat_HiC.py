#!/usr/bin/python3

##==========
## script to extract the contact matrix from the .hic file
## first install the package hic-straw using the following command
## pip install hic-straw
## pip install numpy
##==========

import numpy as np
import hicstraw
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-i", "--hicfile", dest="InpHiCFile", help="Input HiC File")
parser.add_option("-o", "--outfile", dest="OutFile", help="Output File compatible for FitHiChIP")
parser.add_option("-c", "--chromsizefile", dest="chromsizefile", help="Chromosome size file for the reference genome")
parser.add_option("-r", "--resolution", action="store", dest="Resolution", default=5000, help="Target contact matrix resolution. Default = 5000 (5 Kb)")

(options, args) = parser.parse_args()


# ## input HiC file
# InpDir = "/home/sourya"
# InpHiCFile = InpDir + "/GSM2487542_H3K27ac_HiChIP.hic"

# ## output file
# OutFile = InpDir + "/GSM2487542_H3K27ac_HiChIP.bed"

# ## bin size (target resolution)
# TARGETRES = 5000

InpHiCFile = options.InpHiCFile
OutFile = options.OutFile
TARGETRES = int(options.Resolution)

print("Input -  InpHiCFile : " + InpHiCFile)
print("Input -  OutFile : " + OutFile)
print("Input -  TARGETRES : " + str(TARGETRES))

## first read the input .hic file
hic = hicstraw.HiCFile(InpHiCFile)

## print the genome
genomeID = hic.getGenomeID()
print("genomeID : ", genomeID)

## print the resolutions
resolution_values = hic.getResolutions()
print(resolution_values)
if TARGETRES not in resolution_values:
	TARGETRES = min(resolution_values)
	print(" modified TARGETRES : " + str(TARGETRES))

## chromosomes
chromlist = hic.getChromosomes()

## dump the chromosome list and chromosome length
ChromListFile = options.chromsizefile
with open(ChromListFile, 'w') as f:
	for chrom in chromlist:
		if chrom.name != "All":
			currline = str(chrom.name) + "\t" + str(chrom.length)
			f.write(currline)
			f.write('\n')

## now dump the HiC matrix
with open(OutFile, 'w') as f:
	for chrom in chromlist:
		if chrom.name != "All":
			print("processing chromosome : " + str(chrom.name))
			result = hicstraw.straw('observed', 'NONE', InpHiCFile, str(chrom.name), str(chrom.name), 'BP', TARGETRES)
			for i in range(len(result)):
				currline = str(chrom.name) + "\t" + str(int(result[i].binX)) + "\t" + str(int(result[i].binX + TARGETRES)) + "\t" + str(chrom.name) + "\t" + str(int(result[i].binY)) + "\t" + str(int(result[i].binY + TARGETRES)) + "\t" + str(int(result[i].counts))
				# f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(chrom, result[i].binX, (result[i].binX + TARGETRES), chrom, result[i].binY, (result[i].binY + TARGETRES), result[i].counts))    		
				f.write(currline)
				f.write('\n')


