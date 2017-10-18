#!/usr/bin/env python

"""
Created: September 21,2017

This file takes input a HiC-pro pipeline output of validpairs file, Chromosome length file, and a peak detection file
it generates binned intervals from the input chromosomes
and creates coverage vector for individual chromosomes
Dependin on the input peak file, it also marks individual binned segments as either peak or non-peak

Author: Sourya Bhattacharyya
Vijay-Ay lab, LJI
"""

import sys
import os
from optparse import OptionParser

# used to open the gzip compressed text file
import gzip 

#==========================================================
""" 
this class defines a particular chromosome interval
indexed by (chromosome, start, end) - dictionary key 
"""
class ChrmSegment(object):  
    def __init__(self):
        self.is_peak = False
        self.read_count = 0

    def _CheckIfPeak(self):
        return self.is_peak

    def _GetReadCount(self):
        return self.read_count

    def _IncrementRead(self, count=1):
        self.read_count = self.read_count + count

    def _SetPeak(self):
        self.is_peak = True

#===============================================
def main():
    parser = OptionParser() #(usage=usage)
    parser.add_option("-p", "--peakfile", dest="peakfile", help="Peak detection file")
    parser.add_option("-i", "--inpfile", dest="inpfile", help="Input valid pairs file (from HiC-pro pipeline)")
    parser.add_option("-b", "--binsize", dest="binsize", type="int", help="Size of bins employed. DEFAULT 5000 (indicating 5 Kb).")
    parser.add_option("-o", "--outfile", dest="outfile", help="output file storing the coverage of individual genomic bins")
    parser.add_option("-c", "--chrsizefile", dest="chrsizefile", help="File containing chromosome size information")
    parser.set_defaults(peakfile=None, inpfile=None, binsize=5000, outfile=None, chrsizefile=None)
    (options, args) = parser.parse_args()

    # global variables from command line options
    global inpfile
    global peakfile
    global outfile
    global bin_size
    global chrsizefile

    # dictionary structure which stores the instances of the class ChrmSegment
    global SegmentDict
    SegmentDict = dict()

    # two lists which store the chromosome name and their length information
    # with respect to the input chromosome length file
    ChrNameList_Ref = []
    ChrLenList_Ref = []

    # list storing the chromosome names provided in the input validpairs file
    ChrNameList_Current = [] 

    if options.inpfile is not None:
        inpfile = options.inpfile
    else:
        sys.exit("Input validpairs file (output of HiC-pro pipeline) is not provided - quit !!")

    if options.peakfile is not None:
        peakfile = options.peakfile
    else:
        sys.exit("Peak detection file is not provided - quit !!")

    if options.outfile is not None:
        outfile = options.outfile
    else:
        sys.exit("Output file for storing the coverage of segments is not provided - quit !!")

    if options.chrsizefile is not None:
        chrsizefile = options.chrsizefile
    else:
        sys.exit("Chromosome size file is not provided - quit !!")

    # input bin size parameter  
    bin_size = int(options.binsize)

    #=====================================================
    # first read the reference chromosome length file
    # and assign the name and lengths of individual chromosomes
    #=====================================================
    with open(chrsizefile, 'r') as fp:
        for line in fp:
            linecontents = (line.rstrip()).split()
            curr_chr = linecontents[0]
            curr_chr_size = int(linecontents[1])
            ChrNameList_Ref.append(curr_chr)
            ChrLenList_Ref.append(curr_chr_size)

    #=====================================================
    # now read the input valid pairs file
    # and compute the coverage information
    # with respect to 'cis' reads
    #=====================================================
    with gzip.open(inpfile,'r') as fin:    
        for line in fin:
            contents = line.split()
            # check cis interaction
            if (str(contents[1]) == str(contents[4])):
                chrname = str(contents[1])
                read1_pos = int(contents[2])
                read2_pos = int(contents[5])
                
                # check if the chromosome is not in the 
                # already processed list
                # in such a case, include this chromosome
                # and also create the binned intervals
                if chrname not in ChrNameList_Current:
                    # this chromosome is not yet processed
                    ChrNameList_Current.append(chrname)
                    # create python dictionary 
                    # and the bin size specific intervals for this chromosome
                    chrsize = ChrLenList_Ref[ChrNameList_Ref.index(chrname)]
                    if ((chrsize % bin_size) == 0):
                        interval_end = chrsize
                    else:
                        interval_end = (int(chrsize / bin_size) + 1) * bin_size
                    for val in range(0, interval_end, bin_size):
                        curr_key = (chrname, val, (val + bin_size))
                        SegmentDict.setdefault(curr_key, ChrmSegment())
                
                # now process the current paired end read
                # read position 1
                bin_start1 = (int(read1_pos / bin_size)) * bin_size
                # increment the read count in the corresponding dictionary entry of this bin
                curr_key1 = (chrname, bin_start1, (bin_start1 + bin_size))
                SegmentDict[curr_key1]._IncrementRead()
                # read position 2
                bin_start2 = (int(read2_pos / bin_size)) * bin_size
                # increment the read count in the corresponding dictionary entry of this bin
                curr_key2 = (chrname, bin_start2, (bin_start2 + bin_size))
                SegmentDict[curr_key2]._IncrementRead()
    
    #=====================================================
    # scan the peak input file (if provided) and mark the corresponding 
    # chromosome intervals as 1 (having peak)
    #=====================================================
    with open(peakfile, 'r') as fp:
        for line in fp:
            linecontents = (line.rstrip()).split()
            curr_chr = linecontents[0]
            peak_start = int(linecontents[1])
            peak_end = int(linecontents[2])
            # only process those peaks whose chromosome 
            # is present in the current valid pairs file
            # this is possible when the peak information is downloaded from a reference
            if curr_chr in ChrNameList_Current:
                interval_start = (int(peak_start / bin_size)) * bin_size
                if ((peak_end % bin_size) == 0):
                    interval_end = peak_end
                else:
                    interval_end = (int(peak_end / bin_size) + 1) * bin_size
                for val in range(interval_start, interval_end, bin_size):
                    curr_key = (curr_chr, val, (val + bin_size))
                    # mark the corresponding dictionary entry as a peak segment
                    # since it has an overlap with the MACS2 derived peak intervals
                    SegmentDict[curr_key]._SetPeak()
    
    # now write the genomic intervals and corresponding read count + peak information
    fp_cov = open(outfile, 'w')
    fp_cov.write('Chr' + '\t' + 'Start' + '\t' + 'End' + '\t' + 'Coverage' + '\t' + 'IsPeak')
    for i in range(len(ChrNameList_Current)):
        curr_chr = ChrNameList_Current[i]
        curr_chr_size = ChrLenList_Ref[ChrNameList_Ref.index(curr_chr)]
        if ((curr_chr_size % bin_size) == 0):
            interval_end = curr_chr_size
        else:
            interval_end = (int(curr_chr_size / bin_size) + 1) * bin_size
        for val in range(0, interval_end, bin_size):
            curr_key = (curr_chr, val, (val + bin_size))
            fp_cov.write('\n' + str(curr_chr) + '\t' + str(val) + '\t' + str(val + bin_size) + '\t' + str(SegmentDict[curr_key]._GetReadCount()) + '\t' + str(int(SegmentDict[curr_key]._CheckIfPeak())))
    fp_cov.close()

    return

#===============================================
if __name__ == "__main__":
    main()
