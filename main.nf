#!/usr/bin/env nextflow

params.ValidPairs = "none"
params.Interval = "none2"
params.Matrix = "none3"
params.Bed = "none4"
params.PeakFile = "none5"
params.ChrSizeFile = "none6"
params.IntType = 3
params.BINSIZE = 5000
params.LowDistThr = 20000
params.UppDistThr = 2000000
params.UseP2PBackgrnd = 0
params.BiasType = 1
params.MergeInt = 1
params.QVALUE = 0.01
params.PREFIX = "FitHiChIP"
params.OverWrite = 1

vp = file(params.ValidPairs)
interval = file(params.Interval)
matrix = file(params.Matrix)
bed = file(params.Bed)
peakfile = file(params.PeakFile)
chrsize = file(params.ChrSizeFile)


process prepare_config {
    output:
    file nf_conf into config

    script:

    cstr = "OutDir=./fhcout/\nHiCProBasedir=/HiC-Pro-2.11.1/\n"

    if (vp.exists()) 
        cstr = cstr + "ValidPairs=${vp}\n"
    
    if (interval.exists()) 
        cstr = cstr + "Interval=${interval}\n"
    
    if (matrix.exists()) 
        cstr = cstr + "Matrix=${matrix}\n"
    
    if (bed.exists()) 
        cstr = cstr + "Bed=${bed}\n"
    
    if (peakfile.exists()) 
        cstr = cstr + "PeakFile=${peakfile}\n"

    if (chrsize.exists())
        cstr = cstr + "ChrSizeFile=${chrsize}\n"
    else
        cstr = cstr + "ChrSizeFile=/FitHiChIP/TestData/chrom_hg19.sizes\n"
    
    cstr = cstr + "IntType=${params.IntType}\n"
    cstr = cstr + "BINSIZE=${params.BINSIZE}\n"
    cstr = cstr + "LowDistThr=${params.LowDistThr}\n"
    cstr = cstr + "UppDistThr=${params.UppDistThr}\n"
    cstr = cstr + "UseP2PBackgrnd=${params.UseP2PBackgrnd}\n"
    cstr = cstr + "BiasType=${params.BiasType}\n"
    cstr = cstr + "MergeInt=${params.MergeInt}\n"
    cstr = cstr + "QVALUE=${params.QVALUE}\n"
    cstr = cstr + "PREFIX=${params.PREFIX}\n"
    cstr = cstr + "OverWrite=${params.OverWrite}\n"

    """
    echo "$cstr" >> nf_conf
    """

}

process main {

    input:
    file conf from config
    file vpf from vp
    file intf from interval
    file mf from matrix
    file bedf from bed
    file peakfilef from peakfile
    file chrsizef from chrsize

    output:
    file "./fhcout/**/*" into outFiles
    stdout result

    script:

    """
    bash /FitHiChIP/FitHiChIP_HiCPro.sh -C $conf
    """
}

result.println { it.trim() }
