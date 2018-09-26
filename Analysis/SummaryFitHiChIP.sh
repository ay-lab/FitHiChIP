#!/bin/bash

#===============
# script to generate a summary list of output files
# and their links in a single .html file
#===============

# scan the input parameters

OutDir=$1
BIN_SIZE=$2
LowDistThres=$3
UppDistThres=$4
IntType=$5
BiasCorr=$6
BiasType=$7
UseP2PBackgrnd=$8
MergeInteraction=$9
PREFIX=${10}
QVALUE=${11}

# output .html file
html=$OutDir/Summary_results_FitHiChIP.html

echo \<html\> > $html
echo \<head\>\<title\> FitHiChIP summary Report \</title\>\</head\> >>$html

echo \<body\> >> $html
echo \<center\>\<h1\>FitHiChIP summary Report \</h1\>\</center\>  >>$html

echo \<p\> \<strong\>Parameter list for the current execution of FitHiChIP:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/Parameters.txt' \</font\> \</span\> \</p\> >>$html

if [ -f $OutDir'/TimingProfile.txt' ]; then
	echo \<p\> \<strong\>Timing profile for the current execution of FitHiChIP:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/TimingProfile.txt' \</font\> \</span\> \</p\> >>$html
fi

FileIntDistThr=$OutDir'/HiCPro_Matrix_BinSize'$BIN_SIZE'/L_'$LowDistThres'_U'$UppDistThres'/'$PREFIX'.cis.interactions.DistThr.bed'
if [ -f $FileIntDistThr ]; then
	echo \<p\> \<strong\>List of CIS interactions within the specified distance thresholds \(lower threshold: $LowDistThres upper threshold: $UppDistThres \): \</strong\> \<span class=\"tab\"\> \<font color="blue"\> $FileIntDistThr \</font\> \</span\> \</p\> >>$html
	numIntDistThr=`cat $FileIntDistThr | wc -l`
	numIntDistThr=`expr $numIntDistThr - 1`
	echo \<span class=\"tab\"\> \<font color="red"\> Number of interacting \(non zero contacts\) bin pairs \( all to all \) : $numIntDistThr \</font\> \</span\> >>$html
fi

if [ -f $OutDir'/NormFeatures/Coverage_Bias/'$PREFIX'.coverage_Bias.bed' ]; then
	echo \<p\> \<strong\>Bin specific HiChIP coverage, coverage bias, and peak / non-peak boolean values:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/NormFeatures/Coverage_Bias/'$PREFIX'.coverage_Bias.bed' \</font\> \</span\> \</p\> >>$html
	if [ -f $OutDir'/NormFeatures/Coverage_Bias/Plots_Bias/Peak_Bias.pdf' ]; then
		echo \<p\> \<strong\>Variation of coverage bias for peaks:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/NormFeatures/Coverage_Bias/Plots_Bias/Peak_Bias.pdf' \</font\> \</span\> \</p\> >>$html
	fi
	if [ -f $OutDir'/NormFeatures/Coverage_Bias/Plots_Bias/NonPeak_Bias.pdf' ]; then
		echo \<p\> \<strong\>Variation of coverage bias for the non-peaks:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/NormFeatures/Coverage_Bias/Plots_Bias/NonPeak_Bias.pdf' \</font\> \</span\> \</p\> >>$html
	fi
fi

if [ -f $OutDir'/NormFeatures/ICE_Bias/'$PREFIX'.coverage_ICE_Bias.bed' ]; then
	echo \<p\> \<strong\>Bin specific HiChIP coverage, ICE bias, peak / non-peak boolean values:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/NormFeatures/ICE_Bias/'$PREFIX'.coverage_ICE_Bias.bed' \</font\> \</span\> \</p\> >>$html
	echo \<p\> \<strong\>ICE normalized contact matrix:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/NormFeatures/ICE_Bias/'$PREFIX'.norm.Contact.Matrix' \</font\> \</span\> \</p\> >>$html
	if [ -f $OutDir'/NormFeatures/ICE_Bias/Plots_Bias/Peak_Bias.pdf' ]; then
		echo \<p\> \<strong\>Variation of ICE bias for peaks:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/NormFeatures/ICE_Bias/Plots_Bias/Peak_Bias.pdf' \</font\> \</span\> \</p\> >>$html
	fi
	if [ -f $OutDir'/NormFeatures/ICE_Bias/Plots_Bias/NonPeak_Bias.pdf' ]; then
		echo \<p\> \<strong\>Variation of ICE bias for the non-peaks:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $OutDir'/NormFeatures/ICE_Bias/Plots_Bias/NonPeak_Bias.pdf' \</font\> \</span\> \</p\> >>$html
	fi
fi


for looptype in 'Peak2ALL' 'Peak2Peak' 'Peak2NonPeak' 'ALL2ALL'; do

	echo \<br\> ====================================================================================== \</br\> >>$html

	for p2p in 0 1; do

		# only peak to all loop has significance of p2p 0 or 1
		if [[ $p2p == 1 && $looptype != 'Peak2ALL' ]]; then
			continue
		fi

		CurrDir=$OutDir'/FitHiChIP_'$looptype'_b'$BIN_SIZE'_L'$LowDistThres'_U'$UppDistThres
		if [ -d $CurrDir ]; then		
			
			if [ $looptype == 'Peak2ALL' ]; then
				echo \<p\> \<h2\> Statistics for peak to ALL \(both peak and non-peak\) loops \</h2\>\</p\> >>$html
				echo \<p\> \<h3\> FitHiChIP background \(0 means loose, 1 means stringent \) : $p2p \</h3\>\</p\> >>$html
			fi
			if [ $looptype == 'Peak2Peak' ]; then
				echo \<p\> \<h2\> Statistics for peak to peak loops \</h2\>\</p\> >>$html
			fi
			if [ $looptype == 'Peak2NonPeak' ]; then
				echo \<p\> \<h2\> Statistics for peak to non-peak loops \</h2\>\</p\> >>$html
			fi
			if [ $looptype == 'ALL2ALL' ]; then
				if [[ -d $CurrDir'/Coverage_Bias/FitHiC_BiasCorr/' || -d $CurrDir'/ICE_Bias/FitHiC_BiasCorr/' ]]; then
					echo \<p\> \<h2\> Statistics for ALL to ALL \(every possible loops, similar to Hi-C\) \</h3\>\</p\> >>$html
				fi
			fi

			for biastype in 'Coverage_Bias' 'ICE_Bias'; do

				echo \<br\> ==================================================== \</br\> >>$html

				if [ $looptype == 'Peak2ALL' ]; then
					CurrBiasDir=$CurrDir'/P2PBckgr_'$p2p'/'$biastype'/FitHiC_BiasCorr/'
				else
					CurrBiasDir=$CurrDir'/'$biastype'/FitHiC_BiasCorr/'
				fi

				if [ -f $CurrBiasDir$PREFIX'.interactions_FitHiC.bed' ]; then

					if [ $biastype == 'Coverage_Bias' ]; then
						echo \<p\> \<h3\> Statistics for Coverage bias regression \</h3\>\</p\> >>$html
					fi
					if [ $biastype == 'ICE_Bias' ]; then
						echo \<p\> \<h3\> Statistics for ICE bias regression \</h3\>\</p\> >>$html
					fi

					# all loops with FitHiChIP significance
					echo \<p\> \<strong\> All interactions with their FitHiChIP significance values:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $CurrBiasDir$PREFIX'.interactions_FitHiC.bed' \</font\> \</span\> \</p\> >>$html

					numInt=`cat $CurrBiasDir$PREFIX'.interactions_FitHiC.bed' | wc -l`
					numInt=`expr $numInt - 1`
					echo \<span class=\"tab\"\> \<font color="red"\> Number of interactions considered: $numInt \</font\> \</span\> >>$html	

					# only the significant loops
					if [ -f $CurrBiasDir$PREFIX'.interactions_FitHiC_Q'$QVALUE'.bed' ]; then
						echo \<p\> \<strong\> Significant interactions:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $CurrBiasDir$PREFIX'.interactions_FitHiC_Q'$QVALUE'.bed' \</font\> \</span\> \</p\> >>$html

						numInt=`cat $CurrBiasDir$PREFIX'.interactions_FitHiC_Q'$QVALUE'.bed' | wc -l`
						numInt=`expr $numInt - 1`
						echo \<span class=\"tab\"\> \<font color="red"\> Total number of significant interactions \( FDR $QVALUE \) : $numInt \</font\> \</span\>  >>$html
					fi

					# WashU compatible files
					if [ -f $CurrBiasDir$PREFIX'.interactions_FitHiC_Q'$QVALUE'_WashU.bed' ]; then
						echo \<p\> \<strong\> File containing significant loops for visualization in the WashU genome browser: \</strong\> \<span class=\"tab\"\> \<font color="blue"\> $CurrBiasDir$PREFIX'.interactions_FitHiC_Q'$QVALUE'_WashU.bed' \</font\> \</span\> \</p\> >>$html
					fi

					# check if merging adjacent loops are enabled
					if [ -f $CurrBiasDir'Merge_Nearby_Interactions/'$PREFIX'.interactions_FitHiC_Q'$QVALUE'_MergeNearContacts.bed' ]; then
						echo \<p\> \<strong\> After merging adjacent loops \(connected component modeling\) significant interactions:\</strong\> \<span class=\"tab\"\> \<font color="blue"\> $CurrBiasDir'Merge_Nearby_Interactions/'$PREFIX'.interactions_FitHiC_Q'$QVALUE'_MergeNearContacts.bed' \</font\> \</span\> \</p\> >>$html

						numInt=`cat $CurrBiasDir'Merge_Nearby_Interactions/'$PREFIX'.interactions_FitHiC_Q'$QVALUE'_MergeNearContacts.bed' | wc -l`
						numInt=`expr $numInt - 1`
						echo \<span class=\"tab\"\> \<font color="red"\> Total number of significant interactions \( FDR $QVALUE \) after merging adjacent loops: $numInt \</font\> \</span\>  >>$html
					fi

					# check if merging adjacent loops are enabled
					if [ -f $CurrBiasDir'Merge_Nearby_Interactions/'$PREFIX'.interactions_FitHiC_Q'$QVALUE'_MergeNearContacts_WashU.bed' ]; then
						echo \<p\> \<strong\> File containing significant interactions \(After merging adjacent loops\) for visualization in the WashU genome browser: \</strong\> \<span class=\"tab\"\> \<font color="blue"\> $CurrBiasDir'Merge_Nearby_Interactions/'$PREFIX'.interactions_FitHiC_Q'$QVALUE'_MergeNearContacts_WashU.bed' \</font\> \</span\> \</p\> >>$html
					fi	

					echo \<table  border=\"0\"\> >> $html
					echo \<tbody\> >> $html

					# check the spline fitting plot
					# pdffile=$CurrBiasDir'Plots/EqOccBin_SplinePass1.pdf'
					# if [ -f $pdffile ]; then
					# 	pngfile="${pdffile%.*}".png
					# 	if [ ! -f $pngfile ]; then
					# 		convert -verbose -density 500 -resize '800' "${pdffile}" "${pngfile}"
					# 	fi
					# 	imgdata=$( base64 $pngfile )
					# 	echo "<td> <a href=\"$pngfile\"><img src=\"data:image/png;base64,$imgdata\" width=500> </a></td>" >> $html
					# fi
					pngfile=$CurrBiasDir'Plots/EqOccBin_SplinePass1.png'
					if [ -f $pngfile ]; then
						imgdata=$( base64 $pngfile )
						echo "<td> <a href=\"$pngfile\"><img src=\"data:image/png;base64,$imgdata\" width=500> </a></td>" >> $html
					fi

					# check the distance vs significant loop count plot
					# pdffile=$CurrBiasDir$PREFIX'.interactions_FitHiC_Q'$QVALUE'_Dist_CC.pdf'
					# if [ -f $pdffile ]; then
					# 	pngfile="${pdffile%.*}".png
					# 	if [ ! -f $pngfile ]; then
					# 		convert -verbose -density 500 -resize '800' "${pdffile}" "${pngfile}"
					# 	fi
					# 	imgdata=$( base64 $pngfile )
					# 	echo "<td> <a href=\"$pngfile\"><img src=\"data:image/png;base64,$imgdata\" width=500> </a></td>" >> $html
					# fi
					pngfile=$CurrBiasDir$PREFIX'.interactions_FitHiC_Q'$QVALUE'_Dist_CC.png'
					if [ -f $pngfile ]; then
						imgdata=$( base64 $pngfile )
						echo "<td> <a href=\"$pngfile\"><img src=\"data:image/png;base64,$imgdata\" width=500> </a></td>" >> $html
					fi

					# check the number of interactions vs q-value threshold plot
					# pdffile=$CurrBiasDir'Plots/Interaction_vs_qval.pdf'
					# if [ -f $pdffile ]; then
					# 	pngfile="${pdffile%.*}".png
					# 	if [ ! -f $pngfile ]; then
					# 		convert -verbose -density 500 -resize '800' "${pdffile}" "${pngfile}"
					# 	fi
					# 	imgdata=$( base64 $pngfile )
					# 	echo "<td> <a href=\"$pngfile\"><img src=\"data:image/png;base64,$imgdata\" width=500> </a></td>" >> $html
					# fi
					pngfile=$CurrBiasDir'Plots/Interaction_vs_qval.png'
					if [ -f $pngfile ]; then
						imgdata=$( base64 $pngfile )
						echo "<td> <a href=\"$pngfile\"><img src=\"data:image/png;base64,$imgdata\" width=500> </a></td>" >> $html
					fi

					echo \</tbody\>\</table\> >> $html

				fi 	# end file exist condition
			done 	# end bias specific loop

		fi 	# end condition on current interaction type

	done 	# end p2p loop

done 	# end looptype loop


echo \</body\>\</html\> >> $html


