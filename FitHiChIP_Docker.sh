#!/bin/bash
usage(){
cat << EOF
Options:
   	-C  ConfigFile		Name of the configuration file storing the parameters of FitHiChIP.
EOF
}

while getopts "C:" opt;
do
	case "$opt" in
		C) ConfigFile=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done



awk "/^[^#].+=.+/" $ConfigFile | sed "s/=/: /g" > ./nf.yaml

nextflow run main.nf -params-file ./nf.yaml -profile docker


echo "Command completed. It is normal to see some warnings above. You can see if FitHiChIP completed successfully by looking above the warnings."

