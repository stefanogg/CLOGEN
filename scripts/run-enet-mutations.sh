#!/bin/bash

# Shell script to run elastic net model in pyseer using a mutations genotype 
# Run the script in conda pyseer2 environment

# Functions
# Usage
usage () {
	echo "USAGE"
        echo " $0 [options] -d <dir> -g <geno-dir> <pyseer-suffix>"
        echo "OPTIONS:"
        echo " -d : directory with common files"
        echo " -g : directory with genotype files"
        echo " -G : file with all mutations"
        echo " -c : continous phenotype (default: false)"
        echo " -a : alpha (mixing between ridge-lasso regression, default: 1)"
        echo " -p : other pyseer options"
        echo "Please run this script in a conda pyseer2 environment"
}

# Function to print command and then execute
exe() {
        echo "$@"
        eval "$@"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
    	usage
	exit 1
fi

# Parse options
CONT=false
ALPHA=1

while getopts 'g:G:d:ca:p:' option
do
  	case $option in
		g) GENO_DIR=$OPTARG ;;
        G) CORE_G=$OPTARG ;;
		d) PYSEER_DIR=$OPTARG ;;
		c) CONT=true ;;
        a) ALPHA=$OPTARG ;;
  		p) PYOPT=$OPTARG ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))

# Get positional arguments
SUFFIX=$1

# Get full path to files
PYSEER_DIR=$(readlink -e $PYSEER_DIR)
GENO_DIR=$(readlink -e $GENO_DIR)

# Starting statements
echo "Running $0"
echo "The pyseer directory with common files is $PYSEER_DIR"
echo "The genotype directory is $GENO_DIR"
echo "The output files suffix is $SUFFIX"
echo "The alpha parameter (mixing ridge-lasso regression) is $ALPHA"
if [ ! -z "${PYOPT}" ]
then
       	echo "The other pyseer options are: $PYOPT"
fi

# Get phenotype
mkdir phenotype
cp $PYSEER_DIR/phenotype.tab phenotype/phenotype.tab

# Get genotype
mkdir genotype
cp $GENO_DIR/nosyn*.vcf.gz* genotype
if [ ! -z $CORE_G ]
then
    cp $CORE_G genotype
fi

# Create directories for saved variants and saved models
mkdir save-vars
mkdir save-model

if $CONT
then 
	OPT='--continuous'
else
	OPT=''
fi

# Add other pyseer options
OPT="$OPT $PYOPT"

for i in genotype/*gz;
        do NAME=$(basename $i .vcf.gz);
        for j in maf nomaf;
	do 
		if [ "$j" == "maf" ];
                then
                    	echo "Running pyseer on $i with the default -min-af and max-af settings"
			echo "Running:"
                        exe "pyseer $OPT --wg enet --alpha $ALPHA --save-vars save-vars/nomaf.$NAME --save-model save-model/nomaf.$NAME.$SUFFIX --phenotypes phenotype/phenotype.tab --vcf $i > maf.$NAME.$SUFFIX.pyseer-enet.tab 2> >(tee -a maf.$NAME.$SUFFIX.enet.log >&2)"
                else
                    	echo "Running pyseer on $i without maf filter"
			echo "Running:"
                        exe "pyseer $OPT --min-af 0.001 --max-af 0.999 --wg enet --alpha $ALPHA --save-vars save-vars/nomaf.$NAME --save-model save-model/nomaf.$NAME.$SUFFIX --phenotypes phenotype/phenotype.tab --vcf $i > nomaf.$NAME.$SUFFIX.pyseer-enet.tab 2> >(tee -a nomaf.$NAME.$SUFFIX.enet.log >&2)"
                fi;
        done;
done
