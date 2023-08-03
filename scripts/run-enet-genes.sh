#!/bin/bash

# Shell script to run elastic net model in pyseer using a genes genotype 
# Run the script in conda pyseer2 environment

# Functions
# Usage
usage () {
	echo "USAGE"
        echo " $0 [options] -d <dir> -g <geno-dir> -r <pyseer-suffix>"
        echo "OPTIONS:"
        echo " -d : directory with common files"
        echo " -g : directory with genotype files"
        echo " -r : regions file"
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

while getopts 'g:d:r:ca:p:' option
do
  	case $option in
		g) GENO_DIR=$OPTARG ;;
		d) PYSEER_DIR=$OPTARG ;;
        r) REGIONS=$OPTARG ;;
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
REGIONS=$(readlink -e $REGIONS)

# Starting statements
echo "Running $0"
echo "The pyseer directory with common files is $PYSEER_DIR"
echo "The genotype directory is $GENO_DIR"
echo "The files with reference regions is $REGIONS"
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
cp $GENO_DIR/*.vcf.gz* genotype

# Get regions
mkdir regions
cp $REGIONS regions/ref.regions.txt

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
    echo "Running:"
    exe "pyseer $OPT --wg enet --alpha $ALPHA --save-vars save-vars/$NAME --save-model save-model/$NAME.$SUFFIX --phenotypes phenotype/phenotype.tab --vcf $i --burden regions/ref.*regions.txt > $NAME.$SUFFIX.pyseer-enet.tab 2> >(tee -a $NAME.$SUFFIX.enet.log >&2)"  
done
