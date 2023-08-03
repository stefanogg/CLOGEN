#!/bin/bash

# Shell script to run a fixed effects model in pyseer using a gene-burden test
# Run the script in conda pyseer2 environment

# Functions
# Usage
usage () {
	echo "USAGE"
        echo " $0 [options] -d <dir> -g <geno-dir> <pyseer-suffix>"
        echo "OPTIONS:"
        echo " -d : directory with common files"
        echo " -k : distance file"
        echo " -D : number of MDS dimensions to retain (default: 4)"
        echo " -g : directory with genotype files"
        echo " -r : regions file"
        echo " -c : continous phenotype (default: false)"
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
SIMPLE=false
DIM=4

while getopts 'g:d:k:D:r:cp:S' option
do
  	case $option in
		g) GENO_DIR=$OPTARG ;;
		d) PYSEER_DIR=$OPTARG ;;
 		k) DIST_FILE=$OPTARG ;;
        D) DIM=$OPTARG ;;
        r) REGIONS=$OPTARG ;;
		c) CONT=true ;;
  		p) PYOPT=$OPTARG ;;
        S) SIMPLE=true;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))

# Get positional arguments
SUFFIX=$1

# Get full path to files
PYSEER_DIR=$(readlink -e $PYSEER_DIR)
if [ -z $DIST_FILE ]
then
    	DIST_FILE=$(readlink -e $PYSEER_DIR/mash.dist.tab)
else
        DIST_FILE=$(readlink -e $DIST_FILE)
fi
GENO_DIR=$(readlink -e $GENO_DIR)

# Starting statements
echo "Running $0"
echo "The pyseer directory with common files is $PYSEER_DIR"
echo "The distance matrix is $DIST_FILE"
echo "The number of MDS dimensions to retain is : $DIM"
echo "The genotype directory is $GENO_DIR"
echo "The files with reference regions is $REGIONS"
echo "The output files suffix is $SUFFIX"
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

# Get distance matrix
mkdir distance_matrix
cp $DIST_FILE distance_matrix

# Get regions
mkdir regions
cp $REGIONS regions/ref.regions.txt

# Run pyseer
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
    exe "pyseer $OPT --phenotypes phenotype/phenotype.tab --vcf $i --distances distance_matrix/mash.dist.tab --max-dimensions $DIM --burden regions/ref.*regions.txt > $NAME.$SUFFIX.pyseer-fixed.tab"
done
