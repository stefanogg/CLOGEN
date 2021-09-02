#!/bin/bash

# Shell script to run pyseer using a gene-burden genotype
# Run the script in conda pyseer environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <pyseer-suffix>"
	echo "OPTIONS:"
	echo " -d : directory with common files"
	echo " -g : directory with genotype files"
	echo " -r : regions file"
   	echo " -c : continous phenotype (default: false)"
	echo "Please run this script in a conda pyseer environment"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
    	usage
	exit 1
fi

# Parse options
CONT=false

while getopts 'g:d:r:c' option
do
  	case $option in
		g) GENO_DIR=$OPTARG ;;
		d) PYSEER_DIR=$OPTARG ;;
		r) REGIONS=$OPTARG ;;
		c) CONT=true ;;
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

# Get phenotype
mkdir phenotype
cp $PYSEER_DIR/phenotype.tab phenotype/phenotype.tab

# Get genotype
mkdir genotype
cp $GENO_DIR/*.vcf.gz* genotype

# Get kniship matrix
mkdir kinship_matrix
cp $PYSEER_DIR/geno.kinship.tab kinship_matrix

# Get regions
mkdir regions
cp $REGIONS regions/ref.cds.regions.txt

# Run pyseer
if $CONT
then
     	OPT='--continuous'
else
    	OPT=''
fi

for i in genotype/*gz;
	do NAME=$(basename $i .vcf.gz);
	pyseer $OPT --lmm --phenotypes phenotype/phenotype.tab --vcf $i --burden regions/ref.cds.regions.txt --similarity kinship_matrix/geno.kinship.tab > $NAME.$SUFFIX.pyseer.tab;
done
