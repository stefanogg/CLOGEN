#!/bin/bash

# Shell script to run pyseer using a mutations genotype
# Run the script in conda pyseer environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <pyseer-suffix>"
	echo "OPTIONS:"
	echo " -d : directory with common files"
	echo " -g : directory with genotype files"
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

while getopts 'g:d:c' option
do
  	case $option in
		g) GENO_DIR=$OPTARG ;;
		d) PYSEER_DIR=$OPTARG ;;
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

# Starting statements
echo "Running $0"
echo "The pyseer directory with common files is $PYSEER_DIR"
echo "The genotype directory is $GENO_DIR"
echo "The output files suffix is $SUFFIX"

# Get phenotype
mkdir phenotype
cp $PYSEER_DIR/phenotype.tab phenotype/phenotype.tab

# Get genotype
mkdir genotype
cp $GENO_DIR/nosyn*.vcf.gz* genotype

# Get kniship matrix
mkdir kinship_matrix
cp ../geno.kinship.tab kinship_matrix

# Run pyseer
if $CONT
then 
	OPT='--continuous'
else
	OPT=''
fi


for i in genotype/*gz;
        do NAME=$(basename $i .vcf.gz);
        for j in maf nomaf;
                do echo $j;
                if [ "$j" == "maf" ];
                then
                    	echo "Running pyseer on $i with the default -min-af and max-af settings"
                        pyseer $OPT --lmm --phenotypes phenotype/phenotype.tab --vcf $i --similarity kinship_matrix/geno.kinship.tab > maf.$NAME.$SUFFIX.pyseer.tab;
                else
                    	echo "Running pyseer on $i without maf filter"
                        pyseer $OPT --min-af 0.001 --max-af 0.999 --lmm --phenotypes phenotype/phenotype.tab --vcf $i --similarity kinship_matrix/geno.kinship.tab > nomaf.$NAME.$SUFFIX.pyseer.tab;
                fi;
        done;
done
