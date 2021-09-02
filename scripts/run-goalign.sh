#!/bin/bash

# Shell script to run go-align and snippy-coreness
# Run the script in conda base environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <core.full.aln>"
	echo "OPTIONS:"
	echo " -g : minimum percentage of non-gaps to keep as core genome (default: 0.9)"
	echo "Please run this script in conda base environment"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
    	usage
	exit 1
fi

# Parse options
CORE=0.9

while getopts 'g:' option
do
  	case $option in
		g) CORE=$OPTARG ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get positional arguments
ALN=$(readlink -e $1)

# Starting statements
echo "Running $0"
echo "The alignment file is $ALN"
echo "The threshold to keep core genome positions is $CORE non-gaps"

# Create new directory
mkdir goalign
cd goalign

# Clean the alignment (replace all non ACTG with -)
echo "Cleaning the alignment"
cp $ALN core.full.aln
seqkit replace -s -p [NnX] -r - core.full.aln > core.clean.full.aln

# Remove all sites with > 10% gaps
# calculate c option
echo "Removing sites with less than ${CORE} non-gap sites"
OPT=$(echo $CORE | awk '{print 1-$1}')
NAME_TAG=$(echo $CORE | awk '{print $1*100}')
goalign clean sites -c $OPT -i core.clean.full.aln > core${NAME_TAG}.full.aln 2>goalign.c${OPT}.log

# Get variant sites (including sites with < 10%)
# Change conda environment
eval "$(conda shell.bash hook)"
conda activate snippy

snp-sites -o core${NAME_TAG}.aln core${NAME_TAG}.full.aln

# OPTIONAL Get vcf file of variant sites (for pyseer distance matrix)
# snp-sites -c -o core90.vcf core90.full.aln

# Final stats
# Change conda environment
eval "$(conda shell.bash hook)"
conda activate base
seqkit stats *aln > aln.stats.tab
