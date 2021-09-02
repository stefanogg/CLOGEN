#!/bin/bash

# Shell script to process snippy files in order to generate genotype files for pyseer and other GWAS tools (or ML tools)
# Run the script in conda snippy environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <isolates.txt>"
	echo "OPTIONS:"
	echo " -s : path to snippy outout directory"
	echo " -c : path to coreness.tab"
	echo " -p : path to phylogenetic tree"
	echo " -g : minimum percentage of non-gaps to keep as core genome (default: 0.9)"
	echo "Please run this script in conda snippy environment"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
    	usage
	exit 1
fi

# Parse options
CORE=0.9
THREADS=8

while getopts 's:c:p:g:' option
do
  	case $option in
		s) SNIPPY=$OPTARG ;;
		c) CORENESS=$OPTARG ;;
		p) TREE=$OPTARG ;;
		g) CORE=$OPTARG ;;
		t) THREADS=$OPTARG ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get positional arguments
ISO=$(readlink -e $1)
SNIPPY=$(readlink -e $SNIPPY)
CORENESS=$(readlink -e $CORENESS)
TREE=$(readlink -e $TREE)

# Starting statements
echo "Running $0"
echo "The isolates file is $ISO"
echo "The snippy outout directory is $SNIPPY"
echo "The snippy-coreness file is $CORENESS"
echo "The tree file is $TREE"
echo "The threshold to keep core genome positions is $CORE non-gaps"
echo "$0 will use $THREADS threads"

# Create directory
mkdir all_mutations
cd all_mutations

cp $ISO isolates.txt

# Generate a merged vcf and filter core genome mutations
~/perl5/bin/process-mutations.sh -s $SNIPPY -c $CORENESS -g $CORE -t $THREADS isolates.txt
VCF=$(readlink -e merged.core*vcf.gz)

# Run homoplasyFinder
eval "$(conda shell.bash hook)"
conda activate base
mkdir homoplasy_finder
cd homoplasy_finder
NAME=$(basename $TREE)
CLEAN_TREE=clean.$NAME
gotree prune -i $TREE Reference > $CLEAN_TREE
~/perl5/bin/run-homoplasyfinder.sh -p $CLEAN_TREE $VCF
cd ..

# Generate pyseer genotype files
eval "$(conda shell.bash hook)"
conda activate snippy
~/perl5/bin/process-pyseer-genotype.sh -H homoplasy_finder/homoplasy.bed merged.core90.vcf.gz

# Go back to initial directory
cd ..
