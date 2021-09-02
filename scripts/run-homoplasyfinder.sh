#!/bin/bash

# Shell script to preprocess data and run homoplasyFinder in Java
# Run the script in conda base environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <vcf-file>"
	echo "OPTIONS:"
	echo " -p: path to phylogenetic tree"
	echo "Please run this script in conda base environment"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
    	usage
	exit 1
fi

# Parse options
while getopts 'p:' option
do
  	case $option in
                p) TREE=$OPTARG ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get positional arguments
VCF=$1

# Starting statements
echo "Running $0"
echo "The vcf file is $VCF"
echo "The tree file is $TREE"

# Check that number of samples in the vcf and the number of tips are the same (to avoid a homoplasyFinder error)
N_TIPS=$(gotree stats -i $TREE | cut -f3 | sed '1d')

# Change environment to use bcftools
eval "$(conda shell.bash hook)"
conda activate snippy
N_SAMPLES=$(bcftools stats $VCF | grep samples | sed -E 's/.*\t([0-9]+)$/\1/')

echo "The number samples in $VCF is $N_SAMPLES"
echo "The number of tips in $TREE is $N_TIPS"

if [ $N_SAMPLES -ne $N_TIPS ]
then
    	echo "The number of samples and the number of tips are not the same: please prune the tree (using gotree prune) or remove samples"
       	# Command to prune the tree
       	#grep -P -o Refere.* <treefile> # find the reference (here with spelling error, hence use grep, but need to print match only because the newick file has one line)
       	#gotree prune -i <treefile> Refere-ce > clean.<treefile>
        exit 1
fi

# Extract positions
echo "Merging multi-allelic positions"
bcftools norm -m +any -O z $VCF > positions.vcf.gz
bcftools stats positions.vcf.gz > positions.vcf.stats.txt
N_POSITIONS=$(grep records: positions.vcf.stats.txt | sed -E 's/.*\t([0-9]+)$/\1/')
echo "Found $N_POSITIONS mutated positions"

# Generate presence-absence file
echo "Generating presence-absence file"
bcftools view -H positions.vcf.gz | awk -F '\t' 'BEGIN {OFS=FS} {print $1,$2-1,$2,$0}' | cut -f4,5 --complement > positions.vcf.bed
bcftools view -h positions.vcf.gz | grep CHROM | cut -f1-9 --complement | sed 's/^/start\tend\t/' > presence_absence.tab
sed 's/\([0-9]*\)\/[0-9]\S*/\1#/g' positions.vcf.bed | sed 's/[2-9][0-9]*#/1/g' | sed 's/#//g' | cut -f1,4-10 --complement | sed '1d' >> presence_absence.tab

# Change conda environment
eval "$(conda shell.bash hook)"
conda activate base
csvtk tab2csv -t presence_absence.tab > presence_absence.csv

# Run homoplasyFinder in Java
echo "Running homoplasyFinder"
java -jar ~/perl5/bin/HomoplasyFinder.jar --tree $TREE --presenceAbsence presence_absence.csv 

# Generate the bed file to filter the vcf of positions
CHROM=$(cut -f1 positions.vcf.bed | sed '1d' | sort -u)
sed '1d' consistencyIndexReport_*.txt | cut -f1,2 | sed "s/^/$CHROM\t/" > homoplasy.bed
