#!/bin/bash

# Shell script to process vcf files generated by snippy
# Run the script in conda snippy environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <isolates.txt>"
	echo "OPTIONS:"
	echo " -s : path to snippy output directory"
	echo " -c : <coreness.tab>"
	echo " -g : minimum percentage of non-gaps to keep as core genome (default: 0.9)"
	echo " -t : number of threads for bcftools (default 8)"
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

while getopts 's:c:g:t:' option
do
  	case $option in
                s) SNIPPY=$OPTARG ;;
		c) CORENESS=$OPTARG ;;
		g) CORE=$OPTARG ;;
		t) THREADS=$OPTARG ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get positional arguments
ISO=$1

# Starting statements
echo "Running $0"
echo "The snippy output directory is $SNIPPY"
echo "The snippy-coreness file is $CORENESS"
echo "The threshold to keep core genome positions is $CORE non-gaps"

# Get list of strains and generate string of vcf files
ISO=($(cat $ISO))
DIR=(${ISO[@]/#/$SNIPPY/})
VCF=${DIR[@]/%//snps.vcf.gz}

# Check that the string $VCF has 784 words (number of files)
# Count number of snippy vcf files
echo "The number of snps.vcf files that will be merged is"
echo $VCF | wc -w

# Merge vcf files: do not merge multiallelic sites (-m none) and calculate allele counts and allele frequency (note that AC is always double)
bcftools merge -m none -0 -O z --threads $THREADS $VCF | bcftools +fill-tags -O z -- -t AN,AC,AF > merged.vcf.gz
bcftools index merged.vcf.gz

# Check allele count
# bcftools query -f '%POS\t%AN\t%AC\t%AF\n' merged.vcf.gz | head

# Get stats
bcftools stats -s - merged.vcf.gz > merged.vcf.stats.txt

# Remove sites with coreness below threshold
TAG=$(echo $CORE | awk '{print $1*100}')
cp $CORENESS coreness.tab
awk -v core="${CORE}" '$2>=core' coreness.tab > coreness${TAG}.tab
CHROM=$(bcftools query -f '%CHROM\n' merged.vcf.gz | head -n1)
cut -f1 coreness${TAG}.tab | awk -v var=$CHROM -F '\t' 'BEGIN {OFS=FS}{print var,$1-1,$1}' | bedtools merge > core${TAG}.bed
bcftools view -R core${TAG}.bed -O z --threads $THREADS merged.vcf.gz > merged.core${TAG}.vcf.gz
bcftools index merged.core${TAG}.vcf.gz

# Get stats
bcftools stats -s - merged.core${TAG}.vcf.gz > merged.core${TAG}.vcf.stats.txt

# Run homoplasyFinder

# Generate lists of non-synonymous mutations

# Generate lists of non-snynonymous mutations and truncations for the gene-burden test
