#!/bin/bash

# Shell script to process vcf files to generate pyseer genotype files
# Run the script in conda snippy environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <merged.vcf.gz>"
	echo "OPTIONS:"
	echo " -H : bed file with homoplastic sites (if empty skip homoplasy filter)"
	echo " -P : prefix for output files (default core90)"
	echo " -Q : maximum allele frequency threshold for rare mutations (default 0.01)" 
	echo " -m : don't generate gene-burden files (default: false)"
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
MUTATIONS_ONLY=false
PREFIX=core90
THREADS=8
MAX_FRAC=0.01

while getopts 't:mH:P:Q:' option
do
  	case $option in
		t) THREADS=$OPTARG ;;
        m) MUTATIONS_ONLY=true ;;
		H) HOMOPLASY=$OPTARG ;;
		P) PREFIX=$OPTARG ;;
		Q) MAX_FRAC=$OPTARG ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get positional arguments
VCF=$1

# Get full path to files
VCF=$(readlink -e $VCF)


# Starting statements
echo "Running $0"
echo "The merged vcf file is $VCF"
if [ ! -z $HOMOPLASY ]
then
	HOMOPLASY=$(readlink -e $HOMOPLASY)
	echo "The file with homoplastic sites is $HOMOPLASY"
fi
echo "The prefix for output files is $PREFIX"
echo "The maxmimum allele fraction for rare mutations is $MAX_FRAC"

# Generate mutation lists
echo "Generating mutations lists for pyseer"
mkdir mutations
cd mutations

# Mask mutations at non-homoplastic sites. bedtools much faster than bcftools here. Output is the same

#bcftools view -R $HOMOPLASY -O z $VCF > homoplasy.$PREFIX.vcf.gz
if [ ! -z $HOMOPLASY ]
then
	bedtools intersect -u -header -a $VCF -b $HOMOPLASY | bgzip > homoplasy.$PREFIX.vcf.gz
	bcftools view -e 'ANN[0] ~ "[^&]synonymous"' -O z --threads $THREADS homoplasy.$PREFIX.vcf.gz > nosyn.homoplasy.$PREFIX.vcf.gz

	# Get stats
	bcftools stats homoplasy.$PREFIX.vcf.gz > homoplasy.$PREFIX.vcf.stats.txt
	bcftools stats nosyn.homoplasy.$PREFIX.vcf.gz > nosyn.homoplasy.$PREFIX.vcf.stats.gz
fi

# Mask synonymous mutations but NOT homoplastic sites
bcftools view -e 'ANN[0] ~ "[^&]synonymous"' -O z --threads $THREADS $VCF > nosyn.$PREFIX.vcf.gz

# Get stats
bcftools stats nosyn.$PREFIX.vcf.gz > nosyn.$PREFIX.vcf.stats.txt

# Index vcf files
for i in *vcf.gz; do bcftools index $i; done

# Go back to parent directory
cd ..

# Check -m options and exit if true
if $MUTATIONS_ONLY
then
	exit 1
fi

# Generate gene-burden list
echo "Generating gene-burden lists for pyseer"
mkdir genes
cd genes

# Get non synonymous mutations: $PREFIX, no homoplasy filter, no maf filter
cp ../mutations/nosyn.$PREFIX.vcf.gz* .
bcftools stats nosyn.$PREFIX.vcf.gz > nosyn.$PREFIX.vcf.stats.txt

# Get non synonymous mutations: $PREFIX, no homoplasy filter, maf filter
bcftools view -Q $MAX_FRAC -O z --threads $THREADS nosyn.$PREFIX.vcf.gz > rare.nosyn.$PREFIX.vcf.gz
bcftools stats rare.nosyn.$PREFIX.vcf.gz > rare.nosyn.$PREFIX.vcf.stats.txt

# Get truncations: $PREFIX, no homoplasy filter, no maf filter
bcftools view -i 'ANN[0] ~ "frameshift_variant" | ANN[0] ~ "stop_gained" | ANN[0] ~ "stop_lost" | ANN[0] ~ "start_lost"' -O z --threads $THREADS nosyn.$PREFIX.vcf.gz > trunc.$PREFIX.vcf.gz
bcftools stats trunc.$PREFIX.vcf.gz > trunc.$PREFIX.vcf.stats.txt

# Get truncations: $PREFIX, no homoplasy filter, maf filter
bcftools view -Q $MAX_FRAC -O z --threads $THREADS trunc.$PREFIX.vcf.gz > rare.trunc.$PREFIX.vcf.gz
bcftools stats rare.trunc.$PREFIX.vcf.gz > rare.trunc.$PREFIX.vcf.stats.txt

# Check number of mutations
# grep record *stats*

# Index vcf files
for i in *vcf.gz; do bcftools index $i; done

# Go back to parent directory
cd ..
