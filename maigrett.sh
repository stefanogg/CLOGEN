#!/bin/bash

# Shell pipeline to run a GWAS analysis using the snippy output
# maigret performs 3 major steps
# 1) Run snippy-core, extract core genome positions and generate phylogenetic tree (this step is very slow, it can be skipped if needed)
# 2) Extract list of mutations in vcf format from the snippy output
# 3) Run pyseer (mutations and gene burden test)
# Run the script in conda snippy environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <phenotype.tab>"
	echo "OPTIONS:"
	echo " -r : path to reference in gbk format"
	echo " -s : path to snippy outout directory"
	echo " -g : minimum percentage of non-gaps to keep as core genome (default: 0.9)"
	echo " -R : path to file with genes in regions format"
	echo " -d : name of the output directory (default: basename of phenotype.tab files)"
	echo " -S : suffix for pyseer output (default: [mutations|genes].<phenotype>.pyseer.tab)"
	echo " -c : path to snippy-core directory (if empty run process-snippy-core: slow!)"
	echo " -C : continuous phenotype (default: false)"
	echo " -t : number of threads for bcftools (default: 8)"
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
CONT=false

while getopts 'r:s:g:R:d:S:c:Ct:' option
do
  	case $option in
		r) REF=$OPTARG ;;
		s) SNIPPY=$OPTARG ;;
		g) CORE=$OPTARG ;;
		R) REGIONS=$OPTARG ;;
		d) OUTDIR=$OPTARG ;;
		S) SUFFIX=$OPTARG ;;
		c) COREDIR=$OPTARG ;;
		C) CONT=true ;;
		t) THREADS=$OPTARG ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get positional arguments
PHENO=$(readlink -e $1)

# Modify optional arguments
REF=$(readlink -e $REF)
SNIPPY=$(readlink -e $SNIPPY)
REGIONS=$(readlink -e $REGIONS)

if [ ! -z $COREDIR ]
then
	COREDIR=$(readlink -e $COREDIR)
fi

if [ -z $OUTDIR ]
then
	OUTDIR=$(basename $PHENO .tab)
fi

if [ -z $SUFFIX ]
then
	SUFFIX=$(echo $OUTDIR | sed 's/_phenotype//')
fi

# Creating a log file
#LOG='maigret.log'

# Starting statements
echo "Running $0"
echo "The phenotye file is $PHENO"
echo "The reference is $REF"
echo "The snippy outout directory is $SNIPPY"
echo "The threshold to keep core genome positions is $CORE non-gaps"
echo "The regions file is $REGIONS"
echo "The output directory is $OUTDIR"
echo "The suffix for pyseer output files is $SUFFIX"

if $CONT
then
	echo "The phenotype $SUFFIX is continuous"
	CONT='-c'
else
	echo "The phenotype $SUFFIX is binary"
	CONT=''
fi

echo "$0 will use $THREADS threads"
echo ""

if [ ! -z $COREDIR ]
then
	echo "The snippy-core directory is $COREDIR"
	echo "$0 will skip the snippy-core and tree step"
fi

# Create overarching file structrure
if [ -d $OUTDIR ]
then
	echo "$OUTDIR exists. Did you run snippy-core previously?"
else
	mkdir $OUTDIR
fi

cd $OUTDIR
if [ -f phenotype.tab ]
then
	echo "phenotype.tab already present"
else
	cp $PHENO phenotype.tab
	sed '1d' phenotype.tab | cut -f1 > isolates.txt
fi

# Run snippy-core
echo "" 
if [ -z $COREDIR ]
then
	echo "Running snippy-core and inferring ML tree"
	~/perl5/bin/process-snippy-core.sh -r $REF -s $SNIPPY isolates.txt -g $CORE
	COREDIR='snippy-core'
fi

# Extract mutations from snippy
echo "" 
~/perl5/bin/snippy-vcf-to-genotype.sh -s $SNIPPY -c $COREDIR/coreness.tab -p $COREDIR/iqtree/core90.aln.treefile -g $CORE isolates.txt 
# Quick and dirty way of extracting records counts from the vcf files
grep -r records all_mutations | grep SN | cut -f1,4 > all_records.count.tab

# Run pyseer
echo "" 
eval "$(conda shell.bash hook)"
conda activate pyseer
mkdir pyseer
cd pyseer
cp ../isolates.txt .
cp ../phenotype.tab .
VCF=../all_mutations/merged.core90.vcf.gz
python ~/perl5/bin/pyseer/similarity-runner.py --vcf $VCF isolates.txt > geno.kinship.tab

# Use run-pyseer-mutations.sh
mkdir mutations
cd mutations
SUFFIX_M=mutations.$SUFFIX.pyseer.tab
~/perl5/bin/run-pyseer-mutations.sh $CONT -d .. -g ../../all_mutations/mutations $SUFFIX_M
cd ..

# Use run-pyseer-genes.sh
mkdir genes
cd genes
SUFFIX_G=genes.$SUFFIX.pyseer.tab
~/perl5/bin/run-pyseer-genes.sh $CONT -d .. -g ../../all_mutations/genes -r $REGIONS $SUFFIX_G
cd ..

#  Go back to phenotype dir
cd ..

# Go back to parent directory
cd ..
