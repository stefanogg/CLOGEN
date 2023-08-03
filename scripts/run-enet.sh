#!/bin/bash

# Shell pipeline to run the elastic net model in pyseer
# 3 steps
# 1) Prepare file structure
# 2) Run pyseer mutations
# 3) Run pyseer genes
# Run the script in conda pyseer2 environment
# The script writes output to the current directory

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <phenotype.tab>"
	echo "OPTIONS:"
	echo " -g : directory with genotype files"
	echo " -R : path to file with genes in regions format"
	echo " -S : suffix for pyseer output (default: [mutations|genes].<phenotype>.pyseer-enet.tab)"
	echo " -C : continuous phenotype (default: false)"
    echo " -a : alpha (mixing between ridge-lasso regression, default: 1)"
	echo " -t : number of threads (default: 1)"
	echo "Please run this script in conda pyseer2 environment. The script writes output to the current directory"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
    	usage
	exit 1
fi

# Parse options
THREADS=1
CONT=false

while getopts 'g:R:S:a:Ct:' option
do
  	case $option in
		g) GENO=$OPTARG ;;
		R) REGIONS=$OPTARG ;;
		S) SUFFIX=$OPTARG ;;
        a) ALPHA=$OPTARG ;;
		C) CONT=true ;;
		t) THREADS=$OPTARG ;;
        esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get positional arguments
PHENO=$(readlink -e $1)

# Modify optional arguments
GENO=$(readlink -e $GENO)
REGIONS=$(readlink -e $REGIONS)

if [ -z $SUFFIX ]
then
	SUFFIX=$(basename $PHENO .tab | sed 's/_phenotype//')
fi

if [ -z $ALPHA ]
then
    ALPHA=1
fi

# Creating a log file
#LOG='run-pyseer-gemma.log'

# Starting statements
echo "Running $0"
echo "The phenotye file is $PHENO"
echo "The directory with genotype file is $GENO"
echo "The regions file is $REGIONS"
echo "The suffix for pyseer output files is $SUFFIX"
echo "The alpha parameter (mixing ridge-lasso regression) is $ALPHA"

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

# Creating file structure
if [ -f phenotype.tab ]
then
	echo "phenotype.tab already present"
else
	cp $PHENO phenotype.tab
	sed '1d' phenotype.tab | cut -f1 > isolates.txt
fi

# Launch run-enet-mutations.sh
mkdir mutations
cd mutations
SUFFIX_M=mutations.$SUFFIX
run-enet-mutations.sh $CONT -d .. -g $GENO/mutations -G $GENO/merged.core90.vcf.gz -a $ALPHA $SUFFIX_M
cd ..

# Use run-enet-genes.sh
mkdir genes
cd genes
SUFFIX_G=genes.$SUFFIX
run-enet-genes.sh $CONT -d .. -g $GENO/genes -r $REGIONS -a $ALPHA $SUFFIX_G
cd ..

