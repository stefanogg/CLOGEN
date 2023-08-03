#!/bin/bash

# Shell pipeline to run pyseer using the kinship matrix generated by gemma
# 3 steps
# 1) Generate kinship matrix using gemma
# 2) Run pyseer mutations
# 3) Run pyseer genes
# Run the script in conda R environment
# The script writes output to the current directory

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <phenotype.tab>"
	echo "OPTIONS:"
	echo " -g : directory with genotype files"
    echo " -k : kinship matrix (will skip the gemma step if provided)"
	echo " -R : path to file with genes in regions format"
	echo " -S : suffix for pyseer output (default: [mutations|genes].<phenotype>.pyseer.tab)"
	echo " -C : continuous phenotype (default: false)"
    echo " -p : other pyseer options"
	echo " -t : number of threads (default: 1)"
	echo "Please run this script in conda R environment unless the gemma step is skipped. The script writes output to the current directory"
}

# Function to print command and then execute
exe() {
        echo "$@"
        eval "$@"
}

run_gemma () {
    mkdir gemma
    cd gemma    
    VCF=$GENO/merged.core90.vcf.gz
    Rscript ~/R_functions/vcf_to_bimbam.R $VCF
    cut -f2 ../phenotype.tab | sed '1d' > pheno.txt
    gemma -g geno.txt -p pheno.txt -gk 1 -o gemma -notsnp
    mv output/ gemma_gk1
    cd ..
    HEADER=$(sed -z 's/\n/,/g' isolates.txt | sed 's/,$//')
    eval "$(conda shell.bash hook)"
    conda activate base
    paste isolates.txt gemma/gemma_gk1/gemma.cXX.txt | csvtk add-header -t -n ,$HEADER > geno.kinship.tab
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

while getopts 'g:k:R:S:Cp:t:' option
do
  	case $option in
		g) GENO=$OPTARG ;;
        k) DIST_FILE=$OPTARG ;;
		R) REGIONS=$OPTARG ;;
		S) SUFFIX=$OPTARG ;;
		C) CONT=true ;;
        p) PYOPT=$OPTARG ;;
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

# Creating a log file
#LOG='run-pyseer-gemma.log'

# Starting statements
echo "Running $0"
echo "The phenotye file is $PHENO"
echo "The directory with genotype file is $GENO"
echo "The regions file is $REGIONS"
echo "The suffix for pyseer output files is $SUFFIX"

if $CONT
then
	echo "The phenotype $SUFFIX is continuous"
	CONT='-c'
else
	echo "The phenotype $SUFFIX is binary"
	CONT=''
fi

if [ ! -z "${PYOPT}" ]
then
        echo "The other pyseer options are: $PYOPT"
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

if [ -z $DIST_FILE ]
then
    echo "No kinship matrix provided. Running gemma ..."
    run_gemma
else
    cp $DIST_FILE geno.kinship.tab
fi

# Run pyseer
echo "" 
eval "$(conda shell.bash hook)"
conda activate pyseer2

# Fix pyseer options
OPT=$CONT

# Add other pyseer options
if [ ! -z $PYOPT ]
then
    OPT="$OPT -p '${PYOPT}'"
fi

# Use run-pyseer-mutations.sh
mkdir mutations
cd mutations
SUFFIX_M=mutations.$SUFFIX
exe "run-pyseer-mutations.sh $OPT -d .. -g $GENO/mutations $SUFFIX_M"
cd ..

# Use run-pyseer-genes.sh
mkdir genes
cd genes
SUFFIX_G=genes.$SUFFIX
exe "run-pyseer-genes.sh $OPT -d .. -g $GENO/genes -r $REGIONS $SUFFIX_G"
cd ..

