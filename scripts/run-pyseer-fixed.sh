#!/bin/bash

# Shell pipeline to run the fixed effects model in pyseer
# 3 steps
# 1) Generate distance matrix
# 2) Run pyseer mutations
# 3) Run pyseer genes
# Run the script in conda pyseer2 environment
# The script writes output to the current directory

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] -g <dir> -R <regions.txt> <phenotype.tab>"
	echo "OPTIONS:"
	echo " -g : directory with genotype files"
    echo " -m : path to mash sketches (default: ~/GWAS_clinical_outcomes/assemblies/mash/)"
    echo " -k : distance matrix (will skip the mash triangle step if provided)"
    echo " -D : number of MDS dimensions to retain (default: 4)"
	echo " -R : path to file with genes in regions format"
	echo " -S : suffix for pyseer output (default: [mutations|genes].<phenotype>.pyseer.tab)"
	echo " -C : continuous phenotype (default: false)"
    echo " -p : other pyseer options"
	echo " -t : number of threads (default: 1)"
	echo "Please run this script in conda pyseer2. The script writes output to the current directory"
}

# Function to print command and then execute
exe() {
        echo "$@"
        eval "$@"
}

distance () {
    mkdir distance
    cd distance
    for ISO in $(cat ../isolates.txt); do find $MSH -name $ISO.msh ; done > sketches.txt
    mash paste -l sketch sketches.txt
    mash dist sketch.msh sketch.msh | square_mash > mash.dist.tab
    cd ..
    cp distance/mash.dist.tab .
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
DIM=4
MSH=$(readlink -e ~/GWAS_clinical_outcomes/assemblies/mash/)

while getopts 'g:k:m:D:R:S:Cp:t:' option
do
  	case $option in
		g) GENO=$OPTARG ;;
        k) DIST_FILE=$OPTARG ;;
        m) MSH=$OPTARG ;;
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
#MSH=$(readlink -e ${MSH})

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
    echo "No distance matrix provided. Running mash dist ..."
    distance
else
    cp $DIST_FILE mash.dist.tab
fi

# Run pyseer

# Fix pyseer options
OPT=$CONT

# Add other pyseer options
if [ ! -z "${PYOPT}" ]
then
    OPT="$OPT -p '${PYOPT}'"
fi

# Use run-pyseer-mutations.sh
mkdir mutations
cd mutations
SUFFIX_M=mutations.$SUFFIX
exe "run-fixed-mutations.sh $OPT -d .. -g $GENO/mutations -D $DIM $SUFFIX_M"
cd ..

# Use run-pyseer-genes.sh
mkdir genes
cd genes
SUFFIX_G=genes.$SUFFIX
exe "run-fixed-genes.sh $OPT -d .. -g $GENO/genes -r $REGIONS -D $DIM $SUFFIX_G"
cd ..

