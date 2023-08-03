#!/bin/bash

# Shell script to run  snippy-core, snippy-coreness, goalign and iqtree
# Run the script in conda snippy environment

# Functions
usage () {
	echo "USAGE"
	echo " $0 [options] <isolates.txt>"
	echo "OPTIONS:"
	echo " -r : path to reference in gbk format"
	echo " -s : path to snippy outout directory"
	echo " -g : minimum percentage of non-gaps to keep as core genome (default: 0.9)"
	echo " -p : run iqtree (default: false)"
	echo "Please run this script in a conda snippy environment"
}

# Print usage and exit if no arguments
if [ $# -eq 0 ]
then
    	usage
	exit 1
fi

# Parse options
CORE=0.9
TREE=false

while getopts 'r:s:g:p' option
do
  	case $option in
		r) REF=$OPTARG ;;
		s) SNIPPY=$OPTARG ;;
		g) CORE=$OPTARG ;;
        p) TREE=true ;;
		esac
done

# skip to parse command line arguments
shift $((OPTIND-1))


# Get positional arguments
ISO=$(readlink -e $1)
REF=$(readlink -e $REF)
SNIPPY=$(readlink -e $SNIPPY)

# Starting statements
echo "Running $0"
echo "The isolates file is $ISO"
echo "The reference is $REF"
echo "The snippy outout directory is $SNIPPY"
echo "The threshold to keep core genome positions is $CORE non-gaps"
if $TREE
then
	echo "$0 will run IQ-tree"
fi

# Run snippy-core
mkdir snippy-core
cd snippy-core
cp $ISO isolates.txt
mkdir reference
cp $REF reference

# Generate string of snippy directories based on an array
ISO=($(cat isolates.txt))
SNIPPY_DIR=${ISO[@]/#/$SNIPPY/}

# Run snippy-core
snippy-core --ref reference/*.gbk $SNIPPY_DIR 
snippy-coreness core.full.aln > coreness.tab

# Change conda environment
eval "$(conda shell.bash hook)"
conda activate base

# Run goalign
~/perl5/bin/run-goalign.sh -g $CORE core.full.aln
# Get name of final file
TAG=$(echo $CORE | awk '{print $1*100}')
ALN=core${TAG}.aln

# Change conda environment
eval "$(conda shell.bash hook)"
conda activate snippy

if $TREE
then
	# Generate ML tree
	mkdir iqtree
	cd iqtree
	cp ../goalign/$ALN .
	CONST=$(iqtree-calc_const_sites.sh ../core.ref.fa | cut -f2 -d '-' | cut -f2 -d ' ')
	iqtree -fconst $CONST -m GTR+G4 -bb 1000 -alrt 1000 -ntmax 4 -nt AUTO -st DNA -s $ALN
	cd ..
fi

# Go back to initial directory
cd ..
