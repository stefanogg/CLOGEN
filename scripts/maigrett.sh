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
	echo " -p : run IQ-TREE (default: false)"
	echo " -G : run pyseer-genes only (default: false)"
	echo " -C : continuous phenotype (default: false)"
	echo " -t : number of threads for bcftools (default: 8)"
	echo "Please run this script in conda snippy environment"
}

# Function to print command and then execute
exe() {
        echo "$@"
        eval "$@"
}

run_gemma () {
	eval "$(conda shell.bash hook)"
    conda activate R
    mkdir gemma
    cd gemma    
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
CORE=0.9
RUN_TREE=false
GENES_ONLY=false
THREADS=8
CONT=false

while getopts 'r:s:g:R:d:S:c:pGCt:' option
do
  	case $option in
		r) REF=$OPTARG ;;
		s) SNIPPY=$OPTARG ;;
		g) CORE=$OPTARG ;;
		R) REGIONS=$OPTARG ;;
		d) OUTDIR=$OPTARG ;;
		S) SUFFIX=$OPTARG ;;
		c) COREDIR=$OPTARG ;;
		p) RUN_TREE=true ;;
		G) GENES_ONLY=true ;;
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

if $RUN_TREE
then
	echo "$0 will run IQ-tree"
fi

if $GENES_ONLY
then
	echo "$0 will perform gene-burden GWAS only"
fi

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
	if $RUN_TREE
	then
		echo "Running snippy-core and inferring ML tree"
		OPT="-p "
	else
		echo "Running snippy-core"
		OPT=""
	fi
	~/perl5/bin/process-snippy-core.sh $OPT-r $REF -s $SNIPPY isolates.txt -g $CORE
	COREDIR='snippy-core'
fi

# Extract mutations from snippy
echo "" 
if $RUN_TREE
then
	OPT="-p $COREDIR/iqtree/core90.aln.treefile"
else	
	OPT=""
fi
~/perl5/bin/snippy-vcf-to-genotype.sh -s $SNIPPY -c $COREDIR/coreness.tab $OPT-g $CORE isolates.txt 

# Check that snippy-vcf-to-genotype.sh (and process-mutations.sh) worked and exit if not (this is because bcftools merge will raise an error if the number of vcf files exceeds the limit of the system [usually 1,024])
EXIT=$(echo $?)
if [ $EXIT -ne 0 ]
then
    echo "$0 was not able to generate a merged vcf. Exiting"
    exit 1
fi
# Quick and dirty way of extracting records counts from the vcf files
grep -r records all_mutations | grep SN | cut -f1,4 > all_records.count.tab

# Run pyseer
echo "" 
# eval "$(conda shell.bash hook)"
# conda activate pyseer2
mkdir pyseer-gemma
cd pyseer-gemma
cp ../isolates.txt .
cp ../phenotype.tab .
VCF=$(readlink -e ../all_mutations/merged.core90.vcf.gz)
# python ~/perl5/bin/pyseer/similarity-runner.py --vcf $VCF isolates.txt > geno.kinship.tab
run_gemma

eval "$(conda shell.bash hook)"
conda activate pyseer2

if ! $GENES_ONLY
then
	# Use run-pyseer-mutations.sh
	mkdir mutations
	cd mutations
	SUFFIX_M=mutations.$SUFFIX
	~/perl5/bin/run-pyseer-mutations.sh $CONT -d .. -g ../../all_mutations/mutations $SUFFIX_M
	cd ..
fi

# Use run-pyseer-genes.sh
mkdir genes
cd genes
SUFFIX_G=genes.$SUFFIX
~/perl5/bin/run-pyseer-genes.sh $CONT -d .. -g ../../all_mutations/genes -r $REGIONS $SUFFIX_G
cd ..

#  Go back to phenotype dir
cd ..

# Go back to parent directory
cd ..
