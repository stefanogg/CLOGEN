# Script to generate a bimbam file from the merged vcf
# Usage: Rscript vcf_to_bimbam.R <vcf>

# Library
library(vcfR)

# Load file
args <- commandArgs(trailingOnly = T)
file <- args[1]
vcf <- read.vcfR(file)

df_vcf <- vcfR2tidy(vcf, info_only = T)
df_mutations <- df_vcf$fix 
df_mutations$mutation_id <- paste(df_mutations$CHROM, df_mutations$POS, df_mutations$REF, df_mutations$ALT, sep = "_")
df_mutations <- df_mutations[,c("mutation_id", "REF", "ALT")]
rm(df_vcf)

# Extract genotype matrix rows = variants, cols = samples
geno <- extract.gt(vcf, return.alleles = F, as.numeric = T) 
geno <- as.data.frame(geno)
geno <- do.call(cbind, list(df_mutations, geno))

write.table(geno, "geno.txt", sep = ",", col.names = F, row.names = F, quote = F)
