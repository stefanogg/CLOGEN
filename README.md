# CLOGEN

Code to run the statistical genomics analyses included in the manuscript: A statistical genomics framework to trace bacterial genomic predictors of clinical outcomes in Staphylococcus aureus bacteraemia (Clinical Outcomes Genomics, CLOGEN)

The script `maigrett.sh` allows the full bacterial GWAS pipeline using Pyseer. It will need as minimum:
- the output of snippy using a single reference for all isolates included
- a tab-separated phenotype file in the format requested by Pyseer
- the reference in Genbak format
- a regions file (genes, operons, pathways) in the format requested by Pyseer.

The name 'maigrett' was inspired by the famous police detective Jules Maigret created by the Belgian writer Georges Simenon. It also a tribute to the R piping package magrittr. 

The directory `R` contains the scripts used to fit the random forest models of 30 day mortality from S. aureus bacteraemia
