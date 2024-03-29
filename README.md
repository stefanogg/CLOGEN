# CLOGEN

Code to run the statistical genomics analyses included in the manuscript: A statistical genomics framework to trace bacterial genomic predictors of clinical outcomes in Staphylococcus aureus bacteraemia (Clinical Outcomes Genomics, CLOGEN)

The script `maigrett.sh` allows the full bacterial GWAS pipeline using Pyseer. It will need as minimum:
- the output of snippy using a single reference for all isolates included
- a tab-separated phenotype file in the format requested by Pyseer
- the reference in Genbak format
- a regions file (genes, operons, pathways) in the format requested by Pyseer.

(The name 'maigrett' was inspired by the famous police detective Jules Maigret created by the Belgian writer Georges Simenon. It also a tribute to the R piping package magrittr. I don't thinkt that stastical genomics fits into the intuitive, poetic thinking of the 'inspecteur Maigret' who was notoriously suspicious of forensic science, but honestly I was just trying to find a way to pay homage to a writer that I like.)

The directory `R` contains the scripts used to fit the random forest models of 30 day mortality from S. aureus bacteraemia

The script `decompose_variance.py` was modified by an original script by John Lees. It provides an estimate of the contribution of lineage effects to the phenotypic variance.
