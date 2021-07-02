#!/bin/sh


##Workflow

#1. Run Blast
makeblastdb -in GCA_905220365.1_ilVanCard2.1_genomic_chroms.fna -dbtype nucl -parse_seqids -out GCA_905220365.1_ilVanCard2.1_genomic_chroms.fna -title "GCA_905220365.1_ilVanCard2.1_genomic_chroms.fna"  
blastn -db GCA_905220365.1_ilVanCard2.1_genomic_chroms.fna -query pacbio.prim.assembly.fasta -outfmt 6 -out Chrom_scaffolds'.tab'  -evalue 0.010 -word_size 1000 -num_threads 6

#2. sort the output
python blast_order.py > VC_CHROM_1.csv
