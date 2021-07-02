import sys
import string
import numpy as np
import pandas
import math
import re
import random 

#commands that need to be run before you get here
#makeblastdb -in GCA_905220365.1_ilVanCard2.1_genomic_chroms.fna -dbtype nucl -parse_seqids -out GCA_905220365.1_ilVanCard2.1_genomic_chroms.fna -title "GCA_905220365.1_ilVanCard2.1_genomic_chroms.fna"  
#blastn -db GCA_905220365.1_ilVanCard2.1_genomic_chroms.fna -query pacbio.prim.assembly.fasta -outfmt 6 -out Chrom_scaffolds'.tab'  -evalue 0.010 -word_size 1000 -num_threads 6

new_version_chromosomes_4=["LR999924.1","LR999925.1","LR999926.1","LR999927.1","LR999928.1","LR999929.1","LR999930.1","LR999931.1","LR999932.1","LR999933.1","LR999934.1","LR999935.1","LR999936.1","LR999937.1","LR999938.1","LR999939.1","LR999940.1","LR999941.1","LR999942.1","LR999943.1","LR999944.1","LR999945.1","LR999946.1","LR999947.1","LR999948.1","LR999949.1","LR999950.1","LR999951.1","LR999952.1","LR999953.1","LR999954.1","LR999955.1"] # list of scaffolds in the new genome

table=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/Venkat/vanessa_cardui_project/genome_Scans/Chrom_scaffolds.tab",names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"] , sep='\t', ) # Blast Out put 6

table=table.dropna()

fasta_DPlex_new_updates={} #Tracking directory for scaffolds



order=["0","0","Chr_Z","Chr_1","Chr_2","Chr_3","Chr_4","Chr_5","Chr_6","Chr_7","Chr_8","Chr_9","Chr_10","Chr_11","Chr_12","Chr_13","Chr_14","Chr_15","Chr_W","Chr_16","Chr_17","Chr_18","Chr_19","Chr_20","Chr_21","Chr_22","Chr_23","Chr_24","Chr_25","Chr_26","Chr_27","Chr_28","Chr_29","Chr_30"]
#print 'scaffold'+ ',' +'Best_possible_Chrom'+","+",".join([str(i) for i in order]) # output
print 'scaffold'+ ',' +'Best_possible_Chrom'+","+"Proportion_scaffold_of_mapping" # output
for i in set(table['qseqid']): 
	if len(i) > 1:
		table_N=table[table['qseqid']==i]
		Chr_Z=table_N[table_N['sseqid']=="LR999924.1"]['length'].sum()
		Chr_1=table_N[table_N['sseqid']=="LR999925.1"]['length'].sum()
		Chr_2=table_N[table_N['sseqid']=="LR999926.1"]['length'].sum()
		Chr_3=table_N[table_N['sseqid']=="LR999927.1"]['length'].sum()
		Chr_4=table_N[table_N['sseqid']=="LR999928.1"]['length'].sum()
		Chr_5=table_N[table_N['sseqid']=="LR999929.1"]['length'].sum()
		Chr_6=table_N[table_N['sseqid']=="LR999930.1"]['length'].sum()
		Chr_7=table_N[table_N['sseqid']=="LR999931.1"]['length'].sum()
		Chr_8=table_N[table_N['sseqid']=="LR999932.1"]['length'].sum()
		Chr_9=table_N[table_N['sseqid']=="LR999933.1"]['length'].sum()
		Chr_10=table_N[table_N['sseqid']=="LR999934.1"]['length'].sum()
		Chr_11=table_N[table_N['sseqid']=="LR999935.1"]['length'].sum()
		Chr_12=table_N[table_N['sseqid']=="LR999936.1"]['length'].sum()
		Chr_13=table_N[table_N['sseqid']=="LR999937.1"]['length'].sum()
		Chr_14=table_N[table_N['sseqid']=="LR999938.1"]['length'].sum()
		Chr_15=table_N[table_N['sseqid']=="LR999939.1"]['length'].sum()
		Chr_W=table_N[table_N['sseqid']=="LR999940.1"]['length'].sum()
		Chr_16=table_N[table_N['sseqid']=="LR999941.1"]['length'].sum()
		Chr_17=table_N[table_N['sseqid']=="LR999942.1"]['length'].sum()
		Chr_18=table_N[table_N['sseqid']=="LR999943.1"]['length'].sum()
		Chr_19=table_N[table_N['sseqid']=="LR999944.1"]['length'].sum()
		Chr_20=table_N[table_N['sseqid']=="LR999945.1"]['length'].sum()
		Chr_21=table_N[table_N['sseqid']=="LR999946.1"]['length'].sum()
		Chr_22=table_N[table_N['sseqid']=="LR999947.1"]['length'].sum()
		Chr_23=table_N[table_N['sseqid']=="LR999948.1"]['length'].sum()
		Chr_24=table_N[table_N['sseqid']=="LR999949.1"]['length'].sum()
		Chr_25=table_N[table_N['sseqid']=="LR999950.1"]['length'].sum()
		Chr_26=table_N[table_N['sseqid']=="LR999951.1"]['length'].sum()
		Chr_27=table_N[table_N['sseqid']=="LR999952.1"]['length'].sum()
		Chr_28=table_N[table_N['sseqid']=="LR999953.1"]['length'].sum()
		Chr_29=table_N[table_N['sseqid']=="LR999954.1"]['length'].sum()
		Chr_30=table_N[table_N['sseqid']=="LR999955.1"]['length'].sum()
		#print i+ " "+ str(Chr_2)+ " "+str(Chr_3)+ " "+str(Chr_4)+ " "+str(Chr_5)+ " "+str(Chr_6)+ " "+str(Chr_7)+ " "+str(Chr_8)+ " "+str(Chr_9)+ " "+str(Chr_10)+ " "+str(Chr_11)+ " "+str(Chr_12)+ " "+str(Chr_13)+ " "+str(Chr_14)+ " "+str(Chr_15)+ " "+str(Chr_16)+ " "+str(Chr_17)+ " "+str(Chr_18)+ " "+str(Chr_19)+ " "+str(Chr_20)+ " "+str(Chr_21)+ " "+str(Chr_22)+ " "+str(Chr_23)+ " "+str(Chr_24)+ " "+str(Chr_25)+ " "+str(Chr_26)+ " "+str(Chr_27)+ " "+str(Chr_28)+ " "+str(Chr_29)+ " "+str(Chr_30)+ " "+str(Chr_Z)
		#print i+ " "+ str(CDS_Chr_2)+ " "+str(CDS_Chr_3)+ " "+str(CDS_Chr_4)+ " "+str(CDS_Chr_5)+ " "+str(CDS_Chr_6)+ " "+str(CDS_Chr_7)+ " "+str(CDS_Chr_8)+ " "+str(CDS_Chr_9)+ " "+str(CDS_Chr_10)+ " "+str(CDS_Chr_11)+ " "+str(CDS_Chr_12)+ " "+str(CDS_Chr_13)+ " "+str(CDS_Chr_14)+ " "+str(CDS_Chr_15)+ " "+str(CDS_Chr_16)+ " "+str(CDS_Chr_17)+ " "+str(CDS_Chr_18)+ " "+str(CDS_Chr_19)+ " "+str(CDS_Chr_20)+ " "+str(CDS_Chr_21)+ " "+str(CDS_Chr_22)+ " "+str(CDS_Chr_23)+ " "+str(CDS_Chr_24)+ " "+str(CDS_Chr_25)+ " "+str(CDS_Chr_26)+ " "+str(CDS_Chr_27)+ " "+str(CDS_Chr_28)+ " "+str(CDS_Chr_29)+ " "+str(CDS_Chr_30)+ " "+str(CDS_Chr_Z)
		#print " "
		scaffolds=[0,0,Chr_Z,Chr_1,Chr_2,Chr_3,Chr_4,Chr_5,Chr_6,Chr_7,Chr_8,Chr_9,Chr_10,Chr_11,Chr_12,Chr_13,Chr_14,Chr_15,Chr_W,Chr_16,Chr_17,Chr_18,Chr_19,Chr_20,Chr_21,Chr_22,Chr_23,Chr_24,Chr_25,Chr_26,Chr_27,Chr_28,Chr_29,Chr_30]
		scaffolds=[float(n)/sum(scaffolds) for n in scaffolds] #fraction of mapping
		Scaf_max = order[np.argmax(scaffolds)]
		#print i+ ',' +str(Scaf_max)+","+",".join([str(i) for i in scaffolds]) # output
		print i+ ',' +str(Scaf_max)+","+str(scaffolds[np.argmax(scaffolds)]) # output
		if Scaf_max not in fasta_DPlex_new_updates.keys():
			fasta_DPlex_new_updates[Scaf_max]=[i]
		else:
			fasta_DPlex_new_updates[Scaf_max].append(i)


