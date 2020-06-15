import collections
import os, sys, pathlib
import pandas as pd
from Bio import SeqIO
from collections import defaultdict


# The output of this cell is a table compatible with my genoPlotR R pipeline
# It tries to parse out all the individual exons and combine them together to genes

inbred_list = os.listdir(sys.argv[1]+"/exons_inbreds")
coordinate_list = open(sys.argv[1]+"/coordinate_list.tsv", "w")

print("Beginning of step 2.")

import pathlib
pathlib.Path(sys.argv[1]+"/custom_loci").mkdir(parents=True, exist_ok=True)
pathlib.Path(sys.argv[1]+"/custom_genes").mkdir(parents=True, exist_ok=True)

predicted_exons_dict = {}
predicted_genes_dict = defaultdict(str)




for inbred in inbred_list:

    if inbred[-5:] != "fasta":
        continue
    inbred = inbred[6:-6]

    records = SeqIO.parse(sys.argv[1]+"/exons_inbreds/exons_"+inbred+".fasta", "fasta")
    egff = [] # edited_gff to be used in the R portion
    for record in records:
        chromosome = record.name.split(":")[0][3:]
        coords = record.name.split(":")[1].split("|")[0].split("..")
        left, right = int(coords[0]), int(coords[1])
        #if True: #record.name.split(":")[0] in ["chr1","Chr1"] and left<left_lim and left>right_lim: # not relevant for Zx8-10, manually pre filtered
        exon = record.name.split("|")[1]
        strand = [right,left,"-"] if left > right else [left,right,"+"]
        # egff columns are: 'name','feature','start','end','strand'
        egff.append([exon,"exon",strand[0], strand[1], strand[2]])
        predicted_exons_dict[(strand[0],strand[1])] = str(record.seq)

		# This is for creating predicted gene fasta files
		#predicted_exons_dict["bla"] = record.seq()

    egff_df = pd.DataFrame(egff, columns = ['name','feature','start','end','strand'])
    min = egff_df[["start","end"]].min().min()
    #print(inbred,"- Chr1:"+str(egff_df[["start","end"]].min().min())+".."+str(egff_df[["start","end"]].max().max()))
    coordinate_list.write(inbred+"\t"+chromosome+"\t"+str(min)+"\t"+str(egff_df[["start","end"]].max().max())+"\n")

    egff_df[["start","end"]] = egff_df[["start","end"]] - min + 1 # +1 for R indexing
    egff_df = egff_df.sort_values("start").reset_index(drop=True)

	# Doing it in two steps because I'm not sure if the blast fasta hits are sorted
	# First step is getting the coordinates for each fasta and putting it into a table
	# Second step is iterating over the table rows and sticking mRNA when needed (also 2 steps)
    ml = [0]  # mRNA list
    count = 1 # count of mRNAs
    min = min - 1

    for i in range(0,egff_df.shape[0]):
        start = egff_df.loc[i,"start"]
        end = egff_df.loc[i,"end"]
        strand = egff_df.loc[i,"strand"]

        if strand == "+":
            predicted_genes_dict[(count,inbred)] = predicted_genes_dict[(count,inbred)] + predicted_exons_dict[(start+min,end+min)]
        else:
            predicted_genes_dict[(count,inbred)] = predicted_exons_dict[(start+min,end+min)] + predicted_genes_dict[(count,inbred)]

        # Stop the loop after writing the last exon
        if i==egff_df.shape[0]-1:
            break
        next_start = egff_df.loc[i+1,"start"]
        next_strand = egff_df.loc[i+1,"strand"]
        if next_start - end > 1000:
            strand = egff_df.loc[i+1,"strand"]
			# if distance between exons>10000 check that they are not consecutive exon names
            # assuming less than 10 exons [-1] takes the 1-9 exons
            if egff_df.loc[i,"name"] in ["left","right"] or egff_df.loc[i+1,"name"] in ["left","right"]:
                ml.append(i+1)
                count += 1
                continue
            left_exon_num = int(egff_df.loc[i,"name"][-1])
            right_exon_num = int(egff_df.loc[i+1,"name"][-1])
            if strand=="+" and next_strand=="+" and left_exon_num+1==right_exon_num:
                continue
            elif strand=="-" and next_strand=="-" and left_exon_num-1==right_exon_num:
                continue
            else:
                ml.append(i+1)
                count += 1



    ml.append(egff_df.shape[0])

    for i in range(count):
		# for each mRNA add an extra row of exon where mRNA should be and overwrite with mRNA row
        egff_df = pd.concat([egff_df.iloc[:ml[i]+i+1],egff_df.iloc[ml[i]+i:]]).reset_index(drop=True)
        egff_df.loc[ml[i]+i]=["tent_gene",
                              "mRNA",
                              egff_df.loc[ml[i]+i+1,"start"],
                              egff_df.loc[ml[i+1]+i,"end"],
                              egff_df.loc[ml[i]+i+1,"strand"]]
    egff_df.index = list(range(1, len(egff_df) + 1))
    egff_df.to_csv(sys.argv[1]+"/custom_loci/custom_loci_"+inbred+".csv")

"""This part writes the fasta file associated with the predicted custom loci coordinates"""
gene_list = [gene[0] for gene in predicted_genes_dict]
gene_dict = {}
for key, value in predicted_genes_dict.items():
    gene_num = key[0]
    gene_inbred = key[1]
    if gene_num not in gene_dict.keys(): # open a fasta file with the inbred name
        gene_dict[gene_num] = open(sys.argv[1]+"/custom_genes/gene_"+str(gene_num)+".fasta", "w")
    gene_dict[gene_num].write('>'+gene_inbred+"\n")
    gene_dict[gene_num].write(value+"\n")
