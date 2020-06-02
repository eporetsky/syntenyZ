import collections
import pandas as pd
from Bio import SeqIO

# The output of this cell is a table compatible with my genoPlotR R pipeline
# It tries to parse out all the individual exons and combine them together to genes

species_list = ["ZxPI566673", 'B73', 'Mo17', 'W22', 'P39', 'Oh43', 'NC358', 'NC350', 'Mo18W', 'M37W', 
				'Ky21', 'Ki3', 'Ki11', 'Il14H', 'HP301', 'CML69', 'CML52', 'CML333', 'CML277', 
				'CML247', 'CML228', 'B97', 'Tzi8', 'Tx303', 'Oh7B', 'Ms71', 'M162W', 'CML322', 'CML103']


#range_df[range_df['species'].isin(species_list)]

right_lim = 270000000
left_lim = 300000000

# For Teosinte
#right_lim = 119000000
#left_lim = 119900000

print("copy the bellow coordinates to extrat the DNA sequences from fasta files")

for species in species_list:
    records = SeqIO.parse('exons_inbreds/exons_'+species+'.fasta', "fasta")
    egff = [] # edited_gff to be used in the R portion
    for record in records:
        coords = record.name.split(":")[1].split("|")[0].split("..")
        left, right = int(coords[0]), int(coords[1])
        if True:#record.name.split(":")[0] in ["chr1","Chr1"] and left<left_lim and left>right_lim: # not relevant for Zx8-10, manually pre filtered
            exon = record.name.split("|")[1]
            strand = [right,left,"-"] if left > right else [left,right,"+"]
            # egff columns are: 'name','feature','start','end','strand'
            egff.append([exon,"exon",strand[0], strand[1],strand[2]])
    
    egff_df = pd.DataFrame(egff, columns = ['name','feature','start','end','strand'])  
    #print(species,"- Chr1:"+str(egff_df[["start","end"]].min().min())+".."+str(egff_df[["start","end"]].max().max()))
    print(species+"\tChr1\t"+str(egff_df[["start","end"]].min().min())+"\t"+str(egff_df[["start","end"]].max().max()))
    egff_df[["start","end"]] = egff_df[["start","end"]] - egff_df[["start","end"]].min().min() +1
    egff_df = egff_df.sort_values("start").reset_index(drop=True)
    
    
    ml = [0]  # mRNA list
    count = 1 # count of mRNAs
    #print(egff_df)
    for i in range(0,egff_df.shape[0]-1):
        if egff_df.loc[i+1,"start"] - egff_df.loc[i,"end"] > 10000:
            # if distance between exons>1000 check that they are not consecutive exons
            strand = egff_df.loc[i+1,"strand"]
            left_exon_num = int(egff_df.loc[i,"name"][-1])
            right_exon_num = int(egff_df.loc[i+1,"name"][-1])
            if (strand=="+" and left_exon_num+1==right_exon_num) or (strand=="-" and left_exon_num-1==right_exon_num):
                continue
            else:
                ml.append(i+1)
                count += 1
    ml.append(egff_df.shape[0])
    for i in range(count):
        egff_df = pd.concat([egff_df.iloc[:ml[i]+i+1],egff_df.iloc[ml[i]+i:]]).reset_index(drop=True)
        egff_df.loc[ml[i]+i]=["tent_gene",
                              "mRNA",
                              egff_df.loc[ml[i]+i+1,"start"],
                              egff_df.loc[ml[i+1]+i,"end"],
                              egff_df.loc[ml[i]+i+1,"strand"]]
    egff_df.index = list(range(1, len(egff_df) + 1))
    egff_df.to_csv('custom_loci/custom_loci'+species+'.csv')