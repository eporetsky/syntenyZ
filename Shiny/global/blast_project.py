#!/usr/bin/python

import sys
import os
import pandas as pd

project_name = sys.argv[1]
#genotypes = " ".join(sys.argv[2:])

start = 'blastn -db "'
query = " -query projects/"+project_name+"/query.fasta "
#results = "-out projects/"+project_name+"/blast_results.tsv "
formatting = '-num_threads 4 -evalue 1e-4 -outfmt "6 staxid qseqid qstart qend sseqid sstart send"'
output = " >> projects/"+project_name+'/blast_results.tsv'


for db in sys.argv[2:]:
    print(str('blastn -db "BLASTdb/'+ db + '"'+ query + formatting + output))
    os.system('blastn -db "BLASTdb/'+ db + '"'+ query + formatting + output) 

# Returns a dictionary of the taxa ID csv based on the blast DB based on a custom csv file
# For each taxa ID returns a dictionaty key:val like this: '1': 'B73', '2': 'W22' ... 
taxa_dict = {}
title2taxa = open("global/title2taxa.csv","r")
for line in title2taxa:
    line = line.strip().split(",") # strip remove \n from the end
    taxa_dict[int(line[1])] = str(line[0])   # add the lines as a new dicrionary element
# I don't know why it started adding wrong UNICODE encoding to B73 ('\ufeffB73')
taxa_dict[1] = "B73"

df = pd.read_csv("projects/"+project_name+"/blast_results.tsv", sep="\t", header=None)
df[0]= df[0].map(taxa_dict)
df.to_csv("projects/"+project_name+"/blast_results.tsv", sep="\t", header=False, index=False)

#print(str(start+genotypes+quote+query+results+output))
#os.system(str(start+genotypes+quote+query+results+output))
