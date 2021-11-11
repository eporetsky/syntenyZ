#!/usr/bin/python

import sys
import os

project_name = sys.argv[1]
genotypes = " ".join(sys.argv[2:])

start = 'blastn -db "'
quote='" '
query = "-query projects/"+project_name+"/query.fasta "
results = "-out projects/"+project_name+"/blast_results.tsv "
output = '-num_threads 4 -evalue 1e-4 -outfmt "6 staxid qseqid qstart qend sseqid sstart send"'
##print(start+genotypes+quote+query,results,output)
#print(start+genotypes+quote+query,results,output)

os.system(str(start+genotypes+quote+query+results+output))
