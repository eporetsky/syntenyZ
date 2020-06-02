import collections
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# extract a fasta file for each inbred from the raw blast results
# save these as a fasta file with the name of the inbred
# also keep a list of inbreds to iterate over in the next steps
inbred_dict = {}
inbred_list = []
for file in ["exon1",'exon2','exon3','exon4','exon5','exon6','exon7']: # name of blast files
    records = SeqIO.parse('exons_raw/'+str(file)+'.fasta', "fasta")     # fasta files
    for record in records:
        inbred = record.name.split(".")[0] # inbred name, written before "."
        coords = record.name.split(".",1)[1].split("|")[0] # coordinates of blast result
        
        if inbred not in inbred_dict.keys(): # open a fasta file with the inbred name
            inbred_dict[inbred] = open("exons_inbreds/exons_"+inbred+".fasta", "w")
            inbred_list.append(inbred)
        print(file,inbred,coords)
        if int(coords.split(".")[-1]) < 700000000: # not relevant for Zx8-10, manually pre filtered
        	seq = SeqRecord(record.seq, 
        	               id=coords+"|"+file,
            	           description="")
        	SeqIO.write(seq, inbred_dict[inbred], 'fasta')
print(inbred_list)