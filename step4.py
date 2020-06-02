import gzip

import os, sys

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

fasta_dict = {
'ZxPI566673': "Zx-PI566673-REFERENCE-YAN-1.0.fa.gz",
'B73': "Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz",
'Mo17': "Zm-Mo17-REFERENCE-CAU-1.0.fa.gz",
'W22': "Zm-W22-REFERENCE-NRGENE-2.0.fa.gz",
'P39': "Zm-P39-REFERENCE-NAM-1.0.fasta.gz",
'Oh43': "Zm-Oh43-REFERENCE-NAM-1.0.fasta.gz",
'NC358': "Zm-NC358-REFERENCE-NAM-1.0.fasta.gz",
'NC350': "Zm-NC350-REFERENCE-NAM-1.0.fasta.gz",
'Mo18W': "Zm-Mo18W-REFERENCE-NAM-1.0.fasta.gz",
'M37W': "Zm-M37W-REFERENCE-NAM-1.0.fasta.gz",
'Ky21': "Zm-Ky21-REFERENCE-NAM-1.0.fasta.gz",
'Ki3': "Zm-Ki3-REFERENCE-NAM-1.0.fasta.gz",
'Ki11': "Zm-Ki11-REFERENCE-NAM-1.0.fasta.gz",
'Il14H': "Zm-Il14H-REFERENCE-NAM-1.0.fasta.gz",
'HP301': "Zm-HP301-REFERENCE-NAM-1.0.fasta.gz",
'CML69': "Zm-CML69-REFERENCE-NAM-1.0.fasta.gz",
'CML52': "Zm-CML52-REFERENCE-NAM-1.0.fasta.gz",
'CML333': "Zm-CML333-REFERENCE-NAM-1.0.fasta.gz",
'CML277': "Zm-CML277-REFERENCE-NAM-1.0.fasta.gz",
'CML247': "Zm-CML247-REFERENCE-NAM-1.0.fasta.gz",
'CML228': "Zm-CML228-REFERENCE-NAM-1.0.fasta.gz",
'B97': "Zm-B97-REFERENCE-NAM-1.0.fasta.gz",
'Tzi8': "Zm-Tzi8-REFERENCE-NAM-1.0.fasta.gz",
'Tx303': "Zm-Tx303-REFERENCE-NAM-1.0.fasta.gz",
'Oh7B': "Zm-Oh7B-REFERENCE-NAM-1.0.fasta.gz",
'Ms71': "Zm-Ms71-REFERENCE-NAM-1.0.fasta.gz",
'M162W': "Zm-M162W-REFERENCE-NAM-1.0.fasta.gz",
'CML322': "Zm-CML322-REFERENCE-NAM-1.0.fasta.gz",
'CML103': "Zm-CML103-REFERENCE-NAM-1.0.fasta.gz"
}

with open("coordinate_list.tsv","r") as handle:
    for line in handle:
        # ['B73', '1', '283893731', '284129080']
        line = [s.rstrip() for s in line.split("\t")]  # removes \n from all cells in list
        os.system("python3 fasta_retriever.py "+line[0]+" "+fasta_dict[line[0]]+" "+line[1]+" "+line[2]+" "+line[3])
#os.system("cat fasta_loci/*.fasta > fasta_loci/combined.fasta")
