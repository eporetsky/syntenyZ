#!/usr/bin/python

#https://gist.github.com/theJollySin/6eeda4a44db830a35365503178f88788

import gzip

import os, sys

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

import gzip
import io


def get_sequence_from_fasta(inbred_name, fasta_name, chromosome, left, right, padding=0):
    seq = ""
    left=int(left)
    right=int(right)
    with gzip.open("../"+fasta_name, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in ["chr"+chromosome, "Chr"+chromosome, chromosome]:
                print(chromosome, "found in:", fasta_name)
                seq = record.seq[left-1-padding:right+padding]
                with open("fasta_loci/"+inbred_name+".fasta","a") as output:
                        output.write('>'+inbred_name+'\n')
                        output.write(str(seq)+'\n')

                break
    #handle.close()


get_sequence_from_fasta(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
