import os, sys, pathlib
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# extract a fasta file for each inbred from the raw blast results
# save these as a fasta file with the name of the inbred
# also keep a list of inbreds to iterate over in the next steps

print(sys.argv[1])

print("Beginning blastPlotR Step 1")

inbred_dict = {}
inbred_list = []

pathlib.Path(sys.argv[1]+"/exons_inbreds").mkdir(parents=True, exist_ok=True)
file_list = os.listdir(sys.argv[1]+"/exons_raw")

for file in file_list: # name of blast files
    if file.split(".")[1] != "fasta":
        continue
    file = file.split(".")[0]
    records = SeqIO.parse(sys.argv[1]+'/exons_raw/'+str(file)+'.fasta', "fasta")     # fasta files
    for record in records:
        inbred = record.name.split(".")[0] # inbred name, written before "."
        coords = record.name.split(".",1)[1].split("|")[0] # coordinates of blast result

        if inbred not in inbred_dict.keys(): # open a fasta file with the inbred name
            inbred_dict[inbred] = open(sys.argv[1]+"/exons_inbreds/exons_"+inbred+".fasta", "w")
            inbred_list.append(inbred)

        # if int(coords.split(".")[-1]) < 700000000: # not relevant for Zx8-10, manually pre filtered
        seq = SeqRecord(record.seq,
                    id=coords+"|"+file,
                    description="")
        SeqIO.write(seq, inbred_dict[inbred], 'fasta')
print("Finished step 1. List of inbred lines processed: " + ", ".join(inbred_list) + ".")
