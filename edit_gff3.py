import gzip

gff3_dict = {
'ZxPI566673': "Zx-PI566673-REFERENCE-YAN-1.0_Zx00001a.1.gff3.gz",
'B73': "Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz",
'Mo17': "Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3.gz",
'W22': "Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3.gz",
'P39': "Zm-P39-REFERENCE-NAM-1.0_Zm00040a.1.gff3.gz",
'Oh43': "Zm-Oh43-REFERENCE-NAM-1.0_Zm00039a.1.gff3.gz",
'NC358': "Zm-NC358-REFERENCE-NAM-1.0_Zm00037a.1.gff3.gz",
'NC350': "Zm-NC350-REFERENCE-NAM-1.0_Zm00036a.1.gff3.gz",
'Mo18W': "Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034a.1.gff3.gz",
'M37W': "Zm-M37W-REFERENCE-NAM-1.0_Zm00032a.1.gff3.gz",
'Ky21': "Zm-Ky21-REFERENCE-NAM-1.0_Zm00031a.1.gff3.gz",
'Ki3': "Zm-Ki3-REFERENCE-NAM-1.0_Zm00029a.1.gff3.gz",
'Ki11': "Zm-Ki11-REFERENCE-NAM-1.0_Zm00030a.1.gff3.gz",
'Il14H': "Zm-Il14H-REFERENCE-NAM-1.0_Zm00028a.1.gff3.gz",
'HP301': "Zm-HP301-REFERENCE-NAM-1.0_Zm00027a.1.gff3.gz",
'CML69': "Zm-CML69-REFERENCE-NAM-1.0_Zm00020a.1.gff3.gz",
'CML52': "Zm-CML52-REFERENCE-NAM-1.0_Zm00019a.1.gff3.gz",
'CML333': "Zm-CML333-REFERENCE-NAM-1.0_Zm00026a.1.gff3.gz",
'CML277': "Zm-CML277-REFERENCE-NAM-1.0_Zm00024a.1.gff3.gz",
'CML247': "Zm-CML247-REFERENCE-NAM-1.0_Zm00023a.1.gff3.gz",
'CML228': "Zm-CML228-REFERENCE-NAM-1.0_Zm00022a.1.gff3.gz",
'B97': "Zm-B97-REFERENCE-NAM-1.0_Zm00018a.1.gff3.gz",
'Tzi8': "Zm-Tzi8-REFERENCE-NAM-1.0_Zm00042a.1.gff3.gz",
'Tx303': "Zm-Tx303-REFERENCE-NAM-1.0_Zm00041a.1.gff3.gz",
'Oh7B': "Zm-Oh7B-REFERENCE-NAM-1.0_Zm00038a.1.gff3.gz",
'Ms71': "Zm-Ms71-REFERENCE-NAM-1.0_Zm00035a.1.gff3.gz",
'M162W': "Zm-M162W-REFERENCE-NAM-1.0_Zm00033a.1.gff3.gz",
'CML322': "Zm-CML322-REFERENCE-NAM-1.0_Zm00025a.1.gff3.gz",
'CML103': "Zm-CML103-REFERENCE-NAM-1.0_Zm00021a.1.gff3.gz"
}

write_mRNA = False
mRNA = "001"
for key, value in gff3_dict.items():
    with gzip.open("../GFF3/"+value, "rt") as handle:
        #with gzip.open("gff3/"+key+".gff3.gz","wt",compresslevel=9) as fixed:
        with open("gff3/"+key+".gff3","w+") as fixed:
            for row in handle:
                if "\t" not in row or row[0] == "#":
                    continue

                row = row.split("\t")
                desc = row[-1].split(";")

                if row[2] == "gene":
                    gene = desc[0].split("=")[1]
                    if gene.startswith("gene:"):
                        gene = gene[5:]
                    row[-1] = "=".join(["ID",gene])
                    fixed.write("\t".join(row)+"\n")
                    continue

                if row[2] == "mRNA":
                    mRNA = desc[0].split("=")[1]
                    if mRNA.startswith("transcript:"):
                        mRNA = mRNA[11:]
                    row[-1] = "=".join(["ID",mRNA])
                    if write_mRNA:
                        fixed.write("\t".join(row)+"\n")
                    continue

                if mRNA[-3:] != "001":
                    continue

                if row[2] == "exon":
                    row[-1] = "=".join(["Parent",mRNA])
                    fixed.write("\t".join(row)+"\n")
                    continue
