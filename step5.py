import os, sys, pathlib
import pathlib


#coord_list = ['B73', '10', '56446489', '56728327']
with open(sys.argv[1]+"/combined_gff3.tsv", "w+") as combined_gff3:
    with open(sys.argv[1]+"/coordinate_list.tsv","r+") as coordinates:
        for coord_list in coordinates:
            coord_list = [s.rstrip() for s in coord_list.split("\t")] # also removes \n
            coord_list[2] = int(coord_list[2])
            coord_list[3] = int(coord_list[3])

            with open("../GFF3/edited_gff3/"+coord_list[0]+".gff3", "r+") as handle:
                for line in handle:
                    line = [s.rstrip() for s in line.split("\t")]
                    try:
                        line[3], line[4] = int(line[3]), int(line[4])
                    except:
                        pass
                    if line[0] in [coord_list[1], "chr"+coord_list[1], "chr"+coord_list[1]]:
                        if line[4]>coord_list[2] and line[3]<coord_list[3]:

                            gff = line[-1].split("=")[1]

                            line[3] = line[3]-coord_list[2]
                            if line[3]<0:
                                line[3] = 1

                            if line[4] > coord_list[3]:
                                line[4] = coord_list[3] - coord_list[2]
                            else:
                                line[4] = line[4]-coord_list[2]

                            if line[6]=="+":
                                line[6] = "1"
                            else:
                                line[6] = "-1"

                            if line[2]=="gene":
                                line[2] = "mRNA" # I forgot that that's what the R script takes for now
                                combined_gff3.write(coord_list[0]+"\t"+gff+"\t"+line[2]+"\t"+str(line[3])+"\t"+str(line[4])+"\t"+line[6]+"\n")
                            if line[2]=="exon":
                                if gff[-5:] == "_T001":
                                    combined_gff3.write(coord_list[0]+"\t"+gff[:-5]+"\t"+line[2]+"\t"+str(line[3])+"\t"+str(line[4])+"\t"+line[6]+"\n")
