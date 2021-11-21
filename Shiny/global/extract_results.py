import os, sys
from collections import OrderedDict
from collections import defaultdict
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import os
project_folder = "projects/"+sys.argv[1]+"/"

print(sys.argv)

#######################################################################################################################
######################## Part 0 - Generate a list of coordinates using relative start site ############################
#######################################################################################################################
count = 0
identity_dict = {}
chr2identity = open("global/chr2identity.tsv","r")
for line in chr2identity:
    if count==0:
        identity_dict[("title","chr",)] = "identity" # add values to tuple keys
    else:
        line = line.strip().split("\t")            # strip remove \n from the end
        identity_dict[tuple(line[0:2])] = line[2] # add values to tuple keys
    count += 1
identity_dict

#######################################################################################################################
######################## Part 1 - Generate a list of coordinates using relative start site ############################
#######################################################################################################################

blast_df = pd.read_table(project_folder+"blast_results_filtered.tsv", header=None, names=["taxa","ref", "qs","qe","chr","ss","se"])
blast_df['strand'] = blast_df[['ss','se']].apply(lambda x: "+" if x['ss']<x['se'] else "-", axis=1)
blast_tolist = blast_df.sort_values(['taxa', 'chr', 'ss',  ], ascending=[True, True, True]).values.tolist()

overlap_list = [list(blast_tolist[0])]
gene_copies = {overlap_list[-1][1]:1}
list_len = len(blast_tolist)
for row in range(1,list_len):
    last = overlap_list[-1]
    nrow = list(blast_tolist[row]) # the new row list
    ovl, ovr = (last[5],last[6]) if last[-1]=="+" else (last[6],last[5]) #overlap_left/right
    rol, ror = (nrow[5],nrow[6]) if nrow[-1]=="+" else (nrow[6],nrow[5]) #row_left/right

    if nrow[0] != last[0]:
        gene_copies = {nrow[1]:1} # first count of gene_exon copy
    
    if (nrow[0] == last[0]) & (nrow[4] == last[4]) & (nrow[-1] == last[-1]): # if taxa, chrom, strand equal
        if (rol <= ovl) & (ror >= ovr): # if the new blast row spans over then replace with new row completely
            overlap_list[-1] = nrow # if the new row spans previous, replace
        elif (rol < ovl) & (ror >= ovl) & (ror <= ovr):        
            overlap_list[-1][1] = nrow[1]
            if nrow[-1] == "+":
                overlap_list[-1][5] = rol
            else:
                overlap_list[-1][6] = rol
        
        elif (ror > ovr) & (rol >= ovl) & (rol <= ovr):
            overlap_list[-1][1] = nrow[1]
            if nrow[-1] == "+":
                overlap_list[-1][6] = ror
            else:
                overlap_list[-1][5] = ror
        elif (rol >= ovl) & (ror<=ovr): # if new row is within previous row then skip it
            continue
        else:
            overlap_list.append(nrow)
    else:
        overlap_list.append(nrow)

with open(project_folder+"raw_overlap_list.tsv", "w") as overlap_list_out: 
    for row in overlap_list:
        overlap_list_out.write("\t".join(str(e) for e in row)+"\n")
        
#######################################################################################################################
######################## Part 2 - Generate a list of coordinates using relative start site ############################
#######################################################################################################################
left_padding = 0
right_padding = 0

loci_coord = {} 
xlim_coord = {} 
gene_coord = {} 

overlap_list_copy = list(overlap_list.copy())
last = overlap_list_copy[0]

longest_exon_dict = {}
for line in overlap_list_copy:
    if line[1] not in longest_exon_dict.keys():
        longest_exon_dict[line[1]] = line[3]
    else:
        longest_exon_dict[line[1]] = max(longest_exon_dict[line[1]], line[3])

copy_dict = defaultdict(lambda: 1)
for line in overlap_list_copy:
    strand = line[-1]
    left = min(int(line[5]), int(line[6]))
    right = max(int(line[5]), int(line[6]))
    genotype = line[0]
    chromosome = line[4]
    
    do_all=True
    if line[1].split("_")[0] == "GeneID" or do_all==True:
        # Sometimes the splice-sites don't match but it's easy to see it in alignment vs. manually finding the mismatch
        if line[2] < 6:
            if strand == "+":
                left = left - line[2] + 1
            else:
                right = right + line[2] - 1
    # Also check splice-site against exon length, if it is 1 or 2 bp shorter than longest one
    
        if  (longest_exon_dict[line[1]] - line[3]) < 20 and (longest_exon_dict[line[1]] - line[3]) > 0:
            # In this specific instance the stop site in some lines was not within the reference gene range
            padding = 0
            if line[1] == 'TPS31_exon7': 
                padding = 6
            if strand == "+":
                right = right + longest_exon_dict[line[1]] - line[3] + padding
            else:
                left = left - (longest_exon_dict[line[1]] - line[3]) - padding
        
    if len(line[1].split("_")) == 2:
        # Get the gene ID and exon number based on query
        gene = line[1].split("_")[0].split(".")[0]
        exon = line[1].split("_")[1]
    else:
        # Need to check if overlapped blast-hits are renamed to this format
        gene = line[1].split("_")[0] + "_" + line[1].split("_")[1]
        exon = line[1].split("_")[2]
    
    # Because we assume that the overlap list is sorted, first genotype line is upstream
    if genotype not in loci_coord.keys():
        # left_loci is the reference upstream genomic coordinate initated once for each genotype
        left_loci = left - left_padding
        # Start inserting the genomic and relative coordinates to each dictionary
        loci_coord[genotype] = [chromosome, left_loci, right + right_padding]
        xlim_coord[genotype] = [[gene, 1+left_padding, right-left_loci+left_padding]]
        gene_coord[genotype] = [[gene, exon,1+left_padding, right-left_loci+left_padding, strand]]
    else:
        # Every new loci_coord line extend the right loci, because lines are ordered
        loci_coord[genotype][2] = right + right_padding
        # Below assumes the duplicate genes are in order
        if xlim_coord[genotype][-1][0].split(".")[0] == gene: 
            previous_strand = gene_coord[genotype][-1][-1] # I don't think this does anything
            if len(line[1].split("_")) == 2:
                # Need to standarize it, if it has exon written in the exon number or not
                previous_exon = int(gene_coord[genotype][-1][-4].split("exon")[1])
                current_exon = int(exon.split("exon")[1])
            else:
                previous_exon = int(gene_coord[genotype][-1][-4])
                current_exon = int(exon)
            # If exons are consecutive, assume that they are part of the same gene and adjust values
            # if  ((previous_exon+1==current_exon) and strand=="+") or ((previous_exon-1==current_exon) and strand=="-"):
            if  ((previous_exon<current_exon) and strand=="+") or ((previous_exon>current_exon) and strand=="-"):
                xlim_coord[genotype][-1][2] = right-left_loci+left_padding
                gene_coord[genotype][-1].extend([exon, left-left_loci+left_padding, right-left_loci+left_padding, strand])
            # Sometimes minor insertions break exons, so if exon names are the same and arbitrary distance is short keep going
            elif ((left-left_loci-gene_coord[genotype][-1][-2]) < 5000) & (exon==gene_coord[genotype][-1][-4].split("_")[-1]): 
                xlim_coord[genotype][-1][2] = right-left_loci+left_padding
                gene_coord[genotype][-1].extend([exon, left-left_loci+left_padding, right-left_loci+left_padding, strand])
            # Else, assume this is a new duplicate gene
            else: 
                copy_dict[(genotype,gene,)] += 1
                gene = gene+"."+str(copy_dict[(genotype,gene,)])
                xlim_coord[genotype].append([gene, left-left_loci+left_padding, right-left_loci+left_padding])
                gene_coord[genotype].append([gene, exon, left-left_loci+left_padding, right-left_loci+left_padding, strand])            
        else:
            xlim_coord[genotype].append([gene, left-left_loci+left_padding, right-left_loci+left_padding])
            gene_coord[genotype].append([gene, exon, left-left_loci+left_padding, right-left_loci+left_padding, strand])

#######################################################################################################################
############################### Part 3 - Extract the genomic sequence for each line ###################################
#######################################################################################################################

coordinate_list = open(project_folder+"coordinate_list.tsv",'w')
for genotype, coords in loci_coord.items():
    chromosome, left, right = coords[0], coords[1], coords[2]
    coordinate_list.write("\t".join([genotype,chromosome,str(left),str(right)])+"\n")
    os.system('echo ">' + genotype + '" >> '+project_folder+'fasta_loci/combined.fasta')
    os.system('blastdbcmd -db BLASTdb/'+genotype+'.db -entry "'+identity_dict[(genotype,chromosome)]+'" -range "'+str(left)+'-'+str(right)+'" -strand plus -line_length 1000000000| sed 1d  >> ' +project_folder+'fasta_loci/combined.fasta')

#######################################################################################################################
########################### Part 4 - Mask the genomic sequences, keeping only genes   #################################
#######################################################################################################################

# Import the full combined sequences as adictionary
combined_dict = SeqIO.to_dict(SeqIO.parse(open(project_folder+"fasta_loci/combined.fasta"),'fasta'))
for genotype, genes in xlim_coord.items():
    mask_list = []
    for gene in genes:
        mask_list += gene[1:3]
    mask_list = mask_list[1:-1]    
    # make a nested paired nested list of lists
    zipped_coordinates = zip(mask_list[0::2], mask_list[1::2])
    #print(list(zipped_coordinates))
    # Make a masked combined.fasta
    # https://stackoverflow.com/questions/63169902/replace-a-list-of-characters-with-indices-in-a-string-in-python
    for L,R in zipped_coordinates:
        combined_dict[genotype] = combined_dict[genotype][:L] + "N"*(R-L) + combined_dict[genotype][R:]    
        
with open(project_folder+"fasta_loci/combined_masked.fasta", "w") as out_file:
    for name, seq in combined_dict.items():
        SeqIO.write(seq, out_file, "fasta")


#######################################################################################################################
############################### Part 5 - Generate the custom loci and xlim files   ####################################
#######################################################################################################################

# need to make a custom_loci directory

for genotype, genes in gene_coord.items():
    if os.path.exists(project_folder+"custom_loci/custom_loci_"+genotype+".csv"):
        os.remove(project_folder+"custom_loci/custom_loci_"+genotype+".csv")
    count = 1
    with open(project_folder+"custom_loci/custom_loci_"+genotype+".csv", 'a') as custom_loci:
        custom_loci.write(",name,feature,start,end,strand\n")
        for gene in genes:
            name,exon,left,right,strand = gene[0],gene[1],gene[2],gene[-2],gene[4]
            if name == "PROPEP4":
                name = "PROPEP4.1"
            custom_loci.write(",".join([str(count),name,"mRNA",str(left),str(right),strand])+"\n")
            count+=1
            for exon_count in range(int((len(gene)-1)/4)):
                exon,left,right,strand = gene[1+4*exon_count:5+4*exon_count]
                custom_loci.write(",".join([str(count),exon,"exon",str(left),str(right),strand])+"\n")
                count+=1

with open(project_folder+'assigned_xlims.csv', 'w') as delete: 
    pass
delete.close()

with open(project_folder+'assigned_xlims.csv', 'a') as xlims_out:
    xlims_out.write("genotype,left,right\n")
    for genotype, xlims in xlim_coord.items():
        for xlim in xlims:
            xlims_out.write(genotype+","+str(xlim[1])+","+str(xlim[2])+"\n")
            


#######################################################################################################################
################ Part 5 - A fasta parser to save all gene CDS sequences based on coordinates   ########################
#######################################################################################################################

CDS_genes = OrderedDict()
CDS_inbreds = {}
add_id_to_inbred_name = True

#Rename variable. The combined.fasta contains all the full DNA seqs for each inbred
# TODO: Implement parsing using biopython SeqIO instead of python text IO parser
# Assumed fasta sequence for each value is on a single line
coordinate_list = open(project_folder+"fasta_loci/combined.fasta",'r')
for row in coordinate_list:
    # Iterate over fasta file and save inbred name and sequence in temps
    if row[0] == ">":
        temp_inbred = row[1:-1] # skips the linebreak
        continue
    else:
        temp_seq = row[:-1] # skips the linebreak
    
    # Open the speicific inbred file
    with open(project_folder+"custom_loci/custom_loci_"+temp_inbred+".csv", 'r') as loci:
        # Prase the custom loci file for each gene
        skip = True
        temp_CDS = '' 
        for loc in loci:
            if skip == True: skip = False ; continue # skip the first line
            loc = loc.strip().split(",")
            # if the line is mRNA then
            if loc[2]=="mRNA":
                # if the temp_CDS is not empty and mRNA means we are done with parsing
                # the previous gene and that we can add the combined CDS to the dictionary
                if temp_CDS != '':
                    # Save the sequence to the  CDS_genes dict file to write later
                    CDS_genes.update({temp_gene:temp_CDS})
                    # If the gene ID is not in dict, start a new key with avaialble
                    if temp_gene not in CDS_inbreds.keys():
                        CDS_inbreds[temp_gene] = OrderedDict({temp_inbred_id:temp_CDS})
                    # Not sure why, or if, this else statement is needed, could add a print command
                    # Maybe it is for reverse-complement sequences?
                    else:
                        CDS_inbreds[temp_gene].update({temp_inbred_id:temp_CDS})
                temp_gene = loc[1]
                # A boolean of wether to add gene ID to inbred ID. One of the outputs is a fasta
                # file of all CDS sequences for each inbred, and easier to handle if it has built-in ID
                # of both gene and inbred. Alternatively, could be save as just inbred name
                if add_id_to_inbred_name:
                    temp_inbred_id = temp_inbred+"_"+temp_gene
                else:
                    temp_inbred_id = temp_inbred
                    
                temp_gene = loc[1]
                temp_CDS = ''
            
            # When iterating over exons (aka not mRNA) write the DNA sequence to CDS seq
            elif loc[-1]=="+":
                if loc[3]=="1": loc[3]=0 # xlims start with 1 instead of 0, but other coordinates correct
                temp_CDS += temp_seq[int(loc[3]):int(loc[4])+1]
            elif loc[-1]=="-":
                if loc[3]=="1": loc[3]=0 # xlims start with 1 instead of 0, but other coordinates correct
                temp_CDS = str(Seq(temp_seq[int(loc[3]):int(loc[4])+1]).reverse_complement()) + temp_CDS
        
        # Since the last line is not mRNA need to manually add the last gene
        # Copied that from above to add the last gene
        CDS_genes.update({temp_gene:temp_CDS})
        if temp_gene not in CDS_inbreds.keys():
            CDS_inbreds[temp_gene] = OrderedDict({temp_inbred_id:temp_CDS})
        else:
            CDS_inbreds[temp_gene].update({temp_inbred_id:temp_CDS})
    
    with open(project_folder+"custom_CDS/inbreds/"+temp_inbred+".fasta", "w") as inbred_CDS:
        for key, val in CDS_genes.items():
            inbred_CDS.write(">"+key+"\n")
            inbred_CDS.write(val+"\n")

# This has to be written after all inbred sequences were collected            
for gene, gene_dict in CDS_inbreds.items():
    with open(project_folder+"custom_CDS/genes/"+gene+".fasta", "w") as genes_CDS:
        for key, val in gene_dict.items():
            genes_CDS.write(">"+key+"\n")
            genes_CDS.write(val+"\n")

