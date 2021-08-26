library(genoPlotR)
library(RColorBrewer)
`%notin%` <- Negate(`%in%`)
gene_colors <- read.table("assigned_features.csv",sep=",",header=T)

#### Script for using alternating lines ####
all_names <- seg_name
target = "B73"
seg_name <- c()
for (name in all_names){
  if(name!=target){
    seg_name <- c(seg_name,name)
    seg_name <- c(seg_name,target)
  }
}
####   End   ####


#### Script generating the comparative all vs all blast results ####
system("cat fasta_loci/*.fasta > fasta_loci/combined.fasta")
system("makeblastdb -dbtype nucl -in fasta_loci/combined.fasta -out fasta_loci/combined.db")
system("blastn -db fasta_loci/combined.db -query fasta_loci/combined.fasta -out blastall.tsv -num_threads 4 -word_size 100 -evalue 10 -outfmt 6 ")
blastall <- read_comparison_from_blast("blastall.tsv")
####   End   ####

# Order of lines in the final figure
seg_name <- c(
  "PI566673",
  "PH207",
  "EP1",
  "B73",
  "W22",
  "Il14H",
  "NC350",
  "CML247",
  "CML277",
  "Tzi8",
  "Mo17",
  "Ky21",
  "CML69",
  "Mo18W",
  "CML228",
  "CML322",
  "NC358",
  "F7",
  "Oh7B",
  "CML333",
  "B97",
  "Oh43",
  "M37W",
  "Tx303",
  "Ms71",
  "Ki11",
  "Ki3",
  "CML52",
  "HP301",
  "A188",
  "M162W",
  "P39",
  "CML103")


####  Script generating the the genoPlotR figure  ####
#### Load all the libraries below to make it work ####
MygenoPlotR(file_name="test",
            blastall=blastall,
            do_gene_xlim=F,
            plot_height=15, # Don't think it is working
            plot_width=4,
            exon_mode="custom", # gff3 or exon
            gene_annot = "NA",
            annotation_cex = 1,
            annotation_height = 1,
            seg_names=seg_name,
            gene_colors=gene_colors,
            add_intron=F,
            #xlims=xlims_list, #xlims_list. NULL is set for between first and last gff3 element
            annot=T,
            #gene_type = "arrows", 
            arrow_len = 10000,
            text_color = "black")

# Main genoPlotR syntenyZ function that generates the output PDF file
MygenoPlotR <- function(
  file_name,         # choose the file-name to save the PDF figure
  blastall,
  do_gene_xlim=F,
  plot_height=40,
  plot_width=25,
  exon_mode="custom",       # use gff3 or custom coordinate table
  gene_annot="gff3", # support NA, gff3 or custom, blank otherwise
  annotation_cex = 1,
  annotation_height = 1,
  rot = 0,
  seg_names,         # ordered names of the genomes to parse their files
  gene_colors,       # a table of elements and their specified colors
  add_intron=T,
  annot=F,
  #gene_type = "arrows",
  arrow_len=10000,   # change the size of the gene arrows
  text_color="black" # set NULL to remove labels
){
  
  # list of comparisons based on the blast results assuming you ran all-vs-all blast command
  # the coordinate system is not the same as the genomic coordinates. Adjust accordingly
  
  xlims <- xlim_func(seg_names,do_gene_xlim=do_gene_xlim)
  comp_list <- list()
  for (f in 1:(length(seg_names)-1)){
    comp_list[f] <- list(blastall[blastall["name1"]==seg_names[f] & blastall["name2"]==seg_names[f+1],])
    if(do_gene_xlim){comp_list[[f]] <- trim_comp(comp_list[[f]], xlims[[f]], xlims[[f+1]])}
  }
  # parse the genes in the gff or custom coordinate file to generate the dna-seg files
  # that are drown on the plot. There are couple different methods available
  annot_list <- list()
  dna_segs <- list()
  for (f in 1:(length(seg_names))){
    if(exon_mode=="gff3"){
      # takes gff3 and returns 5 column coordinate tabe
      gff3_filename <- ifelse(do_gene_xlim,"combined_xgff3.tsv","combined_gff3.tsv")
      edited_gff <- read.table(gff3_filename, sep="\t",header=F)
      edited_gff <- edited_gff[edited_gff$V1==seg_names[f], c(2:6)]
      colnames(edited_gff) <- c("name","feature","start","end","strand")
      edited_gff <- merge(edited_gff, gene_colors[,c("name","fill")], by = "name", all.x = TRUE)
    }
    if(gene_annot=="gff3"){annot_list[f] <- list(try(as.annotation(gff2annot(edited_gff,color=text_color,rot=0))))}
    custom_gff <- read.table(paste("custom_loci/custom_loci_",seg_name[f],".csv",sep=""), sep=",",header=T,row.names=NULL)
    mutation_gff <- custom_gff[custom_gff$feature  %notin% c("mRNA","exon"),]
    custom_gff <- custom_gff[custom_gff$feature  %in% c("mRNA","exon"),]
    custom_gff <- merge(custom_gff, gene_colors[,c("name","fill")], by = "name", all.x = TRUE)
    custom_gff <- custom_gff[with(custom_gff, order(custom_gff[,"X"])),]
    custom_gff <- subset(custom_gff, select = -X )
    count = 0
    
    # change color of exons if color_exons_as_genes
    #for (i in 1:dim(custom_gff)[1]){
    #  count <- count + 1
    #  if(custom_gff[i,"feature"]=="mRNA"){
    #    gene_color <- custom_gff[i,"fill"]
    #  }
    #  if(custom_gff[i,"feature"]=="exon"){
    #    custom_gff[i,"fill"] <- gene_color
    #  }
    #}
    
    if(gene_annot=="custom"){annot_list[f] <- list(try(as.annotation(gff2annot(custom_gff,color=text_color,rot=0))))}
    if(gene_annot=="NA"){annot_list[f] <- list(NULL)}
    
    if(exon_mode=="gff3"){final_gff <- edited_gff}
    if(exon_mode=="custom"){final_gff <- custom_gff}
    # Add a "color" columnusing a given df, change to black if no color assigned
    if(sum(is.na(final_gff$fill>0))){final_gff[is.na(final_gff$fill),]$fill<-"black"}
    dna_segs[f] <- list(try(as.dna_seg(gff2seg(final_gff,add_intron,feature="exon",
                                               gene_type="arrows",color="black",
                                               lty=1,lwd=1,cex=1))))
    if (dim(mutation_gff)[1]>0){
      mutation_gff["fill"] <- "red"
      print(seg_names[f])
      print(dna_segs[f])
      dna_segs[[f]] <- rbind(dna_segs[f][[1]],gff2seg(mutation_gff,add_intron,feature="frameshift",
                                                      gene_type="starGrob",color="red",
                                                      lty=1,lwd=1,cex=1))
      print(dna_segs[f])
    }
  }
  names(dna_segs)<-seg_names # used to plot the names on the figure
  
  pdf(paste(file_name,".pdf",sep=""),width=plot_width,height=plot_height)
  plot_gene_map(dna_segs=dna_segs,
                seg_plot_height=plot_height,
                #dna_seg_labels = names(tree_barto$leaves),
                annotations=annot_list,
                annotation_cex = annotation_cex,
                annotation_height = annotation_height,
                comparisons=comp_list,
                #global_color_scheme = c("auto", "auto", "blue_red", 0.5),
                #tree=tree_barto,
                #gene_type = gene_type, 
                arrow_head_len=arrow_len,
                main="",
                #gene_type="side_blocks",
                xlims=xlims,
                dna_seg_scale=F,
                scale=F,
                fixed_gap_length = F)
  dev.off()
}

###################################################
### code chunk number 45: ex7_starGrob
###################################################
starGrob <- function(gene, ...){
  ## Coordinates for the star
  x <- sin(((0:5)/2.5)*pi)*(gene$end-gene$start-300)/2 + (gene$end+gene$start+0)/2
  y <- cos(((0:5)/2.5)*pi)*gene$strand*(0.5) + 1
  idx <- c(1, 3, 5, 2, 4, 1)
  ## Attribute line_col only if present in the gene
  line_col <- if (!is.null(gene$line_col)) gene$line_col else gene$col
  ## Having a conditional transparency, depending on a length cut-off
  ## passed via dots
  length_cutoff <- list(...)$length_cutoff
  if (!is.null(length_cutoff)){
    alpha <- if ((gene$end-gene$start) < length_cutoff)  0.3 else  0.8
  } else alpha <- 1
  
  ## Grobs
  g <- polygonGrob(x[idx], y[idx], gp=gpar(fill="black", col="black",
                                           lty=0, lwd=1, alpha=1),
                   default.units="native")
  t <- textGrob(label="***", x=(gene$end+gene$start)/2, y=0.5,
                default.units="native")
  print(gene)
  return(g)
}

# genoPlotR comparisons need to be pre-trimmed when using the xlim function
# This is due to it drawing a comparison between every xlim pair if blast 
# results spand multiple xlims. The trimming makes sure only 1 comparison per xlim pair
trim_comp <- function(comp, xlim1, xlim2){
  comp <- comp[order(comp$aln_len),]
  xlim_ranges <- seq(1, length(xlim1)-1, by=2)
  trimmed_df <- data.frame()
  for(row in row.names(comp)){
    x2_skip <- c()
    for(x1 in xlim_ranges){
      if(xlim1[x1]<=0){break}
      for(x2 in xlim_ranges){
        if(xlim2[x2]<=0){break}
        if(x2 %in% x2_skip){next}
        xs <- c(xlim1[x1],xlim1[x1+1],xlim2[x2],xlim2[x2+1],
                comp[row,"start1"],comp[row,"end1"],comp[row,"start2"],comp[row,"end2"])
        if(
          (xs[1]<=xs[5] & xs[2]<=xs[5]) | #
          (xs[1]>=xs[6] & xs[2]>=xs[6]) | #
          (xs[3]<=xs[7] & xs[4]<=xs[7]) | #
          (xs[3]>=xs[8] & xs[4]>=xs[8]) #
        ){next}
        temp_row <- comp[row,]
        temp_row$start1 <- ifelse(xlim1[x1]>=comp[row,"start1"],xlim1[x1],comp[row,"start1"])
        temp_row$end1 <-   ifelse(xlim1[x1+1]<=comp[row,"end1"],xlim1[x1+1],comp[row,"end1"])
        temp_row$start2 <- ifelse(xlim2[x2]>=comp[row,"start2"],xlim2[x2],comp[row,"start2"])
        temp_row$end2 <-   ifelse(xlim2[x2+1]<=comp[row,"end2"],xlim2[x2+1],comp[row,"end2"])
        trimmed_df <- rbind(trimmed_df,temp_row)
        x2_skip <- seq(1, x2, by=2)
        break
      }
    }
  }
  return(trimmed_df)
}

# generates an xlim list for the selected lines to plot
# do_gene_xlim T returns all xlims and F returns the start and end of loci
xlim_func <- function(seg_name, do_gene_xlim=T){
  xlims_list <- list()  # start a new xlim list
  xlims_df <- read.table("assigned_xlims.csv",header=T,sep=",") # default xlim table
  xlims_df <- xlims_df[xlims_df$genotype %in% seg_name,] # select which inbreds to plot
  freq_table <- table(xlims_df$genotype) # generates a frequency table that counts each gene
  for(i in 1:length(seg_name)){
    filtered_xlims <- xlims_df[xlims_df$genotype==seg_name[i],]
    temp_list <- c()
    for(n in 1:dim(filtered_xlims)[1]){
      #print(filtered_xlims)
      temp_list <- c(temp_list, filtered_xlims$left[n])
      temp_list <- c(temp_list, filtered_xlims$right[n])
    }
    ifelse(do_gene_xlim,
           xlims_list[i] <- list(c(temp_list, rep(c(-2,-1), max(freq_table)-as.numeric(freq_table[seg_name[i]])))),
           xlims_list[i] <- list(c(temp_list[1],temp_list[length(temp_list)]))
    )
  }
  return(xlims_list)
}

# to get the Mo17 gff3 file to work I had to manually replace ";Parent" with whatever
# random string and then replace all "Parent" with "Name". Then it seems to be working.
gff2seg <- function(gff,add_intron=T, feature="exon",gene_type="exons",color="black",lty=1,lwd=1,cex=1){
  gff <- gff[gff$feature==feature,]
  rownames(gff) <- NULL
  gff$X <- NULL
  if(feature %notin% c("mRNA","exon")){
    return(data.frame(name=gff$name[1], start=gff$start[1], end=gff$end[1], strand=gff$strand[1], length=(gff$end[1]-gff$start[1]),
                      feature="CDS", gene_type=gene_type, col=gff$fill[1], fill=gff$fill[1], lty=lty, lwd=lwd, pch=8, cex=cex))
  }
  
  # easier to rbind if dataframe already initialized with the first row
  seg_df <- data.frame(name=gff$name[1], start=gff$start[1], end=gff$end[1], strand=gff$strand[1], length=(gff$end[1]-gff$start[1]),
                       feature="CDS", gene_type=gene_type, col=gff$fill[1], fill=gff$fill[1], lty=lty, lwd=lwd, pch=8, cex=cex)
  
  for(i in 2:dim(gff)[1]){
    n <- dim(seg_df)[1]
    # adds to the dataframe an intron before adding the exon
    if(add_intron & seg_df$name[n]==gff$name[i]){
      seg_df <- rbind(seg_df,
                      data.frame(name=gff$name[i], start=seg_df$end[n]+1, end=gff$start[i]-1, strand=gff$strand[i], length=(gff$start[i]-seg_df$end[n]+2),
                                 feature="CDS_intron", gene_type="introns", col=color, fill="white", lty=lty, lwd=lwd, pch=8, cex=cex)
      )
    }
    seg_df <- rbind(seg_df,
                    data.frame(name=gff$name[i], start=gff$start[i], end=gff$end[i], strand=gff$strand[i], length=(gff$end[i]-gff$start[i]),
                               feature="CDS", gene_type=gene_type, col=gff$fill[i], fill=gff$fill[i], lty=lty, lwd=lwd, pch=8, cex=cex)
    )
  }
  if(feature=="mRNA"){seg_df$feature<-NULL}
  return(seg_df)
}


gff2annot <- function(gff,color="black", rot=0){
  # Create the annotation tables using the more general mRNA rows
  annot <- gff[gff[2]=="mRNA",][c(3,4,1)]
  rownames(annot) <- NULL
  colnames(annot) <- c("x1","x2","text")
  annot$"color" <- color
  annot$"rot" <- rot
  return(annot)
}

fix_exon_name <- function(gff, feature){
  # some of the gff3 files don't repeat the gene name in the exon so we have to fix it
  for(i in 1:dim(gff)[1]){
    if(gff$feature[i]==feature){name<-gff$name[i]}
    gff$name[i] <- name
  }
  return(gff)
}