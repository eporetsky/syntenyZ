library(genoPlotR)
library(RColorBrewer)

gene_colors <- read.table("assigned_exons.csv",sep=",",header=T)
myColors <- brewer.pal(length(unique(gene_colors$symbol)),"Set1")
names(myColors) <- levels(factor(gene_colors$symbol))
gene_colors$fill <- myColors[factor(gene_colors$symbol)]

all_names <- c(
  "NC358","NC350","Ky21","CML69","B73","Oh43","M162W","CML333","CML277","W22",
  "HP301","Tx303","B97","CML247","CML228","CML103","Il14H","Tzi8","Mo18W","M37W",
  "CML322","Oh7B","Mo17","Ms71","CML52","Ki3","Ki11")

seg_name <- c(
  "ZxPI566673",
  "B73",
  "W22",
  "P39",
  "NC358",
  "NC350",
  "Ky21",
  "CML69",
  "Oh43",
  "M162W",
  "CML333",
  "CML277",
  "HP301",
  "Tx303",
  "B97",
  "CML247",
  "Tzi8",
  "Mo18W",
  "Mo17",
  "Il14H",
  "M37W",
  "CML322",
  "Oh7B",
  "Ms71",
  "CML52",
  "Ki3",
  "Ki11",
  "CML228",
  "CML103"
)

all_names <- seg_name

#### script for alternating lines
#target = "ZxPI566673"
target = "W22"
seg_name <- c()
for (name in all_names){
  if(name!=target){
    seg_name <- c(seg_name,name)
    seg_name <- c(seg_name,target)
  }
}
system("cat fasta_loci/*.fasta > fasta_loci/combined.fasta")
system("makeblastdb -dbtype nucl -in fasta_loci/combined.fasta -out fasta_loci/combined.db")
system("blastn -db fasta_loci/combined.db -query fasta_loci/combined.fasta -out blastall.tsv -num_threads 4 -word_size 125 -evalue 10 -outfmt 6 ")
blastall <- read_comparison_from_blast("blastall.tsv")

MygenoPlotR(file_name="test",
            blastall=blastall,
            do_gene_xlim=T,
            plot_height=30, # Don't think it is working
            plot_width=10,
            mode="exons", # gff3 or exon
            seg_names=seg_name,
            gene_colors=gene_colors,
            add_intron=T,
            #xlims=xlims_list, #xlims_list. NULL is set for between first and last gff3 element
            annot=F,
            10000,
            NA)

MygenoPlotR <- function(
  file_name,         # choose the file-name to save the PDF figure
  blastall,
  do_gene_xlim=F,
  plot_height=3,
  plot_width=15,
  mode="exon",       # use gff3 or custom coordinate table
  seg_names,         # ordered names of the genomes to parse their files
  gene_colors,       # a table of elements and their specified colors
  add_intron=T,
  annot=F,
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
    # takes gff3 and returns 5 column coordinate tabe
    if(mode=="gff3"){
      edited_gff <- gff_converter(seg_names[f])
      edited_gff <- merge(edited_gff, gene_colors[,c("name","fill")], by = "name", all.x = TRUE)
    }else{
      edited_gff <- read.table(paste("custom_loci/custom_coordinate_",seg_name[f],".csv",sep=""), sep=",",header=T,row.names=NULL)
      edited_gff <- merge(edited_gff, gene_colors[,c("name","fill")], by = "name", all.x = TRUE)
      edited_gff <- edited_gff[with(edited_gff, order(edited_gff[,"X"])),]
      edited_gff <- subset(edited_gff, select = -X )
      count = 0
      for (i in 1:dim(edited_gff)[1]){
        if(edited_gff[i,"name"]=="tent_gene"){
          count <- count + 1
          edited_gff[i,"name"] <- paste("gene",as.character(count),sep="")
        }
        edited_gff[i,"name"] <- paste("gene",as.character(count),sep="")
      }
    }
    #print(edited_gff)
    #return(list(edited_gff,gene_colors[,c("name","fill")]))
    # Add a "color" columnusing a given df, change to black if no color assigned

    if(sum(is.na(edited_gff$fill>0))){edited_gff[is.na(edited_gff$fill),]$fill<-"black"}
    dna_segs[f] <- list(try(as.dna_seg(gff2seg(edited_gff,add_intron,feature="exon",
                                               gene_type="arrows",color="black",
                                               lty=1,lwd=1,cex=1))))
    ifelse(annot,
           annot_list[f] <- list(try(as.annotation(gff2annot(edited_gff,color=text_color,rot=25)))),
           annot_list[f] <- list(NULL))

    #dna_segs[f] <- list(try(read_dna_seg_from_tab("test_x.csv")))
  }
  #print(annot_list)
  names(dna_segs)<-seg_names # used to plot the names on the figure

  #?plot_gene_map

  pdf(paste(file_name,".pdf",sep=""),width=plot_width,height=plot_height)
  plot_gene_map(dna_segs=dna_segs,
                seg_plot_height=plot_height,
                #dna_seg_labels = names(tree_barto$leaves),
                annotations=annot_list,
                comparisons=comp_list,
                global_color_scheme = c("auto", "auto", "blue_red", 0.5),
                #tree=tree_barto,
                arrow_head_len=arrow_len,
                main="",
                #gene_type="side_blocks",
                xlims=xlims,
                dna_seg_scale=T,
                scale=F,
                fixed_gap_length = T)
  dev.off()
  return(list(dna_segs, annot_list))
}

trim_comp <- function(comp, xlim1, xlim2){
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


xlim_func <- function(seg_name, do_gene_xlim=T){
  # generates an xlim list for the selected lines to plot
  # do_gene_xlim T returns all xlims and F returns the start and end of loci
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
gff_converter <- function(name){
  name_field <- ifelse(name %in% c("W22","B73","Mo17","B73v4"),"Name","Parent")
  gff <- read.table(paste("gff/",name,".gff3",sep=""), header=F,sep="\t")
  gff$V9 <- gsub("transcript:", "", gff$V9)
  gff$V9 <- gsub("gene:", "", gff$V9)
  gff[4:5] <- gff[4:5] - min(gff[4:5]) + 1 # normalize coordinates to first element
  gff <- gff[c(9,3,4,5,7)]                 # reorder columns and filter some of the columns
  gff <- gff[gff$V3 %in% c("mRNA","exon"),]
  parents <- unlist(strsplit(gff$V9,split=";")) # returns a list of all gff3 fields in the last column
  gff$V9 <- parents[grepl(name_field, parents)] # pick only the items in the list starting with Name or Parent
  # Stupid function to add _T001 if it doesn't excist already...
  gff$V9[!grepl("_",gff$V9)] <- unlist(lapply(gff$V9[!grepl("_",gff$V9)],paste,"_T001",sep=""))
  gff$V9 <- unlist(lapply(strsplit(gff$V9,split=";"),`[[`,1))
  gff$V9 <- unlist(lapply(strsplit(gff$V9,split="="),`[[`,2))
  gff$transcript <- unlist(lapply(strsplit(gff$V9,split="_"),`[[`,2))
  gff$transcript <- unlist(lapply(strsplit(gff$transcript,split=".exon"),`[[`,1))
  gff$V9 <- unlist(lapply(strsplit(gff$V9,split="_"),`[[`,1))
  gff <- unique(gff[gff$transcript=="T001",]) # In B73 I get duplicated (pressumably T001 & T002)
  colnames(gff) <- c("name","feature","start","end","strand","transcript")
  gff$strand[gff$strand=="+"] <- 1
  gff$strand[gff$strand=="-"] <- -1
  rownames(gff) <- NULL
  if(name %in% c("W22","Mo17")){gff <- fix_exon_name(gff,"mRNA")}
  #if(name %in% c("Mo17")){gff <- fix_exon_name(gff,"gene")}
  return(gff[c("name","feature","start","end","strand")])
}

gff2seg <- function(gff,add_intron=T, feature="exon",gene_type="exons",color="black",lty=1,lwd=1,cex=1){
  gff <- gff[gff$feature==feature,]
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


### not actually working, but maybe parts will be useful
trim_comp <- function(comp, xlim1, xlim2){
  xlim_ranges <- seq(1, length(xlim1)-2, by=2)
  smallest <- Inf

  for(n in xlim_ranges){
    if(xlim1[n+3]<0) break
    if(xlim2[n+3]<0) break

    temp_smallest <- xlim1[n+2] - xlim1[n+1]
    if(smallest>temp_smallest){smallest<-temp_smallest}
    temp_smallest <- xlim2[n+2] - xlim2[n+1]
    if(smallest>temp_smallest){smallest<-temp_smallest}
  }
  print(smallest)

  trimmed_df <- data.frame()
  for(row in row.names(comp)){
    if(comp[row,"aln_len"]>smallest){
      div_by <- round(ceiling(comp[row,"aln_len"]/smallest))
      comp_row <- comp[row,]
      range1 <- unlist(lapply(seq(comp[row,"start1"],comp[row,"end1"],length.out=div_by),round))
      range2 <- unlist(lapply(seq(comp[row,"start2"],comp[row,"end2"],length.out=div_by),round))
      range1
      for(i in 1:(div_by-1)){
        temp_row <- comp_row
        temp_row$start1 <- range1[i]
        temp_row$end1 <- range1[i+1]-1
        temp_row$start2 <- range2[i]
        temp_row$end2 <- range2[i+1]-1
        trimmed_df <- rbind(trimmed_df,temp_row)
      }
    }else{
      trimmed_df <- rbind(trimmed_df,comp[row,])
    }
  }
  print(trimmed_df)
  return(trimmed_df)
}
