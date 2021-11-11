`%notin%` <- Negate(`%in%`)

syntenyZplotUI <- function(id) {
  ns <- NS(id)
  
  lsf <- list.files(path = paste(getwd(),"custom_CDS/inbreds/", sep="/"), pattern = NULL)
  lsf <- sub('\\.fasta$', '', lsf)[1:10]
  
  navbarPage("Synteny Plot",
  #tabPanel("Synteny Plot",
    sidebarLayout(position="left", fluid=F,
      sidebarPanel(width = 4,
        
        checkboxInput(ns("do_gene_xlim"), "Do gene xlim:", value=F),
        numericInput(ns("plot_height"), "Plot Height:",value=500),
        numericInput(ns("plot_width"), "Plot Width:",value=300),
        numericInput(ns("low_color_percent"), "Select Minimum % Similarity:", value=90),
        actionButton(ns("update_button"), "Update plot"),
        # https://rstudio.github.io/sortable/articles/built_in.html
        bucket_list(
          header = "Drag the items in any desired bucket",
          group_name = "bucket_list_group",
          orientation = "horizontal",
          add_rank_list(
            text = "Drag from here",
            labels = lsf,
            input_id = ns("rank_list_1")
          ),
          add_rank_list(
            text = "to here",
            labels = NULL,
            input_id = ns("rank_list_2")
          )
        )
      ),
  
      mainPanel(width=8,
        tags$p("input$rank_list_1"),
        verbatimTextOutput(ns("results_1")),
        plotOutput(ns("plot") 
                   # width="600px", #paste(ns("input$plot_width"),"px",sep=""), 
                   # height="1200px" #paste(ns("input$plot_height"),"px",sep="")
        )
      )
    )
  )}
 
syntenyZplot <- function(input, output, session) {
  ns <- session$ns
  
  output$results_1 <-
    renderPrint(
      input$rank_list_1 # This matches the input_id of the first rank list
    )
  
  update_plot <- reactive({
    input$update_button
    do_gene_xlim <- isolate(input$do_gene_xlim)
    plot_height <- isolate(input$plot_height)
    plot_width <- isolate(input$plot_width)
    low_color_percent <- isolate(input$low_color_percent)
    # If the selected method is a single reference gene it returns the first gene even if given a gene list
    print(c(do_gene_xlim,plot_height,plot_width,low_color_percent))
    return(c(do_gene_xlim,plot_height,plot_width,low_color_percent))
  })
  
  plot_reactive <- reactive({
    MygenoPlotR(file_name="test_complete",
                blastall=blastall,
                do_gene_xlim=update_plot()[1],
                plot_height=update_plot()[2], # Don't think it is working
                plot_width=update_plot()[3],
                exon_mode="custom", # gff3 or exon
                gene_annot = "NA",
                annotation_cex = 1,
                annotation_height = 1,
                seg_names=input$rank_list_1,
                gene_colors=gene_colors,
                add_intron=F,
                annot=T,
                arrow_len = 10000,
                text_color = "black",
                low_color_percent = update_plot()[4])})
  
  # Call the ggplot function that creates that heatmap
  output$plot <- renderPlot(plot_reactive(),
                            width = function() update_plot()[2],
                            height = function() update_plot()[3])
  
}

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
  arrow_len=10000,   # change the size of the gene arrows
  text_color="black", # set NULL to remove labels,
  low_color_percent = 80
){
  
  # list of comparisons based on the blast results assuming you ran all-vs-all blast command
  # the coordinate system is not the same as the genomic coordinates. Adjust accordingly
  xlims <- xlim_func(seg_names,do_gene_xlim=do_gene_xlim)
  comp_list <- list()
  for (f in 1:(length(seg_names)-1)){
    comp_list[f] <- list(blastall[blastall["name1"]==seg_names[f] & blastall["name2"]==seg_names[f+1],])
    #comp_list[f] <- apply_color_scheme(comp_list[f]$per_id, )
    
    #if(do_gene_xlim){comp_list[[f]] <- trim_comp(comp_list[[f]], xlims[[f]], xlims[[f+1]])}
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
    custom_gff <- read.table(paste("custom_loci/custom_loci_",seg_names[f],".csv",sep=""), sep=",",header=T,row.names=NULL)
    mutation_gff <- custom_gff[custom_gff$feature  %notin% c("mRNA","exon"),]
    custom_gff <- custom_gff[custom_gff$feature  %in% c("mRNA","exon"),]
    custom_gff <- merge(custom_gff, gene_colors[,c("name","fill")], by = "name", all.x = TRUE)
    custom_gff <- custom_gff[with(custom_gff, order(custom_gff[,"X"])),]
    custom_gff <- subset(custom_gff, select = -X )
    count = 0
    
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
      dna_segs[[f]] <- rbind(dna_segs[f][[1]],gff2seg(mutation_gff,add_intron,feature="frameshift",
                                                      gene_type="starGrob",color="red",
                                                      lty=1,lwd=1,cex=1))
    }
  }
  names(dna_segs)<-seg_names # used to plot the names on the figure
  
  color_limit <- rep(1,13)
  color_limit[7] <- low_color_percent
  comp_list[[1]] <- rbind(comp_list[[1]], color_limit)
  
  
  #pdf(paste(file_name,".pdf",sep=""),width=plot_width,height=plot_height)
  return(plot_gene_map(dna_segs=dna_segs,
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
                fixed_gap_length = F))
  #dev.off()
}

gene_colors <- read.table("general/assigned_features.csv",sep=",",header=T)

#system("makeblastdb -dbtype nucl -in fasta_loci/combined.fasta -out fasta_loci/combined.db")
#system("blastn -db fasta_loci/combined.db -query fasta_loci/combined.fasta -out blastall.tsv -num_threads 4 -word_size 100 -evalue 10 -outfmt 6 ")
blastall <- read_comparison_from_blast("blastall.tsv")

xlim_func <- function(seg_name, do_gene_xlim=T){
  # generates an xlim list for the selected lines to plot
  # do_gene_xlim T returns all xlims and F returns the start and end of loci
  xlims_list <- list()  # start a new xlim list
  xlims_df <- read.table("general/assigned_xlims.csv",header=T,sep=",") # default xlim table
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