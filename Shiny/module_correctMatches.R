correctMatchesUI <- function(id) {
  ns <- NS(id)
  tabPanel("Correct Annotation",
     
     sidebarLayout(position="left", fluid=F,
       sidebarPanel(width = 4,
          htmlOutput(ns("genotype_list_select")),
          htmlOutput(ns("select_row")),
          actionButton(ns("add_row"), "Add empty row below selected"),
          actionButton(ns("delete_row"), "Delete selected row"),
          #actionButton(ns("extract_synteny"), "Extract synteny sequences"),
       ),
                   
       mainPanel(width=8,
           fluidRow(
           column(6, style = "",
                  DTOutput(ns('custom_loci_table'))
                  #dataTableOutput(ns("synteny_table"))
            )
          )
       )
))}

correctMatches <- function(input, output, session, project_name) {
  ns <- session$ns
  
  output$current_project_name <- renderText({project_name()})
  
  genotype_list <- reactive(sub('\\.fasta$', '', list.files(path = paste("projects/",project_name(), "/custom_CDS/inbreds/", sep=""), pattern = NULL)))
  # Get the list of available genotypes.  Should be the same as custom_loci 
  output$genotype_list_select <- renderUI({
    selectizeInput(ns("select_genotype"), "Select Genotype:", choices=genotype_list(), multiple=FALSE)
  })
  
  row_list <- reactive(rownames(temp_table()))
  output$select_row <- renderUI({
    selectizeInput(ns("select_row"), "Select Row:", choices=c(0, row_list()), multiple=FALSE)
  })
  
  # Save the custom_loci table edits in this reactiveValues
  temp_file <-  reactiveVal()
  temp_table <- reactiveVal()
  observeEvent(input$select_genotype, {
    temp_file(paste("projects/",project_name(), "/custom_loci/custom_loci_", input$select_genotype, ".csv", sep=""))
    temp_table(read.csv(temp_file()))
  })
  
  # 
  observe({
    output$custom_loci_table = DT::renderDT(DT::datatable(
      temp_table(),
      options = list(
        paging = F
      ),
      rownames = FALSE,
      selection = 'none',
      editable = TRUE
    ))
  })
  #?editData
  # Event for editing specific cells that overwrites the custom_loci file and triggers renderDT observer
  observeEvent(input$custom_loci_table_cell_edit, {
    tmp <- editData(temp_table(), input$custom_loci_table_cell_edit, 'custom_loci_table', rownames=FALSE)
    colnames(tmp)[1] <- ""
    temp_table(tmp)
    write.csv(temp_table(), temp_file(), quote=F, row.names = F)
  })
  
  
  
  observeEvent(input$add_row, {
    # https://stackoverflow.com/questions/11561856/add-new-row-to-dataframe-at-specific-row-index-not-appended
    tmp <- temp_table()
    select_row <- as.numeric(input$select_row)
    tmp <- rbind(tmp[1:select_row,], NA, tmp[-(1:select_row),])
    rownames(tmp) <- NULL # reindexes the row names
    tmp[1] <- rownames(tmp) # reassign the first column new row names
    temp_table(tmp)
    colnames(tmp)[1] <- ""
    temp_table(tmp)
    write.csv(temp_table(), temp_file(), quote=F, row.names = F)
  })
  
  observeEvent(input$delete_row, {
    tmp <- temp_table()
    select_row <- as.numeric(input$select_row)
    tmp <- tmp[-as.numeric(select_row),]
    rownames(tmp) <- NULL # reindexes the row names
    tmp[1] <- rownames(tmp) # reassign the first column new row names
    temp_table(tmp)
    colnames(tmp)[1] <- ""
    temp_table(tmp)
    write.csv(temp_table(), temp_file(), quote=F, row.names = F)
  })
  
  
  
  
  
  
  
  #working_table <- reactiveVal()
  #observeEvent(input$extract_synteny, {
  #  print("Delete previous fasta")
  #  system(paste("rm projects/", project_name(), "/fasta_loci/combined.fasta", sep=""))
  #  system(paste("rm projects/", project_name(), "/fasta_loci/combined_masked.fasta", sep=""))
  #  print("start extract_results")
  #  system(paste("python3 global/extract_results.py", project_name(), sep=" "))
  #  print("read raw_overlap_list")
  #  working_table(read.table(paste("projects/", project_name(),"/raw_overlap_list.tsv", sep=""), sep="\t"))
  #  output$synteny_table <- renderDataTable(working_table(),
  #                                        options = list(pageLength = 100, width="100%", scrollX = TRUE))
  #})
}