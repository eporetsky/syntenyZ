blastSeqsUI <- function(id) {
  ns <- NS(id)
  
  lsf <- list.files(path = paste(getwd(),"BLASTdb/", sep="/"), pattern = "db.nhr$")
  lsf <- sub('\\.db.nhr$', '', lsf)
  
  tabPanel("Data Input",
    sidebarLayout(position="left", fluid=F,
      sidebarPanel(width = 4,
        p("Use to setup blastdb and run queries"),
        selectizeInput(ns("select_genotypes"), "Select Genotypes:",
                       choices=lsf, multiple=TRUE),
        actionButton(ns("save_genotypes"), "Save Genotypes"),
        textInput(ns("project_name"), "Select project name:", value = ""),
        actionButton(ns("project_save"), "Save Genotypes"),
        textOutput({ns("current_project_name")}),
        fileInput(ns("query_fasta"), "Load query fasta:", accept = c(".fasta")),
        p(), hr(), p(),
        actionButton(ns("blast_query"), "Blast query"),
        textOutput({ns("finished_blast_query")})
        
      ),
       
      mainPanel(width=8,
       fluidRow(
         column(6, style = "",
            p("Future main panel"),
         )
       )
     )
))}

blastSeqs <- function(input, output, session) {
  ns <- session$ns
  
  select_genotypes <- reactiveVal()
  observeEvent(input$save_genotypes, {
    #print(input$select_genotypes)
    select_genotypes(input$select_genotypes)
  })
  
  project_name <- reactiveVal()
  observeEvent(input$project_save, {
    system(paste("python3 global/create_project.py", input$project_name, sep=" "))
    output$current_project_name <- renderText(paste("Current project is: ", input$project_name))
    project_name(input$project_name)
  })
  
  #project_name <- reactive(output$current_project_name)
  
  #print(project_name())
  
  observeEvent(input$query_fasta, {
    project_name <- isolate(input$project_name)
    if(project_name==c('')) return()
    if(is.null(input$query_fasta)) return()
    paste("projects/", project_name, sep="")
    file.copy(input$query_fasta$datapath, paste("projects/", project_name, sep=""), overwrite=T)
    file.rename(paste("projects/", project_name, "/0.fasta", sep=""),
                paste("projects/", project_name, "/query.fasta", sep=""))
  })
  
  observeEvent(input$blast_query, {
    project_name <- isolate(input$project_name)
    selected_genotypes <- select_genotypes()
    selected_str <- selected_genotypes[1]
    for(name in selected_genotypes[2:length(selected_genotypes)]){
      selected_str <- paste(selected_str, name, sep=".db ")
    }
    selected_str <- paste(selected_str, ".db", sep="")
    #print(paste("python global/blast_project.py", project_name, selected_str, sep=" "))
    system(paste("python global/blast_project.py", project_name, selected_str, sep=" "))
    output$finished_blast_query <- renderText(paste("Finished blast query ", project_name))
  })
  
  
  return(project_name)
}