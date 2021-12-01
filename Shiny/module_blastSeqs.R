#genotype_list <- reactive(sub('\\.fasta$', '', list.files(path = paste("projects/",project_name(), "/custom_CDS/inbreds/", sep=""), pattern = NULL)))

# genotype_list <- sub('\\.fasta$', '', list.files(path = paste("BLASTdb/",project_name(), 
#                                                             "/custom_CDS/inbreds/", sep=""), pattern = NULL)

blastSeqsUI <- function(id) {
  ns <- NS(id)
  
  
  tabPanel("Data Input",
    sidebarLayout(position="left", fluid=F,
      sidebarPanel(width = 4,
        p("Use to setup blastdb and run queries"),
                   
        htmlOutput(ns("projectList")),
        
        selectizeInput(ns("select_genotypes"), "Select Genotypes:",
                       choices=genotypeList(), multiple=TRUE),
        actionButton(ns("save_genotypes"), "Save Genotypes"),
        
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
  
  project_name <- reactiveVal()
  observeEvent(input$project_save, {
    project_name(input$project_name)
    # Since projectList only updates when Shiny app starts
    projectList(c(projectList(), input$project_name))
  })

  observeEvent(input$project_delete, {
    # Since projectList only updates when Shiny app starts
    projectList(projectList()[!(projectList() %in% input$project_name)])
    if(project_name() %notin% c("PAN1", "Bx1")){system(paste("rm -r ", getwd(), "/projects/", project_name(), sep=""))}
    project_name(projectList()[1])
  })
  
  projectList <- reactiveVal(list.files(path = paste(getwd(),"projects/", sep="/")))
  output$projectList <- renderUI({
    tagList(
      selectizeInput(ns("select_project"), "Select Project:",
                     choices=projectList(), multiple=FALSE),
      textInput(ns("project_name"), "Select project name:", value = ""),
      actionButton(ns("project_save"), "Delete Project"),
      actionButton(ns("project_delete"), "Delete Project"),
      textOutput({ns("current_project_name")})
  )})
  
  observeEvent(input$select_project, {
    project_name(input$select_project)
    output$current_project_name <- renderText(paste("Current project is: ", input$select_project))
    system(paste("python3 global/create_project.py", project_name(), sep=" "))
  })
  
  select_genotypes <- reactiveVal()
  observeEvent(input$save_genotypes, {
    #print(input$select_genotypes)
    select_genotypes(input$select_genotypes)
  })
  
  
  #project_name <- reactive(output$current_project_name)
  
  #print(project_name())
  
  observeEvent(input$query_fasta, {
    project_name <- isolate(project_name())
    if(project_name==c('')) return()
    if(is.null(input$query_fasta)) return()
    paste("projects/", project_name, sep="")
    file.copy(input$query_fasta$datapath, paste("projects/", project_name, sep=""), overwrite=T)
    file.rename(paste("projects/", project_name, "/0.fasta", sep=""),
                paste("projects/", project_name, "/query.fasta", sep=""))
  })
  
  observeEvent(input$blast_query, {
    project_name <- isolate(project_name())
    selected_genotypes <- select_genotypes()
    selected_str <- selected_genotypes[1]
    for(name in selected_genotypes[2:length(selected_genotypes)]){
      selected_str <- paste(selected_str, name, sep=".db ")
    }
    selected_str <- paste(selected_str, ".db", sep="")
    #print(paste("python global/blast_project.py", project_name, selected_str, sep=" "))
    system(paste("rm projects/", project_name, "/blast_results.tsv", sep=""))
    system(paste("python3 global/blast_project.py", project_name(), selected_str, sep=" "))
    output$finished_blast_query <- renderText(paste("Finished blast query ", project_name))
  })
  
  
  return(project_name)
}

genotypeList <- function(){
  # Get the list of genotypes
  BLASTdb_list <- list.files(path = paste(getwd(),"BLASTdb/", sep="/"), pattern = "db.nhr$")
  return(sub('\\.db.nhr$', '', BLASTdb_list))
}