?selectInput
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
        fileInput(ns("query_fasta"), "Load query fasta:", accept = c(".fasta"))
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
  
  observeEvent(input$save_genotypes, {
    print(input$select_genotypes)
  })
  
  observeEvent(input$project_save, {
    system(paste("python Global/create_project.py", input$project_name, sep=" "))
    output$current_project_name <- renderText(paste("Current project is: ", input$project_name))
    #print(input$query_fasta)
  })
  
  observeEvent(input$query_fasta, {
    project_name <- isolate(input$project_name)
    if(project_name==c('')) return()
    if(is.null(input$query_fasta)) return()
    paste("Projects/", project_name, sep="")
    file.copy(input$query_fasta$datapath, paste("Projects/", project_name, sep=""), overwrite=T)
    file.rename(paste("Projects/", project_name, "/0.fasta", sep=""),
                paste("Projects/", project_name, "/query.fasta", sep=""))
  })
}