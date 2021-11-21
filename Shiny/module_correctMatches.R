correctMatchesUI <- function(id) {
  ns <- NS(id)
  tabPanel("Correct Annotation",
     sidebarLayout(position="left", fluid=F,
       sidebarPanel(width = 4,
                    actionButton(ns("extract_synteny"), "Extract synteny sequences"),
       ),
                   
       mainPanel(width=8,
           fluidRow(
           column(6, style = "",
                  dataTableOutput(ns("synteny_table"))
            )
          )
       )
))}

correctMatches <- function(input, output, session, project_name) {
  ns <- session$ns
  
  output$current_project_name <- renderText({project_name()})
  
  working_table <- reactiveVal()
  observeEvent(input$extract_synteny, {
    system(paste("rm projects/", project_name(), "/fasta_loci/combined.fasta", sep=""))
    system(paste("rm projects/", project_name(), "/fasta_loci/combined_masked.fasta", sep=""))
    system(paste("python3 global/extract_results.py", project_name(), sep=" "))
    working_table(read.table(paste("projects/", project_name(),"/raw_overlap_list.tsv", sep=""), sep="\t"))
    output$synteny_table <- renderDataTable(working_table(),
                                          options = list(pageLength = 100, width="100%", scrollX = TRUE))
  })
}