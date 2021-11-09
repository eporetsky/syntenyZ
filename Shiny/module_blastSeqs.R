blastSeqsUI <- function(id) {
  ns <- NS(id)
  navbarMenu("Data Input",
  #tabPanel("Data Input",
    tabPanel("Query"),
    tabPanel("Check BlastDB"),
    tabPanel("Create BlastDB"),
    
    sidebarLayout(position="left", fluid=F,
      sidebarPanel(width = 4,
        p("Use to setup blastdb and run queries"),
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
}