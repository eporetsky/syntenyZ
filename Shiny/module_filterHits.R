filterHitsUI <- function(id) {
  ns <- NS(id)
  tabPanel("Filter Hits",
     sidebarLayout(position="left", fluid=F,
       sidebarPanel(width = 4,
          p("Current project selected:"),
          textOutput({ns("current_project_name")}),
          actionButton(ns("open_raw_results"), "Open raw results"),
          actionButton(ns("open_filtered_results"), "Open filtered results"),
          hr(),
          actionButton(ns("delete_rows"), "Delete selected rows"),
          actionButton(ns("save_filtered_results"), "Save filtered results"),
        ),

  mainPanel(width=8,
    fluidRow(
      column(6, style = "",
        dataTableOutput(ns("blast_table"))
      )
    )
  )
  ))}

filterHits <- function(input, output, session, project_name) {
  ns <- session$ns
  
  output$current_project_name <- renderText({project_name()})
  
  working_table <- reactiveVal()
  observeEvent(input$delete_rows, {
    # https://stackoverflow.com/questions/39136385/delete-row-of-dt-data-table-in-shiny-app
    if(!is.null(input$blast_table_rows_selected)){
      #print(working_table())
      working_table(working_table()[-as.numeric(input$blast_table_rows_selected),])
    }
  })
  
  observeEvent(input$open_raw_results, {
    working_table(read.table(paste("projects/",project_name(),"/blast_results.tsv", sep=""), sep="\t"))
    output$blast_table <- renderDataTable(working_table(),
      options = list(pageLength = 100, width="100%", scrollX = TRUE))
  })
  
  observeEvent(input$open_filtered_results, {
    working_table(read.table(paste("projects/",project_name(),"/blast_results_filtered.tsv", sep=""), sep="\t"))
    output$blast_table <- renderDataTable(working_table(),
                                          options = list(pageLength = 100, width="100%", scrollX = TRUE))
  })
  
  observeEvent(input$save_filtered_results, {
    write.table(working_table(), paste("projects/",project_name(),"/blast_results_filtered.tsv", sep=""), 
                sep="\t", quote=F, row.names=F, col.names=F)
  })
}
