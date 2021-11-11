filterHitsUI <- function(id) {
  ns <- NS(id)
  tabPanel("Filter Hits",
     sidebarLayout(position="left", fluid=F,
       sidebarPanel(width = 4,
          p("Use to setup blastdb and run queries"),
          textOutput({ns("current_project_name")})
        ),

  mainPanel(width=8,
         fluidRow(
           column(6, style = "",
                  p("Future main panel"),
           )
         )
    )
  ))}

filterHits <- function(input, output, session, project_name) {
  ns <- session$ns
  
  output$current_project_name <- renderText({project_name()})

}