######################################################
## syntenyZ Version 0.1
## Written by Elly Poretsky
## Alisa Huffaker Lab
## University of California, San Diego
######################################################
# Install packages that are not yet installed and loads them
# Based on - https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
package.check <- lapply(
  c("shiny","shinythemes", "genoPlotR", "RColorBrewer", "sortable"),
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#source("global.R")                # Loads default_files.csv to specifit which files to load when MutRank starts
source("module_blastSeqs.R")        # Shiny module that handles loading user specified data for the analysis
#source("module_filterHits.R")     # Shiny module that handles the caluclation and printing of the Mutual Rank table 
#source("module_correctMatches.R") # Shiny module that handles the plotting part of the Mutual Rank table as a heatmap
source("module_syntenyZplot.R")   # Shiny module that handles plotting the Mutual Rank table as a network using JS
#source("syntenyZ_functions.R")  # R script file that contains the functions used by the different modules

# Shiny UI function
ui <- fluidPage(
  navbarPage("syntenyZ v0.1",
             theme = shinytheme("flatly"),
             blastSeqsUI("blastSeqsNS"),
             #filterHitsUI("filterHitsNS"),
             #correctMatchesUI("correctMatchesNS"),
             syntenyZplotUI("syntenyZplotNS")
  ),
  tags$head(tags$style('body {color:black; font-size: 175%;}')),
  tags$head(tags$style(HTML('.shiny-bound-input {color:black; font-size: 100%;}')))
)

# Shiny server function
server <- function(input, output, session){
  # Stops the running Shiny app when the app window is closed
  session$onSessionEnded(stopApp)
  
  # Calls module_datainput.R and returns a list containing all the loaded data
  callModule(blastSeqs,"blastSeqsNS")
  callModule(syntenyZplot,"syntenyZplotNS")
  
}

# Create and runs the Shiny app
shinyApp(ui, server)