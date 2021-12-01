######################################################
## syntenyZ Version 0.2
## Written by Elly Poretsky
## Alisa Huffaker Lab
## University of California, San Diego
######################################################
# Install packages that are not yet installed and loads them
# Based on - https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
#system("python3 --version")

package.check <- lapply(
  c("shiny","shinythemes", "genoPlotR", "RColorBrewer", "sortable", "DT"),
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
#source("global.R")                # Loads default_files.csv to specifit which files to load when MutRank starts
source("module_blastSeqs.R")        # Shiny module that handles loading user specified data for the analysis
source("module_filterHits.R")     # Shiny module that handles the caluclation and printing of the Mutual Rank table 
source("module_correctMatches.R") # Shiny module that handles the plotting part of the Mutual Rank table as a heatmap
source("module_syntenyZplot.R")   # Shiny module that handles plotting the Mutual Rank table as a network using JS
#source("syntenyZ_functions.R")  # R script file that contains the functions used by the different modules

# Shiny UI function
ui <- fluidPage(
  navbarPage("syntenyZ v0.2",
             theme = shinytheme("flatly"),
             blastSeqsUI("blastSeqsNS"),
             filterHitsUI("filterHitsNS"),
             correctMatchesUI("correctMatchesNS"),
             syntenyZplotUI("syntenyZplotNS")
  ),
  tags$head(tags$style('body {color:black; font-size: 175%;}')),
  tags$head(tags$style(HTML('.shiny-bound-input {color:black; font-size: 100%;}')))
)

# Shiny server function
server <- function(input, output, session){
  # Stops the running Shiny app when the app window is closed
  session$onSessionEnded(stopApp)
  
  # blastqSeqs prepares BLASTdb files and initiates syntenyZ projects
  project_name <- callModule(blastSeqs,"blastSeqsNS")
  
  # filterHits accesses the project blast results for the selected rows
  callModule(filterHits, "filterHitsNS", project_name)
  
  # correctMatches accesses the projected selected results for final corrections
  callModule(correctMatches, "correctMatchesNS", project_name)
  
  # syntenyZplot accesses the selected project and generates the synteny plot
  callModule(syntenyZplot, "syntenyZplotNS", project_name)
  
}

# Create and runs the Shiny app
shinyApp(ui, server)