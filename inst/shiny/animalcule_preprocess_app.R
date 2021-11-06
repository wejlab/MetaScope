library(shiny)
library(shinythemes)
library(Biostrings)
library(tools)
library(Rsubread)
library(Rsamtools)
library(taxize)
library(rsconnect)

source(file.path("R", "demultiplex.R"),  local = TRUE)
source(file.path("R", "align_target.R"),  local = TRUE)
source(file.path("R", "download_refseq.R"),  local = TRUE)

ui <- navbarPage(
  title = paste("Animalcules Preprocess", sep = ""),
  id = "Animalcules Preprocess",
  fluid = TRUE,
  theme = "bootstrap.min.css",
  source(file.path("inst","shiny","ui", "UI_demultiplex.R"),  local = TRUE)$value,
  source(file.path("inst","shiny","ui", "UI_ref_seq"),  local = TRUE)$value,
  source(file.path("inst","shiny","ui", "UI_align_target"),  local = TRUE)$value,
  footer = includeHTML("www/footer.html")
)

server <- function(input, output, session) {
  source(file.path("inst","shiny","server", "server_demultiplex.R"),  local = TRUE)$value,
  source(file.path("inst","shiny","server", "server_ref_seq"),  local = TRUE)$value,
  source(file.path("inst","shiny","server", "server_align_target"),  local = TRUE)$value
}

shinyApp(ui = ui, server = server)

