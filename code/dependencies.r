load_dependencies <<- function() {
  # Load dependencies
  library(utils)
  if (!require(shiny)) {
    install.packages("shiny")
  }
  library(shiny, quietly = TRUE)
  if (!require(shinyBS)) {
    install.packages("shinyBS")
  }
  if (!require(raster)) {
    install.packages("raster")
  }
  if (!require(rgeos)) {
    install.packages("rgeos")
  }
  if (!require(grDevices)) {
    install.packages("grDevices")
  }
  if (!require(RColorBrewer)) {
    install.packages("RColorBrewer")
  }
  if (!require(devtools)) {
    install.packages("devtools")
  }
}
