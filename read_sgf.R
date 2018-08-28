require(rjson)
require(igraph)
require(plyr)
require(stringr)

readSGF <- function(infile) {
  
  jsongraph <- fromJSON(file=infile)
  
  nodes <- ldply(jsongraph$data$nodes, as.data.frame)
  names(nodes)[1] <- "name"

  edges <- ldply(jsongraph$data$edges, function(edge) {
  data.frame(source=edge$source, target=edge$target)
  })
  
  graph.data.frame(edges, (jsongraph$metadata$directed == "true"), nodes)
}

start <- function() {
  
  files <- list.files("graphs/", full.names = TRUE)
  
  sapply(files, function(file) {convertToGml(file, file)})
}