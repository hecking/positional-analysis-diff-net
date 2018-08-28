#=======================================================================================
katzMPA <- function(g) {
  g
}

#=======================================================================================
removeBotsLatencyBased <- function(g, minLatency) { 
  g %>% 
    delete_edges(which(E(.)$latencyMinutes < minLatency)) %>% 
    delete_vertices(which(degree(.) == 0))
}

#=======================================================================================
weightEdgesTemporally <- function(g) {
  
  endpoints <- ends(g, E(g))
  E(g)$latencyMinutes <- ((as.numeric(V(g)[endpoints[,2]]$timemillis) - 
                             as.numeric(V(g)[endpoints[,1]]$timemillis)) / 60000) + 1
  g
}


#=======================================================================================
weightEdgesTemporallyNormalized <- function(g) {  
  endpoints <- ends(g, E(g))
  latency <- ((as.numeric(V(g)[endpoints[,2]]$timemillis) - 
                 as.numeric(V(g)[endpoints[,1]]$timemillis)) / 60000) + 1 # +1 needed to avoid 0 latency
  E(g)$latencyMinutes <-   2*(-0.5+(1 + exp(-2/median(latency)*latency))^(-1))
  g
}

#=======================================================================================
getKatzMatrix <- function(g, alpha, addSourceSinks=TRUE, weightParam=NULL, reciprocalWeight=TRUE, ...) {
  
  if (addSourceSinks) {
    g <- prepareGraph(g, weightParam)
  }
  if (!is.null(weightParam) && reciprocalWeight) {    
    g <- set.edge.attribute(g, weightParam, value=1-get.edge.attribute(g, weightParam))
  }
  
  if (vcount(g) < 5000) {
    
    I <- diag(nrow=vcount(g))
    solve((I - alpha * get.adjacency(g, attr = weightParam)), ...) - I
  } else {
    
    # truncated katz
    a <- get.adjacency(g, attr = weightParam)
    a + a %*% a + a %*% a %*% a
  }
}

#=======================================================================================
getAccellerationCoefficient <- function(g, alpha, out=TRUE) {
  if (out) {
    rowSums(getKatzMatrix(g, alpha, "latencyMinutes")) - 1 # -1 corresponds to -I in the formula, missing in getKatzMatrix
  } else {
  }
}

#=======================================================================================
katzWeightedSPC <- function(katzMatrix, source, sink) {
  as.matrix(katzMatrix[source,]) %*% t(as.matrix(katzMatrix[,sink]))
}

#=======================================================================================
katzWeightedNPPC <- function(katzMatrix) {  
  as.matrix(colSums(katzMatrix)) %*% t(as.matrix(rowSums(katzMatrix)))
}