require(cluster)
require(NMF)
require(dplyr)
require(igraph)

source("NMF/read_sgf.R")
source("NMF/data_cleaning.R")
source("NMF/acceleration.R")
source("NMF/regularised_network_decomposition.R")
source("KatzNMFUtils.R")

getEdgeBundles <- function(g, vertexAttribute, k) {
  
  nodeDf <- igraph::as_data_frame(g, "vertices")
  
  endpoints <- ends(g, E(g))
  
  edgeDf <- data_frame(v1=nodeDf[endpoints[,1], vertexAttribute], v2=nodeDf[endpoints[,2], vertexAttribute])
  
  #sim <- -as.matrix(dist(edgeDf))
  d <- dist(edgeDf)
  membership <- pam(d, k)$clustering
  
  lapply(unique(membership), function(m) {
    edgeBundle <- which(membership == m)
    
    endpoints <- ends(g, edgeBundle)
    
    list(head=unique(endpoints[,2]), tail=unique(endpoints[,1]))
  })
}

edgeBundleIntitialisation <- function(X, k) {
  
  g <- graph_from_adjacency_matrix(X, weighted = TRUE)
  V(g)$acc <- rowSums(X)
  
  edgeBundles <- getEdgeBundles(g, "acc", k)
  
  W <- sapply(edgeBundles, function(bundle) {
    
    outAffinities <- rep(0, nrow(X))
    outAffinities[bundle$tail] <- sapply(bundle$tail, function(node) {
      
      sqrt(mean(X[node,bundle$head]))
    })
    outAffinities
  })
  
  H <- sapply(edgeBundles, function(bundle) {
    
    inAffinities <- rep(0, nrow(X))
    inAffinities[bundle$head] <- sapply(bundle$head, function(node) {
      
      sqrt(mean(X[bundle$tail, node]))
    })
    inAffinities
  }) %>% t
  
  nmfModel(k, X, W=W, H=H)
}

getNMFModel <- function(X, k, ...) {
  
  nmf(X, k, .options="t", ...)
}

nmfToMembership <- function(W, H) {
  
  outMembership <- sapply(1:nrow(W), function(row) {
    
    which(W[row,] == max(W[row,]))[1]
  })
  
  inMembership <- sapply(1:ncol(H), function(col) {
    
    which(H[,col] == max(H[,col]))[1]
  })
  
  list(inm=inMembership, outm=outMembership)
}

visualiseNodeTypeClusterAssociations <- function(W, H, types) {
  
  cnames <- paste("C", 1:ncol(W), sep="")
  colnames(W) <- paste("C", 1:ncol(W))
  rownames(H) <- paste("C", 1:nrow(H))
  
  p1 <- as.matrix(W) %>% dplyr::as_data_frame() %>%
    mutate(type=types) %>%
    tidyr::gather("Cluster", "Association", -type) %>%
    group_by(Cluster, type) %>%
    summarise_all(sum) %>%
    ggplot(aes(x = Cluster, y = type)) +
    geom_tile(aes(fill = Association))
  
  p2 <- as.matrix(t(H)) %>% dplyr::as_data_frame() %>%
    mutate(type=types) %>%
    tidyr::gather("Cluster", "Association", -type) %>%
    group_by(Cluster, type) %>%
    summarise_all(sum) %>%
    ggplot(aes(x = Cluster, y = type)) +
    geom_tile(aes(fill = Association))
  
  Rmisc::multiplot(p1, p2, cols=2)
}

getClusterDiffusionGraph <- function(W, H, minWeight=1) {
  
  (H %*% W) %>% 
    graph_from_adjacency_matrix(weighted = TRUE) #%>%
    #mst(weights = 1/E(.)$weight)
}

compareClusterings <- function(...) {
  
  clusterings <- list(...)
  
  pairs <- combn(1:length(clusterings), 2)
  
  clSim <- sapply(1:ncol(pairs), function(col) {
    
    mclust::adjustedRandIndex(clusterings[[pairs[1,col]]], clusterings[[pairs[2, col]]])
  })
  
  data_frame(clusterA=pairs[1,], clusterB=pairs[2,], sim=clSim)
}

getArrangedKatzMatrix <- function(km, W, H) {
  require(reshape2)
  require(ggplot2)
  mem <- nmfToMembership(W, H)
  orderOut <- order(mem$outm)
  orderIn <- order(mem$inm)
  arrangedKM <- km[orderOut, orderIn]
  rownames(arrangedKM) <- orderOut
  colnames(arrangedKM) <- orderIn
  
  clSizesOut <- table(mem$out)
  interceptsY <- sapply(1:(length(clSizesOut) - 1), function(i) {
    
    Reduce(sum, clSizesOut[1:i], 0.5)
  })
  
  clSizesIn <- table(mem$inm)
  interceptsX <- sapply(1:(length(clSizesIn) - 1), function(i) {
    
    Reduce(sum, clSizesIn[1:i], 0.5)
  })
    
  arrangedKM %>% 
    as.matrix %>%
    melt(value.name = "coupling") %>% 
    ggplot(aes(x=as.character(Var2), y=as.character(Var1), fill=coupling)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      scale_x_discrete(limits=orderIn) +
      scale_y_discrete(limits=orderOut) +
      geom_vline(xintercept = interceptsX) +
      geom_hline(yintercept = interceptsY)
}

getKatzClusterCorrelation <- function(km, W, H) {
  
  clDiff <- H %*% W
  
  # ++ couplings
  mem <- nmfToMembership(W, H)
  
  blockDensitiesOut <- tidyr::crossing(source=unique(mem$outm), target=unique(mem$outm)) %>%
    rowwise %>%
    do({
      srcs <- mem$outm == .$source
      trgts <- mem$outm == .$target
      data_frame(src=.$source, tar=.$target, density=sum(km[srcs, trgts]) / (length(which(srcs)) * length(which(trgts))))
    })
  
  blockDensitiesIn <- tidyr::crossing(source=unique(mem$outm), target=unique(mem$outm)) %>%
    rowwise %>%
    do({
      srcs <- mem$inm == .$source
      trgts <- mem$inm == .$target
      data_frame(src=.$source, tar=.$target, density=sum(km[srcs, trgts]) / (length(which(srcs)) * length(which(trgts))))
    })
  
  data_frame(
    corrOut=cor.test(as.vector(t(clDiff)), blockDensitiesOut$density)$estimate,
    corrIn=cor.test(as.vector(t(clDiff)), blockDensitiesIn$density)$estimate
  )
}

evaluateK <- function(km, krange) {
  
  res <- nmf(km, krange, method="lee", seed="nndsvd")
  X <- km[-nrow(km), -ncol(km)]
  katzClusterCorr <- data_frame(fit=res$fit) %>%
    rowwise %>%
    do({
    
      getKatzClusterCorrelation(X, basis(.$fit)[-nrow(basis(.$fit)), ], coef(.$fit)[, -ncol(coef(.$fit))])  
    }) %>%
    bind_cols(res$measures)
}

setupDylan <- function() {
  dylan <- readSGF(
    "NMF/Knockin_on_Heaven_s_Door_1480892960751_COMPLETE_DISCOURSE.sgf"
  ) %>%
    cleanData(1.2, TRUE, TRUE) %>% 
    weightEdgesTemporallyNormalized %>%
    extractLargestComponent %>% 
    getKatzMatrix(alpha=0.7, weightParam = "latencyMinutes", addSourceSinks = FALSE)
  
  dylan[dylan < 0] <- 0 # Some numerical inaccuracy can lead so very small negative values.
  
  init <- seed(as.matrix(dylan), 10, method="nndsvd", densify="random")
  coef(init)[1,] <- 0
  coef(init)[,colSums(dylan) == 0] <- c(1, rep(0, 9))
  basis(init)[,10] <- 0
  basis(init)[rowSums(dylan) == 0,] <- 0
  basis(init)[rowSums(dylan) == 0,10] <- 1
  
  nmfResults <- regNetworkDecompositionSparse(dylan, init, lambda=0)
}
