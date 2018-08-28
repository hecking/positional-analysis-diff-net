getKatzMatrix <- function(g, alpha, addSourceSinks=TRUE, weightParam=NULL, reciprocalWeight=TRUE, ...) {
  
  if (addSourceSinks) {
    g <- prepareGraph(g, weightParam)
  }
  if (!is.null(weightParam) && reciprocalWeight) {    
    g <- set.edge.attribute(g, weightParam, value=1-get.edge.attribute(g, weightParam) + 0.001)
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

extractLargestComponent <- function(g) {
  
  # Get components
  comps <- components(g)
  lccNodes <- which(comps$membership == which(comps$csize == max(comps$csize)[1]))
  
  induced_subgraph(g, lccNodes)
}

prepareGraph <- function(g, weightAttr=NULL) {
  
  wAtt <- list()
  if (!is.null(weightAttr)) {
    
    wAtt[[weightAttr]] <- 0
  }
  
  g %>% 
    add_vertices(1) %>%
    add_edges(as.vector(rbind(vcount(.), which(degree(., mode="in") == 0))), attr = wAtt) %>%
    add_vertices(1) %>%
    add_edges(as.vector(rbind(vcount(.), which(degree(., mode="out") == 0))), attr = wAtt)
}

densityMatrixOutIn <- function(A, membershipMatrixOut, membershipMatrixIn) {
  
  memIn <- sapply(1:nrow(A), function(i) {
    which(membershipMatrixIn[,i] == max(membershipMatrixIn[,i]))[1]
  })
  memOut <- sapply(1:nrow(A), function(i) {
    which(membershipMatrixOut[i,] == max(membershipMatrixOut[i,]))[1]
  })
  
  sapply(1:ncol(membershipMatrixOut), function(i) {
    
    sapply(1:nrow(membershipMatrixIn), function(j) {
      # 1 - (length(which(A[mem == i,mem == j] == 0)) / 
      #        length(A[mem == i,mem == j]))
      m <- mean(A[memOut == i,memIn == j])
      ifelse(is.na(m), 0, m)
    })
  }) %>% t
}

nmfToHardClustering <- function(Bout, Bin) {
  
  outMembership <- sapply(1:nrow(Bout), function(row) {
    
    which(Bout[row,] == max(Bout[row,]))[1]
  })
  
  inMembership <- sapply(1:ncol(Bin), function(col) {
    
    which(Bin[,col] == max(Bin[,col]))[1]
  })
  
  lapply(1:ncol(Bout), function(cl) {
    
    list(positive=which(outMembership == cl), negative=which(inMembership == cl))
  })
}