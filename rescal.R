source("NMF/read_sgf.R")
source("NMF/data_cleaning.R")
source("NMF/acceleration.R")
source("NMF/regularised_network_decomposition.R")
source("KatzNMFUtils.R")
require(dplyr)
require(igraph)
require(NMF)

updateR_NN <- function(X, A, R, lambda=0) {
  
  R * ((t(A) %*% X %*% A) / 
    (t(A) %*% A %*% R %*% t(A) %*% A + (lambda * R)))
}

updateR <- function(X, A) {
  
  Z <- kronecker(A, A)
  
  l <-  Matrix::solve(t(Z) %*% Z, sparse=TRUE)
  r <- t(Z) %*% as.vector(X)
  R <- l %*% r
  matrix(R, ncol=ncol(A), nrow=ncol(A))
}

updateA <- function(adj, A, R, eta=0) {
  
  E <- adj %*% A %*% t(R) + t(adj) %*% A %*% R
  B <- R %*% (t(A) %*% A) %*% t(R)
  C <- t(R) %*% (t(A) %*% A) %*% R
  
  BpC <- B + C
  if (eta > 0) {
    
    normalisation <- A %*% BpC + eta
    A <- A * (E / normalisation)
  } else {
    
    tryCatch({
      A <- t(solve(t(BpC), t(E)))
    }, error = function(err) {
      
      print("Matrix is not invertable. Calculate pseudo inverse")
      A <<- E %*% ginv(BpC)
    })
  }
  
  A
}

updateA_NN <- function(X, A, R, lambda=0) {
  
  A * ((X %*% A %*% t(R) + t(X) %*% A %*% R) /
    (A %*% ((R %*% t(A) %*% A %*% t(R) + t(R) %*% t(A) %*% A %*% R) + diag(lambda, nrow(R)))))
}
computeFit <- function(adj, A, R) {
  
  adjApprox <- A %*% R %*% t(A)
  norm((adj - adjApprox), type="F")^2
}

rescal <- function(adj, k, A, R) {
  
  relations <- updateR_NN(adj, A, R)
  
  fit <- computeFit(adj, A, relations)
  
  fits <- sapply(1:150, function(i) {
    
    # optimise A
    A <<- updateA_NN(adj, A, relations)
    
    # update density matrix
    relations <<- updateR_NN(adj, A, relations)
    
    fit <- computeFit(adj, A, relations)
    
    #print(paste("iteration", i, "fit:", fit))
    fit
  })
  
  list(A=A, R=relations, fits=fits)  
}

densityMatrix <- function(A, membershipMatrix) {
  
  mem <- sapply(1:nrow(A), function(i) {
    which(membershipMatrix[i,] == max(membershipMatrix[i,]))[1]
  })
  
  sapply(1:ncol(membershipMatrix), function(i) {
    
    sapply(1:ncol(membershipMatrix), function(j) {
      # 1 - (length(which(A[mem == i,mem == j] == 0)) / 
      #        length(A[mem == i,mem == j]))
      m <- mean(A[mem == i,mem == j])
      ifelse(is.na(m), 0, m)
    })
  }) %>% t
}

startEvaluationDylan <- function() {
  
  dylan <- readSGF(
    "NMF/Knockin_on_Heaven_s_Door_1480892960751_COMPLETE_DISCOURSE.sgf"
  ) %>%
    cleanData(1.2, TRUE, TRUE) %>% 
    weightEdgesTemporallyNormalized %>%
    extractLargestComponent %>% 
    getKatzMatrix(alpha=0.7, weightParam = "latencyMinutes", addSourceSinks = FALSE)
    
  dylan[dylan < 0] <- 0 # Some numerical inaccuracy can lead so very small negative values.
  
  init <- seed(as.matrix(dylan), 10, method="nndsvd", densify="random")
  R <- t(basis(init)) %*% t(coef(init))
  # coef(init)[1,] <- 0
  # coef(init)[,colSums(dylan) == 0] <- c(1, rep(0, 9))
  # basis(init)[,10] <- 0
  # basis(init)[rowSums(dylan) == 0,] <- 0
  # basis(init)[rowSums(dylan) == 0,10] <- 1

  rescalCoefInit <- rescal(dylan, 10, t(coef(init)) + 0.0001, R)
  densRescalCoefInit <- densityMatrix(dylan, rescalCoefInit$A)
  corRescalCoefInit <- cor(as.vector(densRescalCoefInit), as.vector(rescalCoefInit$R))

  rescalBasisInit <- rescal(dylan, 10, basis(init) + 0.0001, R)
  densRescalBasisInit <- densityMatrix(dylan, rescalBasisInit$A)
  corRescalBasisInit <- cor(as.vector(densRescalBasisInit), as.vector(rescalBasisInit$R))
  print(corRescalBasisInit)
  print(corRescalCoefInit)
  # nmfResults <- regNetworkDecompositionSparse(dylan, init, lambda=0)
  # relations <- nmfResults$bin %*% nmfResults$bOut
  # 
  # densNMFBasis <- densityMatrix(dylan, nmfResults$bOut)
  # corNMFBasis <- cor(as.vector(densNMFBasis), as.vector(relations))
  # 
  # densNMFCoef <- densityMatrix(dylan, t(nmfResults$bin))
  # corNMFCoef <- cor(as.vector(densNMFCoef), as.vector(relations))
  # 
  # regNmfResults <- regNetworkDecompositionSparse(dylan, init, lambda = 0.5)
  # relations2 <- regNmfResults$bin %*% regNmfResults$bOut
  # 
  # densNMFBasis2 <- densityMatrix(dylan, regNmfResults$bOut)
  # corNMFBasis2 <- cor(as.vector(densNMFBasis2), as.vector(relations2))
  # 
  # densNMFCoef2 <- densityMatrix(dylan, t(regNmfResults$bin))
  # corNMFCoef2 <- cor(as.vector(densNMFCoef2), as.vector(relations2))
  # 
  # data_frame(rescalCoef=corRescalCoefInit, rescalBasis=corRescalBasisInit, nmfBasis=corNMFBasis, nmfCoef=corNMFCoef,
  #      regNmfBasis=corNMFBasis2, regNmfCoef=corNMFCoef2)
}

startEvaluationSchiaparelli <- function() {
  
  schiaparelle <- readSGF(
    "NMF/Schiaparelli_1480313605884_1_correctedAscii.sgf"
  ) %>%
    cleanData(1.2, TRUE, TRUE) %>% 
    weightEdgesTemporallyNormalized %>%
    extractLargestComponent %>% 
    getKatzMatrix(alpha=0.7, addSourceSinks = FALSE)
  
  schiaparelle[schiaparelle < 0] <- 0 # Some numerical inaccuracy can lead so very small negative values.
  
  init <- seed(as.matrix(schiaparelle), 7, method="nndsvd", densify="random")
  coef(init)[,colSums(schiaparelle) == 0] <- c(1, rep(0, 6))
  basis(init)[rowSums(schiaparelle) == 0,] <- 0
  basis(init)[rowSums(schiaparelle) == 0,7] <- 1
  
  rescalCoefInit <- rescal(schiaparelle, 7, t(coef(init)))
  densRescalCoefInit <- densityMatrix(schiaparelle, rescalCoefInit$A)
  corRescalCoefInit <- cor(as.vector(densRescalCoefInit), as.vector(rescalCoefInit$R))
  
  rescalBasisInit <- rescal(schiaparelle, 7, basis(init))
  densRescalBasisInit <- densityMatrix(schiaparelle, rescalBasisInit$A)
  corRescalBasisInit <- cor(as.vector(densRescalBasisInit), as.vector(rescalBasisInit$R))
  
  nmfResults <- regNetworkDecompositionSparse(schiaparelle, init, lambda=0)
  relations <- nmfResults$bin %*% nmfResults$bOut
  
  densNMFBasis <- densityMatrix(schiaparelle, nmfResults$bOut)
  corNMFBasis <- cor(as.vector(densNMFBasis), as.vector(relations))
  
  densNMFCoef <- densityMatrix(schiaparelle, t(nmfResults$bin))
  corNMFCoef <- cor(as.vector(densNMFCoef), as.vector(relations))
  
  regNmfResults <- regNetworkDecompositionSparse(schiaparelle, init, lambda = 0.5)
  relations2 <- regNmfResults$bin %*% regNmfResults$bOut
  
  densNMFBasis2 <- densityMatrix(schiaparelle, regNmfResults$bOut)
  corNMFBasis2 <- cor(as.vector(densNMFBasis2), as.vector(relations2))
  
  densNMFCoef2 <- densityMatrix(schiaparelle, t(regNmfResults$bin))
  corNMFCoef2 <- cor(as.vector(densNMFCoef2), as.vector(relations2))
  
  list(rescalCoef=corRescalCoefInit, rescalBasis=corRescalBasisInit, nmfBasis=corNMFBasis, nmfCoef=corNMFCoef,
       regNmfBasis=corNMFBasis2, regNmfCoef=corNMFCoef2)
}

getClusterInfo <- function(g, membershipMatrixIn, membershipMatrixOut) {
  
  memIn <- sapply(1:nrow(membershipMatrixIn), function(i) {
    which(membershipMatrixIn[i,] == max(membershipMatrixIn[i,]))[1]
  })
  memOut <- sapply(1:nrow(membershipMatrixOut), function(i) {
    which(membershipMatrixOut[i,] == max(membershipMatrixOut[i,]))[1]
  })
  igraph::as_data_frame(g, what="vertices") %>% 
    mutate(clusterIn=memIn, clusterOut=memOut) %>%
    group_by(clusterIn, clusterOut) %>%
    do({
      types <- table(.$type) %>% dplyr::as_data_frame() %>% spread(key=Var1, value=n)
      medTime <- as.POSIXct(median(as.numeric(.$timemillis)) / 1000, origin="1970-01-01")
      types %>% mutate(median_time=medTime)
    })
}

getGraphSchiaprelli <- function() {
  
  g <- readSGF(
    "NMF/Schiaparelli_1480313605884_1_correctedAscii.sgf"
  ) %>%
    cleanData(1.2, TRUE, TRUE) %>% 
    weightEdgesTemporallyNormalized %>%
    extractLargestComponent 
    
  schiaparelle <- g %>% getKatzMatrix(alpha=0.7, addSourceSinks = FALSE)
  
  schiaparelle[schiaparelle < 0] <- 0 # Some numerical inaccuracy can lead so very small negative values.
  
  init <- seed(as.matrix(schiaparelle), 7, method="nndsvd", densify="random")
  coef(init)[,colSums(schiaparelle) == 0] <- c(1, rep(0, 6))
  basis(init)[rowSums(schiaparelle) == 0,] <- 0
  basis(init)[rowSums(schiaparelle) == 0,7] <- 1
 
  
  nmfResults <- regNetworkDecompositionSparse(schiaparelle, init, lambda=0)
  relations <- nmfResults$bin %*% nmfResults$bOut
  
  
  list(nmf=nmfResults, cinfo=getClusterInfo(g, t(nmfResults$bin), nmfResults$bOut), 
       g=graph_from_adjacency_matrix(as.matrix(relations), weighted = TRUE))
}