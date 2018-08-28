require(NMF)

imbalance <- function(X, seed, lambda) {

  if (lambda > 0) {
    lambda * sum(diag(t(basis(seed)) %*% X %*% t(coef(seed))))  
  } else {
    0
  }
  
  # sapply(1:nrow(basis(seed)), function(i) {
  #   
  #   sapply(1:ncol(coef(seed)), function(j) {
  #     
  #     (basis(seed)[i,] - coef(seed)[,j])^2 * X[i,j]
  #   })
  # }) %>% sum
}

imbalanceSparse <- function(X, Bout, Bin, lambda) {
  
  if (lambda > 0) {
    lambda * sum(diag(t(Bout) %*% X %*% t(Bin))  )
  } else {
    0
  }
} 

regNetworkDecomposition <- function(X, seed, lambda=1) {
  
  eta <- 0.000001 #ifelse(lambda == 0, 0.000001, 0)
  Dr <- diag(rowSums(X))
  Dc <- diag(colSums(X))
  it <- 0

  resNew <- norm(X - basis(seed) %*% coef(seed), type="F")^2 + imbalance(X, seed, lambda)
  repeat({
    it <- it + 1
    resOld <- resNew
    B <- basis(seed)
    C <- coef(seed)
    
    basis(seed) <- basis(seed) * (((1 + lambda) * (X %*% t(coef(seed)) + eta)) / 
                                    ((basis(seed) %*% coef(seed) %*% t(coef(seed)) + (lambda * (Dr %*% basis(seed)))) + eta))
    
    coef(seed) <- coef(seed) * (((1 + lambda) * (t(basis(seed)) %*% X + eta)) / 
                                  ((t(basis(seed)) %*% basis(seed) %*% coef(seed) + (lambda * (coef(seed) %*% Dc))) + eta))
    
    resNew <- norm(X - basis(seed) %*% coef(seed), type="F")^2 + imbalance(X, seed, lambda)
    seed <- trackError(seed, resNew, it, force=TRUE)

    if ((resOld - resNew) < 0.0002) {
      break
    }
  })

  seed
}

# The method is quicker than the normal decomposition because it makes use of sparse matrices, but is does not 
# comply to the interface of the NMF package anymore.
networkDecompositionSparse <- function(X, seed) {
  
  require(Matrix)
  # sparsify
  X <- Matrix(X, sparse=TRUE)
  Bout <- Matrix(basis(seed), sparse=TRUE)
  Bin <- Matrix(coef(seed), sparse=TRUE)
  
  eta <- 0.000001 #ifelse(lambda == 0, 0.000001, 0)
  Dr <- diag(rowSums(X))
  Dc <- diag(colSums(X))
  it <- 0
  
  resNew <- norm(X - Bout %*% Bin, type="F")^2
  repeat({
    
    it <- it + 1
    resOld <- resNew
    
    Bout <- Bout * (((X %*% t(Bin) + eta)) / 
                      ((Bout %*% Bin %*% t(Bin)) + eta))
    
    Bin <- Bin * (((t(Bout) %*% X + eta)) / 
                    ((t(Bout) %*% Bout %*% Bin) + eta))
    
    resNew <- norm(X - Bout %*% Bin, type="F")^2
    
    if ((resOld - resNew) < 0.0002) {
      break
    }
  })
  
  list(bOut=Bout, bin=Bin, iterations=it, fit=resNew)
}

regNetworkDecompositionSparse <- function(X, seed, lambda=1) {
  
  if (lambda == 0) { # avoid unnecessary computations
    
    networkDecompositionSparse(X, seed)
  } else {
    
    require(Matrix)
    # sparsify
    X <- Matrix(X, sparse=TRUE)
    Bout <- Matrix(basis(seed), sparse=TRUE)
    Bin <- Matrix(coef(seed), sparse=TRUE)
    
    eta <- 0.000001 #ifelse(lambda == 0, 0.000001, 0)
    Dr <- diag(rowSums(X))
    Dc <- diag(colSums(X))
    it <- 0
    
    resNew <- norm(X - Bout %*% Bin, type="F")^2 + imbalanceSparse(X, Bout, Bin, lambda)
    repeat({
      
      it <- it + 1
      resOld <- resNew
      
      Bout <- Bout * (((1 + lambda) * (X %*% t(Bin) + eta)) / 
                        ((Bout %*% Bin %*% t(Bin) + (lambda * (Dr %*% Bout))) + eta))
      
      Bin <- Bin * (((1 + lambda) * (t(Bout) %*% X + eta)) / 
                      ((t(Bout) %*% Bout %*% Bin + (lambda * (Bin %*% Dc))) + eta))
      
      resNew <- norm(X - Bout %*% Bin, type="F")^2 + imbalanceSparse(X, Bout, Bin, lambda)
      
      if ((resOld - resNew) < 0.0002) {
        break
      }
    })
    
    list(bOut=Bout, bin=Bin, iterations=it, fit=resNew)  
  }
}