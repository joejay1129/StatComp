## -----------------------------------------------------------------------------
IPW <- function(M, K) {
    n <- nrow(M)
    d <- ncol(M)
    Y <- ifelse(M == 0, 0, 1) 
    
    p_hat <- sum(Y) / (n * d)
    
    Sigma_phat <- crossprod(M) - (1 - p_hat) * diag(colSums(M^2))
    
    eigen_Sigma_phat <- eigen(Sigma_phat)
    
    V <- eigen_Sigma_phat$vectors[, 1:K]
    
    return(V = V)
}

## -----------------------------------------------------------------------------
distribute <- function(matrices_get) {
    n <- nrow(matrices_get[[1]])
    d <- ncol(matrices_get[[1]])
    m <- length(matrices_get)
    
    Sigma_tilde <- matrix(0, d, d)
    V <- vector("list", m)
    
    for (i in 1:m) {
        V[[i]] <- IPW(matrices_get[[i]])
        Sigma_tilde <- Sigma_tilde + (1 / m) * V[[i]] %*% t(V[[i]])
    }
    
    V_res <- eigen(Sigma_tilde)$vec[, 1:2]
    return(V_res)
}

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction(code='double  targetoutofn(NumericVector x, int target, int n) {
    double prob[n+1][target+1];
    prob[0][0] = 1;
    for(int i = 1; i <= n; i++){
        prob[i][0] = prob[i-1][0]*(1-x[i-1]);
    }
    for(int i = 1; i <= n; i++){
        for(int k = 1; k <= i && k <= target; k++){
            prob[i][k] = prob[i-1][k] * (1-x[i-1]) + prob[i-1][k-1] * x[i-1];
        }
    }
    return prob[n][target];
    }')
prob <- c(0.5, 0.5, 0.5, 0.5, 0.5)
target <- 1
n <- 5
targetoutofn(prob, target, n)

