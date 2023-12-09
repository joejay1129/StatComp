#' @title IPW algorithm
#' @description The full name of the IPW algorithm is Inverse Probability Weighting algorithm, which aims to 
#' handle missing matrices in order to obtain their right singular vectors.
#' @param M missing matrices
#' @param K the rank of missing matrices
#' @return the right singular vectors of truly M
#' @export

IPW <- function(M, K) {
    n <- nrow(M)
    d <- ncol(M)
    Y <- ifelse(M == 0, 0, 1)  # 使用向量化操作代替循环
    
    p_hat <- sum(Y) / (n * d)
    
    Sigma_phat <- crossprod(M) - (1 - p_hat) * diag(colSums(M^2))
    
    eigen_Sigma_phat <- eigen(Sigma_phat)
    
    # 得到分解
    V <- eigen_Sigma_phat$vectors[, 1:K]
    
    return(V = V)
}


#' @title distribute algorithm
#' @description Aiming to recover the common right singular vectors of the input matrix columns.
#' @param matrices_get the list of input matrices
#' @return the common right singular vectors of input matrices
#' @export

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