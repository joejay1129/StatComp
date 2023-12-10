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


initial<-function(n, B, seed_1){
    set.seed(seed_1)
    A <- matrix(runif(n * 2, min = -5, max = 5), n, 2)
    M0 <- A%*%t(B)
    
    return(M0)
}


initial_matrices <- function(n, m){
    matrices_initial <- vector("list", m)
    
    d <- round(2 * sqrt(n)) 
    set.seed(1)
    B <- matrix(runif(d * 2, min = -5, max = 5), d, 2)
    
    for (i in 1:m) {
        matrices_initial[[i]] <- initial(n, B = B, seed_1 = i)
    }
    
    return(matrices_initial)
}


SVD_self<-function(M){
    Sigma_hat<-t(M)%*%M
    V<-eigen(Sigma_hat)$vec
}

SVDt_self<-function(M){
    Sigma_hat<-M%*%t(M)
    U<-eigen(Sigma_hat)$vec
}

create <- function(p, M0, seed_2) {
    n <- nrow(M0)
    d <- ncol(M0)
    
    set.seed(seed_2)
    M <- M0 + matrix(rnorm(n * d, 0, 1), n, d)
    P <- matrix(data = rbinom(n * d, 1, p), n, d)
    M_obs <- M * P  # 使用矩阵操作替代循环更新 M_obs
    
    return(M_obs)
}


create_matrices <- function(p, matrices_initial, k = 1) {
    m <- length(matrices_initial)
    matrices_create <- vector("list", m)
    
    for (i in 1:m) {
        matrices_create[[i]] <- create(p = p, M0 = matrices_initial[[i]], seed_2 = k * m + i)
    }
    
    return(matrices_create)
}



rho<-function(U,V){
    norm((U%*%t(U)-V%*%t(V)),type="F")
}


#' @title simulation
#' @description Directly encapsulating the numerical simulation process, you can obtain the distance rho 
#' between the differences in the right singular vectors before and after the response by simply inputting the required parameters.
#' @param m Number of sub-servers
#' @param n data volume on each server
#' @param p missing probability p
#' @param monte number of Monte Carlo iterations
#' @return the common right singular vectors of input matrices
#' @import stats
#' @examples
#' \dontrun{
#' distribute_monte(m = 10, n = 100, p = 0.1, monte = 10)
#' }
#' @export
distribute_monte <- function(m, n, p, monte) {
    rho_ans <- rep(0, monte)
    d <- round(2 * sqrt(n)) 
    
    set.seed(1)
    B <- matrix(runif(d * 2, min = -5, max = 5), d, 2)
    V <- SVDt_self(B)[, 1:2]
    
    matrices_initial <- initial_matrices(n = n, m = m)
    for (i in 1:monte) {
        matrices_get <- create_matrices(p = p, matrices_initial = matrices_initial, k = i)
        V_res <- distribute(matrices_get)
        rho_ans[i] <- rho(V_res, V)
    }
    return(mean(rho_ans))
}