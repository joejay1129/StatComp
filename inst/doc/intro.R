## -----------------------------------------------------------------------------
add <- function(a,b){
    return(a+b)
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

