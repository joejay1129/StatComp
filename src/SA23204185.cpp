#include <Rcpp.h>
using namespace Rcpp;

//' @title A multi-prob calculator for Deliveryman Assignment Problem
//' @description A multi-prob calculator using Rcpp
//' @param x the probability vector of each deliveryman
//' @param target there are excatly target number of deliverymen take order
//' @param n the length of vector x
//' @return the objective probability
//' @examples
//' \dontrun{
//' prob <- c(0.5, 0.5, 0.5, 0.5, 0.5)
//' target <- 1
//' n <- 5
//' DAP(prob, target, n)
//' }
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
NumericVector DAP(NumericVector x, int target, int n) {
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
}






