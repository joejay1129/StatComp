#' @title self coded regression
#' @description apply standard,ridge or locally weighted regression(using Gaussian kernel) on a given data frame(the predicted variable has to be at the last column 
#' of the data frame like 'Boston' in MASS we used in vignettes )
#' @param a 0 to preform standard regression,1 for ridge regression and 2 for locally weighted regression
#' @param b lambda the multiple of identity matrix added to xTx in the process of ridge regression
#' @return the coefficient of the regression whose last column is the value of intercept
#' @examples
#' \dontrun{
#' coef <- add(1,2)
#' }
#' @importFrom Rcpp evalCpp
#' @export

add <-function(a,b){
    return(a+b)
}