#' Two Samples Test for Compositional Data
#' @description
#' \code{cd.test} is used to assess whether there are significant differences between two compositional data samples.
#' @param x1 a numeric vector of compositional data values
#' @param x2 another numeric vector of compositional data values
#' @param paired logical indicator for whether you want a paired test
#'
#' @return the p-value for the test
#' @export
#'
cd.test <- function(x1,x2,paired=FALSE){

  x1 = as.matrix(x1)
  x2 = as.matrix(x2)

  if(ncol(x1)!=ncol(x2)) stop("Error: x1 and x2 must have the same number of columns")
  p = ncol(x1)

  x1 <- log(x1)
  x2 <- log(x2)

  x1 <- x1 - 1/p*rowSums(x1)%*%matrix(1,1,p)
  x2 <- x2 - 1/p*rowSums(x2)%*%matrix(1,1,p)

  if(paired){
    if(dim(x1)[0]!=dim(x2)[0]){
      stop(paste0("Error: for paired test, x1 (", nrow(x1), " rows) and x2 (", nrow(x2), " rows) must have the same number of rows"))
    }
    n <- dim(x1)[1]
    p <- ncol(x1)
    SqStandizeDiff <- ((colSums(x1-x2))/n)^2/((diag(var(x1-x2)))*(n-1)/n^2)
    StatX <- max(SqStandizeDiff)
    pvalue <- 1-exp(-1/sqrt(pi)*exp(-(StatX-(2*log(p)-log(log(p))))/2))

  }
  else {
    n1 <- dim(x1)[1]
    n2 <- dim(x2)[1]
    p <- ncol(x1)
    x1mean <- colSums(x1)/n1
    x2mean <- colSums(x2)/n2
    x1Var <- diag(var(x1))*(n1-1)/n1
    x2Var <- diag(var(x2))*(n2-1)/n2
    xStatVar <- (x1Var*n1 + x2Var*n2)/(n1*n2)
    xStat <- max((x1mean/sqrt(xStatVar) - x2mean/sqrt(xStatVar))^2)
    pvalue <- 1-exp(-1/sqrt(pi)*exp(-(xStat-(2*log(p)-log(log(p))))/2))

  }
  return(pvalue)
}
