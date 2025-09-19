#' Title
#'
#' @param x1 a numeric vector of compositional data values
#' @param x2 another numeric vector of compositional data values
#' @param paired a logical indicating whether you want a paired test
#'
#' @return the p-value for the test
#' @export
#'
cd.test <- function(x1,x2,paired=FALSE){
  if(paired){
    if(dim(x1)[1]!=dim(x2)[1]){
      stop(paste0("Error: x1 (", ncol(x1), " columns) and x2 (", ncol(x2), " columns) must have the same number of columns"))
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
