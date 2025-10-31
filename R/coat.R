
#' Covariance Estimation for Compositional Data
#' @description
#' \code{coat} is used to estimate the covariance matrix of basis data based on the compositional data.
#'
#' @param x n x p composition data matrix (row/column is sample/variable)
#' @param lam tuning parameter for thresholding
#' @param soft indicator for thresholding method: 1 for soft thresholding (default), 0 for hard thresholding
#'
#' @return a list containing following components:
#'    \item{sigma}{covariance matrix }
#'    \item{corr}{correlation matrx derived from the covariance matrix}
#'    \item{lam}{the tuning parameter used}
#'
#' @export
#'
coat <- function(x, lam, nlam=100, soft = 1){
  n <- nrow(x)
  p <- ncol(x)
  clrX <- log(x) - rowSums(log(x)) %*%matrix(1,1,p) / p
  cov <- cov(clrX)*(n-1)/n
  centered.x <- scale(clrX, scale = FALSE)
  theta <- (t(centered.x)^2)%*%(centered.x^2)/n - cov^2
  delta <- cov/(theta^0.5)
  delta <- abs(delta - diag(diag(delta)))

  if(missing(lam)){
    gridInfo <- adaptThresholdRange(clrX)
    lam <- gridInfo$lower + (gridInfo$upper - gridInfo$lower)*rep(1:nlam)/nlam
  }
  nlam = length(lam)
  sigma = array(dim=c(nlam, p, p))
  corr = array(dim=c(nlam, p, p))
  for(i in 1:nlam){
    sigmai <- adaptThreshold(cov, theta, lam[i], soft)
    sigma[i,,] <- sigmai
    corr[i,,] <- diag(diag(sigmai)^(-0.5))%*%sigmai%*%diag(diag(sigmai)^(-0.5))
  }

  return(list(sigma = sigma, corr = corr, lam = lam))
}
#' Select Tunning Parameter for COAT Via Cross-Validation
#' @description
#' \code{cv.coat} is used to choose tuning parameter \code{lam} for \code{coat} by the way of k-fold cross validation.
#' @param x n x p composition data matrix (row/column is sample/variable)
#' @param lam tuning parameters for thresholding. The function will choose the optimized one in the parameter sequence.. If missing, a sequence will be automatically generated based on data.
#' @param nfold number of folds for cross-validation (default: 5)
#' @param soft indicator for thresholding method: 1 for soft thresholding (default), 0 for hard thresholding
#'
#' @return a list containing following components:
#'   \item{sigma}{optimal covariance matrix selected via cross-validation}
#'   \item{corr}{correlation matrix derived from the optimal covariance matrix}
#'   \item{time}{execution time in seconds}
#'   \item{lambda}{optimal tuning parameter selected via cross-validation}
#' @export
#'
cv.coat <- function(x, lam, nlam=100, nfold = 5, foldid=NULL, soft = 1){
  startTime <- proc.time()
  p <- ncol(x)
  clrX <- log(x) - rowSums(log(x)) %*%matrix(1,1,p) / p

  if(missing(lam)){
    coatPred <- adaptThresoldCov(clrX, nFolder=nfold, foldid=foldid, nGrid=nlam, soft = soft, autoGrid = TRUE)
  }
  else{
    coatPred <- adaptThresoldCov(clrX, nFolder=nfold, foldid=foldid, nGrid=nlam, soft = soft, autoGrid = FALSE, grid = lam)
  }
  sigma <- coatPred$sigma
  corr <- coatPred$corr
  lambda_chosen <- coatPred$lambda
  exeTimeClass <- proc.time() - startTime
  exeTime <- as.numeric(exeTimeClass[3])
  return(list(sigma = sigma, corr = corr, time = exeTime, lam = lambda_chosen))
}

#----------------------------------------------------------------------------------------
#  Adaptive thresholding estimation of cov(x)
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#     nFolder ------ number of the foler in cross validation
#        soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#       sigma ------ covariance estimation based on adaptive thresholding
#        corr ------ correlation estimation based on adaptive thresholding
#----------------------------------------------------------------------------------------

adaptThresoldCov <- function(x, nFolder = 5, foldid=NULL, nGrid=100, soft = 1, autoGrid=TRUE, grid=c()){
  n <- nrow(x)
  p <- ncol(x)
  # Set the grid for the choice of tuning parameter

  gridInfo <- adaptThresholdRange(x)
  if(autoGrid){
    grid <- gridInfo$lower + (gridInfo$upper - gridInfo$lower)*rep(1:nGrid)/nGrid
  }
  else{
    grid = grid
  }
  # Multi-folder cross validation
  if(is.null(foldid)) {
    part <- 1 + sample(c(1:n))%%nFolder
  }
  else if(length(foldid) != n){
    warning('length of foldid is different from row numbers of x, foldid will not be used')
    part <- 1 + sample(c(1:n))%%nFolder
  }
  else{
    part <- foldid
  }
  error <- matrix(0, nFolder, nGrid)
  for (i in 1:nFolder){
    xTest <- x[which(part == i),]
    xTrain <- x[which(part != i),]
    gridInfoTrain <- adaptThresholdRange(xTrain)
    covTest <- cov(xTest)*(n-1)/n
    for (j in 1:nGrid){
      sigmaTrain <- adaptThreshold(gridInfoTrain$cov,gridInfoTrain$theta,grid[j],soft)
      error[i,j] <- (norm(sigmaTrain-covTest, "F"))
    }
  }
  errorSum <- colSums(error)
  lambda <- grid[which(errorSum == min(errorSum))][1]
  sigma <- adaptThreshold(gridInfo$cov,gridInfo$theta,lambda,soft)
  corr <- diag(diag(sigma)^(-0.5))%*%sigma%*%diag(diag(sigma)^(-0.5))
  return(list(sigma = sigma, corr = corr, lambda = lambda))
}
#----------------------------------------------------------------------------------------
#  Range of the tuning parameter
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#  Output:
#      A list structure contains:
#       upper ------ upper bound of tuning parameter
#       lower ------ lower bound of tuning parameter
#         cov ------ sample covariance of x
#       theta ------ sample variance of covariance
#----------------------------------------------------------------------------------------

adaptThresholdRange <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  cov <- cov(x)*(n-1)/n
  centered.x <- scale(x, scale = FALSE)
  theta <- (t(centered.x)^2)%*%(centered.x^2)/n - cov^2
  delta <- cov/(theta^0.5)
  delta <- abs(delta - diag(diag(delta)))
  upper <- max(delta)
  lower <- min(delta[which(delta != 0)])
  return(list(upper = upper, lower = lower, theta = theta, cov = cov))
}
#----------------------------------------------------------------------------------------
#  Apply adaptive thresholding to the sample covariance
#  Input:
#           cov ------ p x p covariance matrix
#         theta ------ p x p variance of covariance matrix
#        lambda ------ tuning parameter
#          soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#         sigma ------ p x p matrix, adaptive thresholding result
#----------------------------------------------------------------------------------------

adaptThreshold <- function(cov,theta,lambda,soft){
  covOffDiag <- cov - diag(diag(cov))
  thetaOffDiag <- theta - diag(diag(theta))
  sigmaTmp <- abs(covOffDiag) - lambda*thetaOffDiag^0.5
  sigmaTmp[which(sigmaTmp < 0)] <- 0
  if (soft == 1){
    sigma <- diag(diag(cov)) + sigmaTmp*sign(covOffDiag)
  }else{
    sigma <- cov
    sigma[which(sigmaTmp < 1e-10)] <- 0
    sigma <- sigma + diag(diag(cov))
  }
  return(sigma)
}
#----------------------------------------------------------------------------------------
