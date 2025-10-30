#' @useDynLib compotools cdmm_c
#' @keywords internal
NULL

#' Variable selection in regression with compositional variables
#' @description
#' \code{cdmm} is used to perform variable selection and fit a linear log-contrast model where the predictor variables are compositional data.
#' @param y numeric vector for response variable
#' @param x n x p logged composition data matrix for predictor variable (row/column is sample/variable)
#' @param lam user-specified sequence of tuning parameters for regularization. The function will fit a model for each lam in the sequence. If missing, an automatic sequence is generated based on data (see \code{nlam} and \code{rlam}). Use \code{cv.cdmm} or \code{gic.cdmm} to choose optimised parameter.
#' @param nlam number of lambda values in the automatic sequence (default: 100)
#' @param rlam ratio of the smallest to largest lambda in the automatic sequence (default: 1/nlam)
#' @param mu regularization parameter for the compositional constraint (default: 1)
#' @param std logical indicator for data standardization: TRUE (default) standardizes predictors
#'            and centers the response; FALSE uses raw data
#' @param maxv maximum number of variables allowed in the model (default: 0.4*length(y))
#' @param maxit vector of two integers specifying maximum iterations for outer and inner loops
#'              (default: c(20, 50))
#' @param tol vector of two numeric values specifying convergence tolerances for outer and inner
#'            loops (default: c(1e-4, 1e-7))
#'
#' @return a list containing the following components:
#'   \item{bet}{matrix of coefficient estimates, with columns corresponding to lambda values}
#'   \item{lam}{sequence of lambda values used in the model}
#'   \item{int}{intercept term (all are 0 if std=FALSE)}
#' @export
#'
cdmm <- function(y, x, lam, nlam=100, rlam=1/nlam, mu=1, std=TRUE, maxv=0.4*length(y), maxit=c(20, 50), tol=c(1e-4, 1e-7)) {
	if (std) {
		y <- scale(y, scale=FALSE)
		x <- scale(x, scale=apply(x, 2, sd)*sqrt(nrow(x)-1))
		fac <- 1/attr(x, "scaled:scale")
	} else
		fac <- rep(1, ncol(x))
	if (missing(lam)) {
		lam.max <- max(abs(crossprod(x, y)))
		lam <- lam.max*exp(seq(0, log(rlam), length=nlam))
	}
	res <- .Call("cdmm_c", y, x, fac, lam, mu, maxv, as.integer(maxit), tol)
	names(res) <- c("bet", "lam")
	res$bet <- res$bet*fac
	if (std) res$int <- attr(y, "scaled:center") - drop(crossprod(res$bet, attr(x, "scaled:center")))
	else res$int <- rep(0, length(res$lam))
	res
}

#' CDMM Model Selection By GIC
#' @description
#' \code{gic.cdmm} is used to choose tunning parameter \code{lam} for \code{cdmm} by the way of GIC.
#' @param y numeric vector of the response variable
#' @param x n x p logged composition data matrix of predictor variables (rows = samples, columns = variables)
#' @param lam sequence of tuning parameters for regularization (passed to \code{cdmm} function). The function will choose the optimized one in the parameter sequence. If missing,an automatic sequence is generated based on data (see \code{nlam} and \code{rlam}).
#' @param nlam number of lambda values in the automatic sequence (default: 100)
#' @param rlam ratio of the smallest to largest lambda in the automatic sequence (default: 1/nlam)
#' @param mu regularization parameter for the compositional constraint, passed to \code{cdmm} functionm (default: 1).
#' @param std logical indicator for data standardization: TRUE (default) standardizes predictors and centers the response; FALSE uses raw data
#' @param type character specifying the GIC type: "bic" (default, Bayesian Information Criterion),
#'             "aic" (Akaike Information Criterion), or "ft" (Frequentist T-test criterion)
#'
#' @return a list containing the following components:
#'   \item{bet}{optimal coefficient vector corresponding to the lambda with minimum GIC}
#'   \item{lam}{optimal tuning parameter (lambda) selected by minimum GIC}
#'   \item{int}{intercept term corresponding to the optimal lambda (from CDMM results)}
#' @export
gic.cdmm <- function(y, x, lam, nlam=100, rlam=1/nlam, mu=1, std=TRUE, type="bic") {
  if (missing(lam))
    res <- cdmm(y, x, nlam = nlam, rlam = rlam, mu=mu, std=std)
  else
    res <- cdmm(y, x, lam, mu=mu, std=std)

	res <- cdmm(y, x, lam)

	n <- length(y)
	fit <- log(colMeans((y - matrix(res$int, n, length(res$lam), byrow=TRUE) - x %*% res$bet)^2))
	a <- switch(type, bic=log(n), aic=2, ft=log(log(n))*log(max(ncol(x), n)))/n
	gics <- fit + a*(colSums(res$bet != 0) - 1)
	ilam <- which.min(gics)
	list(bet=res$bet[, ilam], lam=res$lam[ilam], int=res$int[ilam])
}

#' CDMM Model Selection By Cross-Validation
#' @description \code{cv.cdmm} is used to choose tuning parameter \code{lam} for \code{cdmm} by the way of k-fold cross validation.
#'
#' @param y numeric vector of the response variable
#' @param x n x p logged composition data matrix of predictor variables (rows = samples, columns = variables)
#' @param lam sequence of tuning parameters for regularization (passed to \code{cdmm} function). If missing,
#'            an automatic sequence is generated based on data (see \code{nlam} and \code{rlam}).
#' @param nlam number of lambda values in the automatic sequence (default: 100)
#' @param rlam ratio of the smallest to largest lambda in the automatic sequence (default: 1/nlam)
#' @param mu regularization parameter for the compositional constraint, passed to \code{cdmm} functionm (default: 1).
#' @param std logical indicator for data standardization: TRUE (default) standardizes predictors
#'            and centers the response; FALSE uses raw data
#' @param foldid optional integer vector specifying fold assignments for each sample. If missing,
#'               folds are randomly generated (see \code{nfold})
#' @param nfold integer number of folds for cross-validation (default: 5)
#' @param refit logical indicating whether to refit the model on non-validation folds with non-zero
#'              coefficients (default: FALSE, skips refitting)
#' @param type character specifying the optimal lambda selection rule: "min" (default, selects lambda
#'             with minimum CV error) or "1se" (selects the largest lambda within 1 SE of the minimum CV error)
#' @return a list containing the following components:
#'   \item{bet}{optimal coefficient vector corresponding to the selected lambda}
#'   \item{lam}{optimal tuning parameter (lambda) selected by cross-validation}
#'   \item{int}{intercept term corresponding to the optimal lambda (from CDMM results)}
#'   \item{foldid}{integer vector of fold assignments used for cross-validation (for reproducibility)}
#'
#' @export
cv.cdmm <- function(y, x, lam, nlam=100, rlam=1/nlam, mu=1, std=TRUE, foldid, nfold=5, refit=FALSE, type="min") {
  constr = TRUE
  if (missing(lam))
    res <- cdmm(y, x, nlam = nlam, rlam = rlam, mu=mu, std=std)
  else
    res <- cdmm(y, x, lam, mu=mu, std=std)

	if (missing(foldid)) foldid <- sample(rep(1:nfold, length=length(y)))
	pred <- matrix(, nfold, length(res$lam))
	for (i in 1:nfold) {
		yt <- y[foldid != i]; xt <- x[foldid != i, ]
		yv <- y[foldid == i]; xv <- x[foldid == i, ]
		if (constr) {
			fit <- cdmm(yt, xt, res$lam, mu=mu, std=std, maxv=Inf)
			if (refit) for (j in 1:length(res$lam)) {
				supp <- fit$bet[, j] != 0
				if (any(supp)) {
					ans <- cdmm(yt, as.matrix(xt[, supp]), 0, mu=mu, std=std, maxv=Inf)
					fit$bet[supp, j] <- ans$bet; fit$int[j] <- ans$int
				}
			}
		} else {
			fit <- cdmm(yt, xt, res$lam, mu=0, std=std, maxv=Inf, maxit=c(50, 1))
 			if (refit) for (j in 1:length(res$lam)) {
				supp <- fit$sol[, j] != 0
				if (any(supp)) {
					ans <- cdmm(yt, as.matrix(xt[, supp]), 0, mu=0, std=std, maxv=Inf, maxit=c(50, 1))
					fit$bet[supp, j] <- ans$bet; fit$int[j] <- ans$int
				}
			}
		}
		pred[i, ] <- colSums((yv - matrix(fit$int, length(yv), length(res$lam), byrow=TRUE) - xv %*% fit$bet)^2)
	}
	cvm <- colSums(pred)/length(y)
	cvse <- sqrt((colSums(pred^2) - length(y)*cvm^2)/(length(y) - 1))
	imin <- which.min(cvm)
	ilam <- switch(type, min=imin, "1se"=match(TRUE, cvm <= cvm[imin] + cvse[imin]))
	list(bet=res$bet[, ilam], lam=res$lam[ilam], int=res$int[ilam], foldid=foldid)
}

#' Stability Selection for CDMM
#'
#'
#' @param y numeric vector of the response variable
#' @param x n x p composition data matrix of predictor variables (rows = samples, columns = variables)
#' @param lam optional sequence of tuning parameters for regularization. If missing, an automatic
#'            sequence is generated (see \code{nlam} and \code{rlam})
#' @param nlam integer number of lambda values in the automatic sequence (default: 100, only used if \code{lam} is missing)
#' @param rlam ratio of the smallest to largest lambda in the automatic sequence (default: 1/nlam, only used if \code{lam} is missing)
#' @param nsample integer number of subsamples to generate for stability assessment (default: 100)
#' @param nsub integer size of each subsample (default: 50% of total samples, \code{floor(0.5*length(y))})
#' @param constr logical indicating whether to enforce compositional constraints in CDMM:
#'               TRUE (default, uses standard CDMM) or FALSE (disables constraints via \code{mu=0})
#' @param seed Integer random seed for reproducible subsample generation (default: 213)
#'
#' @return a list containing the following components:
#'   \item{path}{p x nlam matrix of selection frequencies: each entry (i,j) is the proportion of subsamples
#'               where predictor i was selected at lambda j}
#'   \item{prob}{numeric vector of maximum selection probabilities: each entry is the highest frequency
#'               of selection for a predictor across all lambda values}
#'
stab.cdmm <- function(y, x, lam, nlam=100, rlam=1/nlam, nsample=100, nsub=floor(0.5*length(y)), constr=TRUE, seed=213) {
	set.seed(seed)
	if (missing(lam)) {
		y <- scale(y, scale=FALSE)
		x <- scale(x, scale=apply(x, 2, sd)*sqrt(nrow(x)-1))
		lam.max <- max(abs(crossprod(x, y)))
		lam <- lam.max*exp(seq(0, log(rlam), length=nlam))
	}
	path <- matrix(0, ncol(x), nlam)
	for (i in 1:nsample) {
		isub <- sample(1:length(y), nsub)
		ysub <- y[isub]; xsub <- x[isub, ]
		if (constr)
			sol <- cdmm(ysub, xsub, lam, maxv=Inf)$sol
		else
			sol <- cdmm(ysub, xsub, lam, mu=0, maxv=Inf, maxit=c(50, 1))$sol
		path <- path + (sol != 0)
	}
	path <- path/nsample
	prob <- apply(path, 1, max)
	list(path=path, prob=prob)
}
