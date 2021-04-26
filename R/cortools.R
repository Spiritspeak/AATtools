#' @name correlation-tools
#' @title Correlation tools
#' @description Helper functions to compute important statistics from correlation coefficients.
#' @param r,r1,r2 a correlation value
#' @param z a Z-score
#' @param n,n1,n2 sample sizes
#' @param alpha the significance level to use
#' @examples
#' z <- r2z(.5)
#' r <- z2r(z)
#' t<-r2t(r,30)
#' r2p(r,30)
#' print(rconfint(r,30))
#' print(compcorr(.5,.7,20,20))
NULL

#' @export
#' @describeIn correlation-tools converts correlation coefficients to z-scores
r2z<-function(r){ .5 * (log(1+r) - log(1-r)) }
#' @export
#' @describeIn correlation-tools converts z-scores to correlation coefficients
z2r<-function(z){ (exp(2*z)-1)/(exp(2*z)+1) }
#' @export
#' @describeIn correlation-tools Converts correlation coefficients to t-scores
r2t<-function(r,n){ (r*sqrt(n-2))/(1-r^2) }
#' @export
#' @describeIn correlation-tools Computes the p-value for a given correlation
r2p<-function(r,n){ 2*pt(abs(r2t(r,n)),n-2,lower.tail=FALSE) }
#' @export
#' @describeIn correlation-tools Computes confidence intervals for a given correlation coefficient
rconfint<-function(r,n,alpha=.05){
  z<-r2z(r)
  zint<-qnorm(1-alpha/2) * sqrt(1/(n-3))
  confints<-c(z2r(z-zint),z2r(z+zint))
}
#' @export
#' @describeIn correlation-tools computes the significance of the difference between two correlation coefficients
compcorr<-function(r1,r2,n1,n2){
  zval<-abs(r2z(r1)-r2z(r2)) / sqrt((1/(n1-3)) + (1/(n2-3)))
  pval<-min(1,pnorm(abs(zval),lower.tail=F)*2)
  return(structure(list(zscore=zval,pvalue=pval),class="compcorr"))
}

print.compcorr<-function(x,...){
  cat("Two-tailed Z-test for the difference between two correlation coefficients.",
      "\nZ =",x$zscore,"\np =",x$pvalue,"\n")
}


#' Compute a minimally average correlation
#'
#' This function computes a minimally biased average of correlation values.
#' This is needed because simple averaging of correlations is negatively biased,
#' and the often used z-transformation method of averaging correlations is positively biased.
#' The algorithm was developed by Olkin & Pratt (1958) and implemented by Jan Seifert.
#'
#' @param r a vector containing correlation values
#' @param n a single value or vector containing sample sizes
#' @param na.rm Logical. Should missing values be removed?
#'
#' @return An average correlation.
#' @export
#'
#' @references
#' Olkin, I., & Pratt, J. (1958). Unbiased estimation of certain correlation coefficients.
#' The Annals of Mathematical Statistics, 29. https://doi.org/10.1214/aoms/1177706717
#'
#' https://github.com/SigurdJanson/AveragingCorrelations/blob/master/CorrAggBias.R
#'
#' https://medium.com/@jan.seifert/averaging-correlations-part-ii-9143b546860b
#'
#' @examples
#' cormean(c(0,.3,.5),c(30,30,60))
cormean <- function(r, n, na.rm=F) {
  if(any(n <= 4)) stop("Sample size must be at least 5")
  df <- n-1
  k <- (9*sqrt(2)-7)/2
  G <- r * (1 + ((1-r^2) / (2 * (df - k))))

  NaCount <- ifelse(na.rm == TRUE, sum(is.na(G)), 0)
  if(length(n) < length(G)) n <- rep(n, length.out = length(G))

  Num <- sum(df * G, na.rm)
  Den <- sum(n, na.rm) - (length(n)-NaCount)
  if(any(Den == 0)) Den[Den == 0] <- NA
  return(Num / Den)
}
