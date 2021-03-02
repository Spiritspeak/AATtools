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
  pval<-pnorm(abs(zval)/2,lower.tail=F)
  return(structure(list(zscore=zval,pvalue=pval),class="compcorr"))
}

print.compcorr<-function(x,...){
  cat("Two-tailed Z-test for the difference between two correlation coefficients.",
      "\nZ =",x$zscore,"\np =",x$pvalue,"\n")
}

