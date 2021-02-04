
# utils ####
#' @name reliability-coefficients
#' @title Correct a correlation coefficient for being based on only a subset of the data
NULL

#' @describeIn reliability-coefficients Perform a Spearman-Brown correction on the provided correlation score.
#'
#' @param corr To-be-corrected correlation coefficient
#' @param ntests An integer indicating how many times larger the full test is, for which the corrected correlation coefficient is being computed.
#' When \code{ntests=2}, the formula will compute what the correlation coefficient would be if the test were twice as long.
#' @param fix.negative Determines how to deal with a negative value. "nullify" sets it to zero,
#' "bilateral" applies the correction as if it were a positive number, and then sets it to negative.
#' "none" gives the raw value. It should be noted that negative values are not supposed to occur,
#' and there is no commonly accepted way to deal with them when they do occur.
#' @return Spearman-Brown-corrected correlation coefficient.
#' @export
#'
#' @examples
#'
#' SpearmanBrown(.5)
SpearmanBrown<-function(corr,ntests=2,fix.negative=c("nullify","bilateral","none")){
  fix.negative<-match.arg(fix.negative)
  if(fix.negative=="bilateral"){
    s<-sign(corr)
    corr<-abs(corr)
    sb<-ntests*corr / (1+(ntests-1)*corr)
    return(s*sb)
  }else{
    sb<-ntests*corr / (1+(ntests-1)*corr)
    if(fix.negative=="nullify"){
      return(ifelse(sb<0,0,sb))
    }else{
      return(sb)
    }
  }
}

#' @describeIn reliability-coefficients Compute the true reliability using the Flanagan-Rulon formula,
#' which takes into account inequal variances between split halves
#' @param x1 scores from half 1
#' @param x2 scores from half 2
#' @param x scores from the full sample (more accurate if provided)
#' @export
#'
#' @examples
#' FlanaganRulon(a<-rnorm(50),rnorm(50)+a*.5,fix.negative="bilateral")
FlanaganRulon<-function(x1,x2,x=NULL,fix.negative=c("nullify","bilateral","none")){
  fix.negative<-match.arg(fix.negative)
  d<-var(x1-x2)

  if(missing(x)){
    k<-var(x1+x2)
  }else{
    k<-var(x)
  }

  if(fix.negative=="none"){
    return(1-d/k)
  }else if(fix.negative=="bilateral"){
    fr<-(1-d/k)
    #fr<-ifelse(fr>0,fr,fr / (1-fr))
    fr<-fr/max(1, 1-fr)
    return(fr)
  }else if(fix.negative=="nullify"){
    fr<-1-d/k
    return(ifelse(fr>0,fr,0))
  }
}

#' @describeIn reliability-coefficients Compute split-half reliability using the Raju formula,
#' which takes into account unequal split-halves and variances.
#'
#' @param prop Proportion of the first half to the complete sample
#'
#' @export
#'
#' @examples
#' a<-rnorm(50)
#' b<-rnorm(50)+a*.5
#' RajuCoefficient(a,b,prop=.4,fix.negative="bilateral")
RajuCoefficient<-function(x1,x2,prop,fix.negative=c("nullify","bilateral","none")){
  fix.negative<-match.arg(fix.negative)
  covar<-cov(x1,x2)
  if(fix.negative=="bilateral"){
    sumvar<-var(x1)+var(x2)+2*abs(covar)
  }else{
    sumvar<-var(x1)+var(x2)+2*covar
  }

  raju<-covar / (prop * (1-prop) * sumvar)
  return(ifelse(fix.negative=="nullify" & raju<0,0,raju))
}
