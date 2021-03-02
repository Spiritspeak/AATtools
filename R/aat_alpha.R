
#' Compute reliability using Cronbach's Alpha
#'
#' @description \code{aat_alpha} computes approach-avoidance scores for each stimulus
#' within each participant. It then computes Cronbach's alpha by treating each stimulus
#' as an item on a questionnaire and measuring how much these scores correlate within
#' each participant.
#'
#' \code{aat_alpha_jackknife} additionally leaves out one stimulus or participant
#' at a time and computes Cronbach's alpha for each such exclusion. This gives a
#' glimpse of which stimuli may have especially confused participants or otherwise
#' induced responses that are unlike those to other stimuli.
#'
#'
#' @param ds a data.frame
#' @param subjvar character naming the column with participant IDs
#' @param stimvar character naming the column with stimulus IDs
#' @param pullvar character naming the column for movement direction
#' @param rtvar character indicating the column with reaction times
#' @param algorithm character name of the algorithm used to compute
#' per-stimulus approach bias scores
#' @param delete.missing character denoting the strategy to be taken when missing
#' values are encountered, which may occur when one participant saw a stimulus which
#' another participant did not see, or when all approach trials for a certain stimulus
#' were excluded for one participant. Defaults to estimating the covariance matrix.
#' @param iters Numeric. When there are missing values, the covariance matrix is estimated
#' instead of computed; how many iterations should be used?
#'
#' @return \code{aat_alpha} returns an \code{aat_alpha} S3 object,
#' \code{aat_alpha_jackknife} returns an \code{aat_alpha_jackknife} S3 object,
#' both with print and plot methods.
#' @export
#' @name aat_alpha
#' @author Sercan Kahveci
#' @details When there are missing values, the covariance matrix is estimated using
#' 100 iterations of an expectation-maximization algorithm developed by Beale and Little (1975)
#' and implemented in R by StackExchange user Eric at URL
#' https://stats.stackexchange.com/questions/182718/ml-covariance-estimation-from-expectation-maximization-with-missing-data
#' @references Cousijn, J., Goudriaan, A. E., & Wiers, R. W. (2011).
#' Reaching out towards cannabis: Approach‐bias in heavy cannabis users predicts
#' changes in cannabis use. Addiction, 106(9), 1667-1674.
#'
#' Beale, E. M. L., & Little, R. J. A.. (1975). Missing Values in Multivariate Analysis.
#' Journal of the Royal Statistical Society. Series B (methodological), 37(1), 129–145.
#'
#' @examples
#' data(erotica)
#' #We artificially reduce the number of stimuli here because the original
#' #erotica dataset is not suitable for computing Cronbach's alpha.
#' erotica$stimulus<- substr(as.character(erotica$stimulus),5,5)
#'
#' myalpha<-aat_alpha(erotica,"subject","stimulus","is_pull","RT")
#' print(myalpha)
#' plot(myalpha)
#'
#' myalpha2<-aat_alpha_jackknife(erotica,"subject","stimulus","is_pull","RT")
#' print(myalpha2)
#' plot(myalpha2)
aat_alpha<-function(ds,subjvar,stimvar,pullvar,rtvar,
                    algorithm=c("aat_singlemeandiff","aat_singlemediandiff"),
                    delete.missing=c("estimate","subjects","stimuli","both","none"),
                    iters=100){
  algorithm<-match.arg(algorithm)
  delete.missing<-match.arg(delete.missing)

  pps<-unique(ds[[subjvar]])
  stims<-unique(ds[[stimvar]])
  biases<-list()
  for(u in seq_along(pps)){
    biases[[u]]<-
    do.call(algorithm,list(ds=ds[ds[[subjvar]]==pps[u],],
                                    subjvar=stimvar,pullvar=pullvar,rtvar=rtvar)) %>%
      setNames(c(stimvar,paste0("subject-",pps[u])))
  }
  biasset<-Reduce(function(x,y){merge(x,y,by=stimvar,all=T)},x=biases[-1],init=biases[1])
  biasmat<-t(as.matrix(biasset[,-1]))
  colnames(biasmat)<-biasset[[stimvar]]

  deletions<-list(subjects=numeric(0),stimuli=numeric(0))
  if(delete.missing=="estimate"){
    covmat<-covEM(biasmat,iters)
  }else{
    if(delete.missing %in% c("subjects","both")){
      deletions$subjects<-biasmat %>% apply(1,function(x){any(is.na(x))}) %>% which()
    }
    if(delete.missing %in% c("stimuli","both")){
      deletions$stimuli<-biasmat %>% apply(2,function(x){any(is.na(x))}) %>% which()
    }
    if(length(deletions$subjects)>0){
      biasmat<-biasmat[-deletions$subjects,]
      warning("Deleted ",length(deletions$subjects)," subjects due to missing values.")
    }
    if(length(deletions$stimuli)>0){
      biasmat<-biasmat[,-deletions$stimuli]
      warning("Deleted ",length(deletions$stimuli)," stimuli due to missing values.")
    }
    covmat<-cov(biasmat)
  }

  alpha<-calpha(covmat)

  out<-structure(list(alpha=alpha,biasmat=biasmat,covmat=covmat,
                      deletions=list(subjects=pps[deletions$subjects],
                                     stimuli=stims[deletions$stimuli])),
                 class="aat_alpha")
  return(out)
}

#' @rdname aat_alpha
#' @export
aat_alpha_jackknife<-function(ds,subjvar,stimvar,pullvar,rtvar,
                              algorithm=c("aat_singlemeandiff","aat_singlemediandiff"),
                              delete.missing=c("estimate","subjects","stimuli","both","none")){
  delete.missing<-match.arg(delete.missing)
  simple_alpha<-aat_alpha(ds,subjvar,stimvar,pullvar,rtvar,algorithm,delete.missing)

  #stimulus exclusions
  stim_outvec<-setNames(numeric(ncol(simple_alpha$biasmat)),colnames(simple_alpha$biasmat))
  for(i in seq_along(stim_outvec)){
    stim_outvec[i]<-calpha(cov(simple_alpha$biasmat[,-i]))
  }

  #participant exclusions
  pp_outvec<-setNames(numeric(nrow(simple_alpha$biasmat)),rownames(simple_alpha$biasmat))
  for(i in seq_along(pp_outvec)){
    pp_outvec[i]<-calpha(cov(simple_alpha$biasmat[-i,]))
  }

  out<-structure(class="aat_alpha_jackknife",.Data=c(simple_alpha,
                                                     list(stimulus_holdout=stim_outvec,participant_holdout=pp_outvec)))
}

#' @rdname aat_alpha
#' @param x an \code{aat_alpha} or \code{aat_alpha_jackknife} object
#' @param ... Ignored.
# @describeIn aat_alpha A print method for \code{aat_alpha}
print.aat_alpha<-function(x,...){
  cat(sep="","Alpha = ",format(x$alpha,digits=2,scientific=F),
      "\nBased on ",ncol(x$biasmat)," valid stimuli and ",
      nrow(x$biasmat)," valid participants.")
}

#' @rdname aat_alpha
# @describeIn aat_alpha A plot method for \code{aat_alpha}
plot.aat_alpha<-function(x,...){
  axis=1
  rmean<-apply(x$biasmat,axis,mean)
  rvar<-apply(x$biasmat,axis,var)

  # circle drawing code. Uncomment once it's figured out
  # how to meaningfully scale the circle width and height

  # if(axis==1){  covmat<-cov(t(x$biasmat))
  # }else{  covmat<-cov(x$biasmat)  }
  # invar<-sum(diag(covmat))/sum(covmat)
  # scalar<-sd(rmean) + mean(sqrt(rvar))
  # correction<-nrow(covmat)/(nrow(covmat)-1)
  # wv<-(invar)*scalar /correction
  # bv<-(1-invar)*scalar*correction

  # bv<-sd(rmean)
  # wv<-mean(sqrt(rvar))

  # lineset<-data.frame(x=mean(x$biasmat) +
  #                       cos(0:100 / 100 * 2*pi)*bv * 1/2*sqrt(2) -
  #                       sin(0:100 / 100 * 2*pi)*wv * 1/2*sqrt(2),
  #                     y=mean(x$biasmat) +
  #                       cos(0:100 / 100 * 2*pi)*bv * 1/2*sqrt(2) +
  #                       sin(0:100 / 100 * 2*pi)*wv * 1/2*sqrt(2))
  dispval<-0#(bv+wv)/100
  plotset<-data.frame(xstart=c(rmean+dispval,rmean-dispval),
                      ystart=c(rmean-dispval,rmean+dispval),
                      xend=c(rmean+sqrt(rvar) *1/2*sqrt(2)/2,
                             rmean-sqrt(rvar) *1/2*sqrt(2)/2),
                      yend=c(rmean-sqrt(rvar) *1/2*sqrt(2)/2,
                             rmean+sqrt(rvar) *1/2*sqrt(2)/2))

  lims<-c(min(c(plotset$xstart,plotset$xend)),
          max(c(plotset$xstart,plotset$xend)))
  plot(1,type="n",xlim=lims,ylim=lims,
       main=paste0("Reliability\n","alpha = ",round(x$alpha,digits=2)),
       xlab="Participants' scores",ylab="Participants' scores")
  #circle
  #lines(lineset$x,lineset$y)
  #mean bias scores
  points(rmean,rmean)
  #variance lines
  segments(plotset$xstart,plotset$ystart,plotset$xend,plotset$yend)
}

#' @rdname aat_alpha
# @describeIn aat_alpha A print method for \code{aat_alpha_jackknife}
print.aat_alpha_jackknife<-function(x,...){
  cat(sep="","Alpha = ",format(x$alpha,digits=2,scientific=F),
      "\nBased on ",ncol(x$biasmat)," valid stimuli and ",
      nrow(x$biasmat)," valid participants.",
      "\nLargest alpha achieveable through single stimulus omission = ",
      format(max(x$stimulus_holdout),digits=2,scientific=F),
      " (",names(x$stimulus_holdout)[which.max(x$stimulus_holdout)],")",
      "\nLargest alpha achieveable through single participant omission = ",
      format(max(x$participant_holdout),digits=2,scientific=F),
      " (",names(x$participant_holdout)[which.max(x$participant_holdout)],")")
}

#' @rdname aat_alpha
#' @param exclusions Character. Should the function display Cronbach's alpha for
#' individually excluded participants, stimuli, or both?
#  @describeIn aat_alpha A plot method for \code{aat_alpha_jackknife}
plot.aat_alpha_jackknife<-function(x,exclusions=c("both","stimulus","participant"),...){
  exclusions<-match.arg(exclusions)

  if(exclusions=="both"){
    par(mfrow=c(1,2))
  }

  if(exclusions %in% c("both","stimulus")){
    plot(sort(x$stimulus_holdout),ylab="Alpha",xlab="",cex.main=.9,
         main="Alpha when single\nstimuli are excluded")
    abline(h=x$alpha)
  }

  if(exclusions %in% c("both","participant")){
    plot(sort(x$participant_holdout),ylab="Alpha",xlab="",cex.main=.9,
         main="Alpha when individual\nparticipants are excluded")
    abline(h=x$alpha)
  }
  par(mfrow=c(1,1))
}

#' @name covmat-tools
#' @title Covariance Matrix Tools
#' These functions allow for the computation of covariance matrices with missing data
#' and for the computation of Cronbach's alpha from them.
NULL

#' @describeIn covmat-tools Compute Cronbach's alpha from a covariance matrix
#' @param covmat a covariance matrix
#' @export
calpha<-function(covmat){
  (nrow(covmat)/(nrow(covmat)-1))*(1 - sum(diag(covmat))/sum(covmat))
}

#' @describeIn covmat-tools This function computes a covariance matrix from data with some values missing at random.
#' The code was written by Eric from StackExchange. https://stats.stackexchange.com/questions/182718/ml-covariance-estimation-from-expectation-maximization-with-missing-data
#' @param dat_missing a matrix with missing values
#' @param iters the number of iterations to perform to estimate missing values
#' @references Beale, E. M. L., & Little, R. J. A.. (1975). Missing Values in Multivariate Analysis. Journal of the Royal Statistical Society. Series B (methodological), 37(1), 129–145.
#' @export
#' @examples
#' # make data with missing values
#' missing_mtcars <- mtcars
#' for(i in 1:20){
#'   missing_mtcars[sample(1:nrow(mtcars),1),sample(1:ncol(mtcars),1)]<-NA
#' }
#' covmat<-covEM(as.matrix(missing_mtcars))
#' calpha(covmat)
covEM<-function(dat_missing,iters=1000){
  if(!anyNA(dat_missing)){
    return(cov(dat_missing))
  }

  n <- nrow(dat_missing)
  nvar <- ncol(dat_missing)
  is_na <- apply(dat_missing,2,is.na) # index if NAs

  dat_impute <- dat_missing # data matrix for imputation
  # set initial estimates to means from available data
  for(i in 1:ncol(dat_impute)){
    dat_impute[is_na[,i],i] <- colMeans(dat_missing,na.rm = TRUE)[i]
  }

  # starting values for EM
  means <- colMeans(dat_impute)
  # NOTE: multiplying by (nrow-1)/(nrow) to get ML estimate
  sigma <- cov(dat_impute)*(nrow(dat_impute)-1)/nrow(dat_impute)

  # carry out EM over 100 iterations
  for(j in 1:iters){
    bias <- matrix(0,nvar,nvar)
    for(i in 1:n){
      row_dat <- dat_missing[i,]
      miss <- which(is.na(row_dat))
      if(length(miss)>0){
        bias[miss,miss] <- bias[miss,miss] + sigma[miss,miss] -
          sigma[miss,-miss] %*% solve(sigma[-miss,-miss]) %*% sigma[-miss,miss]
        dat_impute[i,miss] <- means[miss] +
          (sigma[miss,-miss] %*% solve(sigma[-miss,-miss])) %*%
          (row_dat[-miss]-means[-miss])
      }
    }
    # get updated means and covariance matrix
    means <- colMeans(dat_impute)
    biased_sigma <- cov(dat_impute)*(n-1)/n
    # correct for bias in covariance matrix
    sigma <- biased_sigma + bias/n
  }
  return(sigma)
}
