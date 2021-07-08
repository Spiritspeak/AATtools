

#' Compute stimulus-specific bias scores
#' Computes mean single-difference scores (push - pull) for each stimulus.
#'
#' @param ds the \code{data.frame} to use
#' @param subjvar Name of the subject-identifying variable
#' @param stimvar Name of the stimulus-identifying variable
#' @param pullvar Name of the movement-direction identifying variable
#' @param targetvar Optional. Name of the stimulus-category identifying variable
#' @param rtvar Name of the reaction-time identifying variable
#' @param iters If there are missing values (which is almost inevitable) then
#' multiple imputation will be used to complete the covariance matrix - this argument sets
#' the number of multiple imputations to be used.
#'
#' @return Exports a \code{list} containing
#' a \code{data.frame} with stimulus-specific bias scores, indicated in the column names,
#' a covariance matrix of that same data, and
#' a \code{data.frame} indicating to which stimulus category each stimulus belongs.
#' @export
#'
#' @examples
#' ds<-aat_simulate(biasfx_jitter=40,nstims=16)
#' ds$stim<-paste0(ds$stim,"-",ds$is_target)
#' aat_stimulusscores(ds,"subj","stim","is_pull","is_target","rt")
aat_stimulusscores<-function(ds,subjvar,stimvar,pullvar,targetvar=NULL,rtvar,iters=5){
  ds<-aat_preparedata(ds,subjvar=subjvar,pullvar=pullvar,stimvar=stimvar,targetvar=targetvar,rtvar=rtvar)

  pps<-unique(ds[[subjvar]])
  stims<-unique(ds[[stimvar]])

  if(!is.null(targetvar)){
    stimcats<-distinct(ds[c(stimvar,targetvar)]) %>% setNames(c("stim","cat"))
  }else{
    stimcats<-data.frame(stim=stims,cat=0,stringsAsFactors=F)
  }

  biases<-list()
  for(u in seq_along(pps)){
    biases[[u]]<-
      do.call(aat_singlemeandiff,list(ds=ds[ds[[subjvar]]==pps[u],],
                                      subjvar=stimvar,pullvar=pullvar,rtvar=rtvar)) %>%
      setNames(c(stimvar,paste0("",pps[u]))) #subject-
  }
  biasset<-Reduce(function(x,y){merge(x,y,by=stimvar,all=T)},x=biases[-1],init=biases[1])
  biasmat<-t(as.matrix(biasset[,-1]))
  colnames(biasmat)<-biasset[[1]]

  unmissing<-covEM(biasmat,iters)
  covmat<-unmissing$sigma
  rownames(covmat)<-colnames(covmat)
  dataset<-unmissing$data

  out<-list(data=dataset,covmat=covmat,stimcats=stimcats)
  return(out)
}


#' Compute a dataset's reliability from its covariance matrix
#'
#' This function computes mean single-difference scores (push minus pull) for individual stimuli,
#' and computes the reliability from that information.
#' Missing values are dealt with using multiple imputation.
#'
#' When only one stimulus category is indicated, one of the commonly known reliability algorithms
#' provided with the \code{algorithm} argument is used.
#' When two stimulus categories are indicated, this function uses Lord's (1963) algorithm to
#' compute the reliability of a double mean difference score, using the algorithms in \code{algorithm}
#' to estimate the reliability of indiviau lstimulus categories.
#'
#' When one wants to compute the reliability of a double median difference score or D-score,
#' \code{aat_splithalf()} is recommended instead.
#'
#' @param ds the \code{data.frame} to use
#' @param subjvar Name of the subject-identifying variable
#' @param stimvar Name of the stimulus-identifying variable
#' @param pullvar Name of the movement-direction identifying variable
#' @param targetvar Optional. Name of the stimulus-category identifying variable
#' @param rtvar Name of the reaction-time identifying variable
#' @param algorithm The reliability formula to use. Defaults to Cronbach's alpha, but Guttman's Lambda-2 is recommended instead.
#' @param iters If there are missing values (which is almost inevitable) then
#' multiple imputation will be used to complete the covariance matrix - this option sets
#' the number of multiple imputations to be used.
#'
#' @return Returns an \code{aat_covreliability} object containing the reliability value
#' as well as the dataset and covariance matrix with replaced missing values. When
#' the argument \code{targetvar} is provided, the output also contains the reliability of the
#' individual stimulus categories and their intercorrelation.
#'
#' @export
#'
#' @references
#' Lord, F.Y. (1963), "Elementary Models for Measuring Change",
#' in Problems in Measuring Change, C.W. Harris, ed.. Madison. Wisconsin:
#' University of Wisconsin.
#'
#' @examples
#' #We generate a dataset with 16 stimuli in each category
#' ds<-aat_simulate(biasfx_jitter=40,nstims=16)
#' ds$stim<-paste0(ds$stim,"-",ds$is_target)
#'
#' # If Lord's formula and
#' # bootstrapped splithalf measure something similar,
#' # then the outcomes should be close to each other.
#' aat_covreliability(ds=ds,subjvar="subj",stimvar="stim",pullvar="is_pull",
#'                            targetvar="is_target",rtvar="rt")
#' aat_splithalf(ds=ds,subjvar="subj",pullvar="is_pull",targetvar="is_target",rtvar="rt",
#'               algorithm="aat_doublemeandiff",iters=100,plot=FALSE)
#'
#' #Testing reliability for single-difference scores
#' ds<-ds[ds$is_target==1,]
#' aat_covreliability(ds=ds,subjvar="subj",stimvar="stim",pullvar="is_pull",rtvar="rt")
aat_covreliability<-function(ds,subjvar,stimvar,pullvar,targetvar=NULL,rtvar,algorithm=c("calpha","lambda2","lambda4"),iters=5){
  algorithm<-match.arg(algorithm)
  sc<-aat_stimulusscores(ds,subjvar=subjvar,stimvar=stimvar,pullvar=pullvar,targetvar=targetvar,
                         rtvar=rtvar,iters=iters)

  if(!is.null(targetvar)){
    dia<-diag(sc$covmat)
    firstcat <-which(names(dia) %in% sc$stimcats$stim[sc$stimcats$cat==0])
    secondcat<-which(names(dia) %in% sc$stimcats$stim[sc$stimcats$cat==1])
    n1<-length(firstcat )
    n2<-length(secondcat)

    r11<-do.call(algorithm,list(covmat=sc$covmat[firstcat, firstcat ]))
    r22<-do.call(algorithm,list(covmat=sc$covmat[secondcat,secondcat]))
    # r12<-cor(x=rowSums(sc$dataset[,firstcat]),
    #          y=rowSums(sc$dataset[,secondcat]))
    r12<-sum(sc$covmat[firstcat,secondcat])/sqrt(sum(sc$covmat[firstcat,firstcat])*sum(sc$covmat[secondcat,secondcat]))
    s1<-sqrt(sum(sc$covmat[firstcat,firstcat]))/n1
    s2<-sqrt(sum(sc$covmat[secondcat,secondcat]))/n2
    rel<-(s1^2*r11+s2^2*r22-2*s1*s2*r12)/
      (s1^2+s2^2-2*s1*s2*r12)
  }else{
    rel<-do.call(algorithm,list(covmat=sc$covmat))
  }

  out<-structure(list(rel=rel,data=sc$data,covmat=sc$covmat,algorithm=algorithm),
                 class="aat_covreliability")
  if(!is.null(targetvar)){
    out$components<-list(r11=r11,r22=r22,r12=r12,n1=n1,n2=n2,s1=s1,s2=s2)
  }
  return(out)
}

print.aat_covreliability<-function(x,...){
  cat(sep="","r = ",mf(x$rel),
      "\nBased on ",ncol(x$data)," valid stimuli, ",
      nrow(x$data)," valid participants, and the ",
      x$algorithm," algorithm.\n")
  if(any("components"==names(x))){
    cat(sep="",
        "Reliability of stimulus category 1: r = ",mf(x$components$r11),", n = ",x$components$n1,", sd = ",mf(x$components$s1),"\n",
        "Reliability of stimulus category 2: r = ",mf(x$components$r22),", n = ",x$components$n2,", sd = ",mf(x$components$s2),"\n",
        "Category intercorrelation: r = ",mf(x$components$r12),"\n")
  }
}
