# frame<-data.frame(jitter=rep(c(1:10*10,0:9),2),intercat=NA,split=NA)
# dslist<-list()
# for(i in seq_len(nrow(frame))){
#   message(i)
#   ds<-aat_simulate(biasfx_jitter=frame$jitter[i],nstims=16)
#   ds$stim<-paste0(ds$stim,"-",ds$is_target)
#   ds<-ds%>%filter(!(stim %in% paste0(1:12,"-",1)))
#   dslist[[i]]<-ds
#
#   rel<-aat_splithalf(ds=ds,subjvar="subj",pullvar="is_pull",targetvar="is_target",rtvar="rt",
#                      algorithm="aat_doublemeandiff",iters=1000,plot=F)
#   frame$split[i]<-rel$uncorrected$r%>%SpearmanBrown(fix.negative = "none")
# }


#' Compute the reliability of a double difference score
#'
#' This function uses Lord's (1963) algorithm to quickly compute the reliability
#' of a double mean difference score. Missing values are dealt with using multiple imputation.
#'
#' When the study in question uses double median difference scores or D-scores,
#' \code{aat_splithalf()} is recommended instead; and if only single mean difference scores
#' are used, \code{aat_alpha} suffices.
#'
#' @param ds the \code{data.frame} to use
#' @param subjvar Name of the subject-identifying variable
#' @param stimvar Name of the stimulus-identifying variable
#' @param pullvar Name of the movement-direction identifying variable
#' @param targetvar Name of the stimulus-category identifying variable
#' @param rtvar Name of the reaction-time identifying variable
#' @param iters If there are missing values (which is almost inevitable) then
#' multiple imputation will be used to complete the covariance matrix - this option sets
#' the number of multiple imputations to be used.
#'
#' @return Returns an \code{aat_intercategorical_alpha} object containing the alpha value
#' as well as the dataset and covariance matrix with replaced missing values.
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
#' #If the intercategorical alpha and
#' # bootstrapped splithalf measure something similar,
#' # then the outcomes should be close to each other.
#' aat_intercategorical_alpha(ds=ds,subjvar="subj",stimvar="stim",pullvar="is_pull",
#'                            targetvar="is_target",rtvar="rt")
#' aat_splithalf(ds=ds,subjvar="subj",pullvar="is_pull",targetvar="is_target",rtvar="rt",
#'               algorithm="aat_doublemeandiff",iters=100,plot=FALSE)
aat_intercategorical_alpha<-function(ds,subjvar,stimvar,pullvar,targetvar,rtvar,iters=5){
  pps<-unique(ds[[subjvar]])
  stims<-unique(ds[[stimvar]])
  stimcats<-distinct(ds[c(stimvar,targetvar)])
  biases<-list()
  for(u in seq_along(pps)){
    biases[[u]]<-
      do.call(aat_singlemeandiff,list(ds=ds[ds[[subjvar]]==pps[u],],
                             subjvar=stimvar,pullvar=pullvar,rtvar=rtvar)) %>%
      setNames(c(stimvar,paste0("subject-",pps[u])))
  }
  biasset<-Reduce(function(x,y){merge(x,y,by=stimvar,all=T)},x=biases[-1],init=biases[1])
  biasmat<-t(as.matrix(biasset[,-1]))
  colnames(biasmat)<-stimcats[[targetvar]][match(biasset[[stimvar]],stimcats$stim)]

  unmissing<-covEM(biasmat,iters)
  covmat<-unmissing$sigma
  rownames(covmat)<-colnames(covmat)
  dataset<-unmissing$data

  dia<-diag(covmat)
  firstcat<-which(names(dia) == unique(stimcats[[targetvar]])[1])
  secondcat<-which(names(dia) == unique(stimcats[[targetvar]])[2])
  n1<-length(firstcat)
  n2<-length(secondcat)

  # r11<-calpha(covmat[firstcat,firstcat])
  # r22<-calpha(covmat[secondcat,secondcat])
  # r12<-cor(x=rowSums(dataset[,firstcat]),
  #          y=rowSums(dataset[,secondcat]))
  # s1<-sqrt(sum(covmat[firstcat,firstcat]))/n1
  # s2<-sqrt(sum(covmat[secondcat,secondcat]))/n2
  # alpha<-(s1^2*r11+s2^2*r22-2*s1*s2*r12)/
  #   (s1^2+s2^2-2*s1*s2*r12)

  co1<-2*sum(covmat[firstcat,firstcat][upper.tri(covmat[firstcat,firstcat])])/n1/n1
  co2<-2*sum(covmat[secondcat,secondcat][upper.tri(covmat[secondcat,secondcat])])/n2/n2
  co12<-sum(covmat[firstcat,secondcat])/n1/n2
  v1<-sum(covmat[firstcat,firstcat])/n1/n1
  v2<-sum(covmat[secondcat,secondcat])/n2/n2

  alpha<-(co1*n1/(n1-1) + co2*n2/(n2-1) - 2*co12)/ (v1+v2-2*co12)

  out<-structure(list(alpha=alpha,data=dataset,covmat=covmat),
                 class="aat_intercategorical_alpha")
  return(out)
}

# for(i in seq_len(nrow(frame))){
#   frame$intercat[i]<-
#     aat_intercategorical_alpha(ds=dslist[[i]],subjvar="subj",stimvar="stim",
#                                pullvar="is_pull",targetvar="is_target",rtvar="rt",
#                                algorithm="aat_singlemeandiff")$alpha
# }
# frame %>% filter(split>0) %$% plot(split,intercat,main=mean((split-intercat)))
# lines(-2:2,-2:2)
# lm(split~intercat,data=frame)

#' @describeIn aat_intercategorical_alpha
print.aat_intercategorical_alpha<-function(x,...){
  cat("alpha = ",mf(x$alpha))
}



