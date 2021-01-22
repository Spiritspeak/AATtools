
# Score computation algorithms ####

#' @title AAT score computation algorithms
#' @name Algorithms
#' @description
#' \itemize{
#' \item \code{aat_doublemeandiff} computes a mean-based double-difference score:
#'
#' \code{(mean(push_target) - mean(pull_target)) - (mean(push_control) - mean(pull_control))}
#' \item \code{aat_doublemediandiff} computes a median-based double-difference score:
#'
#' \code{(median(push_target) - median(pull_target)) - (median(push_control) - median(pull_control))}
#' \item \code{aat_dscore} computes D-scores for a 2-block design (see Greenwald, Nosek, and Banaji, 2003):
#'
#' \code{((mean(push_target) - mean(pull_target)) - (mean(push_control) - mean(pull_control))) / sd(participant_reaction_times)}
#' \item \code{aat_dscore_multiblock} computes D-scores for pairs of sequential blocks
#' and averages the resulting score (see Greenwald, Nosek, and Banaji, 2003).
#' Requires extra \code{blockvar} argument, indicating the name of the block variable.
#' \item \code{aat_regression} and \code{aat_standardregression} fit regression models to participants' reaction times and extract a term that serves as AAT score.
#' \code{aat_regression} extracts the raw coefficient, equivalent to a mean difference score.
#' \code{aat_standardregression} extracts the t-score of the coefficient, standardized on the basis of the variability of the participant's reaction times.
#' These algorithms can be used to regress nuisance variables out of the data before computing AAT scores.
#' When using these functions, additional arguments must be provided:
#' \itemize{
#' \item \code{formula} - a formula to fit to the data
#' \item \code{aatterm} - the term within the formula that indicates the approach bias; this is usually the interaction of the pull and target terms.
#' }
#' \item \code{aat_doublemeanquotient} and \code{aat_doublemedianquotient} compute a log-transformed ratio of approach to avoidance for both stimulus categories and subtract these ratios:
#'
#' \code{log(mean(pull_target) / mean(push_target)) - log(mean(pull_control) / mean(push_control))}
#' \item \code{aat_singlemeandiff} and \code{aat_singlemediandiff} subtract the mean or median approach reaction time from the mean or median avoidance reaction time.
#' These algorithms are only sensible if the supplied data contain a single stimulus category.
#' }
#' @param ds A long-format data.frame
#' @param subjvar Column name of the participant identifier variable
#' @param pullvar Column name of the movement variable (0: avoid; 1: approach)
#' @param targetvar Column name of the stimulus category variable (0: control stimulus; 1: target stimulus)
#' @param rtvar Column name of the reaction time variable
#' @param ... Other arguments passed on by functions (ignored)
#'
#' @return A data.frame containing participant number and computed AAT score.
NULL

#' @export
#' @rdname Algorithms
aat_doublemeandiff<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  group_by(ds,!!sym(subjvar)) %>%
    summarise(ab=(mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=TRUE) -
                    mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=TRUE)) -
                (mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=TRUE) -
                   mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=TRUE)))
}

#' @export
#' @rdname Algorithms
aat_doublemediandiff<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  group_by(ds,!!sym(subjvar)) %>%
    summarise(ab=(median(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=TRUE) -
                    median(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=TRUE)) -
                (median(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=TRUE) -
                   median(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=TRUE)))
}

#' @export
#' @rdname Algorithms
aat_dscore<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  group_by(ds,!!sym(subjvar)) %>%
    summarise(ab=((mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=TRUE) -
                     mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=TRUE)) -
                    (mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=TRUE) -
                       mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=TRUE))) /
                sd(!!sym(rtvar),na.rm=TRUE))
}

#' @param blockvar name of the variable indicating block number
#' @export
#' @rdname Algorithms
#note: this matches sequential blocks with one another.
aat_dscore_multiblock<-function(ds,subjvar,pullvar,targetvar,rtvar,blockvar,...){
  ds %>% mutate(.blockset = floor((!!sym(blockvar) - min(!!sym(blockvar)))/2) ) %>%
    group_by(!!sym(subjvar),.data$.blockset) %>%
    summarise(ab=((mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=TRUE) -
                     mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=TRUE)) -
                    (mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=TRUE) -
                       mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=TRUE))) /
                sd(!!sym(rtvar),na.rm=TRUE)) %>%
    group_by(!!sym(subjvar)) %>% summarise(ab=mean(ab,na.rm=TRUE))
}

#' @param formula A regression formula to fit to the data to compute an AAT score
#' @param aatterm A character naming the formula term representing the approach bias.
#' Usually this is the interaction of the movement-direction and stimulus-category terms.
#' @export
#' @rdname Algorithms
aat_regression<-function(ds,subjvar,formula,aatterm,...){
  output<-data.frame(pp=unique(ds[[subjvar]]),ab=NA,var=NA)
  for(i in seq_len(nrow(output))){
    mod<-coef(summary(lm(formula,data=ds[ds[[subjvar]]==output[i,"pp"],])))
    if(aatterm %in% rownames(mod)){
      output[i,"ab"]<- -mod[rownames(mod)==aatterm,1]
      output[i,"var"]<- mod[rownames(mod)==aatterm,2]
    }
  }
  colnames(output)[colnames(output)=="pp"]<-subjvar
  return(output)
}

#' @export
#' @rdname Algorithms
aat_standardregression<-function(ds,subjvar,formula,aatterm,...){
  output<-data.frame(pp=unique(ds[[subjvar]]),ab=NA,var=NA)
  for(i in seq_len(nrow(output))){
    mod<-coef(summary(lm(formula,data=ds[ds[[subjvar]]==output[i,"pp"],])))
    if(aatterm %in% rownames(mod)){
      output[i,"ab"]<- -mod[rownames(mod)==aatterm,1]
      output[i,"var"]<- mod[rownames(mod)==aatterm,2]
    }
  }
  colnames(output)[colnames(output)=="pp"]<-subjvar
  output$ab<-output$ab/output$var
  return(output)
}

#' @export
#' @rdname Algorithms
aat_doublemedianquotient<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  group_by(ds,!!sym(subjvar)) %>%
    summarise(ab=log(median(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=TRUE) /
                       median(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=TRUE)) -
                log(median(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=TRUE) /
                      median(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=TRUE)))
}

#' @export
#' @rdname Algorithms
aat_doublemeanquotient<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  group_by(ds,!!sym(subjvar)) %>%
    summarise(ab=log(mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=TRUE) /
                       mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=TRUE)) -
                log(mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=TRUE) /
                      mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=TRUE)))
}

#' @export
#' @rdname Algorithms
aat_singlemeandiff<-function(ds,subjvar,pullvar,rtvar,...){
  group_by(ds,!!sym(subjvar))%>%
    summarise(ab=mean(subset(!!sym(rtvar),!!sym(pullvar)==1)) -
                mean(subset(!!sym(rtvar),!!sym(pullvar)==0)))
}

#' @export
#' @rdname Algorithms
aat_singlemediandiff<-function(ds,subjvar,pullvar,rtvar,...){
  group_by(ds,!!sym(subjvar))%>%
    summarise(ab=median(subset(!!sym(rtvar),!!sym(pullvar)==1)) -
                median(subset(!!sym(rtvar),!!sym(pullvar)==0)))
}
