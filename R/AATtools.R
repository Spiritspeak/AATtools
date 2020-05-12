#AATtools

.onLoad<-function(libname, pkgname){
  registerS3method("print",class="aat_splithalf",method=print.aat_splithalf)
  registerS3method("plot",class="aat_splithalf",method=plot.aat_splithalf)
  registerS3method("print",class="aat_bootstrap",method=print.aat_bootstrap)
  registerS3method("plot",class="aat_bootstrap",method=plot.aat_bootstrap)
  registerS3method("print",class="qreliability",method=print.qreliability)
  registerS3method("plot",class="qreliability",method=plot.qreliability)
  packageStartupMessage("Thank you for loading AATtools v0.0.1")
}

# splithalf engine ####

#multicore splithalf
#' @title Compute the bootstrapped split-half reliability for approach-avoidance task data
#' @description Compute bootstrapped split-half reliability for approach-avoidance task data.
#' @param ds a longformat data.frame
#' @param subjvar Quoted name of the participant identifier column
#' @param pullvar Quoted name of the column indicating pull trials.
#' Pull trials should either be represented by 1, or by the second level of a factor.
#' @param targetvar Name of the column indicating trials featuring the target stimulus.
#' Target stimuli should either be represented by 1, or by the second level of a factor.
#' @param rtvar Name of the reaction time column.
#' @param iters Total number of desired iterations. At least 200 are recommended for reasonable confidence intervals;
#' If you want to see plots of your data, 1 iteration is enough.
#' @param plot Create a scatterplot of the AAT scores computed from each half of the data from the last iteration.
#' This is highly recommended, as it helps to identify outliers that can inflate or diminish the reliability.
#' @param include.raw logical indicating whether raw split-half data should be included in the output object.
#' @param cluster pre-specified registered multi-core DoParallel cluster that can be used to speed up computations if multiple calls to aat_splithalf are made.
#' If no cluster is provided, aat_splithalf will start up a cluster each time it is called, which takes some extra time.
#' @param algorithm Function (without brackets or quotes) to be used to compute AAT scores. See \link{Algorithms} for a list of usable algorithms.
#' @param trialdropfunc Function (without brackets or quotes) to be used to exclude outlying trials in each half.
#' The way you handle outliers for the reliability computation should mimic the way you do it in your regular analyses.
#' It is recommended to exclude outlying trials when computing AAT scores using the mean double-dfference scores and regression scoring approaches,
#' but not when using d-scores or median double-difference scores.
#' \itemize{
#' \item \code{prune_nothing} excludes no trials (default)
#' \item \code{trial_prune_3SD} excludes trials deviating more than 3SD from the mean per participant.
#' \item \code{trial_prune_SD_dropcases} removes trials deviating more than a specific number of standard deviations from the participant's mean,
#' and removes participants with an excessive percentage of outliers.
#' Required arguments:
#' \itemize{
#' \item \code{trialsd} - trials deviating more than \code{trialsd} standard deviations from the participant's mean are excluded (optional; default is 3)
#' \item \code{maxoutliers} - participants with a higher percentage of outliers are removed from the data. (optional; default is .15)
#' }
#' \item \code{trial_recode_SD} recodes outlying reaction times to the nearest non-outlying value,
#' with outliers defined as reaction times deviating more than a certain number of standard deviations from the participant's mean. Required argument:
#' \itemize{
#' \item \code{trialsd} - trials deviating more than this many standard deviations from the mean are classified as outliers.
#' }
#' \item \code{trial_prune_percent_subject} and \code{trial_prune_percent_sample} remove trials below and/or above certain percentiles,
#' on a subject-by-subject basis or sample-wide, respectively. The following arguments are available:
#' \itemize{
#' \item \code{lowerpercent} and \code{uppperpercent} (optional; defaults are .01 and .99).
#' }
#' }
#' @param errortrialfunc Function (without brackets or quotes) to apply to an error trial.
#'
#' \itemize{
#' \item \code{prune_nothing} removes no errors (default).
#' \item \code{error_replace_blockmeanplus} replaces error trial reaction times with the block mean, plus an arbitrary extra quantity.
#' If used, the following additional arguments are required:
#' \itemize{
#' \item \code{blockvar} - Quoted name of the block variable (mandatory)
#' \item \code{errorvar} - Quoted name of the error variable, where errors are 1 or TRUE and correct trials are 0 or FALSE (mandatory)
#' \item \code{errorbonus} - Amount to add to the reaction time of error trials. Default is 0.6 (recommended by \code{Greenwald, Nosek, & Banaji, 2003})
#' }
#' \item \code{error_prune_dropcases} removes errors and drops participants if they have more errors than a given percentage. The following arguments are available:
#' \itemize{
#' \item \code{errorvar} - Quoted name of the error variable, where errors are 1 or TRUE and correct trials are 0 or FALSE (mandatory)
#' \item \code{maxerrors} - participants with a higher percentage of errors are excluded from the dataset. Default is .15.
#' }
#' }
#' @param casedropfunc Function (without brackets or quotes) to be used to exclude outlying participant scores in each half.
#' The way you handle outliers here should mimic the way you do it in your regular analyses.
#' \itemize{
#' \item \code{prune_nothing} excludes no participants (default)
#' \item \code{case_prune_3SD} excludes participants deviating more than 3SD from the sample mean.
#' }
#' @param ... Other arguments, to be passed on to the algorithm or outlier rejection functions (see arguments above)
#'
#' @return A list, containing the mean bootstrapped split-half reliability, bootstrapped 95% confidence intervals,
#' a list of data.frames used over each iteration, and a vector containing the split-half reliability of each iteration.
#'
#' @author Sercan Kahveci
#' @seealso \link{q_reliability}
#' @examples #Not Run
#' aat_splithalf(ds=ds2,subjvar="subjectid",pullvar="is_pull",targetvar="is_food",
#'               rtvar="rt",iters=1000,trialdropfunc=trial_prune_3SD,
#'               casedropfunc=case_prune_3SD,plot=T,algorithm=aat_dscore)
#' #Mean reliability: 0.521959
#' #Spearman-Brown-corrected r: 0.6859041
#' #95%CI: [0.4167018, 0.6172474]
#'
#' #Regression Splithalf
#' aat_splithalf(ds=ds2,subjvar="subjectid",pullvar="is_pull",targetvar="is_food",
#'               rtvar="rt",iters=100,trialdropfunc=trial_prune_3SD,
#'               casedropfunc=case_prune_3SD,plot=T,algorithm=aat_regression,
#'               formula = "rt ~ is_pull * is_food",
#'               aatterm = "is_pull:is_food")
#' #Mean reliability: 0.5313939
#' #Spearman-Brown-corrected r: 0.6940003
#' #95%CI: [0.2687186, 0.6749176]
#' @export
aat_splithalf<-function(ds,subjvar,pullvar,targetvar=NULL,rtvar,iters,plot=T,include.raw=F,cluster=NULL,
                        algorithm=c("aat_doublemeandiff","aat_doublemediandiff",
                                    "aat_dscore","aat_dscore_multiblock",
                                    "aat_regression","aat_standardregression",
                                    "aat_doublemedianquotient","aat_doublemeanquotient",
                                    "aat_singlemeandiff","aat_singlemediandiff"),
                        trialdropfunc=c("prune_nothing","trial_prune_3SD","trial_prune_SD_dropcases","trial_recode_SD",
                                        "trial_prune_percent_subject","trial_prune_percent_sample"),
                        errortrialfunc=c("prune_nothing","error_replace_blockmeanplus","error_prune_dropcases"),
                        casedropfunc=c("prune_nothing","case_prune_3SD"),
                        ...){
  packs<-c("magrittr","dplyr","AATtools")

  #Handle arguments
  args<-list(...)
  if(!(algorithm %in% c("aat_singlemeandiff","aat_singlemediandiff","aat_regression","aat_standardregression")) & is.null(targetvar)){
    stop("Argument targetvar missing but required for algorithm!")
  }
  algorithm<-ifelse(is.function(algorithm),deparse(substitute(algorithm)),match.arg(algorithm))
  trialdropfunc<-ifelse(is.function(trialdropfunc),deparse(substitute(trialdropfunc)),match.arg(trialdropfunc))
  casedropfunc<-ifelse(is.function(casedropfunc),deparse(substitute(casedropfunc)),match.arg(casedropfunc))
  errortrialfunc<-ifelse(is.function(errortrialfunc),deparse(substitute(errortrialfunc)),match.arg(errortrialfunc))
  if(errortrialfunc=="error_replace_blockmeanplus"){
    stopifnot(!is.null(args$blockvar),!is.null(args$errorvar))
    if(is.null(args$errorbonus)){ args$errorbonus<- 0.6 }
    if(is.null(args$blockvar)){ args$blockvar<- 0 }
    if(is.null(args$errorvar)){ args$errorvar<- 0 }
  }
  stopifnot(!(algorithm=="aat_dscore_multiblock" & is.null(args$blockvar)))
  if(algorithm %in% c("aat_regression","aat_standardregression")){
    if(!any(c("formula","aatterm") %in% names(args))){
      args$formula<-paste0(rtvar,"~",pullvar,"*",targetvar)
      args$aatterm<-paste0(pullvar,":",targetvar)
      warning("No regression formula or AAT-term provided. Defaulting to formula ",
              args$formula," and AAT-term ",args$aatterm)
    }
  }
  ds%<>%aat_preparedata(subjvar,pullvar,targetvar,rtvar,...)

  #splithalf loop
  if(is.null(cluster)){
    cluster<-makeCluster(detectCores()-1)
    registerDoParallel(cluster)
    custom.cluster<-F
  }else{
    custom.cluster<-T
  }
  results<-
    foreach(iter = seq_len(iters), .packages=packs) %dopar% {
      #Split data
      if(is.null(targetvar)){
        iterds<-ds%>%group_by(!! sym(subjvar), !! sym(pullvar))%>%
          mutate(key=sample(n())%%2)%>%ungroup()
      }else{
        iterds<-ds%>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
          mutate(key=sample(n())%%2)%>%ungroup()
      }
      #Handle outlying trials
      iterds<-do.call(trialdropfunc,c(args,list(ds=iterds,subjvar=subjvar,rtvar=rtvar)))
      #Handle error trials
      iterds<-do.call(errortrialfunc,c(args,list(ds=iterds,subjvar=subjvar,rtvar=rtvar)))
      #intermediate prune of empty cases
      iterds<-drop_empty_cases(iterds,subjvar)

      # abds<-do.call(algorithm,c(list(iterds=iterds,subjvar=subjvar,pullvar=pullvar,
      #                                targetvar=targetvar,rtvar=rtvar),args))

      #Compute AB
      half0set<-iterds%>%filter(key==0)
      half1set<-iterds%>%filter(key==1)
      abds<-merge(
        do.call(algorithm,c(list(ds=half0set,subjvar=subjvar,pullvar=pullvar,
                                 targetvar=targetvar,rtvar=rtvar),args)),
        do.call(algorithm,c(list(ds=half1set,subjvar=subjvar,pullvar=pullvar,
                                 targetvar=targetvar,rtvar=rtvar),args)),
        by=subjvar,suffixes=c("half0","half1"))
      #Remove outlying participants
      abds<-do.call(casedropfunc,list(ds=abds))
      #Compute reliability
      currcorr<-cor(abds$abhalf0,abds$abhalf1,use="complete.obs")
      list(corr=currcorr,abds=abds,rawdata=iterds)
    }
  if(!custom.cluster){
    stopCluster(cluster)
  }

  cors<-sapply(results,FUN=function(x){x$corr})
  ordering<-order(cors)
  cors<-cors[ordering]

  #cat(scan("splithalfmessages.txt",what=character(),quiet=TRUE))
  output<-list(rsplithalf=mean(cors),
               lowerci=cors[round(iters*0.025)],
               upperci=cors[round(iters*0.975)],
               rSB=SpearmanBrown(mean(cors)),
               parameters=c(list(ds=ds,
                                 subjvar=subjvar,
                                 pullvar=pullvar,
                                 targetvar=targetvar,
                                 rtvar=rtvar,
                                 iters=iters,
                                 algorithm=algorithm,
                                 trialdropfunc=trialdropfunc,
                                 errortrialfunc=errortrialfunc,
                                 casedropfunc=casedropfunc),
                            args),
               itercors=cors,
               iterdata=lapply(results,function(x){ x$abds })[ordering]) %>%
    structure(class = "aat_splithalf")
  if(include.raw){
    output$rawiterdata<-lapply(results,function(x){ x$rawdata })[ordering]
  }

  if(plot){ plot(output) }
  return(output)
}

#' @export
#' @rdname aat_splithalf
print.aat_splithalf<-function(x){
  cat("\nr = ",format(x$rsplithalf, digits=2),
      "\nSpearman-Brown-corrected r = ",format(x$rSB,digits=2),
      "\n95%CI = [", format(x$lowerci,digits=2), ", ", format(x$upperci,digits=2),"]\n",
      sep="")
}

#' @title Plot split-half scatterplots
#'
#' @param x an \code{aat_splithalf} object
#' @param type Character argument indicating which iteration should be chosen. Must be an abbreviation of
#' \code{"median"} (default), \code{"minimum"}, \code{"maximum"}, or \code{"random"}.
#'
#' @export
#' @rdname aat_splithalf
#' @examples
#' #Coming soon
plot.aat_splithalf<-function(x,type=c("median","minimum","maximum","random")){
  type<-match.arg(type)
  if(type=="median"){
    title<-"Split-half Scatterplot for Iteration with Median Reliability"
    idx<-ceiling(x$parameters$iters/2)
  }else if(type=="minimum"){
    title<-"Split-half Scatterplot for Iteration with the Lowest Reliability"
    idx<-1
  }else if(type=="maximum"){
    title<-"Split-half Scatterplot for Iteration with the Highest Reliability"
    idx<-x$parameters$iters
  }else if(type=="random"){
    title<-"Split-half Scatterplot for Random Iteration"
    idx<-sample(1:x$parameters$iters,1)
  }
  abds<-x$iterdata[[idx]]
  plot(abds$abhalf0,abds$abhalf1,pch=20,main=
         paste0(title,
                "\n(r = ", round(x$itercors[idx],digits=2),")"),
       xlab="Half 1 computed bias",ylab="Half 2 computed bias")
  text(abds$abhalf0,abds$abhalf1,abds[,1],cex= 0.7, pos=3, offset=0.3)
}

# Outlier removing algorithms ####

#' @export
prune_nothing<-function(ds,...){
  ds
}

#' @export
trial_prune_percent_subject<-function(ds,subjvar,rtvar,lowerpercent=.01,upperpercent=.99,...){
  ds %>% group_by(!!sym(subjvar)) %>% mutate(percentile=rank(!!sym(rtvar))/n()) %>%
    filter(percentile > lowerpercent & percentile< upperpercent) %>% ungroup()
}

#' @export
trial_prune_percent_sample<-function(ds,rtvar,lowerpercent=.01,upperpercent=.99,...){
  ds %>% mutate(percentile=rank(!!sym(rtvar))/n()) %>%
    filter(percentile > lowerpercent & percentile< upperpercent)
}

#' @export
trial_prune_3SD<-function(ds,subjvar,rtvar,...){
  ds %>% group_by(!!sym(subjvar)) %>% filter(abs(scale(!!sym(rtvar))) <3) %>% ungroup()
}

#' @export
trial_prune_SD_dropcases<-function(ds,subjvar,rtvar,trialsd=3,maxoutliers=.15,...){
  ds %>% group_by(!!sym(subjvar)) %>% mutate(is.ol=as.numeric(abs(scale(!!sym(rtvar))) >=3),avg.ol=mean(is.ol)) %>%
    ungroup %>% filter(is.ol==0 & avg.ol<maxoutliers)
}

#' @export
trial_recode_SD<-function(ds,subjvar,rtvar,trialsd=3,...){
  dsa<- ds %>% group_by(!!sym(subjvar)) %>% mutate(ol.z.score=scale(!!sym(rtvar)),
                                      ol.type=(ol.z.score >= trialsd) - (ol.z.score <= -trialsd),
                                      is.ol=abs(ol.type),
                                      ol.max.rt=mean(!!sym(rtvar))+sd(!!sym(rtvar))*trialsd,
                                      ol.min.rt=mean(!!sym(rtvar))-sd(!!sym(rtvar))*trialsd)
  dsa[dsa$is.ol!=0,]<-ifelse(dsa$ol.type==1,dsa$ol.max.rt,dsa$ol.min.rt)
  dsa %>% dplyr::select(-ol.type,-ol.max.rt,-ol.min.rt,-ol.z.score)
}

#' @export
case_prune_3SD<-function(ds,...){
  dplyr::filter(ds,(abhalf0 < mean(abhalf0,na.rm=T)+3*sd(abhalf0,na.rm=T) &
                      abhalf0 > mean(abhalf0,na.rm=T)-3*sd(abhalf0,na.rm=T)) &
                  (abhalf1 < mean(abhalf1,na.rm=T)+3*sd(abhalf1,na.rm=T) &
                     abhalf1 > mean(abhalf1,na.rm=T)-3*sd(abhalf1,na.rm=T)))
}

#Replace error trial latencies with correct block mean RT + 600

#' @export
error_replace_blockmeanplus<-function(ds,subjvar,rtvar,blockvar,errorvar,errorbonus, ...){
  if(!("is.ol" %in% colnames(ds))){ ds$is.ol<-0; }
  ds%<>%group_by(!!sym(subjvar),!!sym(blockvar), key)%>%
    mutate(newrt=mean((!!sym(rtvar))[!(!!sym(errorvar)) & is.ol==0])+errorbonus)%>%ungroup()
  ds[ds[,errorvar]==1,rtvar]<-ds[ds[,errorvar]==1,]$newrt
  dplyr::select(ds,-newrt)
}

#' @export
error_prune_dropcases<-function(ds,subjvar, errorvar, maxerrors = .15, ...){
  ds%>%group_by(!!sym(subjvar), key)%>%
    filter(mean(!!sym(errorvar))<maxerrors & !!sym(errorvar) == F)
}

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
#' \item \code{formula} - a quoted formula to fit to the data;
#' \item \code{aatterm} the quoted random effect within the subject variable that indicates the approach bias; this is usually the interaction of the pull and target terms.
#' }
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
    summarise(ab=(mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=T) -
                    mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=T)) -
                (mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=T) -
                   mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=T)))
}

#' @export
#' @rdname Algorithms
aat_doublemediandiff<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  group_by(ds,!!sym(subjvar)) %>%
    summarise(ab=(median(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=T) -
                    median(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=T)) -
                (median(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=T) -
                   median(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=T)))
}

#' @export
#' @rdname Algorithms
aat_dscore<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  group_by(ds,!!sym(subjvar)) %>%
    summarise(ab=((mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=T) -
                     mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=T)) -
                    (mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=T) -
                       mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=T))) /
                sd(!!sym(rtvar),na.rm=T))
}

#' @param blockvar name of the variable indicating block number
#' @export
#' @rdname Algorithms
#note: this matches sequential blocks with one another.
aat_dscore_multiblock<-function(ds,subjvar,pullvar,targetvar,rtvar,blockvar,...){
  ds %>% mutate(.blockset = floor((!!sym(blockvar) - min(!!sym(blockvar)))/2) ) %>%
    group_by(!!sym(subjvar),.blockset) %>%
    summarise(ab=((mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=T) -
                     mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=T)) -
                    (mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=T) -
                       mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=T))) /
                sd(!!sym(rtvar),na.rm=T)) %>%
    group_by(!!sym(subjvar)) %>% summarise(ab=mean(ab,na.rm=T))
}

#' @param formula A character string containing a regression formula to fit to the data to compute an AAT score
#' @param aatterm The formula term representing the approach bias. Usually this is the interaction of the movement-direction and stimulus-category terms.
#' @export
#' @rdname Algorithms
aat_regression<-function(ds,subjvar,formula,aatterm,...){
  output<-data.frame(pp=unique(ds[[subjvar]]),ab=NA,var=NA)
  for(i in seq_len(nrow(output))){
    mod<-coef(summary(lm(as.formula(formula),data=ds[ds[[subjvar]]==output[i,"pp"],])))
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
    mod<-coef(summary(lm(as.formula(formula),data=ds[ds[[subjvar]]==output[i,"pp"],])))
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
    summarise(ab=(median(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=T) /
                    median(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=T)) -
                (median(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=T) /
                   median(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=T)))
}

#' @export
#' @rdname Algorithms
aat_doublemeanquotient<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  group_by(ds,!!sym(subjvar)) %>%
    summarise(ab=(mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 1),na.rm=T) /
                    mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 1),na.rm=T)) -
                (mean(subset(!!sym(rtvar),!!sym(pullvar)==0 & !!sym(targetvar) == 0),na.rm=T) /
                   mean(subset(!!sym(rtvar),!!sym(pullvar)==1 & !!sym(targetvar) == 0),na.rm=T)))
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

#############################
# bootstrapped bias score computation
#' @title Compute bootstrapped approach-bias scores
#' @description Compute bootstrapped approach-bias scores with confidence intervals.
#' @param ds a longformat data.frame
#' @param subjvar Quoted name of the participant identifier column
#' @param pullvar Quoted name of the column indicating pull trials.
#' Pull trials should either be represented by 1, or by the second level of a factor.
#' @param targetvar Name of the column indicating trials featuring the target stimulus.
#' Target stimuli should either be represented by 1, or by the second level of a factor.
#' @param rtvar Name of the reaction time column.
#' @param iters Total number of desired iterations. At least 200 are required to get confidence intervals that make sense.
#' @param plot Plot the bias scores and their confidence intervals after computation is complete. This gives a good overview of the data.
#' @param algorithm Function (without brackets or quotes) to be used to compute AAT scores. See \link{aat_doublemeandiff} for a list of usable algorithms.
#' @param trialdropfunc Function (without brackets or quotes) to be used to exclude outlying trials in each half.
#' The way you handle outliers for the reliability computation should mimic the way you do it in your regular analyses.
#' It is recommended to exclude outlying trials when computing AAT scores using the mean double-dfference scores and regression scoring approaches,
#' but not when using d-scores or median double-difference scores.
#' \itemize{
#' \item \code{prune_nothing} excludes no trials (default)
#' \item \code{trial_prune_3SD} excludes trials deviating more than 3SD from the mean per participant.
#' \item \code{trial_prune_SD_dropcases} removes trials deviating more than a specific number of standard deviations from the participant's mean,
#' and removes participants with an excessive percentage of outliers.
#' Required arguments:
#' \itemize{
#' \item \code{trialsd} - trials deviating more than \code{trialsd} standard deviations from the participant's mean are excluded (optional; default is 3)
#' \item \code{maxoutliers} - participants with a higher percentage of outliers are removed from the data. (optional; default is .15)
#' }
#' \item \code{trial_recode_SD} recodes outlying reaction times to the nearest non-outlying value,
#' with outliers defined as reaction times deviating more than a certain number of standard deviations from the participant's mean. Required argument:
#' \itemize{
#' \item \code{trialsd} - trials deviating more than this many standard deviations from the mean are classified as outliers.
#' }
#' \item \code{trial_prune_percent_subject} and \code{trial_prune_percent_sample} remove trials below and/or above certain percentiles,
#' on a subject-by-subject basis or sample-wide, respectively. The following arguments are available:
#' \itemize{
#' \item \code{lowerpercent} and \code{uppperpercent} (optional; defaults are .01 and .99).
#' }
#' }
#' @param errortrialfunc Function (without brackets or quotes) to apply to an error trial.
#'
#' \itemize{
#' \item \code{prune_nothing} removes no errors (default).
#' \item \code{error_replace_blockmeanplus} replaces error trial reaction times with the block mean, plus an arbitrary extra quantity.
#' If used, the following additional arguments are required:
#' \itemize{
#' \item \code{blockvar} - Quoted name of the block variable (mandatory)
#' \item \code{errorvar} - Quoted name of the error variable, where errors are 1 or TRUE and correct trials are 0 or FALSE (mandatory)
#' \item \code{errorbonus} - Amount to add to the reaction time of error trials. Default is 0.6 (recommended by \code{Greenwald, Nosek, & Banaji, 2003})
#' }
#' \item \code{error_prune_dropcases} removes errors and drops participants if they have more errors than a given percentage. The following arguments are available:
#' \itemize{
#' \item \code{errorvar} - Quoted name of the error variable, where errors are 1 or TRUE and correct trials are 0 or FALSE (mandatory)
#' \item \code{maxerrors} - participants with a higher percentage of errors are excluded from the dataset. Default is .15.
#' }
#' }
#' @param ... Other arguments, to be passed on to the algorithm or outlier rejection functions (see arguments above)
#'
#'
#' @return A list, containing bootstrapped bias scores, their variance, bootstrapped 95 percent confidence intervals,
#' the number of iterations, and a matrix of bias scores for each iteration.
#'
#' @author Sercan Kahveci
#' @examples
#' @export
aat_bootstrap<-function(ds,subjvar,pullvar,targetvar=NULL,rtvar,iters,plot=T,include.raw=F,
                        algorithm=c("aat_doublemeandiff","aat_doublemediandiff",
                                    "aat_dscore","aat_dscore_multiblock",
                                    "aat_regression","aat_standardregression",
                                    "aat_doublemeanquotient","aat_doublemedianquotient",
                                    "aat_singlemeandiff","aat_singlemediandiff"),
                        trialdropfunc=c("prune_nothing","trial_prune_3SD","trial_prune_SD_dropcases","trial_recode_SD",
                                        "trial_prune_percent_subject","trial_prune_percent_sample"),
                        errortrialfunc=c("prune_nothing","error_replace_blockmeanplus","error_prune_dropcases"),
                        ...){
  packs<-c("magrittr","dplyr","AATtools")

  #Handle arguments
  args<-list(...)
  if(!(algorithm %in% c("aat_singlemeandiff","aat_singlemediandiff","aat_regression","aat_standardregression")) & is.null(targetvar)){
    stop("Argument targetvar missing but required for algorithm!")
  }
  algorithm<-ifelse(is.function(algorithm),deparse(substitute(algorithm)),match.arg(algorithm))
  trialdropfunc<-ifelse(is.function(trialdropfunc),deparse(substitute(trialdropfunc)),match.arg(trialdropfunc))
  errortrialfunc<-ifelse(is.function(errortrialfunc),deparse(substitute(errortrialfunc)),match.arg(errortrialfunc))
  if(errortrialfunc=="error_replace_blockmeanplus"){
    stopifnot(!is.null(args$blockvar),!is.null(args$errorvar))
    if(is.null(args$errorbonus)){ args$errorbonus<- 0.6 }
    if(is.null(args$blockvar)){ args$blockvar<- 0 }
    if(is.null(args$errorvar)){ args$errorvar<- 0 }
  }
  stopifnot(!(algorithm=="aat_dscore_multiblock" & is.null(args$blockvar)))
  if(algorithm %in% c("aat_regression","aat_standardregression")){
    if(!any(c("formula","aatterm") %in% names(args))){
      args$formula<-paste0(rtvar,"~",pullvar,"*",targetvar)
      args$aatterm<-paste0(pullvar,":",targetvar)
      warning("No formula and/or AAT-term provided. Defaulting to formula ",
              args$formula," and AAT-term ",args$aatterm)
    }
  }
  ds %<>% aat_preparedata(subjvar,pullvar,targetvar,rtvar,...) %>% mutate(key=1)

  #bootstrap loop
  cluster<-makeCluster(detectCores()-1)
  registerDoParallel(cluster)
  results<-
    foreach(iter = seq_len(iters), .packages=packs, .combine=cbind) %dopar% {
      #Split data
      iterds<-ds %>% group_by(!!sym(subjvar), !!sym(pullvar), !!sym(targetvar)) %>%
        sample_n(size=n(),replace=T) %>% ungroup()

      #Handle outlying trials
      iterds<-do.call(trialdropfunc,list(ds=iterds,subjvar=subjvar,rtvar=rtvar))
      #Handle error trials
      iterds<-do.call(errortrialfunc,list(ds=iterds,subjvar=subjvar,rtvar=rtvar,
                                          blockvar=args$blockvar,errorvar=args$errorvar,
                                          errorbonus=args$errorbonus))

      abds<-do.call(algorithm,c(list(ds=iterds,subjvar=subjvar,pullvar=pullvar,
                                     targetvar=targetvar,rtvar=rtvar),args))

      #colnames(abds)<-c(subjvar,paste0("iter", formatC(iter, width = nchar(iters), format = "d", flag = "0")))
      outvar<-abds$ab
      names(outvar)<-abds[[subjvar]]
      outvar
    }
  stopCluster(cluster)

  statset<-data.frame(ppidx=rownames(results),
                      bias=rowMeans(results,na.rm=T),
                      var=apply(results,MARGIN = 1,FUN=var,na.rm=T),
                      lowerci=apply(results,MARGIN=1,FUN=function(x){quantile(x,0.025,na.rm=T)}),
                      upperci=apply(results,MARGIN=1,FUN=function(x){quantile(x,0.975,na.rm=T)}))
  statset$ci<-statset$upperci-statset$lowerci

  #q-reliability
  bv<-var(statset$bias,na.rm=T)
  wv<-mean(statset$var,na.rm=T)
  q<-(bv-wv)/(bv+wv)

  output<-list(bias=statset,
               reliability=q,
               parameters=c(list(ds=ds,
                                 subjvar=subjvar,
                                 pullvar=pullvar,
                                 targetvar=targetvar,
                                 rtvar=rtvar,
                                 iters=iters,
                                 algorithm=algorithm,
                                 trialdropfunc=trialdropfunc,
                                 errortrialfunc=errortrialfunc),args)) %>%
    structure(class = "aat_bootstrap")
  if(include.raw){
    output$iterdata<-results
  }
  if(plot){ plot(output) }
  return(output)
}

#' @export
#' @rdname aat_bootstrap
print.aat_bootstrap<-function(x){
  cat("Bootstrapped bias scores and confidence intervals",
      "Mean bias score: ", mean(x$bias$bias,na.rm=T),
      "Mean confidence interval: ",mean(x$bias$ci,na.rm=T),
      "reliability: q = ",mean(x$bias$bias,na.rm=T),
      "Number of iterations: ",xparameters$iters,sep="")
}

#' @export
#' @rdname aat_bootstrap
plot.aat_bootstrap <- function(x){
  statset<-x$bias
  statset<-statset[!is.na(statset$bias) & !is.na(statset$upperci) & !is.na(statset$lowerci),]
  rank<-rank(statset$bias)
  wideness<-max(statset$upperci) - min(statset$lowerci)
  plot(x=statset$bias,y=rank,xlim=c(min(statset$lowerci)-0.01*wideness,max(statset$upperci)+0.01*wideness),
       xlab="Bias score",main=paste0("Individual bias scores with 95%CI",
                                     "\nEstimated reliability: q = ",x$reliability))
  segments(x0=statset$lowerci,x1=statset$bias-0.005*wideness,y0=rank,y1=rank)
  segments(x0=statset$bias+0.005*wideness,x1=statset$upperci,y0=rank,y1=rank)
  abline(v=0)
  #text(x=statset$bias,y=statset$rownr,labels=statset$ppidx,cex=0.5)
}

#' @title Compute simple AAT scores
#' @description Compute simple AAT scores, with optional outlier exclusion and error trial recoding.
#' @param ds a long-format data.frame
#' @param subjvar column name of subject variable
#' @param pullvar column name of pull/push indicator variable, must be numeric or logical (where pull is 1 or TRUE)
#' @param targetvar column name of target stimulus indicator, must be numeric or logical (where target is 1 or TRUE)
#' @param rtvar column name of reaction time variable
#' @param trialdropfunc Function (without brackets or quotes) to be used to exclude outlying trials in each half.
#' The way you handle outliers for the reliability computation should mimic the way you do it in your regular analyses.
#' It is recommended to exclude outlying trials when computing AAT scores using the mean double-dfference scores and regression scoring approaches,
#' but not when using d-scores or median double-difference scores.
#' \itemize{
#' \item \code{prune_nothing} excludes no trials (default)
#' \item \code{trial_prune_3SD} excludes trials deviating more than 3SD from the mean per participant.
#' \item \code{trial_prune_SD_dropcases} removes trials deviating more than a specific number of standard deviations from the participant's mean,
#' and removes participants with an excessive percentage of outliers.
#' Required arguments:
#' \itemize{
#' \item \code{trialsd} - trials deviating more than \code{trialsd} standard deviations from the participant's mean are excluded (optional; default is 3)
#' \item \code{maxoutliers} - participants with a higher percentage of outliers are removed from the data. (optional; default is .15)
#' }
#' \item \code{trial_recode_SD} recodes outlying reaction times to the nearest non-outlying value,
#' with outliers defined as reaction times deviating more than a certain number of standard deviations from the participant's mean. Required argument:
#' \itemize{
#' \item \code{trialsd} - trials deviating more than this many standard deviations from the mean are classified as outliers.
#' }
#' \item \code{trial_prune_percent_subject} and \code{trial_prune_percent_sample} remove trials below and/or above certain percentiles,
#' on a subject-by-subject basis or sample-wide, respectively. The following arguments are available:
#' \itemize{
#' \item \code{lowerpercent} and \code{uppperpercent} (optional; defaults are .01 and .99).
#' }
#' }
#' @param errortrialfunc Function (without brackets or quotes) to apply to an error trial.
#'
#' \itemize{
#' \item \code{prune_nothing} removes no errors (default).
#' \item \code{error_replace_blockmeanplus} replaces error trial reaction times with the block mean, plus an arbitrary extra quantity.
#' If used, the following additional arguments are required:
#' \itemize{
#' \item \code{blockvar} - Quoted name of the block variable (mandatory)
#' \item \code{errorvar} - Quoted name of the error variable, where errors are 1 or TRUE and correct trials are 0 or FALSE (mandatory)
#' \item \code{errorbonus} - Amount to add to the reaction time of error trials. Default is 0.6 (recommended by \code{Greenwald, Nosek, & Banaji, 2003})
#' }
#' \item \code{error_prune_dropcases} removes errors and drops participants if they have more errors than a given percentage. The following arguments are available:
#' \itemize{
#' \item \code{errorvar} - Quoted name of the error variable, where errors are 1 or TRUE and correct trials are 0 or FALSE (mandatory)
#' \item \code{maxerrors} - participants with a higher percentage of errors are excluded from the dataset. Default is .15.
#' }
#' }
#' @param ... Other arguments, to be passed on to the algorithm or outlier rejection functions (see arguments above)
#'
#'
#' @export
aat_compute<-function(ds,subjvar,pullvar,targetvar=NULL,rtvar,
                        algorithm=c("aat_doublemeandiff","aat_doublemediandiff",
                                    "aat_dscore","aat_dscore_multiblock",
                                    "aat_regression","aat_standardregression",
                                    "aat_doublemeanquotient","aat_doublemedianquotient",
                                    "aat_singlemeandiff","aat_singlemediandiff"),
                        trialdropfunc=c("prune_nothing","trial_prune_3SD","trial_prune_SD_dropcases","trial_recode_SD",
                                        "trial_prune_percent_subject","trial_prune_percent_sample"),
                        errortrialfunc=c("prune_nothing","error_replace_blockmeanplus","error_prune_dropcases"),
                        ...){
  #Handle arguments
  args<-list(...)
  if(!(algorithm %in% c("aat_singlemeandiff","aat_singlemediandiff","aat_regression","aat_standardregression")) & is.null(targetvar)){
    stop("Argument targetvar missing but required for algorithm!")
  }
  algorithm<-ifelse(is.function(algorithm),deparse(substitute(algorithm)),match.arg(algorithm))
  trialdropfunc<-ifelse(is.function(trialdropfunc),deparse(substitute(trialdropfunc)),match.arg(trialdropfunc))
  errortrialfunc<-ifelse(is.function(errortrialfunc),deparse(substitute(errortrialfunc)),match.arg(errortrialfunc))
  if(errortrialfunc=="error_replace_blockmeanplus"){
    stopifnot(!is.null(args$blockvar),!is.null(args$errorvar))
    if(is.null(args$errorbonus)){ args$errorbonus<- 0.6 }
    if(is.null(args$blockvar)){ args$blockvar<- 0 }
    if(is.null(args$errorvar)){ args$errorvar<- 0 }
  }
  stopifnot(!(algorithm=="aat_dscore_multiblock" & is.null(args$blockvar)))
  if(algorithm %in% c("aat_regression","aat_standardregression")){
    if(!any(c("formula","aatterm") %in% names(args))){
      args$formula<-paste0(rtvar,"~",pullvar,"*",targetvar)
      args$aatterm<-paste0(pullvar,":",targetvar)
      warning("No regression formula or AAT-term provided. Defaulting to formula ",
              args$formula," and AAT-term ",args$aatterm)
    }
  }
  ds %<>% aat_preparedata(subjvar,pullvar,targetvar,rtvar,...) %>% mutate(key=1)

  #Handle outlying trials
  iterds<-do.call(trialdropfunc,list(ds=ds,subjvar=subjvar,rtvar=rtvar))
  #Handle error trials
  iterds<-do.call(errortrialfunc,list(ds=ds,subjvar=subjvar,rtvar=rtvar,
                                      blockvar=args$blockvar,errorvar=args$errorvar,
                                      errorbonus=args$errorbonus))

  abds<-do.call(algorithm,c(list(ds=ds,subjvar=subjvar,pullvar=pullvar,
                                 targetvar=targetvar,rtvar=rtvar),args))

  abds <- merge(x=abds,by=subjvar,all=T,y=iterds %>% group_by(!!sym(subjvar)) %>% summarise(trials=n()))
  return(abds)
}

# utils ####
#' Correct a correlation coefficient for being based on only a subset of the data.
#' @title Spearman-Brown corrections for Correlation Coefficients
#' @description Perform a Spearman-Brown correction on the provided correlation score.
#'
#' @param corr To-be-corrected correlation coefficient
#' @param ntests An integer indicating how many times larger the full test is, for which the corrected correlation coefficient is being computed.
#' When \code{ntests=2}, the formula will compute what the correlation coefficient would be if the test were twice as long.
#' @param fix.negative Determines how to deal with a negative value. "nullify" sets it to zero,
#' "bilateral" applies the correction as if it were a positive number, and then sets it to negative. "none" gives the raw value.
#' @return Spearman-Brown-corrected correlation coefficient.
#' @export
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

aat_preparedata<-function(ds,subjvar,pullvar,targetvar=NULL,rtvar,...){
  args<-list(...)
  stopifnot(all(c(subjvar,pullvar,targetvar,rtvar) %in% colnames(ds)))
  ds[[subjvar]]%<>%as.factor()
  if(is.logical(ds[,pullvar])){
    warning("Recoded ",pullvar," from logical to numeric. Please make sure that FALSE ",
            "represents push trials and TRUE represents pull trials")
    ds[,pullvar]%<>%as.numeric()
  }
  if(is.factor(ds[,pullvar])){
    warning("Recoded ",pullvar," from factor to numeric. Please make sure that ",
            levels(ds[,pullvar])[1], " represents push trials and ",levels(ds[,pullvar])[2],
            " represents pull trials")
    ds[,pullvar]<-as.numeric(ds[,pullvar])-1
  }
  if(!is.null(targetvar)){
    if(is.logical(ds[,targetvar])){
      warning("Recoded ",targetvar," from logical to numeric. Please make sure that FALSE ",
              "represents control/neutral stimuli and TRUE represents target stimuli")
      ds[,targetvar]%<>%as.numeric()
    }
    if(is.factor(ds[,targetvar])){
      warning("Recoded ",targetvar," from factor to numeric. Please make sure that ",
              levels(ds[,targetvar])[1], " represents control/neutral stimuli and ",levels(ds[,targetvar])[2],
              " represents target stimuli")
      ds[,targetvar]<-as.numeric(ds[,targetvar])-1
    }
  }
  rmindices<- which(is.na(ds[[subjvar]]) & is.na(ds[[pullvar]]) & is.na(ds[[targetvar]]) & is.na(ds[[rtvar]]))
  if(length(rmindices)>0){
    ds<-ds[-rmindices,]
    warning("Removed ",length(rmindices)," rows due to presence of NA in critical variable(s)")
  }

  cols<-c(subjvar,pullvar,targetvar,rtvar,args$errorvar,args$blockvar)
  if("formula" %in% names(args)){
    formterms <- terms(args$formula) %>% attr("variables") %>% as.character()
    formterms <- formterms[-1]
    if(any(!(formterms %in% colnames(ds)))){
      stop("Formula term(s) ",paste(formterms[!(formterms %in% colnames(ds))],collapse=", ")," missing from dataset")
    }
    cols <- c(cols,formterms)
  }
  ds<-ds[,cols]
  return(ds)
}

val_between<-function(x,lb,ub){x>lb & x<ub}

drop_empty_cases<-function(iterds,subjvar){
  iterds%>%group_by(!!sym(subjvar))%>%filter(val_between(mean(key),0,1))
}


#' Compute psychological experiment reliability
#' @description This function can be used to compute an exact reliability score for a psychological task whose results involve a difference score.
#' The resulting q coefficient is equivalent to the average all possible split-half reliability scores.
#' @param ds a long-format data.frame
#' @param subjvar name of the subject variable
#' @param formula a formula predicting the participant's reaction time using trial-level variables such as movement direction and stimulus category
#' @param aatterm a string denoting the term in the formula that contains the participant's approach bias
#'
#' @return a qreliability object, containing the reliability coefficient, and a data.frame with participants' bias scores and score variance.
#' @export
#' @author Sercan Kahveci
#' @examples
#' q_reliability(ds=dataset,subjvar="subjectid",formula= rt ~ is_pull * is_food, aatterm = "is_pull:is_food")
#'
q_reliability<-function(ds,subjvar,formula,aatterm=NA){
  # argument checks
  cols<-c(subjvar,rownames(attr(terms(formula),"factors")))
  stopifnot(all(cols %in% colnames(ds)))
  ds<-ds[apply(!is.na(ds[,cols]),MARGIN=1,FUN=all),]

  # functional part
  coefs<-data.frame(pp=unique(ds[[subjvar]]),ab=NA,var=NA)
  for(u in 1:nrow(coefs)){
    iterset<-ds[ds[[subjvar]]==coefs[u,]$pp,]
    mod<-lm(formula,data=iterset)
    coefs[u,]$ab <- -coef(mod)[ifelse(is.na(aatterm),length(coef(mod)),aatterm)]
    coefs[u,]$var <- (diag(vcov(mod)))[ifelse(is.na(aatterm),length(coef(mod)),aatterm)] # squared standard error
  }

  bv<-var(coefs$ab,na.rm=T)
  wv<-mean(coefs$var,na.rm=T)
  q<-(bv-wv)/(bv+wv)
  # alternative: (2*wv)/(bv*wv) -1

  return(structure(list(q=q,coefs=coefs),class="qreliability"))
}

#' @export
print.qreliability<-function(x){
  cat("q = ",x$q,"\n",sep="")
}

#' @export
#' @rdname q_reliability
plot.qreliability<-function(x){
  bv<-var(x$coefs$ab,na.rm=T) / nrow(x$coefs)*1.96 *2
  wv<-mean(x$coefs$var,na.rm=T) / nrow(x$coefs)*1.96 *2
  plotset<-data.frame(x=mean(x$coefs$ab) + cos(0:100 / 100 * 2*pi)*bv * 1/2*sqrt(2) - sin(0:100 / 100 * 2*pi)*wv * 1/2*sqrt(2),
                      y=mean(x$coefs$ab) + cos(0:100 / 100 * 2*pi)*bv * 1/2*sqrt(2) + sin(0:100 / 100 * 2*pi)*wv * 1/2*sqrt(2))
  plot(plotset$x,plotset$y,type="l",main=paste0("Reliability\n","q = ",round(x$q,digits=2)),xlab="Participants' scores",ylab="Participants' scores")
  points(x$coefs$ab,x$coefs$ab)
  dispval<-(bv+wv)/100
  plotset<-data.frame(xstart=c(x$coefs$ab+dispval,x$coefs$ab-dispval),
                      ystart=c(x$coefs$ab-dispval,x$coefs$ab+dispval),
                      xend=c(x$coefs$ab+sqrt(x$coefs$var)*1.96 *1/2*sqrt(2),
                             x$coefs$ab-sqrt(x$coefs$var)*1.96 *1/2*sqrt(2)),
                      yend=c(x$coefs$ab-sqrt(x$coefs$var)*1.96 *1/2*sqrt(2),
                             x$coefs$ab+sqrt(x$coefs$var)*1.96 *1/2*sqrt(2)))
  segments(plotset$xstart,plotset$ystart,plotset$xend,plotset$yend)
}









