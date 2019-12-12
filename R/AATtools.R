#AATtools

.onLoad<-function(libname, pkgname){
  registerS3method("print",class="aat_splithalf",method=print.aat_splithalf)
  registerS3method("plot",class="aat_splithalf",method=plot.aat_splithalf)
  registerS3method("plot",class="aat_bootstrap",method=plot.aat_bootstrap)
  packageStartupMessage("Thank you for loading AATtools v0.0.1")
}

# splithalf engine ####

#multicore splithalf
#' @title Compute the bootstrapped split-half reliability for approach-avoidance task data
#' @description Compute bootstrapped split-half reliability for approach-avoidance task data.
#' \code{aat_splithalf()} uses multicore computation, which is fast, but provides no clear output when there are errors.
#' \code{aat_splithalf_singlecore()} is much slower, but more easily debugged.
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
#' @param algorithm Function (without brackets or quotes) to be used to compute AAT scores. See \link{Algorithms} for a list of usable algorithms.
#' @param trialdropfunc Function (without brackets or quotes) to be used to exclude outlying trials in each half.
#' The way you handle outliers for the reliability computation should mimic the way you do it in your regular analyses.
#' It is recommended to exclude outlying trials when computing AAT scores using the mean double-dfference scores and multilevel scoring approaches,
#' but not when using d-scores or median double-difference scores.
#' \code{prune_nothing} excludes no trials, \code{trial_prune_3SD} excludes trials deviating more than 3SD from the mean per participant.
#' \code{trial_prune_SD_dropcases} allows you to set the maximum standard deviation to include using argument \code{trialsd} (default is 3)
#' and prune participants altogether if they have more than a certain proportion of outliers using argument \code{maxoutliers} (default is .15)
#' @param errortrialfunc Function (without brackets or quotes) to apply to an error trial.
#'
#' \code{error_replace_blockmeanplus} replaces error trial reaction times with the block mean plus an arbitrary extra amount of time.
#' If used, the following additional arguments are required:
#' \itemize{
#' \item \code{blockvar} - Quoted name of the block variable
#' \item \code{errorvar} - Quoted name of the error variable, where errors are 1 or TRUE and correct trials are 0 or FALSE
#' \item \code{errorbonus} - Amount to add to the reaction time of error trials. Default is 0.6 (recommended by \code{Greenwald, Nosek, & Banaji, 2003})
#' }
#' \code{error_prune_dropcases} removes errors and drops participants if they have a larger proportion of errors
#'  than given by argument \code{maxerrors}, default is .15
#' @param casedropfunc Function (without brackets or quotes) to be used to exclude outlying participant scores in each half.
#' The way you handle outliers here should mimic the way you do it in your regular analyses.
#' \code{prune_nothing} excludes no participants, while \code{case_prune_3SD} excludes participants deviating more than 3SD from the sample mean.
#' @param ... Other arguments, to be passed on to the algorithm functions (see \code{algorithm} above)
#'
#' @return A list, containing the mean bootstrapped split-half reliability, bootstrapped 95% confidence intervals,
#' a list of data.frames used over each iteration, and a vector containing the split-half reliability of each iteration.
#'
#' @author Sercan Kahveci
#' @seealso \link{plot.aat_splithalf}
#' @examples #Not Run
#' aat_splithalf(ds=ds2,subjvar="subjectid",pullvar="is_pull",targetvar="is_food",
#'               rtvar="rt",iters=1000,trialdropfunc=trial_prune_3SD,
#'               casedropfunc=case_prune_3SD,plot=T,algorithm=aat_dscore)
#' #Mean reliability: 0.521959
#' #Spearman-Brown-corrected r: 0.6859041
#' #95%CI: [0.4167018, 0.6172474]
#'
#' #Multilevel Splithalf
#' aat_splithalf(ds=ds2,subjvar="subjectid",pullvar="is_pull",targetvar="is_food",
#'               rtvar="rt",iters=100,trialdropfunc=trial_prune_3SD,
#'               casedropfunc=case_prune_3SD,plot=T,algorithm=aat_multilevelscore,
#'               formula = "rt ~ is_pull * is_food + (is_pull * is_food | subjectid)",
#'               aatterm = "is_pull:is_food")
#' #Mean reliability: 0.5313939
#' #Spearman-Brown-corrected r: 0.6940003
#' #95%CI: [0.2687186, 0.6749176]
#' @export
aat_splithalf<-function(ds,subjvar,pullvar,targetvar,rtvar,iters,plot=T,
                        algorithm=c("aat_doublemeandiff","aat_doublemediandiff","aat_dscore",
                                    "aat_dscore_multiblock","aat_multilevelscore"),
                        trialdropfunc=c("prune_nothing","trial_prune_3SD","trial_prune_SD_dropcases"),
                        errortrialfunc=c("prune_nothing","error_replace_blockmeanplus","error_prune_SD_dropcases"),
                        casedropfunc=c("prune_nothing","case_prune_3SD"),
                        ...){
  packs<-c("magrittr","dplyr","AATtools")

  #Handle arguments
  args<-list(...)
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
  if(algorithm=="aat_multilevelscore"){
    packs<-c(packs,"lme4")
    if(!any(c("formula","aatterm") %in% names(args))){
      args$formula<-paste0(rtvar,"~",1,"+(",pullvar,"*",targetvar,"|",subjvar,")")
      args$aatterm<-paste0(pullvar,":",targetvar)
      warning("No multilevel formula or AAT-term provided. Defaulting to formula ",
              args$formula," and AAT-term ",args$aatterm)
    }
  }
  ds%<>%aat_preparedata(subjvar,pullvar,targetvar,rtvar,...)

  #splithalf loop
  cluster<-makeCluster(detectCores()-1)#,outfile="splithalfmessages.txt")
  registerDoParallel(cluster)
  results<-
    foreach(iter = seq_len(iters), .packages=packs) %dopar% {
      #Split data
      iterds<-ds%>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
        mutate(key=sample(n())%%2)%>%ungroup()
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
  stopCluster(cluster)

  cors<-sapply(results,FUN=function(x){x$corr})
  ordering<-order(cors)
  cors<-cors[ordering]

  #cat(scan("splithalfmessages.txt",what=character(),quiet=TRUE))
  output<-list(rsplithalf=mean(cors),
               lowerci=cors[round(iters*0.025)],
               upperci=cors[round(iters*0.975)],
               rSB=SpearmanBrown(mean(cors)),
               iters=iters,
               itercors=cors,
               iterdata=lapply(results,function(x){ x$abds })[ordering],
               rawiterdata=lapply(results,function(x){ x$rawdata })[ordering]) %>%
    structure(class = "aat_splithalf")
  if(plot){ plot(output) }
  return(output)
}

#Singlecore splithalf (slower but produces output)
#' @rdname aat_splithalf
#' @export
aat_splithalf_singlecore<-function(ds,subjvar,pullvar,targetvar,rtvar,iters,plot=T,
                                   algorithm=c("aat_doublemeandiff","aat_doublemediandiff","aat_dscore",
                                               "aat_dscore_multiblock","aat_multilevelscore"),
                                   trialdropfunc=c("prune_nothing","trial_prune_3SD"),
                                   errortrialfunc=c("prune_nothing","error_replace_blockmeanplus"),
                                   casedropfunc=c("prune_nothing","case_prune_3SD"),
                                   ...){

  #Handle arguments
  args<-list(...)
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
  if(algorithm=="aat_multilevelscore" & !any(c("formula","aatterm") %in% names(args))){
    args$formula<-paste0(rtvar,"~",1,"+(",pullvar,"*",targetvar,"|",subjvar,")")
    args$aatterm<-paste0(pullvar,":",targetvar)
    warning("No multilevel formula or AAT-term provided. Defaulting to formula ",
            args$formula," and AAT-term ",args$aatterm)
  }
  ds%<>%aat_preparedata(subjvar,pullvar,targetvar,rtvar)

  #splithalf loop
  results<-list()
  for(iter in seq_len(iters)){
    #Split data
    iterds<-ds%>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
      mutate(key=sample(n())%%2)%>%ungroup()
    #Handle outlying trials
    iterds<-do.call(trialdropfunc,list(ds=iterds,subjvar=subjvar,rtvar=rtvar))
    #Handle error trials
    iterds<-do.call(errortrialfunc,list(ds=iterds,subjvar=subjvar,rtvar=rtvar,
                                        blockvar=args$blockvar,errorvar=args$errorvar,
                                        errorbonus=args$errorbonus))

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
    results[[iter]]<-list(corr=currcorr,abds=abds)
  }
  cors<-sapply(results,FUN=function(x){x$corr}) %>% sort
  #cat(scan("splithalfmessages.txt",what=character(),quiet=TRUE))
  output<-list(rsplithalf=mean(cors),
               lowerci=cors[round(iters*0.025)],
               upperci=cors[round(iters*0.975)],
               rSB=SpearmanBrown(mean(cors)),
               iters=iters,
               itercors=sapply(results,function(x){ x$corr }),
               iterdata=lapply(results,function(x){ x$abds })) %>%
    structure(class = "aat_splithalf")
  if(plot){ plot(output) }
  return(output)
}


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
#'
#' @examples
#' #Coming soon
plot.aat_splithalf<-function(x,type=c("median","minimum","maximum","random")){
  type<-match.arg(type)
  if(type=="median"){
    title<-"Split-half Scatterplot for Iteration with Median Reliability"
    idx<-ceiling(x$iters/2)
  }else if(type=="minimum"){
    title<-"Split-half Scatterplot for Iteration with the Lowest Reliability"
    idx<-1
  }else if(type=="maximum"){
    title<-"Split-half Scatterplot for Iteration with the Highest Reliability"
    idx<-x$iters
  }else if(type=="random"){
    title<-"Split-half Scatterplot for Random Iteration"
    idx<-sample(1:x$iters,1)
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
trial_prune_3SD<-function(ds,subjvar,rtvar){
  ds$ol.rt.var<-ds[[rtvar]]
  ds%>%group_by(!! sym(subjvar),key)%>%
    mutate( prune.ol.mean=mean(ol.rt.var,na.rm=T),prune.ol.sd=sd(ol.rt.var,na.rm=T))%>%
    dplyr::filter((ol.rt.var<prune.ol.mean+3*prune.ol.sd) &
                    (ol.rt.var>prune.ol.mean-3*prune.ol.sd) ) %>%
    dplyr::select(-prune.ol.mean,-prune.ol.sd,-ol.rt.var)
}

#' @export
trial_prune_SD_dropcases<-function(ds,subjvar,rtvar,trialsd=3,maxoutliers=.15){
  ds$ol.rt.var<-ds[[rtvar]]
  ds%>%group_by(!! sym(subjvar),key)%>%
    mutate( prune.ol.mean=mean(ol.rt.var,na.rm=T),prune.ol.sd=sd(ol.rt.var,na.rm=T))%>%
    dplyr::filter(abs(scale(ol.rt.var)) < trialsd & (mean(abs(scale(ol.rt.var)) > trialsd)) < maxoutliers) %>%
    dplyr::select(-prune.ol.mean,-prune.ol.sd,-ol.rt.var)
}

#' @export
case_prune_3SD<-function(ds){
  dplyr::filter(ds,(abhalf0 < mean(abhalf0,na.rm=T)+3*sd(abhalf0,na.rm=T) &
                      abhalf0 > mean(abhalf0,na.rm=T)-3*sd(abhalf0,na.rm=T)) &
                  (abhalf1 < mean(abhalf1,na.rm=T)+3*sd(abhalf1,na.rm=T) &
                     abhalf1 > mean(abhalf1,na.rm=T)-3*sd(abhalf1,na.rm=T)))
}

#Replace error trial latencies with correct block mean RT + 600

#' @export
error_replace_blockmeanplus<-function(ds,subjvar,rtvar,blockvar,errorvar,errorbonus, ...){
  ds%<>%group_by(!!sym(subjvar),!!sym(blockvar), key)%>%
    mutate(newrt=mean((!!sym(rtvar))[!(!!sym(errorvar))])+errorbonus)%>%ungroup()
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
#' \item \code{aat_multilevelscore} fits a multilevel model using lme4 and extracts a random effect serving as AAT score. When using this function, additional arguments must be provided:
#' \itemize{
#' \item \code{formula} - a quoted formula to fit to the data;
#' \item \code{aatterm} the quoted random effect within the subject variable that indicates the approach bias; this is usually the interaction of the pull and target terms.
#' }
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
  ds%<>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
    summarise(means =mean(!! sym(rtvar),na.rm=T))%>%group_by()
  abds<-ds%>%mutate(groupcol=paste0("pull",!!sym(pullvar),"target",!!sym(targetvar)))%>%
    dplyr::select(!! subjvar, groupcol,means)%>%tidyr::spread("groupcol",value="means",drop=F)
  abds%<>%
    mutate(ab=(pull0target1-pull1target1)-(pull0target0-pull1target0))%>%
    dplyr::select(!!sym(subjvar),ab)
  return(abds)
}

#' @export
#' @rdname Algorithms
aat_doublemediandiff<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  ds%<>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
    summarise(medians =median(!! sym(rtvar),na.rm=T))%>%group_by()
  abds<-ds%>%mutate(groupcol=paste0("pull",!!sym(pullvar),"target",!!sym(targetvar)))%>%
    dplyr::select(!! subjvar, groupcol,medians)%>%tidyr::spread("groupcol",value="medians",drop=F)
  abds%<>%
    mutate(ab=(pull0target1-pull1target1)-(pull0target0-pull1target0))%>%
    dplyr::select(!!sym(subjvar),ab)
  return(abds)
}

#' @export
#' @rdname Algorithms
aat_dscore<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  setmeans <- ds %>% group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar)) %>%
    summarise(means =mean(!! sym(rtvar),na.rm=T)) %>% group_by() %>%
    mutate(groupcol=paste0("pull",!! sym(pullvar),"target",!! sym(targetvar))) %>%
    dplyr::select(-!!sym(pullvar),-!!sym(targetvar)) %>% tidyr::spread(key=groupcol,value=means)

  setsds <- ds%>%group_by(!! sym(subjvar)) %>%
    summarise(sds =sd(!! sym(rtvar),na.rm=T))

  abds <- merge(setmeans,setsds,by=subjvar) %>%
    mutate(ab=((pull0target1-pull1target1)-(pull0target0-pull1target0))/sds) %>%
    dplyr::select(!!sym(subjvar),ab)
  return(abds)
}

#' @param blockvar name of the variable indicating block number
#' @export
#' @rdname Algorithms
#note: this matches sequential blocks with one another.
aat_dscore_multiblock<-function(ds,subjvar,pullvar,targetvar,rtvar,blockvar,...){
  setmeans <- ds %>% group_by(!!sym(subjvar), !!sym(pullvar), !!sym(targetvar), !!sym(blockvar)) %>%
    summarise(means=mean(!!sym(rtvar),na.rm=T)) %>% group_by() %>%
    mutate(groupcol=paste0("pull", !!sym(pullvar),"target", !!sym(targetvar)),
           blockset=LETTERS[1+floor((!!sym(blockvar) - min(!!sym(blockvar)))/2)]) %>%
    dplyr::select(-!!sym(pullvar),-!!sym(targetvar),-!!sym(blockvar)) %>%
    tidyr::spread(key=groupcol,value=means) %>% mutate(mergekey=paste0(!!sym(subjvar),blockset))
  setsds <- ds%>%group_by(!!sym(subjvar),
                          blockset=LETTERS[1+floor((!!sym(blockvar)-min(!!sym(blockvar)))/2)]) %>%
    summarise(sds =sd(!!sym(rtvar),na.rm=T)) %>% mutate(mergekey=paste0(!!sym(subjvar),blockset))
  abds <- merge(setmeans,setsds,by="mergekey",suffixes=c("",".sd")) %>%
    mutate(ab=((pull0target1-pull1target1)-(pull0target0-pull1target0))/sds) %>%
    group_by(!!sym(subjvar)) %>% summarise(ab=mean(ab)) %>% dplyr::select(!!sym(subjvar),ab)
  return(abds)
}

#' @param formula A character string containing a formula to fit to the data and derive multilevel scores from
#' @param aatterm The random term, grouped under the subject variable, which represents the approach bias. Usually this is the interaction of the pull and target terms.
#' @export
#' @rdname Algorithms
aat_multilevelscore<-function(ds,subjvar,formula,aatterm,...){
  fit<- lme4::lmer(as.formula(formula),data=ds,control=
                     lme4::lmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 2e6),calc.derivs=F))
  output<-data.frame(ab=-lme4::ranef(fit,whichel=subjvar)[[subjvar]][[aatterm]])
  testset<<-fit
  output[[subjvar]]<-fit@flist[[subjvar]]%>%levels
  return(output)
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
#' \code{prune_nothing} excludes no trials, while \code{trial_prune_3SD} excludes trials deviating more than 3SD from the mean per participant.
#' @param errortrialfunc Function (without brackets or quotes) to apply to an error trial.
#' \code{error_replace_blockmeanplus} replaces error trial reaction times with the block mean plus an arbitrary extra amount of time.
#' If used, the following additional arguments are required:
#' \itemize{
#' \item \code{blockvar} - Quoted name of the block variable
#' \item \code{errorvar} - Quoted name of the error variable, where errors are 1 or TRUE and correct trials are 0 or FALSE
#' \item \code{errorbonus} - Amount to add to the reaction time of error trials. Default is 0.6 (recommended by \code{Greenwald, Nosek, & Banaji, 2003})
#' }
#' @param ... Other arguments, to be passed on to the algorithm functions (see \code{algorithm} above)
#'
#' @return A list, containing bootstrapped bias scores, a data frame with bootstrapped 95% confidence intervals,
#' the number of iterations, and a matrix of bias scores for each iteration.
#'
#' @author Sercan Kahveci
#' @examples
#' @export
aat_bootstrap<-function(ds,subjvar,pullvar,targetvar,rtvar,iters,plot=T,
                        algorithm=c("aat_doublemeandiff","aat_doublemediandiff","aat_dscore",
                                    "aat_dscore_multiblock","aat_multilevelscore"),
                        trialdropfunc=c("prune_nothing","trial_prune_3SD"),
                        errortrialfunc=c("prune_nothing","error_replace_blockmeanplus"),
                        ...){
  packs<-c("magrittr","dplyr","AATtools")

  #Handle arguments
  args<-list(...)
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
  if(algorithm=="aat_multilevelscore"){
    packs<-c(packs,"lme4")
    if(!any(c("formula","aatterm") %in% names(args))){
      args$formula<-paste0(rtvar,"~",1,"+(",pullvar,"*",targetvar,"|",subjvar,")")
      args$aatterm<-paste0(pullvar,":",targetvar)
      warning("No multilevel formula or AAT-term provided. Defaulting to formula ",
              args$formula," and AAT-term ",args$aatterm)
    }
  }
  ds %<>% aat_preparedata(subjvar,pullvar,targetvar,rtvar) %>% mutate(key=1)

  #bootstrap loop
  cluster<-makeCluster(detectCores()-1)
  registerDoParallel(cluster)
  results<-
    foreach(iter = seq_len(iters), .packages=c("magrittr","dplyr","lme4"), .combine=cbind) %dopar% {
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
                      bias=rowMeans(results),
                      lowerci=apply(results,MARGIN=1,FUN=function(x){quantile(x,0.025)}),
                      upperci=apply(results,MARGIN=1,FUN=function(x){quantile(x,0.975)}))
  statset$ci<-statset$upperci-statset$lowerci

  output<-list(bias=statset,iters=iters,iterdata=results) %>%
    structure(class = "aat_bootstrap")
  if(plot){ plot(output) }
  return(output)
}

plot.aat_bootstrap <- function(x){
  statset<-x$bias
  rank<-order(statset$bias)
  wideness<-max(statset$upperci) - min(statset$lowerci)
  plot(x=statset$bias,y=rank,xlim=c(min(statset$lowerci)-0.01*wideness,max(statset$upperci)+0.01*wideness),
       xlab="Bias score",main=paste0("Individual bias scores with 95%CI",
                                     "\nMean confidence interval: ",round(mean(statset$ci),digits=2)))
  segments(x0=statset$lowerci,x1=statset$bias-0.005*wideness,y0=rank,y1=rank)
  segments(x0=statset$bias+0.005*wideness,x1=statset$upperci,y0=rank,y1=rank)
  #text(x=statset$bias,y=statset$rownr,labels=statset$ppidx,cex=0.5)
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



aat_preparedata<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  args<-list(...)
  stopifnot(all(c(subjvar,pullvar,targetvar,rtvar) %in% colnames(ds)))
  ds[,subjvar]%<>%as.factor()
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
  rmindices<- which(is.na(ds[[subjvar]]) & is.na(ds[[pullvar]]) & is.na(ds[[targetvar]]) & is.na(ds[[rtvar]]))
  if(length(rmindices)>0){
    ds<-ds[-rmindices,]
    warning("Removed ",length(rmindices)," rows due to presence of NA in critical variable(s)")
  }
  ds<-ds[,c(subjvar,pullvar,targetvar,rtvar,args$errorvar,args$blockvar)]
  return(ds)
}

val_between<-function(x,lb,ub){x>lb & x<ub}

drop_empty_cases<-function(iterds,subjvar){
  iterds%>%group_by(!!sym(subjvar))%>%filter(val_between(mean(key),0,1))
}




