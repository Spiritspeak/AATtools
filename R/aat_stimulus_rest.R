
subtraction.matrix<-function(avec,bvec){
  na<-length(avec)
  nb<-length(bvec)
  out<-matrix(NA,nrow=na,ncol=nb)
  for(i in seq_len(na)){
    out[i,]<-avec[i]-bvec
  }
  return(out)
}

meanpercentile<-function(sample,population){
  sample %>% sapply(function(x) mean(x<population,na.rm=T)) %>% mean(na.rm=T)
}

#' Compute stimulus-rest correlations of double-difference scores
#' This function provides a statistic that can give an indication of how deviant
#' the responses to specific stimuli are, in comparison to the rest of the stimulus set.
#' The algorithm computes stimulus-rest correlations of stimulus-specific double-difference scores.
#' It takes single-difference approach-avoidance scores for each stimulus, and computes
#' every possible subtraction between individual stimuli from both stimulus categories.
#' It then computes correlations between every such subtraction of stimuli on one hand, and
#' the mean double difference score of all other stimuli. Stimulus-rest correlations are then
#' computed by averaging every such subtraction-rest correlation involving a specific stimulus.
#'
#' @param ds a \code{data.frame}
#' @param subjvar the label of the participant identifier variable
#' @param stimvar the label of the stimulus identifier variable
#' @param pullvar the label of the movement direction identifier variable
#' @param targetvar the label of the stimulus category identifier variable
#' @param rtvar the label of the reaction time variable
#'
#' @return Returns a \code{aat_stimulus_rest} object containing statistics for each stimulus.
#' Stats include the average stimulus-rest correlation (mcor); the standard deviation of
#' dyad-rest correlations for this stimulus (sdcor);
#' the number of valid correlations involved in these statistic (n);
#' the average percentile of dyad-rest correlations involving the stimulus within
#' the distribution of all other dyad-rest correlations (restpercentile);
#' as well as z-scores (zpercentile) and p-values for this percentile (pval).
#'
#' @export
#'
#' @examples
#'
#' ds<-aat_simulate()
#' stimrest<-aat_stimulus_rest(ds,subjvar="subj",stimvar="stim",pullvar="is_pull",
#'                      targetvar="is_target",rtvar="rt")
#' plot(stimrest)
#' print(stimrest)
aat_stimulus_rest<-function(ds,subjvar,stimvar,pullvar,targetvar,rtvar){
  # check data
  ds<-aat_preparedata(ds,subjvar,pullvar,targetvar,rtvar,stimvar=stimvar)

  #compute single-difference scores
  biasset<-ds%>%group_by(!!sym(subjvar),!!sym(stimvar),!!sym(targetvar))%>%
    summarise(bias=mean(subset(!!sym(rtvar),!!sym(pullvar)==0),na.rm=T)-
                   mean(subset(!!sym(rtvar),!!sym(pullvar)==1),na.rm=T))

  #compute dataset with all subtractions
  pp<-unique(biasset[[subjvar]])
  dl<-as.list(numeric(length(pp)))
  for(i in seq_along(pp)){
    mm<-subtraction.matrix(biasset[biasset[[subjvar]]==pp[i] & biasset[[targetvar]]==0,][["bias"]],
                           biasset[biasset[[subjvar]]==pp[i] & biasset[[targetvar]]==1,][["bias"]])
    cols<-unique(biasset[biasset[[subjvar]]==pp[i] & biasset[[targetvar]]==0,][[stimvar]])
    rows<-unique(biasset[biasset[[subjvar]]==pp[i] & biasset[[targetvar]]==1,][[stimvar]])

    dl[[i]]<-do.call(data.frame,
                list(as.vector(mm),rep(cols,each=nrow(mm)),rep(rows,n=ncol(mm)),F) %>%
                setNames(c(paste0("subject_",pp[i]),"first","second","stringsAsFactors")))
  }
  diffset<-Reduce(function(x,y){merge(x,y,by=c("first","second"),all=T,sort=F,
                                      no.dups=F)},x=dl[-1],init=dl[1])

  #compute correlations between specific subtraction and mean of all subtractions
  # not involving the stimuli involved in this particular subtraction
  diffmat<-as.matrix(diffset[,-1:-2])
  corvec<-numeric(nrow(diffset))
  counts<-numeric(nrow(diffset))
  for(i in seq_along(corvec)){
    relrows<-which(diffset$first[i] != diffset$first & diffset$second[i] != diffset$second)
    corvec[i]<-cor(diffmat[i,],cs<-colSums(diffmat[relrows,],na.rm=T),use="complete.obs")
    counts[i]<-sum(!is.na(diffmat[i,]) & !is.na(cs))
  }

  #compute stats based on correlations
  difflist<-data.frame(cor=corvec,
                       count=counts,
                       firstimg=diffset$first,
                       secondimg=diffset$second,
                       rownum=seq_along(corvec))
  stimstats<-rbind(difflist%>%rename(img=.data$firstimg)%>%group_by(.data$img)%>%
                   summarise(mcor=cormean(.data$cor,.data$count,na.rm=T),
                             sdcor=sd(.data$cor,na.rm=T),
                             n=sum(!is.na(cor)),
                             restpercentile=meanpercentile(.data$cor,corvec[-.data$rownum]),
                             lowerci=quantile(.data$cor,.025,na.rm=T),
                             upperci=quantile(.data$cor,.975,na.rm=T)),
                 difflist%>%rename(img=.data$secondimg)%>%group_by(.data$img)%>%
                   summarise(mcor=cormean(.data$cor,.data$count,na.rm=T),
                             sdcor=sd(.data$cor,na.rm=T),
                             n=sum(!is.na(.data$cor)),
                             restpercentile=meanpercentile(.data$cor,corvec[-.data$rownum]),
                             lowerci=quantile(.data$cor,.025,na.rm=T),
                             upperci=quantile(.data$cor,.975,na.rm=T)))
  stimstats$zpercentile<-qnorm(stimstats$restpercentile)
  stimstats$pval<-1-abs(.5-stimstats$restpercentile)*2
  # selfpercentile=mean(cor>mean(corvec))
  # stimstats$percentile<-stimstats$mcor %>% sapply(function(x) mean(x < corvec))
  # stimstats$compcorrz<-(r2z(stimstats$mcor)-r2z(mean(corvec))) / sqrt((1/(stimstats$n-3)) + (1/(length(corvec)-3)))
  # stimstats$compcorrp<-pnorm(abs(stimstats$compcorrz)/2,lower.tail=F)
  # stimstats$zcor<-(stimstats$mcor-mean(corvec))/sd(corvec) #full sample-perspective P-value
  # stimstats$zcor2<-(stimstats$mcor-mean(corvec))/stimstats$sdcor #stimulus-perspective P-value

  return(structure(list(stimstats=stimstats,difflist=difflist),class="aat_stimulus_rest"))
}

#' @rdname aat_stimulus_rest
#' @param x an \code{aat_stimulus_rest} object
#' @param ... Ignored.
#' @export
plot.aat_stimulus_rest<-function(x,...){
  statset<-x$stimstats
  statset<-statset[!is.na(statset$mcor) & !is.na(statset$upperci) & !is.na(statset$lowerci),]
  ranks<-rank(statset$mcor)
  wideness<-max(statset$upperci) - min(statset$lowerci)

  plot(x=statset$mcor,y=ranks,
       xlim=c(min(statset$lowerci)-0.01*wideness,max(statset$upperci)+0.01*wideness),
       xlab="Stimulus-rest correlation",main=paste0("Stimulus-rest correlations with 95%CI"),
       yaxt="n")
  segments(x0=statset$lowerci,x1=statset$mcor-0.005*wideness,y0=ranks,y1=ranks)
  segments(x0=statset$mcor+0.005*wideness,x1=statset$upperci,y0=ranks,y1=ranks)
  abline(v=mean(statset$mcor))

  axis(2, labels=statset$img,at=ranks,las=1,cex.axis=.5)
}

#' @rdname aat_stimulus_rest
#' @param x an \code{aat_stimulus_rest} object
#' @param ... Ignored.
#' @export
print.aat_stimulus_rest<-function(x,...){
  print(x$stimstats[order(x$stimstats$restpercentile),
                    c("img","mcor","restpercentile","zpercentile","pval")],n=nrow(x$stimstats))
}
