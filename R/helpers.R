FlanaganRulonBilateral<-function(x1,x2){
  key<-!is.na(x1) & !is.na(x2)
  x1<-x1[key]
  x2<-x2[key]
  fr<-(1-var(x1-x2)/var(x1+x2))
  return(fr/max(1, 1-fr))
}

RajuBilateral<-function(x1,x2,prop){
  covar<-cov(x1,x2)
  sumvar<-var(x1)+var(x2)+2*abs(covar)
  raju<-covar / (prop * (1-prop) * sumvar)
  return(raju)
}

val_between<-function(x,lb,ub){x>lb & x<ub}

drop_empty_cases<-function(iterds,subjvar){
  iterds%>%group_by(!!sym(subjvar))%>%filter(val_between(mean(key),0,1))
}

form2char<-function(x){
  if(is.character(x)){ return(x) }
  fs<-as.character(x)
  fs<-paste(fs[2],fs[1],fs[3])
  return(fs)
}

is.formula <- function(x){
  inherits(x,"formula")
}

mf<-function(x,digits=2){
  s<-format(round(x,digits=digits),
            digits=digits,drop0trailing=T,scientific=F,nsmall=digits)
  s<-gsub("^0\\.","\\.",s)
  return(s)
}

r_check_limit_cores <- function() {
  Rcheck <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  return((nchar(Rcheck[1]) > 0) & (Rcheck != "false"))
}

unregisterDoParallel <- function(cluster) {
  stopCluster(cluster)
  registerDoSEQ()
  #env <- foreach:::.foreachGlobals
  #rm(list=ls(name=env), pos=env)
}

aat_preparedata<-function(ds,subjvar,pullvar,targetvar=NULL,rtvar,...){
  args<-list(...)

  cols<-c(subjvar,pullvar,targetvar,rtvar,args$errorvar,args$blockvar,args$stimvar)
  if("formula" %in% names(args)){
    formterms <- args$formula %>% as.formula() %>% terms() %>%
      attr("variables") %>% as.character()
    formterms <- formterms[-1]
    if(any(!(formterms %in% colnames(ds)))){
      stop("Formula term(s) ",paste(formterms[!(formterms %in% colnames(ds))],collapse=", ")," missing from dataset")
    }
    cols <- c(cols,formterms)
  }
  stopifnot(all(cols %in% colnames(ds)))
  ds<-ds[,cols]

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
              levels(ds[,targetvar])[1], " represents control/neutral stimuli and ",
              levels(ds[,targetvar])[2], " represents target stimuli")
      ds[,targetvar]<-as.numeric(ds[,targetvar])-1
    }
  }

  rmindices <- ds[,cols] %>% lapply(FUN=is.na) %>% as.data.frame %>%
    apply(MARGIN=1,FUN=any) %>% which

  if(length(rmindices)>0){
    ds<-ds[-rmindices,]
    warning("Removed ",length(rmindices),
            " rows due to presence of NA in critical variable(s)")
  }
  return(ds)
}
