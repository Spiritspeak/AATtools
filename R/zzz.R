#' @import dplyr
#' @import magrittr
#' @import doParallel
#' @import foreach
#' @importFrom magrittr %>% %<>% %$%
#' @importFrom dplyr group_by ungroup mutate summarise sample_n n filter select
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach getDoParRegistered registerDoSEQ
#' @importFrom stats var median sd lm vcov terms as.formula coef cor cov setNames quantile pt
#' rnorm rgamma pnorm qnorm
#' @importFrom graphics abline points segments text plot par
.onLoad<-function(libname, pkgname){
  #avoid CRAN errors
  utils::globalVariables(c("abhalf0","abhalf1","ab","key","."),"AATtools")

  #register generic functions
  registerS3method("print",class="aat_splithalf",method=print.aat_splithalf)
  registerS3method("plot",class="aat_splithalf",method=plot.aat_splithalf)
  registerS3method("print",class="aat_bootstrap",method=print.aat_bootstrap)
  registerS3method("plot",class="aat_bootstrap",method=plot.aat_bootstrap)
  registerS3method("print",class="qreliability",method=print.qreliability)
  registerS3method("plot",class="qreliability",method=plot.qreliability)
  registerS3method("print",class="aat_alpha",method=print.aat_alpha)
  registerS3method("plot",class="aat_alpha",method=plot.aat_alpha)
  registerS3method("print",class="aat_alpha_jackknife",method=print.aat_alpha_jackknife)
  registerS3method("plot",class="aat_alpha_jackknife",method=plot.aat_alpha_jackknife)

  #set max number of cores to use
  if (r_check_limit_cores()) {
    num_workers <- 2L
  } else {
    num_workers <- max(parallel::detectCores()-1,1)
  }
  options(AATtools.workers=num_workers)

  #greet user
  #packageStartupMessage("Thank you for loading AATtools v0.0.1")
}
