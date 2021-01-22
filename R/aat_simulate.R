#' Simulate AAT datasets and predict parameters
#' \code{aat_simulate()} generates approach-avoidance task datasets.
#'
#' @param npps Number of participants
#' @param nstims Number of stimuli
#' @param stimreps Number of repetitions of each stimulus within each group
#' (i.e. within approach target, avoid target, approach control, avoid control)
#' @param meanrt Mean sample reaction time
#' @param meanrt_jitter Extent by which participants' mean RTs
#' deviate from mean sample RT.
#' @param sdrt Standard deviation of samplewide RTs,
#' ignoring effects of movement, stimulus, and approach bias.
#' In essence, this represents the amount of pure noise present in the data.
#' @param sdrt_jitter Extent by which standard deviations of individual participants' RTs
#' are larger or smaller than the samplewide SD.
#' @param pullfx size of the effect of approach-versus-avoidance, in milliseconds
#' @param pullfx_jitter Individual variation in the effect of approach-versus-avoidance
#' @param stimfx size of the effect of stimulus category, in milliseconds
#' @param stimfx_jitter Individual variation in the effect of stimulus category
#' @param biasfx Size of the approach bias effect, in milliseconds
#' @param biasfx_jitter Individual variation in the approach bias effect
#'
#' @return \code{aat_simulate()} returns a \code{data.frame} with the following columns:
#' subj (participant ID), stim (stimulus number), rep (stimulus repetition number),
#' is_pull (0 = avoid, 1 = approach), is_target (0 = control stimulus, 1 = target stimulus),
#' meanrt (participant's mean RT), sdrt (participant's residual standard deviation),
#' pullfx (participant approach-avoidance effect size in ms),
#' stimfx (participant stimulus category effect size in ms),
#' biasfx (participant approach bias effect size in ms),
#' and rt (trial reaction time).
#' Additionally, the data.frame has the attribute \code{population_reliability} which represents
#' the expected reliability of the data given the provided parameters.
#' @details Defaults of \code{aat_simulate()} are based on
#' Kahveci, Van Alebeek, Berking, & Blechert (2021).
#' @export
#'
#' @examples
#' ts<- aat_simulate(pullfx = 50, stimfx = 10, biasfx = 100)
#' mod<-lm(rt~is_pull*is_target,data=ts)
#' coef(mod) #these should be somewhat close to the provided coefficients
#' print(attr(ts,"population_reliabilty"))
#' print(q_reliability(ts,"subj",rt~is_pull*is_target,"is_pull:is_target"))
#' #these two should be not too far apart,
#' # and should converge when the process is repeated a bunch
#'
#'
#' # Here's how to derive the parameters used in this function from a real dataset
#' \dontrun{
#' mod<-lmer(decisiontime ~ is_pull * is_food + (is_pull * is_food | subjectid),data=dsa)
#' fixef(mod) # from here, all the fx and mean RTs are derived
#' ranef(mod)$subjectid %>% apply(2,sd) #from here, all the fx jitters are derived
#' dsa %>% group_by(subjectid) %>% summarise(sd=sd(resid)) %>%
#' summarise(m=mean(sd),s=sd(sd)) # from here, sdrt_jitter is derived
#' }
aat_simulate<-function(npps=40,nstims=32,stimreps=2,
                  meanrt=743,meanrt_jitter=66,
                  sdrt=133,sdrt_jitter=38,
                  pullfx=25,pullfx_jitter=40,
                  stimfx=10,stimfx_jitter=35,
                  biasfx=35,biasfx_jitter=75){

  #set properties
  subjprops<-data.frame(subj=1:npps,
                        meanrt=meanrt+meanrt_jitter*rnorm(npps),
                        sdrt=sdrt+sdrt_jitter*rnorm(npps),
                        pullfx=pullfx+pullfx_jitter*rnorm(npps),
                        stimfx=stimfx+stimfx_jitter*rnorm(npps),
                        biasfx=biasfx+biasfx_jitter*rnorm(npps))

  #initialize dataset
  ds<-expand.grid(subj=1:npps,stim=1:nstims,rep=1:stimreps,is_pull=0:1,is_target=0:1)
  ds<-merge(ds,subjprops,by="subj",all.x=T)

  #Generate RTs
  gshape<-3
  gscale<-1
  ds$rt<-(rgamma(n=nrow(ds),shape=gshape,scale=gscale)-gshape*gscale) *
    ds$sdrt/(sqrt(gshape)*gscale) + ds$meanrt +
    (ds$is_pull-.5)*ds$pullfx + (ds$is_target-.5)*ds$stimfx +
    (ds$is_pull*ds$is_target-.25)*ds$biasfx

  #compute true "population" reliability (Kahveci's Q)
  alt_q <- (biasfx_jitter^2)/(biasfx_jitter^2 + sdrt^2 /(nstims*stimreps) *4)
  attr(ds,"population_reliability")<-alt_q

  #output
  return(ds)
}


