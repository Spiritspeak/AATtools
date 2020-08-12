
<!-- badges: start -->
![CRAN status](https://www.r-pkg.org/badges/version/AATtools)
[![Travis build status](https://travis-ci.org/Spiritspeak/AATtools.svg?branch=master)](https://travis-ci.org/Spiritspeak/AATtools)
<!-- badges: end -->

# Introduction

AATtools provides tools to deal with data from implicit psychological tasks. It provides methods to compute reliability and aggregate data into bias scores. 

Install it this way:
```{r install, eval=FALSE}
install.packages("AATtools")
```

# Rationale

Reliability scores are typically not computed for psychological tasks that produce a sum score. This has led to a literature full of inconsistent results and methodological decisions that have been guided by intuition rather than empirical decisionmaking. AATtools tries to solve these problems by providing multiple methods of computing the reliability of the Approach-Avoidance task as well as other implicit psychological tasks. 
Importantly, it enables researchers to compute the reliability of their entire data processing pipeline, factoring in the influence of decisions to remove or keep outliers in the final reliability score. This gives the researcher a clear overview over how reliable the data are that have actually been used in the study's analyses and enables them to explore the best ways to deal with non-normality, outliers and error trials. 

# Getting Started

## Getting your data in the right format

AATtools works with long-format `data.frames` that follow a specific format. Your `data.frame` should contain one row per trial, and it should contain variables that designate the relevant conditions (approach/avoidance, control/target stimuli) with 0's and 1's. 

## Compute reliability

If your data is in the right format, you can get started right away.

```
dataset <- erotica[erotica$is_irrelevant==0,]
split <- aat_splithalf(ds=dataset, #our dataset
                       subjvar="subject", #the column designating participant IDs
                       pullvar="is_pull", #the column designating approach (1) and avoidance (0) trials
                       targetvar="is_target", #the column designating target (1) and control (0) stimulus trials
                       rtvar="RT", #the column designating reaction times
                       iters=10000, #Set the number of bootstrapping iterations (more is better)
                       trialdropfunc="trial_prune_3SD", #Indicate whether outliers should be removed, and if so, how
                       casedropfunc="case_prune_3SD", #Indicate whether outlying approach bias scores should be removed
                       algorithm="aat_dscore" #The algorithm by which approach bias scores should be computed
                       )
                       
```


