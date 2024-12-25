# ICHe9r1
Inference for Cumulative Incidences and Treatment Effects in Randomized Controlled Trials with Time-to-Event Outcomes

## Estimate the counterfactual cumulative incidences
    fit = surv.ICH(A,Time,cstatus,strategy='composite',cov1=NULL,method='np',weights=NULL,subset=NULL)
    # A: treatment
    # Time: time to event
    # cstatus: event indicator, 1 for primary event, 2 for intercurrent event, 0 for censoring
    # strategy: analysis strategy, "treatment", "composite", "natural", "removed", "whileon", "principal"
    # cov1: baseline covariates to adjust (for propensity score weighting)
    # method: estimation method, nonparametric (np) or efficient influence function based
    # weights: weights for each subject (applicable for nonparametric estimation)
    # subset: subset of individuals used for analysis

## Plot the estimated countarfactual cumulative incidences
    plot.inc(fit,decrease=FALSE,conf.int=.95,nboot=0,seed=0,xlab='Time',xlim=NULL,ylim=c(0,1),
             legend=c('Treated','Controlled'),cex=0.9,...)
    # fit: model
    # decrease: FALSE for cumulative incidences, TRUE for survival probabilities
    # conf.int: level of confidence interval
    # nboot: resampling times for bootstrap to construct confidence interval

## Plot the estimated treatment effect
    plot.ate(fit,decrease=FALSE,conf.int=.95,nboot=0,seed=0,xlab='Time',xlim=NULL,ylim=c(-1,1),...)
    # fit: model
    # decrease: FALSE for diff in cumulative incidences, TRUE for diff in survival probabilities
    # conf.int: level of confidence interval
    # nboot: resampling times for bootstrap to construct confidence interval
