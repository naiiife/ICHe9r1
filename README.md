# ICHe9r1
Inference for Cumulative Incidences and Treatment Effects in Randomized Controlled Trials with Time-to-Event Outcomes

## Estimate the counterfactual cumulative incidences
    fit = surv.ICH(A,Time,cstatus,strategy='composite',cov1=NULL,weights=NULL,subset=NULL)
        A: treatment
## Plot the estimated countarfactual cumulative incidences
    plot.inc(fit,decrease=FALSE,conf.int=.95,nboot=0,seed=0,xlab='Time',xlim=NULL,ylim=c(0,1),legend=c('Treated','Controlled'),cex=0.8,...)

## Plot the estimated treatment effect
    plot.ate(fit,nboot=0,seed=0,decrease=FALSE,conf.int=.95,xlab='Time',xlim=NULL,ylim=c(-1,1),...)
