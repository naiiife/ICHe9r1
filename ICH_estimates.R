## Treatment Effect Estimation for Time-to-Event Data under ICH E9(R1)
## Competing risks data structure: 
## Only one of the primary event or intercurrent event is observed
## Methodology: Inverse probability weighting (without variance output)

library(survival)

matchy <- function(yvec,xvec,newx,z){
  tm = c(0,xvec)
  s = c(z,yvec)
  res = NULL
  for (x in newx){
    ti = max(which(tm<=x))
    res = append(res,s[ti])
  }
  return(res)
}

matcht <- function(xvec,newx){
  tm = c(0,xvec)
  res = NULL
  for (x in newx){
    ti = max(which(tm<=x))
    res = append(res,ti)
  }
  return(res)
}

plot.inc <- function(fit,decrease=FALSE,ylim=c(0,1),legend=c('Treated','Controlled'),...){
  tm = fit$time
  cif1 = fit$cif1
  cif0 = fit$cif0
  x = 'topleft'
  ylab = 'Cumulative incidence'
  if (decrease==TRUE){
    cif1 = 1-cif1
    cif0 = 1-cif0
    x = 'bottomleft'
    ylab = 'Survival probability'
  }
  plot(tm,cif1,type='s',col='brown',lwd=2,
       xlab='Time',ylab=ylab,ylim=ylim,...)
  points(tm,cif0,type='s',col='darkcyan',lwd=2)
  legend(x,cex=0.8,col=c('brown','darkcyan'),lwd=c(2,2),
         legend=legend)
}

plot.boot <- function(fit,decrease=FALSE,conf.int=.95,ylim=c(-1,1),...){
  tm = fit$time
  dcif = fit$cif1 - fit$cif0
  ciu = dcif - qnorm((1-conf.int)/2)*fit$sdt
  cil = dcif + qnorm((1-conf.int)/2)*fit$sdt
  ylab = 'Difference in cumulative incidences'
  if (decrease==TRUE){
    dcif = -dcif
    ciu = -ciu
    cil = -cil
    ylab = 'Difference in survival probabilities'
  }
  plot(tm,dcif,type='s',ylim=ylim,xlab='Time',ylab=ylab,lwd=2,...)
  abline(h=0,lty=2)
  points(tm,ciu,type='s',lty=5,lwd=1.5)
  points(tm,cil,type='s',lty=5,lwd=1.5)
}


surv.composite_gkm <- function(A,Time,cstatus,weights,tm=NULL,subset=NULL){
  if (!is.null(subset)){
    A = A[subset]
    Time = Time[subset]
    cstatus = cstatus[subset]
    weights = weights[subset]
  }
  fit.co1 = survfit(Surv(Time,cstatus!=0)~1, weights, subset=(A==1))
  fit.co0 = survfit(Surv(Time,cstatus!=0)~1, weights, subset=(A==0))
  if (is.null(tm)) tm = sort(unique(c(fit.co1$time,fit.co0$time)))
  cumhazco1 = matchy(fit.co1$cumhaz,fit.co1$time,tm,0)
  cumhazco0 = matchy(fit.co0$cumhaz,fit.co0$time,tm,0)
  dcif1 = exp(-cumhazco1)*diff(c(0,cumhazco1))
  dcif0 = exp(-cumhazco0)*diff(c(0,cumhazco0))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  return(list(time=tm,cif1=cif1,cif0=cif0))
}

surv.composite <- function(A,Time,cstatus,weights=NULL,
                           bootstrap=FALSE,nboot=100,seeds=0){
  N = length(A)
  if (is.null(weights)) weights = rep(1,N)
  result = surv.composite_gkm(A,Time,cstatus,weights,tm=NULL)
    if (bootstrap) {
      set.seed(seeds)
      tm = result$time
      ate = matrix(0,nboot,length(tm))
      for (b in 1:nboot){
        bt = surv.composite_gkm(A,Time,cstatus,weights,tm=tm,subset=sample(N,replace=TRUE))
        ate[b,] = bt$cif1-bt$cif0
      }
      sdt = apply(ate,2,sd)
      result = c(result,list(sdt=sdt))
    }
  return(result)
}


surv.natural_gkm <- function(A,Time,cstatus,weights,tm=NULL,subset=NULL){
  if (!is.null(subset)){
    A = A[subset]
    Time = Time[subset]
    cstatus = cstatus[subset]
    weights = weights[subset]
  }
  fit.11 = survfit(Surv(Time,cstatus==1)~1, weights, subset=(A==1))
  fit.10 = survfit(Surv(Time,cstatus==1)~1, weights, subset=(A==0))
  fit.21 = survfit(Surv(Time,cstatus==2)~1, weights, subset=(A==1))
  fit.20 = survfit(Surv(Time,cstatus==2)~1, weights, subset=(A==0))
  if (is.null(tm)) tm = sort(unique(c(fit.11$time,fit.10$time,fit.21$time,fit.20$time)))
  cumhaz11 = matchy(fit.11$cumhaz,fit.11$time,tm,0)
  cumhaz10 = matchy(fit.10$cumhaz,fit.10$time,tm,0)
  cumhaz21 = matchy(fit.21$cumhaz,fit.21$time,tm,0)
  cumhaz20 = matchy(fit.20$cumhaz,fit.20$time,tm,0)
  dcif1 = exp(-cumhaz11-cumhaz20)*diff(c(0,cumhaz11))
  dcif0 = exp(-cumhaz10-cumhaz20)*diff(c(0,cumhaz10))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  return(list(time=tm,cif1=cif1,cif0=cif0))
}

surv.natural <- function(A,Time,cstatus,weights=NULL,
                         bootstrap=FALSE,nboot=100,seeds=0){
  N = length(A)
  if (is.null(weights)) weights = rep(1,N)
  result = surv.natural_gkm(A,Time,cstatus,weights,tm=NULL)
    if (bootstrap) {
      set.seed(seeds)
      tm = result$time
      ate = matrix(0,nboot,length(tm))
      for (b in 1:nboot){
        bt = surv.natural_gkm(A,Time,cstatus,weights,tm=tm,subset=sample(N,replace=TRUE))
        ate[b,] = bt$cif1-bt$cif0
      }
      sdt = apply(ate,2,sd)
      result = c(result,list(sdt=sdt))
    }
  return(result)
}


surv.removed_gkm <- function(A,Time,cstatus,weights,tm=NULL,subset=NULL){
  if (!is.null(subset)){
    A = A[subset]
    Time = Time[subset]
    cstatus = cstatus[subset]
    weights = weights[subset]
  }
  fit.11 = survfit(Surv(Time,cstatus==1)~1, weights, subset=(A==1))
  fit.10 = survfit(Surv(Time,cstatus==1)~1, weights, subset=(A==0))
  if (is.null(tm)) tm = sort(unique(c(fit.11$time,fit.10$time)))
  cumhaz11 = matchy(fit.11$cumhaz,fit.11$time,tm,0)
  cumhaz10 = matchy(fit.10$cumhaz,fit.10$time,tm,0)
  cif1 = 1-exp(-cumhaz11)
  cif0 = 1-exp(-cumhaz10)
  return(list(time=tm,cif1=cif1,cif0=cif0))
}

surv.removed <- function(A,Time,cstatus,weights=NULL,
                         bootstrap=FALSE,nboot=100,seeds=0){
  N = length(A)
  result = surv.removed_gkm(A,Time,cstatus,weights,tm=NULL)
    if (bootstrap) {
      set.seed(seeds)
      tm = result$time
      ate = matrix(0,nboot,length(tm))
      for (b in 1:nboot){
        bt = surv.removed_gkm(A,Time,cstatus,weights,tm=tm,subset=sample(N,replace=TRUE))
        ate[b,] = bt$cif1-bt$cif0
      }
      sdt = apply(ate,2,sd)
      result = c(result,list(sdt=sdt))
    }
  return(result)
}


surv.whileon_gkm <- function(A,Time,cstatus,weights,tm=NULL,subset=NULL){
  if (!is.null(subset)){
    A = A[subset]
    Time = Time[subset]
    cstatus = cstatus[subset]
    weights = weights[subset]
  }
  fit.11 = survfit(Surv(Time,cstatus==1)~1, weights, subset=(A==1))
  fit.10 = survfit(Surv(Time,cstatus==1)~1, weights, subset=(A==0))
  fit.21 = survfit(Surv(Time,cstatus==2)~1, weights, subset=(A==1))
  fit.20 = survfit(Surv(Time,cstatus==2)~1, weights, subset=(A==0))
  if (is.null(tm)) tm = sort(unique(c(fit.11$time,fit.10$time,fit.21$time,fit.20$time)))
  cumhaz11 = matchy(fit.11$cumhaz,fit.11$time,tm,0)
  cumhaz10 = matchy(fit.10$cumhaz,fit.10$time,tm,0)
  cumhaz21 = matchy(fit.21$cumhaz,fit.21$time,tm,0)
  cumhaz20 = matchy(fit.20$cumhaz,fit.20$time,tm,0)
  dcif1 = exp(-cumhaz11-cumhaz21)*diff(c(0,cumhaz11))
  dcif0 = exp(-cumhaz10-cumhaz20)*diff(c(0,cumhaz10))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  return(list(time=tm,cif1=cif1,cif0=cif0))
}

surv.whileon <- function(A,Time,cstatus,weights=NULL,
                         bootstrap=FALSE,nboot=100,seeds=0){
  N = length(A)
  result = surv.whileon_gkm(A,Time,cstatus,weights,tm=NULL)
    if (bootstrap) {
      set.seed(seeds)
      tm = result$time
      ate = matrix(0,nboot,length(tm))
      for (b in 1:nboot){
        bt = surv.whileon_gkm(A,Time,cstatus,weights,tm=tm,subset=sample(N,replace=TRUE))
        ate[b,] = bt$cif1-bt$cif0
      }
      sdt = apply(ate,2,sd)
      result = c(result,list(sdt=sdt))
    }
  return(result)
}

surv.principal_gkm <- function(A,Time,cstatus,weights,principalscore='regression',
                               tm=NULL,subset=NULL){
  if (!is.null(subset)){
    A = A[subset]
    Time = Time[subset]
    cstatus = cstatus[subset]
    weights = weights[subset]
  }
  fit.11 = survfit(Surv(Time,cstatus==1)~1, weights, subset=(A==1))
  fit.10 = survfit(Surv(Time,cstatus==1)~1, weights, subset=(A==0))
  fit.21 = survfit(Surv(Time,cstatus==2)~1, weights, subset=(A==1))
  fit.20 = survfit(Surv(Time,cstatus==2)~1, weights, subset=(A==0))
  if (is.null(tm)) tm = sort(unique(c(fit.11$time,fit.10$time,fit.21$time,fit.20$time)))
  cumhaz11 = matchy(fit.11$cumhaz,fit.11$time,tm,0)
  cumhaz10 = matchy(fit.10$cumhaz,fit.10$time,tm,0)
  cumhaz21 = matchy(fit.21$cumhaz,fit.21$time,tm,0)
  cumhaz20 = matchy(fit.20$cumhaz,fit.20$time,tm,0)
  dcif1 = exp(-cumhaz11-cumhaz21)*diff(c(0,cumhaz11))
  dcif0 = exp(-cumhaz10-cumhaz20)*diff(c(0,cumhaz10))
  if (principalscore=='regression'){
    dcif1nps = exp(-cumhaz11-cumhaz21)*diff(c(0,cumhaz21))
    dcif0nps = exp(-cumhaz10-cumhaz20)*diff(c(0,cumhaz20))
    PI1 = 1 - max(cumsum(dcif1nps))
    PI0 = 1 - max(cumsum(dcif0nps))
  } else {
    fit.C1 = survfit(Surv(Time,cstatus==0)~1, weights, subset=(A==1))
    fit.C0 = survfit(Surv(Time,cstatus==0)~1, weights, subset=(A==0))
    SC1 = matchy(fit.C1$surv,fit.C1$time,Time,1)
    SC0 = matchy(fit.C0$surv,fit.C0$time,Time,1)
    SC1[SC1<=0] = 0.001
    SC0[SC0<=0] = 0.001
    PI1 = 1-mean(((cstatus==2)/SC1)[A==0&SC1>0.001])
    PI0 = 1-mean(((cstatus==2)/SC0)[A==1&SC0>0.001])
  }
  cif1 = cumsum(dcif1)/PI1
  cif0 = cumsum(dcif0)/PI0
  return(list(time=tm,cif1=cif1,cif0=cif0))
}

surv.principal <- function(A,Time,cstatus,weights=NULL,
                           principalscore='regression',
                           bootstrap=FALSE,nboot=100,seeds=0){
  N = length(A)
  if (is.null(weights)) weights = rep(1,N)
    result = surv.principal_gkm(A,Time,cstatus,weights,principalscore,tm=NULL)
    if (bootstrap) {
      set.seed(seeds)
      tm = result$time
      ate = matrix(0,nboot,length(tm))
      for (b in 1:nboot){
        bt = surv.principal_gkm(A,Time,cstatus,weights,principalscore,
                                tm=tm,subset=sample(N,replace=TRUE))
        ate[b,] = bt$cif1-bt$cif0
      }
      sdt = apply(ate,2,sd)
      result = c(result,list(sdt=sdt))
    }
  return(result)
}

treatmentpolicy <- function(A,T,dT,weights=NULL,tm=NULL,subset=NULL){
  if (!is.null(subset)){
    A = A[subset]
    T = T[subset]
    dT = dT[subset]
    weights = weights[subset]
  }
  fit.1 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=(A==1))
  fit.0 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=(A==0))
  if (is.null(tm)) tm = sort(unique(c(fit.1$time,fit.0$time)))
  cumhaz1 = matchy(fit.1$cumhaz,fit.1$time,tm,0)
  cumhaz0 = matchy(fit.0$cumhaz,fit.0$time,tm,0)
  cif1 = 1-exp(-cumhaz1)
  cif0 = 1-exp(-cumhaz0)
  return(list(time=tm,cif1=cif1,cif0=cif0))
}

surv.treatment <- function(A,Time,cstatus,weights=NULL,
                           bootstrap=FALSE,nboot=100,seeds=0){
  N = length(A)
  if (is.null(weights)) weights = rep(1,N)
  result = treatmentpolicy(A,Time,cstatus,weights)
    if (bootstrap) {
      set.seed(seeds)
      tm = result$time
      ate = matrix(0,nboot,length(tm))
      for (b in 1:nboot){
        bt = treatmentpolicy(A,Time,cstatus,weights,tm,subset=sample(N,replace=TRUE))
        ate[b,] = bt$cif1-bt$cif0
      }
      sdt = apply(ate,2,sd)
      result = c(result,list(sdt=sdt))
    }
  return(result)
}


### Example

setwd('D:/Data/白血病')
dat = read.csv('leukemiaPKU.csv')
A = 1-dat$TRANSPLANT
X = cbind(dat$MRD,dat$CR,dat$ALL)
cstatus = dat$RELAPSE + 2*dat$TRM
Time = dat$RELAPSET
## Input: 
# <A>, treatment
# <Time>, survival time (outcome or intercurrent event)
# <cstatus>, event type (1 if outcome, 2 if intercurrent event, 0 if censoring)
# Strategies:
# <treatment>, treatment policy strategy, which is also known as the intention-to-treat analysis
# <composite>, composite variable strategy, which combines intercurrent event and outcome as a single event
# <natural>, hypothetical strategy, which envisions a scenario where the hazards of intercurrent event are identical between groups
# <removed>, hypothetical strategy, which envisions a scenario where intercurrent events are removed
# <whileon>, while on strategy, or competing risk, which assumes outcome would never happen if intercurrent event happens
# <principal>, principal stratum strategy, which restricts to a subpopulation who never experience intercurrent events regardless of treatments
fit = surv.treatment(A,Time,cstatus,bootstrap=TRUE)
fit = surv.composite(A,Time,cstatus,bootstrap=TRUE)
fit = surv.natural(A,Time,cstatus,bootstrap=TRUE)
fit = surv.removed(A,Time,cstatus,bootstrap=TRUE)
fit = surv.whileon(A,Time,cstatus,bootstrap=TRUE)
fit = surv.principal(A,Time,cstatus,bootstrap=TRUE)
plot.inc(fit, decrease=FALSE)
plot.boot(fit)
