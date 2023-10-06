## Causal inference for time-to-event data with intercurrent events
## main analysis function: surv.ICH(A,Time,cstatus,strategy,weights,subset)
## A: treatment
## Time: Event time
## cstatus: Event type, 1 for primary outcome, 2 for intercurrent event, 0 for censoring
## cov1: baseline covariates to adjust
## weights: weights (default by NULL)
## subset: which subset of data to use, logical (default by ALL)
## draw cumulative incidence curv: plot.inc
## Calculate confidence interval: surv.boot(fit,nboot)
## nboot: runs of bootstrap (default by 0, which gives analytical form)
## Draw treatment effect: plot.boot(fit.b)


library(survival)

matchy <- function(yvec,xvec,newx){
  ivec = sapply(newx, function(x) max(which(xvec<=x)))
  if (is.vector(yvec)) {
    return(yvec[ivec])
  } else {
    return(yvec[ivec,])
  }
}

ipscore <- function(A,covA){
  fps = glm(A~covA,family='binomial')
  ps = predict(fps,type='response')
  ips = A/ps + (1-A)/(1-ps)
  return(ips)
}

surv.treatment <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = rep(TRUE,length(A))
  fit1 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==1))
  fit0 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==0))
  cif1 = c(0, 1 - exp(-fit1$cumhaz))
  cif0 = c(0, 1 - exp(-fit0$cumhaz))
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  return(list(time1=c(0,fit1$time),time0=c(0,fit0$time),cif1=cif1,cif0=cif0,
              se1=se1,se0=se0))
}

surv.composite <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = rep(TRUE,length(A))
  fit1 = survfit(Surv(Time,cstatus!=0)~1, weights=weights, subset=subset&(A==1))
  fit0 = survfit(Surv(Time,cstatus!=0)~1, weights=weights, subset=subset&(A==0))
  cif1 = c(0, 1 - exp(-fit1$cumhaz))
  cif0 = c(0, 1 - exp(-fit0$cumhaz))
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  return(list(time1=c(0,fit1$time),time0=c(0,fit0$time),cif1=cif1,cif0=cif0,
              se1=se1,se0=se0))
}

surv.removed <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = rep(TRUE,length(A))
  fit1 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==1))
  fit0 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==0))
  cif1 = c(0, 1 - exp(-fit1$cumhaz))
  cif0 = c(0, 1 - exp(-fit0$cumhaz))
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  return(list(time1=c(0,fit1$time),time0=c(0,fit0$time),cif1=cif1,cif0=cif0,
              se1=se1,se0=se0))
}

surv.natural <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = rep(TRUE,length(A))
  fit11 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==1))
  fit10 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==0))
  fit21 = survfit(Surv(Time,cstatus==2)~1, weights=weights, subset=subset&(A==1))
  fit20 = survfit(Surv(Time,cstatus==2)~1, weights=weights, subset=subset&(A==0))
  Time = c(0, unique(sort(c(fit11$time,fit21$time,fit10$time,fit20$time))))
  fit11 = matchy(rbind(0,cbind(fit11$cumhaz,fit11$std.err)),c(0,fit11$time),Time)
  fit10 = matchy(rbind(0,cbind(fit10$cumhaz,fit10$std.err)),c(0,fit10$time),Time)
  fit21 = matchy(rbind(0,cbind(fit21$cumhaz,fit21$std.err)),c(0,fit21$time),Time)
  fit20 = matchy(rbind(0,cbind(fit20$cumhaz,fit20$std.err)),c(0,fit20$time),Time)
  dcif1 = exp(-fit11[,1]-fit20[,1])*diff(c(0,fit11[,1]))
  dcif0 = exp(-fit10[,1]-fit20[,1])*diff(c(0,fit10[,1]))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  M1 = diff(c(0,fit11[,2]^2))
  M0 = diff(c(0,fit20[,2]^2))
  G11 = cumsum(M1)*cif1^2 + cumsum(M1*(exp(-fit11[,1]-fit20[,1])+cif1)^2) -
    2*cif1*cumsum(M1*(exp(-fit11[,1]-fit20[,1])+cif1))
  G01 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  se1 = sqrt(G11+G01)
  M1 = diff(c(0,fit10[,2]^2))
  G10 = cumsum(M1)*cif0^2 + cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0)^2) -
    2*cif0*cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0))
  G00 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  se0 = sqrt(G10+G00)
  dcif = cif1-cif0
  G0 = cumsum(M0)*dcif^2 + cumsum(M0*dcif^2) - 2*dcif*cumsum(M0*dcif)
  se = sqrt(G11+G10+G0)
  return(list(time1=Time,time0=Time,cif1=cif1,cif0=cif0,se1=se1,se0=se0,se=se))
}

surv.whileon <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = rep(TRUE,length(A))
  fit11 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==1))
  fit10 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==0))
  fit21 = survfit(Surv(Time,cstatus==2)~1, weights=weights, subset=subset&(A==1))
  fit20 = survfit(Surv(Time,cstatus==2)~1, weights=weights, subset=subset&(A==0))
  time1 = c(0, fit11$time)
  time0 = c(0, fit10$time)
  fit11 = rbind(0,cbind(fit11$cumhaz,fit11$std.err))
  fit10 = rbind(0,cbind(fit10$cumhaz,fit10$std.err))
  fit21 = rbind(0,cbind(fit21$cumhaz,fit21$std.err))
  fit20 = rbind(0,cbind(fit20$cumhaz,fit20$std.err))
  dcif1 = exp(-fit11[,1]-fit21[,1])*diff(c(0,fit11[,1]))
  dcif0 = exp(-fit10[,1]-fit20[,1])*diff(c(0,fit10[,1]))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  M1 = diff(c(0,fit11[,2]^2))
  M0 = diff(c(0,fit21[,2]^2))
  G1 = cumsum(M1)*cif1^2 + cumsum(M1*(exp(-fit11[,1]-fit21[,1])+cif1)^2) -
    2*cif1*cumsum(M1*(exp(-fit11[,1]-fit21[,1])+cif1))
  G0 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  se1 = sqrt(G1+G0)
  M1 = diff(c(0,fit10[,2]^2))
  M0 = diff(c(0,fit20[,2]^2))
  G1 = cumsum(M1)*cif0^2 + cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0)^2) -
    2*cif0*cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0))
  G0 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  se0 = sqrt(G1+G0)
  return(list(time1=time1,time0=time0,cif1=cif1,cif0=cif0,se1=se1,se0=se0))
}

surv.principal <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = rep(TRUE,length(A))
  fit11 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==1))
  fit10 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset&(A==0))
  fit21 = survfit(Surv(Time,cstatus==2)~1, weights=weights, subset=subset&(A==1))
  fit20 = survfit(Surv(Time,cstatus==2)~1, weights=weights, subset=subset&(A==0))
  time1 = c(0, fit11$time)
  time0 = c(0, fit10$time)
  fit11 = rbind(0,cbind(fit11$cumhaz,fit11$std.err))
  fit10 = rbind(0,cbind(fit10$cumhaz,fit10$std.err))
  fit21 = rbind(0,cbind(fit21$cumhaz,fit21$std.err))
  fit20 = rbind(0,cbind(fit20$cumhaz,fit20$std.err))
  dcif1 = exp(-fit11[,1]-fit21[,1])*diff(c(0,fit11[,1]))
  dcif0 = exp(-fit10[,1]-fit20[,1])*diff(c(0,fit10[,1]))
  dcif1 = exp(-fit11[,1]-fit21[,1])*diff(c(0,fit11[,1]))
  dcif0 = exp(-fit10[,1]-fit20[,1])*diff(c(0,fit10[,1]))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  PR1 = 1 - sum(exp(-fit11[,1]-fit21[,1])*diff(c(0,fit21[,1])))
  PR0 = 1 - sum(exp(-fit10[,1]-fit20[,1])*diff(c(0,fit20[,1])))
  M1 = diff(c(0,fit11[,2]^2))
  M0 = diff(c(0,fit21[,2]^2))
  G1 = cumsum(M1)*cif1^2 + cumsum(M1*(exp(-fit11[,1]-fit21[,1])+cif1)^2) -
    2*cif1*cumsum(M1*(exp(-fit11[,1]-fit21[,1])+cif1))
  G0 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  se1 = sqrt(G1+G0)/PR1
  M1 = diff(c(0,fit10[,2]^2))
  M0 = diff(c(0,fit20[,2]^2))
  G1 = cumsum(M1)*cif0^2 + cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0)^2) -
    2*cif0*cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0))
  G0 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  se0 = sqrt(G1+G0)/PR0
  return(list(time1=time1,time0=time0,cif1=cif1/PR1,cif0=cif0/PR0,se1=se1,se0=se0))
}

surv.ICH <- function(A,Time,cstatus,strategy='composite',cov1=NULL,
                     weights=NULL,subset=NULL){
  N = length(A)
  if (is.null(weights)) weights = rep(1,N)
  if (!is.null(cov1)) weights = weights*ipscore(A,cov1)
  if (strategy=='treatment') fit = surv.treatment(A,Time,cstatus,weights,subset)
  if (strategy=='composite') fit = surv.composite(A,Time,cstatus,weights,subset)
  if (strategy=='natural') fit = surv.natural(A,Time,cstatus,weights,subset)
  if (strategy=='removed') fit = surv.removed(A,Time,cstatus,weights,subset)
  if (strategy=='whileon') fit = surv.whileon(A,Time,cstatus,weights,subset)
  if (strategy=='principal') fit = surv.principal(A,Time,cstatus,weights,subset)
  return(c(fit,list(A=A,Time=Time,cstatus=cstatus,strategy=strategy,cov1=cov1,
                    weights=weights,subset=subset)))
}

surv.boot <- function(fit,nboot=0,seed=0){
  N = length(fit$A)
  Time = sort(unique(c(0,fit$Time[fit$cstatus>0],max(fit$Time))))
  cif1 = matchy(fit$cif1,fit$time1,Time)
  cif0 = matchy(fit$cif0,fit$time0,Time)
  se1 = matchy(fit$se1,fit$time1,Time)
  se0 = matchy(fit$se0,fit$time0,Time)
  ate = cif1-cif0
  se = sqrt(se1^2+se0^2)
  if (fit$strategy=='natural') se = matchy(fit$se,fit$time1,Time)
  if (nboot>1){
  cif1l = cif0l = te = NULL
  for(b in 1:nboot){
    wt = as.vector(rmultinom(1,N,rep(1/N,N)))
    fitb = surv.ICH(fit$A,fit$Time,fit$cstatus,fit$strategy,fit$cov1,fit$weights*wt,fit$subset)
    cifb1 = matchy(fitb$cif1,fitb$time1,Time)
    cifb0 = matchy(fitb$cif0,fitb$time0,Time)
    cif1l = rbind(cif1l, cifb1)
    cif0l = rbind(cif0l, cifb0)
    te = rbind(te, cifb1-cifb0)
  }
  se1 = apply(cif1l,2,sd)
  se0 = apply(cif0l,2,sd)
  se = apply(te,2,sd)
  }
  return(list(Time=Time,ate=ate,cif1=cif1,cif0=cif0,se1=se1,se0=se0,se=se))
}

plot.inc <- function(fit,decrease=FALSE,conf.int=.95,nboot=0,seed=0,xlab='Time',xlim=NULL,
                     ylim=c(0,1),legend=c('Treated','Controlled'),cex=0.8,...){
  if (decrease==TRUE){
    cif1 = 1-fit$cif1
    cif0 = 1-fit$cif0
    x = 'bottomleft'
    ylab = 'Survival probability'
  } else {
    cif1 = fit$cif1
    cif0 = fit$cif0
    x = 'topleft'
    ylab = 'Cumulative incidence'
  }
  if (!is.null(xlim)){
  i1 = fit$time1<=xlim[2]
  i0 = fit$time0<=xlim[2]
  t1 = fit$time1[i1]
  t0 = fit$time0[i0]
  cif1 = cif1[i1]
  cif0 = cif0[i0]
  } else {
  t1 = fit$time1
  t0 = fit$time0
  i1 = rep(TRUE,length(t1))
  i0 = rep(TRUE,length(t0))
  }
  plot(t1,cif1,type='s',col='brown',lwd=2,
       xlab=xlab,ylab=ylab,ylim=ylim,...)
  points(t0,cif0,type='s',col='darkcyan',lwd=2)
  if (!is.null(conf.int)){
    fit.b = surv.boot(fit,nboot=0,seed=0)
    z = -qnorm((1-conf.int)/2)
    points(t1,cif1+fit.b$se1[i1]*z,type='s',lty=2,lwd=1.5,col='brown')
    points(t1,cif1-fit.b$se1[i1]*z,type='s',lty=2,lwd=1.5,col='brown')
    points(t0,cif0+fit.b$se0[i0]*z,type='s',lty=2,lwd=1.5,col='darkcyan')
    points(t0,cif0-fit.b$se0[i0]*z,type='s',lty=2,lwd=1.5,col='darkcyan')
  }
  legend(x,cex=cex,col=c('brown','darkcyan'),lwd=c(2,2),legend=legend)
}

plot.ate <- function(fit,nboot=0,seed=0,decrease=FALSE,conf.int=.95,xlab='Time',
                      xlim=NULL,ylim=c(-1,1),...){
  fit.b = surv.boot(fit,nboot=nboot,seed=seed)
  tm = fit.b$Time
  dcif = fit.b$ate
  se = fit.b$se
  ciu = dcif - qnorm((1-conf.int)/2)*se
  cil = dcif + qnorm((1-conf.int)/2)*se
  ylab = 'Diff in cumulative incidences'
  if (decrease==TRUE){
    dcif = -dcif
    ciu = -ciu
    cil = -cil
    ylab = 'Diff in survival probabilities'
  }
  if (!is.null(xlim)) {
    id = tm<=xlim[2]
    tm = tm[id]
    dcif = dcif[id]
    ciu = ciu[id]
    cil = cil[id]
  }
  plot(tm,dcif,type='s',ylim=ylim,xlab=xlab,ylab=ylab,lwd=2,...)
  abline(h=0,lty=2)
  points(tm,ciu,type='s',lty=5,lwd=1.5,col='brown')
  points(tm,cil,type='s',lty=5,lwd=1.5,col='brown')
}
