## Causal inference for time-to-event data with intercurrent events
## surv.ICH(A,Time,cstatus,strategy,cov1,weights,subset): Calculate cumulative incidences
#### A: treatment
#### Time: Event time
#### cstatus: Event type, 1 for primary outcome, 2 for intercurrent event, 0 for censoring
#### strategy: which analysis strategy is used
#### cov1: baseline covariates to adjust
#### method: estimation method, np for IPW, otherwise EIF based
#### weights: weights (default by NULL)
#### subset: which subset of data to use, logical (default by using all)
## surv.boot(fit): Calculate treatment effect and confidence interval
#### nboot: 0 using analytical form, positive integer using bootstrap
## plot.inc(fit): Plot cumulative incidence curve with confidence interval
## plot.ate(fit): Plot treatment effect with confidence interval


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
  if (is.null(subset)) subset = 1:length(A)
  fit1 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==1])
  fit0 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==0])
  cif1 = c(0, 1 - exp(-fit1$cumhaz))
  cif0 = c(0, 1 - exp(-fit0$cumhaz))
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  surv_diff = survdiff(Surv(Time,cstatus==1)~A,subset=subset)
  p = 1 - pchisq(surv_diff$chisq, length(surv_diff$n)-1)
  tt1 = c(0,fit1$time)
  tt0 = c(0,fit0$time)
  tt = sort(unique(c(tt1,tt0)))
  ate = matchy(cif1,tt1,tt)-matchy(cif0,tt0,tt)
  se = sqrt(matchy(se1,tt1,tt)^2+matchy(se0,tt0,tt)^2)
  return(list(time1=tt1,time0=tt0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.composite <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  fit1 = survfit(Surv(Time,cstatus>0)~1, weights=weights, subset=subset[A[subset]==1])
  fit0 = survfit(Surv(Time,cstatus>0)~1, weights=weights, subset=subset[A[subset]==0])
  cif1 = c(0, 1 - exp(-fit1$cumhaz))
  cif0 = c(0, 1 - exp(-fit0$cumhaz))
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  surv_diff = survdiff(Surv(Time,cstatus>0)~A,subset=subset)
  p = 1 - pchisq(surv_diff$chisq, length(surv_diff$n)-1)
  tt1 = c(0,fit1$time)
  tt0 = c(0,fit0$time)
  tt = sort(unique(c(tt1,tt0)))
  ate = matchy(cif1,tt1,tt)-matchy(cif0,tt0,tt)
  se = sqrt(matchy(se1,tt1,tt)^2+matchy(se0,tt0,tt)^2)
  return(list(time1=tt1,time0=tt0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.removed <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  fit1 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==1])
  fit0 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==0])
  cif1 = c(0, 1 - exp(-fit1$cumhaz))
  cif0 = c(0, 1 - exp(-fit0$cumhaz))
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  surv_diff = survdiff(Surv(Time,cstatus==1)~A,subset=subset)
  p = 1 - pchisq(surv_diff$chisq, length(surv_diff$n)-1)
  tt1 = c(0,fit1$time)
  tt0 = c(0,fit0$time)
  tt = sort(unique(c(tt1,tt0)))
  ate = matchy(cif1,tt1,tt)-matchy(cif0,tt0,tt)
  se = sqrt(matchy(se1,tt1,tt)^2+matchy(se0,tt0,tt)^2)
  return(list(time1=tt1,time0=tt0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.natural <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  fit11 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==1])
  fit10 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==0])
  fit21 = survfit(Surv(Time,cstatus>1)~1, weights=weights, subset=subset[A[subset]==1])
  fit20 = survfit(Surv(Time,cstatus>1)~1, weights=weights, subset=subset[A[subset]==0])
  time = c(0, unique(sort(c(fit11$time,fit21$time,fit10$time,fit20$time))))
  fit11 = matchy(rbind(0,cbind(fit11$cumhaz,fit11$std.err)),c(0,fit11$time),time)
  fit10 = matchy(rbind(0,cbind(fit10$cumhaz,fit10$std.err)),c(0,fit10$time),time)
  fit21 = matchy(rbind(0,cbind(fit21$cumhaz,fit21$std.err)),c(0,fit21$time),time)
  fit20 = matchy(rbind(0,cbind(fit20$cumhaz,fit20$std.err)),c(0,fit20$time),time)
  dcif1 = exp(-fit11[,1]-fit20[,1])*diff(c(0,fit11[,1]))
  dcif0 = exp(-fit10[,1]-fit20[,1])*diff(c(0,fit10[,1]))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  V1 = fit11[,2]^2
  V0 = fit20[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  G11 = cumsum(M1)*cif1^2 + cumsum(M1*(exp(-fit11[,1]-fit20[,1])+cif1)^2) -
    2*cif1*cumsum(M1*(exp(-fit11[,1]-fit20[,1])+cif1))
  G01 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  se1 = sqrt(G11+G01)
  V1 = fit10[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  M1 = diff(c(0,V1))
  G10 = cumsum(M1)*cif0^2 + cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0)^2) -
    2*cif0*cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0))
  G00 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  se0 = sqrt(G10+G00)
  dcif = cif1-cif0
  G0 = cumsum(M0)*dcif^2 + cumsum(M0*dcif^2) - 2*dcif*cumsum(M0*dcif)
  se = sqrt(G11+G10+G0)
  surv_diff = survdiff(Surv(Time,cstatus==1)~A,subset=subset)
  p = 1 - pchisq(surv_diff$chisq, length(surv_diff$n)-1)
  ate = cif1-cif0
  return(list(time1=time,time0=time,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=time,ate=ate,se=se,p.val=p))
}

surv.whileon <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  fit11 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==1])
  fit10 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==0])
  fit21 = survfit(Surv(Time,cstatus>1)~1, weights=weights, subset=subset[A[subset]==1])
  fit20 = survfit(Surv(Time,cstatus>1)~1, weights=weights, subset=subset[A[subset]==0])
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
  V1 = fit11[,2]^2
  V0 = fit21[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  G1 = cumsum(M1)*cif1^2 + cumsum(M1*(exp(-fit11[,1]-fit21[,1])+cif1)^2) -
    2*cif1*cumsum(M1*(exp(-fit11[,1]-fit21[,1])+cif1))
  G0 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  se1 = sqrt(G1+G0)
  V1 = fit10[,2]^2
  V0 = fit20[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  M1[is.infinite(M1)] = 0
  M0[is.infinite(M0)] = 0
  G1 = cumsum(M1)*cif0^2 + cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0)^2) -
    2*cif0*cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0))
  G0 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  se0 = sqrt(G1+G0)
  tt = sort(unique(c(time1,time0)))
  ate = matchy(cif1,time1,tt)-matchy(cif0,time0,tt)
  se = sqrt(matchy(se1,time1,tt)^2+matchy(se0,time0,tt)^2)
  return(list(time1=time1,time0=time0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=NULL))
}

surv.principal <- function(A,Time,cstatus,weights=rep(1,length(A)),subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  fit11 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==1])
  fit10 = survfit(Surv(Time,cstatus==1)~1, weights=weights, subset=subset[A[subset]==0])
  fit21 = survfit(Surv(Time,cstatus>1)~1, weights=weights, subset=subset[A[subset]==1])
  fit20 = survfit(Surv(Time,cstatus>1)~1, weights=weights, subset=subset[A[subset]==0])
  time1 = c(0, fit11$time)
  time0 = c(0, fit10$time)
  fit11 = rbind(0,cbind(fit11$cumhaz,fit11$std.err))
  fit10 = rbind(0,cbind(fit10$cumhaz,fit10$std.err))
  fit21 = rbind(0,cbind(fit21$cumhaz,fit21$std.err))
  fit20 = rbind(0,cbind(fit20$cumhaz,fit20$std.err))
  S1 = exp(-fit11[,1]-fit21[,1])
  S0 = exp(-fit10[,1]-fit20[,1])
  dcif1 = S1*diff(c(0,fit11[,1]))
  dcif0 = S0*diff(c(0,fit10[,1]))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  PR1 = min(S1 + cif1)
  PR0 = min(S0 + cif0)
  #PR1 = 1 - sum(S1*diff(c(0,fit21[,1])))
  #PR0 = 1 - sum(S0*diff(c(0,fit20[,1])))
  V1 = fit11[,2]^2
  V0 = fit21[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  G1 = cumsum(M1)*cif1^2 + cumsum(M1*(S1+cif1)^2) -
                 2*cif1*cumsum(M1*(S1+cif1))
  G0 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  G3 = cif1^2/PR1^2*sum(M1*(S1+cif1-PR1)^2)
  G2 = cif1^2/PR1^2*sum(M0*(cif1-PR1)^2)
  G5 = 2*cif1/PR1*(cumsum(M1*(S1+cif1)^2) + cumsum(M1)*cif1*PR1 -
                 cumsum(M1*(S1+cif1))*(PR1+cif1))
  G4 = 2*cif1/PR1*(cumsum(M0*cif1^2) + cumsum(M0)*cif1*PR1 -
                 cumsum(M0*cif1)*(PR1+cif1))
  se1 = sqrt(G1+G0+G3+G2-G5-G4)/PR1
  se1[is.nan(se1)] = rev(na.omit(se1))[1]
  
  V1 = fit10[,2]^2
  V0 = fit20[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  G1 = cumsum(M1)*cif0^2 + cumsum(M1*(S0+cif0)^2) -
    2*cif0*cumsum(M1*(S0+cif0))
  G0 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  G3 = cif0^2/PR0^2*sum(M1*(S0+cif0-PR0)^2) 
  G2 = cif0^2/PR0^2*sum(M0*(cif0-PR0)^2)
  G5 = 2*cif0/PR0*(cumsum(M1*(S0+cif0)^2) + cumsum(M1)*cif0*PR0 -
                 cumsum(M1*(S0+cif0))*(PR0+cif0))
  G4 = 2*cif0/PR0*(cumsum(M0*cif0^2) + cumsum(M0)*cif0*PR0 -
                 cumsum(M0*cif0)*(PR0+cif0))
  se0 = sqrt(G1+G0+G3+G2-G5-G4)/PR0
  se0[is.nan(se0)] = rev(na.omit(se0))[1]
  cif1 = cif1/PR1
  cif0 = cif0/PR0
  tt = sort(unique(c(time1,time0)))
  ate = matchy(cif1,time1,tt)-matchy(cif0,time0,tt)
  se = sqrt(matchy(se1,time1,tt)^2+matchy(se0,time0,tt)^2)
  return(list(time1=time1,time0=time0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=NULL))
}

surv.treatment.eff <- function(A,Time,cstatus,X=NULL,subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  n = length(A[subset])
  if (is.null(X)){
    psfit = glm(A~NULL, family='binomial', subset=subset)
    fit1 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==1])
    fit0 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==0])
  } else {
    psfit = glm(A~X, family='binomial', subset=subset)
    fit1 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==1])
    fit0 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==0])
  }
  ps = predict(psfit, type='response')
  tt1 = c(0,basehaz(fit1)$time)
  tt0 = c(0,basehaz(fit0)$time)
  tt = sort(unique(c(tt1,tt0)))
  if (!is.null(X)){
    X = as.matrix(X)
    Xb1 = as.numeric(X[subset,]%*%fit1$coefficients)
    Xb0 = as.numeric(X[subset,]%*%fit0$coefficients)
    Xb1c = as.numeric(X[subset,]%*%fit1c$coefficients)
    Xb0c = as.numeric(X[subset,]%*%fit0c$coefficients)
  } else {
    Xb1 = Xb0 = Xb1c = Xb0c = rep(0,n)
  }
  cumhaz1 = matchy(c(0,basehaz(fit1)$hazard),tt1,tt)
  cumhaz1 = exp(Xb1)%*%t(cumhaz1)
  cumhaz0 = matchy(c(0,basehaz(fit0)$hazard),tt0,tt)
  cumhaz0 = exp(Xb0)%*%t(cumhaz0)
  cumhaz1c = matchy(c(0,basehaz(fit1c)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = matchy(c(0,basehaz(fit0c)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  dN = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]==1))
  Y = sapply(tt, function(l) as.numeric(Time[subset]>=l))
  lam1 = cbind(0,t(apply(cumhaz1,1,diff)))
  dMP1 = (dN-Y*lam1)/exp(-cumhaz1-cumhaz1c)
  lam0 = cbind(0,t(apply(cumhaz0,1,diff)))
  dMP0 = (dN-Y*lam0)/exp(-cumhaz0-cumhaz0c)
  cif1x = A[subset]/ps*exp(-cumhaz1)*t(apply(dMP1,1,cumsum))+1-exp(-cumhaz1)
  cif0x = (1-A[subset])/(1-ps)*exp(-cumhaz0)*t(apply(dMP0,1,cumsum))+1-exp(-cumhaz0)
  cif1 = colMeans(cif1x)
  cif0 = colMeans(cif0x)
  se1 = apply(cif1x,2,sd)/sqrt(n)
  se0 = apply(cif0x,2,sd)/sqrt(n)
  ate = cif1-cif0
  se = sqrt(se1^2+se0^2)
  eif1 = t(t(cif1x)-cif1)
  eif0 = t(t(cif0x)-cif0)
  Tt = sum(((cif1-cif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt)])
  V1 = apply((t(eif1-eif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt),],2,sum)
  V2 = apply(((cif1-cif0)*apply(cbind(0,eif1+eif0),1,diff))[tt<0.99*max(tt),],2,sum)
  p = 2*pnorm(-abs(Tt)/sd(V1+V2)*sqrt(n))
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.composite.eff <- function(A,Time,cstatus,X=NULL,subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  n = length(A[subset])
  if (is.null(X)){
    psfit = glm(A~NULL, family='binomial', subset=subset)
    fit1 = coxph(Surv(Time,cstatus>0)~NULL, subset=subset[A[subset]==1])
    fit0 = coxph(Surv(Time,cstatus>0)~NULL, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==0])
  } else {
    psfit = glm(A~X, family='binomial', subset=subset)
    fit1 = coxph(Surv(Time,cstatus>0)~X, subset=subset[A[subset]==1])
    fit0 = coxph(Surv(Time,cstatus>0)~X, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==0])
  }
  ps = predict(psfit, type='response')
  tt1 = c(0,basehaz(fit1)$time)
  tt0 = c(0,basehaz(fit0)$time)
  tt = sort(unique(c(tt1,tt0)))
  if (!is.null(X)){
    X = as.matrix(X)
    Xb1 = as.numeric(X[subset,]%*%fit1$coefficients)
    Xb0 = as.numeric(X[subset,]%*%fit0$coefficients)
    Xb1c = as.numeric(X[subset,]%*%fit1c$coefficients)
    Xb0c = as.numeric(X[subset,]%*%fit0c$coefficients)
  } else {
    Xb1 = Xb0 = Xb1c = Xb0c = rep(0,n)
  }
  cumhaz1 = matchy(c(0,basehaz(fit1)$hazard),tt1,tt)
  cumhaz1 = exp(Xb1)%*%t(cumhaz1)
  cumhaz0 = matchy(c(0,basehaz(fit0)$hazard),tt0,tt)
  cumhaz0 = exp(Xb0)%*%t(cumhaz0)
  cumhaz1c = matchy(c(0,basehaz(fit1c)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = matchy(c(0,basehaz(fit0c)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  dN = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]>0))
  Y = sapply(tt, function(l) as.numeric(Time[subset]>=l))
  lam1 = cbind(0,t(apply(cumhaz1,1,diff)))
  dMP1 = (dN-Y*lam1)/exp(-cumhaz1-cumhaz1c)
  lam0 = cbind(0,t(apply(cumhaz0,1,diff)))
  dMP0 = (dN-Y*lam0)/exp(-cumhaz0-cumhaz0c)
  cif1x = A[subset]/ps*exp(-cumhaz1)*t(apply(dMP1,1,cumsum))+1-exp(-cumhaz1)
  cif0x = (1-A[subset])/(1-ps)*exp(-cumhaz0)*t(apply(dMP0,1,cumsum))+1-exp(-cumhaz0)
  cif1 = colMeans(cif1x)
  cif0 = colMeans(cif0x)
  se1 = apply(cif1x,2,sd)/sqrt(n)
  se0 = apply(cif0x,2,sd)/sqrt(n)
  ate = cif1-cif0
  se = sqrt(se1^2+se0^2)
  eif1 = t(t(cif1x)-cif1)
  eif0 = t(t(cif0x)-cif0)
  Tt = sum(((cif1-cif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt)])
  V1 = apply((t(eif1-eif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt),],2,sum)
  V2 = apply(((cif1-cif0)*apply(cbind(0,eif1+eif0),1,diff))[tt<0.99*max(tt),],2,sum)
  p = 2*pnorm(-abs(Tt)/sd(V1+V2)*sqrt(n))
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.removed.eff <- function(A,Time,cstatus,X=NULL,subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  n = length(A[subset])
  if (is.null(X)){
    psfit = glm(A~NULL, family='binomial', subset=subset)
    fit11 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==1])
    fit10 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==0])
    fit21 = coxph(Surv(Time,cstatus>1)~NULL, subset=subset[A[subset]==1])
    fit20 = coxph(Surv(Time,cstatus>1)~NULL, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==0])
  } else {
    psfit = glm(A~X, family='binomial', subset=subset)
    fit11 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==1])
    fit10 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==0])
    fit21 = coxph(Surv(Time,cstatus>1)~X, subset=subset[A[subset]==1])
    fit20 = coxph(Surv(Time,cstatus>1)~X, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==0])
  }
  ps = predict(psfit, type='response')
  tt1 = c(0,basehaz(fit11)$time)
  tt0 = c(0,basehaz(fit10)$time)
  tt = sort(unique(c(tt1,tt0)))
  if (!is.null(X)){
    X = as.matrix(X)
    Xb11 = as.numeric(X[subset,]%*%fit11$coefficients)
    Xb10 = as.numeric(X[subset,]%*%fit10$coefficients)
    Xb21 = as.numeric(X[subset,]%*%fit21$coefficients)
    Xb20 = as.numeric(X[subset,]%*%fit20$coefficients)
    Xb1c = as.numeric(X[subset,]%*%fit1c$coefficients)
    Xb0c = as.numeric(X[subset,]%*%fit0c$coefficients)
  } else {
    Xb11 = Xb10 = Xb21 = Xb20 = Xb1c = Xb0c = rep(0,n)
  }
  cumhaz11 = matchy(c(0,basehaz(fit11)$hazard),tt1,tt)
  cumhaz11 = exp(Xb11)%*%t(cumhaz11)
  cumhaz10 = matchy(c(0,basehaz(fit10)$hazard),tt0,tt)
  cumhaz10 = exp(Xb10)%*%t(cumhaz10)
  cumhaz21 = matchy(c(0,basehaz(fit21)$hazard),c(0,basehaz(fit21)$time),tt)
  cumhaz21 = exp(Xb21)%*%t(cumhaz21)
  cumhaz20 = matchy(c(0,basehaz(fit20)$hazard),c(0,basehaz(fit20)$time),tt)
  cumhaz20 = exp(Xb20)%*%t(cumhaz20)
  cumhaz1c = matchy(c(0,basehaz(fit1c)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = matchy(c(0,basehaz(fit0c)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  dN = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]==1))
  Y = sapply(tt, function(l) as.numeric(Time[subset]>=l))
  lam1 = cbind(0,t(apply(cumhaz11,1,diff)))
  dMP1 = (dN-Y*lam1)/exp(-cumhaz11-cumhaz21-cumhaz1c)
  lam0 = cbind(0,t(apply(cumhaz10,1,diff)))
  dMP0 = (dN-Y*lam0)/exp(-cumhaz10-cumhaz20-cumhaz0c)
  cif1x = A[subset]/ps*exp(-cumhaz11)*t(apply(dMP1,1,cumsum))+1-exp(-cumhaz11)
  cif0x = (1-A[subset])/(1-ps)*exp(-cumhaz10)*t(apply(dMP0,1,cumsum))+1-exp(-cumhaz10)
  cif1 = colMeans(cif1x)
  cif0 = colMeans(cif0x)
  se1 = apply(cif1x,2,sd)/sqrt(n)
  se0 = apply(cif0x,2,sd)/sqrt(n)
  ate = cif1-cif0
  se = sqrt(se1^2+se0^2)
  eif1 = t(t(cif1x)-cif1)
  eif0 = t(t(cif0x)-cif0)
  Tt = sum(((cif1-cif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt)])
  V1 = apply((t(eif1-eif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt),],2,sum)
  V2 = apply(((cif1-cif0)*apply(cbind(0,eif1+eif0),1,diff))[tt<0.99*max(tt),],2,sum)
  p = 2*pnorm(-abs(Tt)/sd(V1+V2)*sqrt(n))
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.natural.eff <- function(A,Time,cstatus,X=NULL,subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  n = length(A[subset])
  if (is.null(X)){
    psfit = glm(A~NULL, family='binomial', subset=subset)
    fit11 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==1])
    fit10 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==0])
    fit21 = coxph(Surv(Time,cstatus>1)~NULL, subset=subset[A[subset]==1])
    fit20 = coxph(Surv(Time,cstatus>1)~NULL, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==0])
  } else {
    psfit = glm(A~X, family='binomial', subset=subset)
    fit11 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==1])
    fit10 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==0])
    fit21 = coxph(Surv(Time,cstatus>1)~X, subset=subset[A[subset]==1])
    fit20 = coxph(Surv(Time,cstatus>1)~X, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==0])
  }
  ps = predict(psfit, type='response')
  tt11 = c(0,basehaz(fit11)$time)
  tt10 = c(0,basehaz(fit10)$time)
  tt21 = c(0,basehaz(fit21)$time)
  tt20 = c(0,basehaz(fit20)$time)
  tt = sort(unique(c(tt11,tt10,tt21,tt20)))
  if (!is.null(X)){
    X = as.matrix(X)
    Xb11 = as.numeric(X[subset,]%*%fit11$coefficients)
    Xb10 = as.numeric(X[subset,]%*%fit10$coefficients)
    Xb21 = as.numeric(X[subset,]%*%fit21$coefficients)
    Xb20 = as.numeric(X[subset,]%*%fit20$coefficients)
    Xb1c = as.numeric(X[subset,]%*%fit1c$coefficients)
    Xb0c = as.numeric(X[subset,]%*%fit0c$coefficients)
  } else {
    Xb11 = Xb10 = Xb21 = Xb20 = Xb1c = Xb0c = rep(0,n)
  }
  cumhaz11 = matchy(c(0,basehaz(fit11)$hazard),tt11,tt)
  cumhaz11 = exp(Xb11)%*%t(cumhaz11)
  cumhaz10 = matchy(c(0,basehaz(fit10)$hazard),tt10,tt)
  cumhaz10 = exp(Xb10)%*%t(cumhaz10)
  cumhaz21 = matchy(c(0,basehaz(fit21)$hazard),tt21,tt)
  cumhaz21 = exp(Xb21)%*%t(cumhaz21)
  cumhaz20 = matchy(c(0,basehaz(fit20)$hazard),tt20,tt)
  cumhaz20 = exp(Xb20)%*%t(cumhaz20)
  cumhaz1c = matchy(c(0,basehaz(fit1c)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = matchy(c(0,basehaz(fit0c)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  cumhaz1 = cumhaz11+cumhaz21
  cumhaz0 = cumhaz10+cumhaz20
  dN1 = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]==1))
  dN2 = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]>1))
  Y = sapply(tt, function(l) as.numeric(Time[subset]>=l))
  lam11 = cbind(0,t(apply(cumhaz11,1,diff)))
  lam10 = cbind(0,t(apply(cumhaz10,1,diff)))
  lam21 = cbind(0,t(apply(cumhaz21,1,diff)))
  lam20 = cbind(0,t(apply(cumhaz20,1,diff)))
  dMP11 = (dN1-Y*lam11)/exp(-cumhaz1-cumhaz1c)
  dMP21 = (dN2-Y*lam21)/exp(-cumhaz1-cumhaz1c)
  dMP10 = (dN1-Y*lam10)/exp(-cumhaz0-cumhaz0c)
  dMP20 = (dN2-Y*lam20)/exp(-cumhaz0-cumhaz0c)
  cif1 = t(apply(exp(-cumhaz11-cumhaz20)*t(apply(cbind(0,cumhaz11),1,diff)),1,cumsum))
  cif0 = t(apply(exp(-cumhaz0)*t(apply(cbind(0,cumhaz10),1,diff)),1,cumsum))
  cif1x = A[subset]/ps*t(apply((exp(-cumhaz11-cumhaz20)+cif1)*dMP11,1,cumsum))-
    A[subset]/ps*cif1*t(apply(dMP11,1,cumsum))+
    (1-A[subset])/(1-ps)*t(apply(cif1*dMP20,1,cumsum))-
    (1-A[subset])/(1-ps)*cif1*t(apply(dMP20,1,cumsum))+cif1
  cif0x = (1-A[subset])/(1-ps)*t(apply(exp(-cumhaz0)*dMP10,1,cumsum))-
    (1-A[subset])/(1-ps)*cif0*t(apply(dMP10+dMP20,1,cumsum))+
    (1-A[subset])/(1-ps)*t(apply(cif0*(dMP10+dMP20),1,cumsum))+cif0
  cif1 = colMeans(cif1x)
  cif0 = colMeans(cif0x)
  se1 = apply(cif1x,2,sd)/sqrt(n)
  se0 = apply(cif0x,2,sd)/sqrt(n)
  ate = cif1-cif0
  se = apply(cif1x-cif0x,2,sd)/sqrt(n)
  eif1 = t(t(cif1x)-cif1)
  eif0 = t(t(cif0x)-cif0)
  Tt = sum(((cif1-cif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt)])
  V1 = apply((t(eif1-eif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt),],2,sum)
  V2 = apply(((cif1-cif0)*apply(cbind(0,eif1+eif0),1,diff))[tt<0.99*max(tt),],2,sum)
  p = 2*pnorm(-abs(Tt)/sd(V1+V2)*sqrt(n))
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.whileon.eff <- function(A,Time,cstatus,X=NULL,subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  n = length(A[subset])
  if (is.null(X)){
    psfit = glm(A~NULL, family='binomial', subset=subset)
    fit11 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==1])
    fit10 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==0])
    fit21 = coxph(Surv(Time,cstatus>1)~NULL, subset=subset[A[subset]==1])
    fit20 = coxph(Surv(Time,cstatus>1)~NULL, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==0])
  } else {
    psfit = glm(A~X, family='binomial', subset=subset)
    fit11 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==1])
    fit10 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==0])
    fit21 = coxph(Surv(Time,cstatus>1)~X, subset=subset[A[subset]==1])
    fit20 = coxph(Surv(Time,cstatus>1)~X, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==0])
  }
  ps = predict(psfit, type='response')
  tt11 = c(0,basehaz(fit11)$time)
  tt10 = c(0,basehaz(fit10)$time)
  tt21 = c(0,basehaz(fit21)$time)
  tt20 = c(0,basehaz(fit20)$time)
  tt = sort(unique(c(tt11,tt10,tt21,tt20)))
  if (!is.null(X)){
    X = as.matrix(X)
    Xb11 = as.numeric(X[subset,]%*%fit11$coefficients)
    Xb10 = as.numeric(X[subset,]%*%fit10$coefficients)
    Xb21 = as.numeric(X[subset,]%*%fit21$coefficients)
    Xb20 = as.numeric(X[subset,]%*%fit20$coefficients)
    Xb1c = as.numeric(X[subset,]%*%fit1c$coefficients)
    Xb0c = as.numeric(X[subset,]%*%fit0c$coefficients)
  } else {
    Xb11 = Xb10 = Xb21 = Xb20 = Xb1c = Xb0c = rep(0,n)
  }
  cumhaz11 = matchy(c(0,basehaz(fit11)$hazard),tt11,tt)
  cumhaz11 = exp(Xb11)%*%t(cumhaz11)
  cumhaz10 = matchy(c(0,basehaz(fit10)$hazard),tt10,tt)
  cumhaz10 = exp(Xb10)%*%t(cumhaz10)
  cumhaz21 = matchy(c(0,basehaz(fit21)$hazard),tt21,tt)
  cumhaz21 = exp(Xb21)%*%t(cumhaz21)
  cumhaz20 = matchy(c(0,basehaz(fit20)$hazard),tt20,tt)
  cumhaz20 = exp(Xb20)%*%t(cumhaz20)
  cumhaz1c = matchy(c(0,basehaz(fit1c)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = matchy(c(0,basehaz(fit0c)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  cumhaz1 = cumhaz11+cumhaz21
  cumhaz0 = cumhaz10+cumhaz20
  dN1 = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]==1))
  dN2 = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]>1))
  Y = sapply(tt, function(l) as.numeric(Time[subset]>=l))
  lam11 = cbind(0,t(apply(cumhaz11,1,diff)))
  lam10 = cbind(0,t(apply(cumhaz10,1,diff)))
  lam21 = cbind(0,t(apply(cumhaz21,1,diff)))
  lam20 = cbind(0,t(apply(cumhaz20,1,diff)))
  dMP11 = (dN1-Y*lam11)/exp(-cumhaz1-cumhaz1c)
  dMP21 = (dN2-Y*lam21)/exp(-cumhaz1-cumhaz1c)
  dMP10 = (dN1-Y*lam10)/exp(-cumhaz0-cumhaz0c)
  dMP20 = (dN2-Y*lam20)/exp(-cumhaz0-cumhaz0c)
  cif1 = t(apply(exp(-cumhaz1)*t(apply(cbind(0,cumhaz11),1,diff)),1,cumsum))
  cif0 = t(apply(exp(-cumhaz0)*t(apply(cbind(0,cumhaz10),1,diff)),1,cumsum))
  cif1x = A[subset]/ps*t(apply(exp(-cumhaz1)*dMP11,1,cumsum))-
    A[subset]/ps*cif1*t(apply(dMP11+dMP21,1,cumsum))+
    A[subset]/ps*t(apply(cif1*(dMP11+dMP21),1,cumsum))+cif1
  cif0x = (1-A[subset])/(1-ps)*t(apply(exp(-cumhaz0)*dMP10,1,cumsum))-
    (1-A[subset])/(1-ps)*cif0*t(apply(dMP10+dMP20,1,cumsum))+
    (1-A[subset])/(1-ps)*t(apply(cif0*(dMP10+dMP20),1,cumsum))+cif0
  cif1 = colMeans(cif1x)
  cif0 = colMeans(cif0x)
  se1 = apply(cif1x,2,sd)/sqrt(n)
  se0 = apply(cif0x,2,sd)/sqrt(n)
  ate = cif1-cif0
  se = apply(cif1x-cif0x,2,sd)/sqrt(n)
  eif1 = t(t(cif1x)-cif1)
  eif0 = t(t(cif0x)-cif0)
  Tt = sum(((cif1-cif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt)])
  V1 = apply((t(eif1-eif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt),],2,sum)
  V2 = apply(((cif1-cif0)*apply(cbind(0,eif1+eif0),1,diff))[tt<0.99*max(tt),],2,sum)
  p = 2*pnorm(-abs(Tt)/sd(V1+V2)*sqrt(n))
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.principal.eff <- function(A,Time,cstatus,X=NULL,subset=NULL){
  if (is.null(subset)) subset = 1:length(A)
  n = length(A[subset])
  if (is.null(X)){
    psfit = glm(A~NULL, family='binomial', subset=subset)
    fit11 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==1])
    fit10 = coxph(Surv(Time,cstatus==1)~NULL, subset=subset[A[subset]==0])
    fit21 = coxph(Surv(Time,cstatus>1)~NULL, subset=subset[A[subset]==1])
    fit20 = coxph(Surv(Time,cstatus>1)~NULL, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~NULL, subset=subset[A[subset]==0])
  } else {
    psfit = glm(A~X, family='binomial', subset=subset)
    fit11 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==1])
    fit10 = coxph(Surv(Time,cstatus==1)~X, subset=subset[A[subset]==0])
    fit21 = coxph(Surv(Time,cstatus>1)~X, subset=subset[A[subset]==1])
    fit20 = coxph(Surv(Time,cstatus>1)~X, subset=subset[A[subset]==0])
    fit1c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==1])
    fit0c = coxph(Surv(Time,cstatus==0)~X, subset=subset[A[subset]==0])
  }
  ps = predict(psfit, type='response')
  tt11 = c(0,basehaz(fit11)$time)
  tt10 = c(0,basehaz(fit10)$time)
  tt21 = c(0,basehaz(fit21)$time)
  tt20 = c(0,basehaz(fit20)$time)
  tt = sort(unique(c(tt11,tt10,tt21,tt20)))
  if (!is.null(X)){
    X = as.matrix(X)
    Xb11 = as.numeric(X[subset,]%*%fit11$coefficients)
    Xb10 = as.numeric(X[subset,]%*%fit10$coefficients)
    Xb21 = as.numeric(X[subset,]%*%fit21$coefficients)
    Xb20 = as.numeric(X[subset,]%*%fit20$coefficients)
    Xb1c = as.numeric(X[subset,]%*%fit1c$coefficients)
    Xb0c = as.numeric(X[subset,]%*%fit0c$coefficients)
  } else {
    Xb11 = Xb10 = Xb21 = Xb20 = Xb1c = Xb0c = rep(0,n)
  }
  cumhaz11 = matchy(c(0,basehaz(fit11)$hazard),tt11,tt)
  cumhaz11 = exp(Xb11)%*%t(cumhaz11)
  cumhaz10 = matchy(c(0,basehaz(fit10)$hazard),tt10,tt)
  cumhaz10 = exp(Xb10)%*%t(cumhaz10)
  cumhaz21 = matchy(c(0,basehaz(fit21)$hazard),tt21,tt)
  cumhaz21 = exp(Xb21)%*%t(cumhaz21)
  cumhaz20 = matchy(c(0,basehaz(fit20)$hazard),tt20,tt)
  cumhaz20 = exp(Xb20)%*%t(cumhaz20)
  cumhaz1c = matchy(c(0,basehaz(fit1c)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = matchy(c(0,basehaz(fit0c)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  cumhaz1 = cumhaz11+cumhaz21
  cumhaz0 = cumhaz10+cumhaz20
  dN1 = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]==1))
  dN2 = sapply(tt, function(l) (Time[subset]==l)*(cstatus[subset]>1))
  Y = sapply(tt, function(l) as.numeric(Time[subset]>=l))
  lam11 = cbind(0,t(apply(cumhaz11,1,diff)))
  lam10 = cbind(0,t(apply(cumhaz10,1,diff)))
  lam21 = cbind(0,t(apply(cumhaz21,1,diff)))
  lam20 = cbind(0,t(apply(cumhaz20,1,diff)))
  dMP11 = (dN1-Y*lam11)/exp(-cumhaz1-cumhaz1c)
  dMP21 = (dN2-Y*lam21)/exp(-cumhaz1-cumhaz1c)
  dMP10 = (dN1-Y*lam10)/exp(-cumhaz0-cumhaz0c)
  dMP20 = (dN2-Y*lam20)/exp(-cumhaz0-cumhaz0c)
  cif11 = t(apply(exp(-cumhaz1)*t(apply(cbind(0,cumhaz11),1,diff)),1,cumsum))
  cif10 = t(apply(exp(-cumhaz0)*t(apply(cbind(0,cumhaz10),1,diff)),1,cumsum))
  cif21 = t(apply(exp(-cumhaz1)*t(apply(cbind(0,cumhaz21),1,diff)),1,cumsum))
  cif20 = t(apply(exp(-cumhaz0)*t(apply(cbind(0,cumhaz20),1,diff)),1,cumsum))
  cif1x = A[subset]/ps*t(apply(exp(-cumhaz1)*dMP11,1,cumsum))-
    A[subset]/ps*cif11*t(apply(dMP11+dMP21,1,cumsum))+
    A[subset]/ps*t(apply(cif11*(dMP11+dMP21),1,cumsum))+cif11
  cif0x = (1-A[subset])/(1-ps)*t(apply(exp(-cumhaz0)*dMP10,1,cumsum))-
    (1-A[subset])/(1-ps)*cif10*t(apply(dMP10+dMP20,1,cumsum))+
    (1-A[subset])/(1-ps)*t(apply(cif10*(dMP10+dMP20),1,cumsum))+cif10
  cif.wo1 = colMeans(cif1x)
  cif.wo0 = colMeans(cif0x)
  if.wo1 = t(t(cif1x)-cif.wo1)
  if.wo0 = t(t(cif0x)-cif.wo0)
  cif1x = A[subset]/ps*exp(-cumhaz1)*t(apply(dMP11+dMP21,1,cumsum))+1-exp(-cumhaz1)
  cif0x = (1-A[subset])/(1-ps)*exp(-cumhaz0)*t(apply(dMP10+dMP20,1,cumsum))+1-exp(-cumhaz0)
  cif.cv1 = colMeans(cif1x)
  cif.cv0 = colMeans(cif0x)
  if.cv1 = t(t(cif1x)-cif.cv1)
  if.cv0 = t(t(cif0x)-cif.cv0)
  cif1x = t((t(if.wo1)+cif.wo1)/min(1-cif.cv1+cif.wo1))+
    ((if.cv1-if.wo1)[,ncol(if.cv1)]/min(1-cif.cv1+cif.wo1)^2)%*%t(cif.wo1)
  cif0x = t((t(if.wo0)+cif.wo0)/min(1-cif.cv0+cif.wo0))+
    ((if.cv0-if.wo0)[,ncol(if.cv1)]/min(1-cif.cv0+cif.wo0)^2)%*%t(cif.wo0)
  cif1 = colMeans(cif1x)
  cif0 = colMeans(cif0x)
  se1 = apply(cif1x,2,sd)/sqrt(n)
  se0 = apply(cif0x,2,sd)/sqrt(n)
  ate = cif1-cif0
  se = apply(cif1x-cif0x,2,sd)/sqrt(n)
  eif1 = t(t(cif1x)-cif1)
  eif0 = t(t(cif0x)-cif0)
  Tt = sum(((cif1-cif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt)])
  V1 = apply((t(eif1-eif0)*diff(c(0,cif1+cif0)))[tt<0.99*max(tt),],2,sum)
  V2 = apply(((cif1-cif0)*apply(cbind(0,eif1+eif0),1,diff))[tt<0.99*max(tt),],2,sum)
  p = 2*pnorm(-abs(Tt)/sd(V1+V2)*sqrt(n))
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}

surv.ICH <- function(A,Time,cstatus,strategy='composite',cov1=NULL,method='np',
                     weights=NULL,subset=NULL){
  N = length(A)
  if (is.null(weights)) weights = rep(1,N)
  if (!is.null(cov1)) weights = weights*ipscore(A,cov1)
  if (method=='np'){
    if (strategy=='treatment') fit = surv.treatment(A,Time,cstatus,weights,subset)
    if (strategy=='composite') fit = surv.composite(A,Time,cstatus,weights,subset)
    if (strategy=='natural') fit = surv.natural(A,Time,cstatus,weights,subset)
    if (strategy=='removed') fit = surv.removed(A,Time,cstatus,weights,subset)
    if (strategy=='whileon') fit = surv.whileon(A,Time,cstatus,weights,subset)
    if (strategy=='principal') fit = surv.principal(A,Time,cstatus,weights,subset)
  } else {
    if (strategy=='treatment') fit = surv.treatment.eff(A,Time,cstatus,cov1,subset)
    if (strategy=='composite') fit = surv.composite.eff(A,Time,cstatus,cov1,subset)
    if (strategy=='natural') fit = surv.natural.eff(A,Time,cstatus,cov1,subset)
    if (strategy=='removed') fit = surv.removed.eff(A,Time,cstatus,cov1,subset)
    if (strategy=='whileon') fit = surv.whileon.eff(A,Time,cstatus,cov1,subset)
    if (strategy=='principal') fit = surv.principal.eff(A,Time,cstatus,cov1,subset)
  }
  return(c(fit,list(A=A,Time=Time,cstatus=cstatus,strategy=strategy,cov1=cov1,
                    method=method,weights=weights,subset=subset)))
}

surv.boot <- function(fit,nboot=0,seed=0){
  N = length(fit$A)
  time1 = fit$time1
  time0 = fit$time0
  Time = sort(unique(c(0,fit$Time[fit$cstatus>0],max(fit$Time))))
  cif1 = fit$cif1
  cif0 = fit$cif0
  se1 = fit$se1
  se0 = fit$se0
  ate = matchy(fit$ate,fit$time,Time)
  se = matchy(fit$se,fit$time,Time)
  if (nboot>1){
    cif1l = cif0l = te = NULL
    set.seed(seed)
    for(b in 1:nboot){
      #wt = as.vector(rmultinom(1,N,rep(1/N,N)))
      subset = sample(1:N,replace=TRUE)
      fitb = surv.ICH(fit$A,fit$Time,fit$cstatus,fit$strategy,fit$cov1,fit$method,
                      fit$weights,subset)
      cifb1 = matchy(fitb$cif1,fitb$time1,Time)
      cifb0 = matchy(fitb$cif0,fitb$time0,Time)
      cif1l = rbind(cif1l, cifb1)
      cif0l = rbind(cif0l, cifb0)
      te = rbind(te, cifb1-cifb0)
    }
    se1 = apply(cif1l,2,sd)
    se0 = apply(cif0l,2,sd)
    se1 = matchy(se1,Time,time1)
    se0 = matchy(se0,Time,time0)
    se = apply(te,2,sd)
  }
  return(list(time1=time1,time0=time0,cif1=fit$cif1,cif0=fit$cif0,
              se1=se1,se0=se0,time=Time,ate=ate,se=se,
              strategy=fit$strategy,method=fit$method))
}

plot.inc <- function(fit,decrease=FALSE,conf.int=.95,nboot=0,seed=0,xlab='Time',xlim=NULL,
                     ylim=c(0,1),legend=c('Treated','Control'),cex=0.9,...){
  if (fit$strategy=='treatment') stname = 'Treatment policy'
  if (fit$strategy=='composite') stname = 'Composite variable'
  if (fit$strategy=='natural') stname = 'Hypothetical I (natural)'
  if (fit$strategy=='removed') stname = 'Hypothetical II (removed)'
  if (fit$strategy=='whileon') stname = 'While on treatment'
  if (fit$strategy=='principal') stname = 'Principal stratum'
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
  plot(t1,cif1,type='s',col='brown',lwd=2,main=stname,
       xlab=xlab,ylab=ylab,ylim=ylim,...)
  points(t0,cif0,type='s',col='darkcyan',lwd=2)
  if (!is.null(conf.int)){
    fit.b = surv.boot(fit,nboot,seed)
    z = -qnorm((1-conf.int)/2)
    points(t1,cif1+fit.b$se1[i1]*z,type='s',lty=2,lwd=1.5,col='brown')
    points(t1,cif1-fit.b$se1[i1]*z,type='s',lty=2,lwd=1.5,col='brown')
    points(t0,cif0+fit.b$se0[i0]*z,type='s',lty=2,lwd=1.5,col='darkcyan')
    points(t0,cif0-fit.b$se0[i0]*z,type='s',lty=2,lwd=1.5,col='darkcyan')
  }
  legend(x,cex=cex,col=c('brown','darkcyan'),lwd=c(2,2),legend=legend)
}

plot.ate <- function(fit,decrease=FALSE,conf.int=.95,nboot=0,seed=0,xlab='Time',
                     xlim=NULL,ylim=c(-1,1),...){
  if (fit$strategy=='treatment') stname = 'Treatment policy'
  if (fit$strategy=='composite') stname = 'Composite variable'
  if (fit$strategy=='natural') stname = 'Hypothetical I (natural)'
  if (fit$strategy=='removed') stname = 'Hypothetical II (removed)'
  if (fit$strategy=='whileon') stname = 'While on treatment'
  if (fit$strategy=='principal') stname = 'Principal stratum'
  fit.b = surv.boot(fit,nboot=nboot,seed=seed)
  tm = fit.b$time
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
  plot(tm,dcif,type='s',main=stname,ylim=ylim,xlab=xlab,ylab=ylab,lwd=2,...)
  abline(h=0,lty=2)
  points(tm,ciu,type='s',lty=5,lwd=1.5,col='brown')
  points(tm,cil,type='s',lty=5,lwd=1.5,col='brown')
}
