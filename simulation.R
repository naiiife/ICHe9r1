lid = Sys.getenv("SLURM_ARRAY_TASK_ID")
if (lid=='') lid=0
lid = as.numeric(lid)

a1 = 0.05
a0 = 0.03
c1 = 0.04
c0 = 0.05
N = 500
generatedata <- function(N){
  A = rbinom(N,1,0.5)
  T1 = rweibull(N, 2, sqrt(2/a1))
  T0 = rweibull(N, 2, sqrt(2/a0))
  R1 = rexp(N, c1)
  R0 = rexp(N, c0)
  C = runif(N,4,8)
  C[C>7] = 7
  T = T1*A + T0*(1-A)
  R = R1*A + R0*(1-A)
  R[R>=T] = 99
  dT = as.numeric(T <= C)
  dR = as.numeric(R <= C)
  T = T*dT + C*(1-dT)
  R = R*dR + C*(1-dR)
  Time = (T+R-abs(T-R))/2
  cstatus = dT + 2*dR
  cstatus[cstatus>2] = 2
  return(list(A=A,T=T,R=R,dT=dT,dR=dR,Time=Time,cstatus=cstatus))
}

source('ICH_functions.R')
B = 5000
x = 1:6
mt = 7
ate.tp = (1 - exp(-a1*x^2/2)) - (1 - exp(-a0*x^2/2))
ate.cv = (1 - exp(-a1*x^2/2-c1*x)) - (1 - exp(-a0*x^2/2-c0*x))
ate.wo = (1 - exp(-a1*x^2/2-c1*x) - exp(c1^2/2/a1)*sqrt(2*pi*c1^2/a1)*(
  pnorm(sqrt(a1)*(x+c1/a1)) - pnorm(c1/sqrt(a1)))) -
  (1 - exp(-a0*x^2/2-c0*x) - exp(c0^2/2/a0)*sqrt(2*pi*c0^2/a0)*(
    pnorm(sqrt(a0)*(x+c0/a0)) - pnorm(c0/sqrt(a0))))
ate.h1 = (1 - exp(-a1*x^2/2-c0*x) - exp(c0^2/2/a1)*sqrt(2*pi*c0^2/a1)*(
  pnorm(sqrt(a1)*(x+c0/a1)) - pnorm(c0/sqrt(a1)))) -
  (1 - exp(-a0*x^2/2-c0*x) - exp(c0^2/2/a0)*sqrt(2*pi*c0^2/a0)*(
    pnorm(sqrt(a0)*(x+c0/a0)) - pnorm(c0/sqrt(a0))))
ate.h2 = (1 - exp(-a1*x^2/2)) - (1 - exp(-a0*x^2/2))
ate.ps = (1 - exp(-a1*x^2/2-c1*x) - exp(c1^2/2/a1)*sqrt(2*pi*c1^2/a1)*(
  pnorm(sqrt(a1)*(x+c1/a1)) - pnorm(c1/sqrt(a1))))/
  (1 - exp(c1^2/2/a1)*sqrt(2*pi*c1^2/a1)*(
    pnorm(sqrt(a1)*(mt+c1/a1)) - pnorm(c1/sqrt(a1)))) - 
  (1 - exp(-a0*x^2/2-c0*x) - exp(c0^2/2/a0)*sqrt(2*pi*c0^2/a0)*(
    pnorm(sqrt(a0)*(x+c0/a0)) - pnorm(c0/sqrt(a0))))/
  (1 - exp(c0^2/2/a0)*sqrt(2*pi*c0^2/a0)*(
    pnorm(sqrt(a0)*(mt+c0/a0)) - pnorm(c0/sqrt(a0))))

set.seed(2024+lid)
dat = generatedata(N)
fit = surv.ICH(dat$A,dat$T,dat$dT,'treatment')
hat = matchy(fit$ate, fit$time, x)
se = matchy(fit$se, fit$time, x)
hat.tp = hat
cr.tp = (hat+1.96*se >= ate.tp)*(hat-1.96*se <= ate.tp)
wd.tp = 2*1.96*se
fit = surv.boot(fit, nboot=100)
se = matchy(fit$se, fit$time, x)
crb.tp = (hat+1.96*se >= ate.tp)*(hat-1.96*se <= ate.tp)
wdb.tp = 2*1.96*se

fit = surv.ICH(dat$A,dat$Time,dat$cstatus,'composite')
hat = matchy(fit$ate, fit$time, x)
se = matchy(fit$se, fit$time, x)
hat.cv = hat
cr.cv = (hat+1.96*se >= ate.cv)*(hat-1.96*se <= ate.cv)
wd.cv = 2*1.96*se
fit = surv.boot(fit, nboot=100)
se = matchy(fit$se, fit$time, x)
crb.cv = (hat+1.96*se >= ate.cv)*(hat-1.96*se <= ate.cv)
wdb.cv = 2*1.96*se

fit = surv.ICH(dat$A,dat$Time,dat$cstatus,'whileon')
hat = matchy(fit$ate, fit$time, x)
se = matchy(fit$se, fit$time, x)
hat.wo = hat
cr.wo = (hat+1.96*se >= ate.wo)*(hat-1.96*se <= ate.wo)
wd.wo = 2*1.96*se
fit = surv.boot(fit, nboot=100)
se = matchy(fit$se, fit$time, x)
crb.wo = (hat+1.96*se >= ate.wo)*(hat-1.96*se <= ate.wo)
wdb.wo = 2*1.96*se

fit = surv.ICH(dat$A,dat$Time,dat$cstatus,'natural')
hat = matchy(fit$ate, fit$time, x)
se = matchy(fit$se, fit$time, x)
hat.h1 = hat
cr.h1 = (hat+1.96*se >= ate.h1)*(hat-1.96*se <= ate.h1)
wd.h1 = 2*1.96*se
fit = surv.boot(fit, nboot=100)
se = matchy(fit$se, fit$time, x)
crb.h1 = (hat+1.96*se >= ate.h1)*(hat-1.96*se <= ate.h1)
wdb.h1 = 2*1.96*se

fit = surv.ICH(dat$A,dat$Time,dat$cstatus,'removed')
hat = matchy(fit$ate, fit$time, x)
se = matchy(fit$se, fit$time, x)
hat.h2 = hat
cr.h2 = (hat+1.96*se >= ate.h2)*(hat-1.96*se <= ate.h2)
wd.h2 = 2*1.96*se
fit = surv.boot(fit, nboot=100)
se = matchy(fit$se, fit$time, x)
crb.h2 = (hat+1.96*se >= ate.h2)*(hat-1.96*se <= ate.h2)
wdb.h2 = 2*1.96*se

fit = surv.ICH(dat$A,dat$Time,dat$cstatus,'principal')
hat = matchy(fit$ate, fit$time, x)
se = matchy(fit$se, fit$time, x)
hat.ps = hat
cr.ps = (hat+1.96*se >= ate.ps)*(hat-1.96*se <= ate.ps)
wd.ps = 2*1.96*se
fit = surv.boot(fit, nboot=100)
se = matchy(fit$se, fit$time, x)
crb.ps = (hat+1.96*se >= ate.ps)*(hat-1.96*se <= ate.ps)
wdb.ps = 2*1.96*se

res = rbind(hat.tp,hat.cv,hat.wo,hat.h1,hat.h2,hat.ps,
            cr.tp,cr.cv,cr.wo,cr.h1,cr.h2,cr.ps,
            wd.tp,wd.cv,wd.wo,wd.h1,wd.h2,wd.ps,
            crb.tp,crb.cv,crb.wo,crb.h1,crb.h2,crb.ps,
            wdb.tp,wdb.cv,wdb.wo,wdb.h1,wdb.h2,wdb.ps)
save(res, file=paste0('res',lid,'.Rdata'))


#setwd(...)
#true = rbind(ate.tp,ate.cv,ate.wo,ate.h1,ate.h2,ate.ps)
#bias = crate = width = crateb = widthb = 0
#for (b in 1:B){
#  res = try(load(paste0('res',lid,'.Rdata')))
#  if ('try-error' %in% class(res)) next
#  bias = bias + (res[1:6,]-true)
#  crate = crate + res[7:12,]
#  width = width + res[13:18,]
#  crateb = crateb + res[19:24,]
#  widthb = widthb + res[25:30,]
#  count = count+1
#}
#bias = bias/count
#crate = crate/count
#width = width/count
#crateb = crateb/count
#widthb = widthb/count
#row.names(bias) = row.names(crate) = row.names(width) = 
#  row.names(crateb) = row.names(widthb) =
#  c('Treatment policy','Composite variable','While on treatment',
#    'Hypothetical I (natural)','Hypothetical II (removed)','Principal stratum')
#colnames(bias) = colnames(crate) = colnames(widthb) = colnames(crateb) = 
#  colnames(width) = 1:6

#library(xtable)
#xtable(crate, digits=3, caption='Coverage rate of confidence intervals')
#xtable(width, digits=3, caption='Width of confidence intervals')