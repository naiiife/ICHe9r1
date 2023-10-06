library(haven)
library(survival)
library(survminer)
library(cmprsk)
library(nnet)

setwd('D:/Data/PKU- Beijing University/PKU- Beijing University/LEADER EX2211-3748/ADaM')
dat = read_sas('adtte.sas7bdat')
#dat.adadj = read_sas('adadj.sas7bdat')
dat.mace = dat[dat$PARAMCD=='MACEEVTM',]
dat.noncvdeath = dat[dat$PARAMCD=='NONCVTM',]
dat.cvdeath = dat[dat$PARAMCD=='MCECVDTM',]
T.mace = dat.mace$AVAL
C.mace = 1-dat.mace$CNSR
T.noncvdeath = dat.noncvdeath$AVAL
C.noncvdeath = 1-dat.noncvdeath$CNSR
T.cvdeath = dat.cvdeath$AVAL
C.cvdeath = 1-dat.cvdeath$CNSR
C.mace[T.mace>54] = 0
C.noncvdeath[T.noncvdeath>54] = 0
C.cvdeath[T.cvdeath>54] = 0
T.mace[T.mace>54] = 54
T.noncvdeath[T.noncvdeath>54] = 54
T.cvdeath[T.cvdeath>54] = 54
T.death = (T.cvdeath+T.noncvdeath-abs(T.cvdeath-T.noncvdeath))/2
C.death = C.cvdeath+C.noncvdeath
Z = as.numeric(dat.cvdeath$TRTA=='Lira')
A = Z

## Analysis 1
setwd('D:/PAPER/ICH_time_to_event/Figures1')
## Primary outcome: MACE
## Intercurrent event: NCV death
cstatus = C.mace + 2*C.noncvdeath
cstatus[cstatus>2&C.mace==C.noncvdeath] = 1
cstatus[cstatus>2] = 2
Time = (T.mace+T.noncvdeath-abs(T.mace-T.noncvdeath))/2
table(Z,cstatus)
# Composite variable
fit.cv = surv.ICH(Z,Time,cstatus,'composite')
plot.inc(fit.cv, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='Composite Variable', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
text(6,0.2, 'P = 0.0222')
fit.cvb = surv.boot(fit.cv)
plot.boot(fit.cvb, ylim=c(-0.1,0.05), main='Composite Variable', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# Hypothetical
fit.h1 = surv.ICH(Z,Time,cstatus,'natural')
plot.inc(fit.h1, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='Hypothetical Scenario I', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
text(6,0.2, 'P = 0.0090')
fit.h1b = surv.boot(fit.h1)
plot.boot(fit.h1b, ylim=c(-0.1,0.05), main='Hypothetical Scenario I', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
fit.h2 = surv.ICH(Z,Time,cstatus,'removed')
plot.inc(fit.h2, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='Hypothetical Scenario II', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
text(6,0.2, 'P = 0.0090')
fit.h2b = surv.boot(fit.h2)
plot.boot(fit.h2b, ylim=c(-0.1,0.05), main='Hypothetical Scenario II', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# While on
fit.wo = surv.ICH(Z,Time,cstatus,'whileon')
plot.inc(fit.wo, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='While On', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
fit.wob = surv.boot(fit.wo)
plot.boot(fit.wob, ylim=c(-0.1,0.05), main='While On', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# Principal stratum
fit.ps = surv.ICH(Z,Time,cstatus,'principal')
plot.inc(fit.ps, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='Principal Stratum', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
fit.psb = surv.boot(fit.ps)
plot.boot(fit.psb, ylim=c(-0.1,0.05), main='Principal Stratum', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# All
plot(fit.cvb$Time,fit.cvb$ate, ylim=c(-0.05,0.03), lwd=1.3, type='s', col=3,
     main='Treatment Effect (Analysis I)', xlab='Month', ylab='Diff in cumulative incidences', xaxt='n')
axis(1,seq(0,54,9))
points(fit.h1b$Time,fit.h1b$ate, type='s', lwd=1.3, col=4)
points(fit.h2b$Time,fit.h2b$ate, type='s', lwd=1.3, col=5)
points(fit.wob$Time,fit.wob$ate, type='s', lwd=1.3, col=6)
points(fit.psb$Time,fit.psb$ate, type='s', lwd=1.3, col=7)
abline(h=0, lty=2)
legend('bottomleft',cex=0.8,col=3:7,lwd=rep(1.3,5),
       legend=c('Composite variable','Hypothetical scenario I',
                'Hypothetical scenario II','While on','Principal stratum'))

## Test
sdf = survdiff(Surv(Time,cstatus>=1)~Z)
1 - pchisq(sdf$chisq, length(sdf$n) - 1) #0.0222
sdf = survdiff(Surv(Time,cstatus==1)~Z)
1 - pchisq(sdf$chisq, length(sdf$n) - 1) #0.0090

## Analysis 2
setwd('D:/PAPER/ICH_time_to_event/Figures2')
## Primary outcome: all cause death
## Intercurrent event: MACE
cstatus = C.death + 2*C.mace
cstatus[cstatus>2&C.death==C.mace] = 1
cstatus[cstatus>2] = 2
Time = (T.death+T.mace-abs(T.death-T.mace))/2
table(Z,cstatus)
# Treatment policy
fit.tp = surv.ICH(Z,T.death,C.death,'treatment')
plot.inc(fit.tp, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='Treatment Policy', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
text(6,0.2, 'P = 0.0170')
fit.tpb = surv.boot(fit.tp)
plot.boot(fit.tpb, ylim=c(-0.1,0.05), main='Treatment Policy', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# Composite variable
fit.cv = surv.ICH(Z,Time,cstatus,'composite')
plot.inc(fit.cv, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='Composite Variable', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
text(6,0.2, 'P = 0.0222')
fit.cvb = surv.boot(fit.cv)
plot.boot(fit.cvb, ylim=c(-0.1,0.05), main='Composite Variable', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# Hypothetical
fit.h1 = surv.ICH(Z,Time,cstatus,'natural')
plot.inc(fit.h1, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='Hypothetical Scenario I', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
text(6,0.2, 'P = 0.0145')
fit.h1b = surv.boot(fit.h1)
plot.boot(fit.h1b, ylim=c(-0.1,0.05), main='Hypothetical Scenario I', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
fit.h2 = surv.ICH(Z,Time,cstatus,'removed')
plot.inc(fit.h2, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'),
          main='Hypothetical Scenario II', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
text(6,0.2, 'P = 0.0145')
fit.h2b = surv.boot(fit.h2)
plot.boot(fit.h2b, ylim=c(-0.1,0.05), main='Hypothetical Scenario II', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# While on
fit.wo = surv.ICH(Z,Time,cstatus,'whileon')
plot.inc(fit.wo, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='While On', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
fit.wob = surv.boot(fit.wo)
plot.boot(fit.wob, ylim=c(-0.1,0.05), main='While On', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# Principal stratum
fit.ps = surv.ICH(Z,Time,cstatus,'principal')
plot.inc(fit.ps, ylim=c(0,0.25), legend=c('Liraglutide','Placebo'), 
         main='Principal Stratum', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
fit.psb = surv.boot(fit.ps)
plot.boot(fit.psb, ylim=c(-0.1,0.05), main='Principal Stratum', xlab='Month', xaxt='n')
axis(1,seq(0,54,9))
# All
plot(fit.tpb$Time,fit.tpb$ate, ylim=c(-0.05,0.03), lwd=1.3, type='s', col=2,
     main='Treatment Effect (Analysis II)', xlab='Month', ylab='Diff in cumulative incidences', xaxt='n')
axis(1,seq(0,54,9))
points(fit.cvb$Time,fit.cvb$ate, type='s', lwd=1.3, col=3)
points(fit.h1b$Time,fit.h1b$ate, type='s', lwd=1.3, col=4)
points(fit.h2b$Time,fit.h2b$ate, type='s', lwd=1.3, col=5)
points(fit.wob$Time,fit.wob$ate, type='s', lwd=1.3, col=6)
points(fit.psb$Time,fit.psb$ate, type='s', lwd=1.3, col=7)
abline(h=0, lty=2)
legend('bottomleft',cex=0.8,col=2:7,lwd=rep(1.3,6),
       legend=c('Treatment policy','Composite variable','Hypothetical scenario I',
                'Hypothetical scenario II','While on','Principal stratum'))

## Test
sdf = survdiff(Surv(T.death,C.death)~Z)
1 - pchisq(sdf$chisq, length(sdf$n) - 1) #0.0170
sdf = survdiff(Surv(Time,cstatus>=1)~Z)
1 - pchisq(sdf$chisq, length(sdf$n) - 1) #0.0222
sdf = survdiff(Surv(Time,cstatus==1)~Z)
1 - pchisq(sdf$chisq, length(sdf$n) - 1) #0.0145

