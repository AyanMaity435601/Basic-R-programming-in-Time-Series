# Soumyajyoti Das B2130046

dev.off()

data <- read.csv("OvarianCancerData.csv")
data

Status = data$status
Group = data$treatment

#Survival Analysis needs to be performed on the data

#time -> time to event
#status -> censoring
#treatment, age -> time independent co-variate

library(survival)
library(MASS)
library(fitdistrplus)

head(data)        # Have a glimpse of data
names(data)       # the column names of the data
str(data)         # types of variables
summary(data) 

# Step 1. Preparing Right Censor Data
#-------------------------------------------------------------------------
Lifetime = data$time                  # Right Censor Survival Data
# Data preparation for AFT parametric fit
CensorData = data.frame(cbind(Lifetime,data$status))                 # Censor Data      
SurvObject = Surv(data$time,(data$status)==1)    # Survival Object
# Preparing right-censored data for parametric fit
left = Lifetime; right = Lifetime; right[which(data$status==0)] = NA # left=right for failed and NA for censored 
LRdata = data.frame(cbind(left,right))

# Kaplan Meier estimate of Survival function
Result = survfit(Surv(Lifetime,data$status)~1,
                 CensorData,stype=1)  # stype=1: direct method (i.e. KM) to estimate S(t)
plot(Result,xlab = "time", 
     ylab="Survival probability",ylim=c(0,1),col=1) # KM estimator and confidence intervals
title("KM Estimate for Survival Function")
legend('bottomleft',c("KM Estimate","95% CI"),lty=c(1,2),col=c(1,1))

dev.off()

# Nelson Aalen estimate of Cumulative Hazard function
ResultH = survfit(Surv(Lifetime,data$status)~1,CensorData,conf.type="plain",
                  stype=2,ctype=1)  # stype=2,ctype=1: indirect method to estimate S(t) from H(t)
plot(ResultH$time,ResultH$cumhaz,xlab = "time", 
     ylab="Cummulative hazard",ylim=c(0,3),type="l",col=1)
lines(ResultH$time,ResultH$cumhaz+1.96*ResultH$std.err,lty=2,col=1)
lines(ResultH$time,ResultH$cumhaz-1.96*ResultH$std.err,lty=2,col=1)
title("NA Estimate for Cumulative Hazard Function")
legend('topleft',c("NA Estimate for H(t)","95% CI"),lty=c(1,2),col=c(1,1))

dev.off()

plot(ResultH,xlab = "time", 
     ylab="Survival probability",ylim=c(0,1),col=1)
title("NA Estimate for Survival Function")
legend('bottomleft',c("NA Estimate for S(t)","95% CI"),lty=c(1,2),col=c(1,1))

dev.off()

# Kernel estimate of hazard function for right-censored data

library(muhaz)

mhFit=muhaz(Lifetime,Status,kern="epanechnikov")
plot(mhFit$est.grid,mhFit$haz.est,type="l",xlab="time",ylab="Hazard rate")
title("Kernel Estimate of Hazard Function by the Epanechnikov Kernel")

dev.off()

# Comparison of survival/Hazard functions between groups [KM/NA]
# Grouping on Treatment Type (Type=1 & Type=2)

head(cbind(Lifetime,Status,Group))
par(mfrow=c(1,2))

# Right censored data: graphical comparison of two groups

KM_G=survfit(Surv(Lifetime,Status)~Group)
plot(KM_G,col=c(1,2),xlab="time",ylab="survival probability",
     xlim=c(0,max(KM_G$time)+1),ylim=c(0,1.05))
legend('topright',c("Tr:Type 1","Tr:Type 2"),lty=c(1,1),col=c(1,2),cex=.6)
title(main='KM Estimate for Survival of 2 Treatment Types')

NA_G=survfit(Surv(Lifetime,Status)~Group,stype=2,ctype=1)
plot(NA_G$time[seq(1,13)],NA_G$cumhaz[seq(1,13)],col=c(1),
     xlab="time",ylab="cumulative hazard",type='s',
     xlim=c(0,max(NA_G$time)+1),ylim=c(0,ceiling(max(NA_G$cumhaz))+.1))
lines(NA_G$time[seq(14,26)],NA_G$cumhaz[seq(14,26)],col=c(2),type='s')
legend('topleft',c("Tr: Type 1","Tr: Type 2"),lty=c(1,1),col=c(1,2),cex=.6)
title(main='NA Estimate for Cumulative Hazard of 2 Treatment Types')

dev.off()

# Right censored data: statistical comparison of two groups
survdiff(Surv(Lifetime,Status)~Group,rho = 0)           # logrank test
survdiff(Surv(Lifetime,Status)~Group,rho = 1)           # Gehan-Wilcoxon test

# Building Cox proportional hazard models Semi-Parametric

ResultSP = coxph(SurvObject ~ treatment+age, data=data)     # Build the coxph model
summary(ResultSP)
BaselineHazard=basehaz(ResultSP)
plot(BaselineHazard$time,BaselineHazard$hazard,type='s',xlab='Time',ylab='Baseline Hazard')

ResultSP = coxph(SurvObject ~ age, data=data)     # Build the coxph model
summary(ResultSP)
BaselineHazard=basehaz(ResultSP)
lines(BaselineHazard$time,BaselineHazard$hazard,type='s',lty=2)
legend('topleft',c("H0(t), 2 covariates","H0(t), 1 covariate"),lty=c(1,2),col=c(1,1))
title("Cox Proportional Model")

dev.off()

# Compare predictions for X = mu+sd and X = mu-sd (on Age at mean)
newdata <- data.frame(rep(mean(data$age)+sd(data$age),2))
colnames(newdata)=c('age')

plot(survfit(ResultSP, newdata, se.fit=F), xlab = "Time", ylab="Survival probability")
ts = seq(min(Lifetime),max(Lifetime))
newdata[,'age']=rep(mean(data$age)-sd(data$age),2)
lines(survfit(ResultSP, newdata, se.fit=F), lty = 2)
legend('topright',c("Cox, X = mu+sd","Cox, X = mu-sd"),lty=c(1,2),col=c(1,1))
title("Prediction using the Age Co-variate")

# dev.off()
