library(survival)


mydata=read.table("GBMLGG.txt", sep = "\t", header = TRUE)
mydata$T1=mydata$PFI.time
mydata$T2=mydata$OS.time
mydata$censor1=mydata$PFI
mydata$censor2=mydata$OS
km_fit1 <- survfit(Surv(mydata$T1,mydata$censor1) ~ 1,data=mydata)
plot(km_fit1,col="blue",lty=1)
