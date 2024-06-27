mydata=read.table("GBMLGG.txt", sep = "\t", header = TRUE)
mydata$T1=mydata$PFI.time
mydata$T2=mydata$OS.time
mydata$censor1=mydata$PFI
mydata$censor2=mydata$OS
plot(mydata$T1[mydata$censor1==0 | mydata$censor2==0],  mydata$T2[mydata$censor1==0 | mydata$censor2==0],type="n",xlab="",ylab="",xaxt="n",yaxt="n")
points(mydata$T1[mydata$censor1==1 & mydata$censor2==1],  mydata$T2[mydata$censor1==1 & mydata$censor2==1],pch=20,cex=1,lty=4)
points(mydata$T1[mydata$censor1==0 | mydata$censor2==0],  mydata$T2[mydata$censor1==0 | mydata$censor2==0],pch=1,cex=1,lty=4)
abline(0, 1,lty=2)
mtext(seq(0,4000,1000),side=1,las=1,at=seq(0,4000,1000),cex=1,font=1,line=0.4)
axis(side=1,seq(0,4000,1000),tcl=-0.2,labels=FALSE)
mtext(seq(0,6500,1000),side=2,las=3,at=seq(0,6500,1000),cex=1,font=1,line=0.4)
axis(side=2,seq(0,6500,1000),tcl=-0.1,labels=FALSE)
mtext(expression(paste("TTP ","(","day",")")),side=1,las=1,at=2000,cex=1,font=1,line=1.3)
mtext(expression(paste("OS ","(","day",")")),side=2,las=3,at=3250,cex=1,font=1,line=1.2)

