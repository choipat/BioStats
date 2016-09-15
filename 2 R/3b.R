library(ALL); data(ALL)
B1<-exprs(ALL[,ALL$BT=="B1"])
mean.B1<-apply(B1,1,mean)
mean.B1