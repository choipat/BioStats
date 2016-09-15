library(ALL); data(ALL)
B1<-exprs(ALL[,ALL$BT=="B1"])
mean.B1<-apply(B1,1,mean)
order.B1<-order(mean.B1, decreasing=TRUE)
mean.B1[order.B1[1:3]]