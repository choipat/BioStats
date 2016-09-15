library(multtest); data(golub)
sd.all<-apply(golub,1,sd)
sd.all
sum(sd.all>1)