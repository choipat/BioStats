library(multtest); data(golub)
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL", "AML"))
golub.AML<-golub[,gol.fac=="AML"]
AML5<-golub.AML[1:5,]
getwd()
write.csv(AML5,file="AML5.csv")