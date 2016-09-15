library(multtest); data(golub)
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL", "AML"))
golub.ALL<-golub[,gol.fac=="ALL"]
ALL5<-golub.ALL[1:5,]
getwd()
write.table(ALL5,file="ALL5.txt")