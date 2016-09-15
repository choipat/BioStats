library(multtest); data(golub)
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL", "AML"))
meanALL <- apply(golub[,gol.fac=="ALL"], 1, mean)
meanALL