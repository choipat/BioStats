library(multtest); data(golub)
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL", "AML"))
meanAML <- apply(golub[,gol.fac=="AML"], 1, mean)
meanAML