library(multtest); data(golub)
x<-golub[101,]
y<-golub[102,]
plot(x,y,xlab=golub.gnames[101,2],ylab=golub.gnames[102,2])