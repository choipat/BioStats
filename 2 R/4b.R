str(trees)
plot(trees$Height~trees$Girth,xlim=c(8,21),ylim=c(0,90),col="blue",type="o",xlab="Girth", ylab="Height and Volume", pch="+")
points(trees$Girth,trees$Volume,col="red",type="o", pch="o")
legend("bottomright",c("Height","Volume"),col=c("blue","red"),lty=c(1,1))