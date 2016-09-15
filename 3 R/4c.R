Y_range<-(0:3)
EY<-sum(Y_range*dbinom(Y_range,size=3,p=0.25))
VarY<-sum((Y_range-EY)^2*dbinom(Y_range,size=3,p=0.25))
VarY



n<-3
p<-0.25
n*p*(1-p)