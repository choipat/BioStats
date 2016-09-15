fx<-function(x) (2^x)*exp(-2)/factorial(x)
f_x<-function(x) fx(x)*(x>=0)
integrate(f_x, lower=-2, upper=4)
