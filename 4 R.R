# 1
# Finding sd
3/sqrt(5)


# 2
# Find the probability that Y is great than 15.
n <- 20
p <- 0.7
Mean <- n*p # Calculating Mean
Mean # Display Mean
Var <- n*p*(1-p) # Calculating Variance
Var # Display Variance
1-pnorm(15, mean=14, sd=sqrt(4.2)/sqrt(100)) # Dividing with sqrt(100) as the number of microRNAs = 100


# 4
# Find the mean of Y.
x1<-rchisq(10000, df = 10) # chisq
x2<-rgamma(10000, 1, 2) # Gamma
x3<-rt(10000, 3) # t-distribution
y<- (sqrt(x1)*x2) + (4*(x3^2)) # Defining y
mean(y) # Calculating Mean


#5
# A plot and brief explanation
# Used few websites like http://www.statmethods.net/graphs/density.html to solve this problem
an <- sqrt(2*log(n)) - 0.5*(log(log(n))+log(4*pi))*(2*log(n))^(-1/2)
bn <- (2*log(n))^(-1/2)
f<-function(x){exp(-x)*exp(-exp(-x))}
nsim <- 10000 # 10000 Simulations
for (i in 1:1000) {e[i] <- (max(rnorm(nsim))-an)/bn} # Subtract an and divide by bn, load extreme values into e
plot(density(e), col = "red", ylim = c(0,0.45), main = "Extreme data density, normal and extreme values distribution") # To generate the density of extreme data
curve(f, add=TRUE, col = "blue", ylim=c(0,0.45), xlim=c(-4,10)) # Extreme values ; add: add to existing plot
curve(dnorm, add=TRUE, col = "yellow") # Normal values
legend("topright",c("Density", "Extreme value", "Normal distribution"),col=c("red", "blue", "yellow"),lty=c(1,1,1))


#3
require(mvtnorm)
data<-rmvnorm(50,mean=c(9,10),sigma=matrix(c(3,2,2,5),nrow=2))
meanxy<-apply(data,1,mean)
varxy<-apply(data,1,var)
mean(meanxy) + c(-1,1)*1.96*sqrt(var(meanxy)/50) #95%CI for mean
mean(varxy) + c(-1,1)*1.96*sqrt(var(varxy)/50) #95%CI for variance
sqrt(1.912846)
sqrt(4.173590)
rnorm(10000, mean=c(9.27048, 10.40593), sd=c(1.383057, 4.173590))