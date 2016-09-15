################## 1b Analytical and numerical MLE ############### 1b

#Analytical MLE
data <- c(1.636, 0.374, 0.534, 3.015, 0.932, 0.179)
mean(data)
1/mean(data)

#Numerical MLE
numMLE<-function(x) -sum(log(dexp(c(1.636, 0.374, 0.534, 3.015, 0.932, 0.179),x)))
optim(1, numMLE)

######## 2b One-sided 90% lower confidence interval of m ########## 2b

Mean<-100.8
sd<-12.4
lc<-Mean+(qt(0.1,52)*(12.4/sqrt(53)))
lc

############### 3a: ALL: Bootstrap 95% CIs for mean ############### 3a

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="ALL"]
n<-length(Zyxin)
nboot<-1000
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin[sample(1:n,replace=TRUE)]
  boot.xbar[i]<-mean(data.star)
}
ci.ALL<-quantile(boot.xbar,c(0.025,0.975))
ci.ALL

############### 3a: AML: Bootstrap 95% CIs for mean ############### 3a

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="AML"]
n<-length(Zyxin)
nboot<-1000
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin[sample(1:n,replace=TRUE)]
  boot.xbar[i]<-mean(data.star)
}
ci.AML<-quantile(boot.xbar,c(0.025,0.975))
ci.AML

############### 3a: ALL: Bootstrap 95% CIs for variance ############### 3a

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="ALL"]
n<-length(Zyxin)
nboot<-1000
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin[sample(1:n,replace=TRUE)]
  boot.xbar[i]<-var(data.star)
}
ci.ALLVar<-quantile(boot.xbar,c(0.025,0.975))
ci.ALLVar

############### 3a: AML: Bootstrap 95% CIs for variance ############### 3a

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-Zyxin<-golub[2124,gol.fac=="AML"]
n<-length(Zyxin)
nboot<-1000
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin[sample(1:n,replace=TRUE)]
  boot.xbar[i]<-var(data.star)
}
ci.AMLVar<-quantile(boot.xbar,c(0.025,0.975))
ci.AMLVar

############### 3b: ALL: Parametric 95% CIs for mean ############### 3b

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="ALL"]
n<-length(Zyxin)
ci.meanALL<-mean(Zyxin)+qt(c(0.025,0.975),df=n-1)*sd(Zyxin)/sqrt(n)
ci.meanALL

############### 3b: AML: Parametric 95% CIs for mean ############### 3b

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="AML"]
n<-length(Zyxin)
ci.meanAML<-mean(Zyxin)+qt(c(0.025,0.975),df=n-1)*sd(Zyxin)/sqrt(n)
ci.meanAML

############### 3b: ALL: Parametric 95% CIs for variance ############### 3b

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="ALL"]
n<-length(Zyxin)
ci.varALL<-(n-1)*sd(Zyxin)^2/qchisq(c(0.975,0.025),df=n-1)
ci.varALL
############### 3b: AML: Parametric 95% CIs for variance ############### 3b

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="AML"]
n<-length(Zyxin)
ci.varAML<-(n-1)*sd(Zyxin)^2/qchisq(c(0.975,0.025),df=n-1)
ci.varAML

############### 3c: ALL: Bootstrap 95% CIs for median ############### 3c

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="ALL"]
n<-length(Zyxin)
nboot<-1000
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin[sample(1:n,replace=TRUE)]
  boot.xbar[i]<-median(data.star)
}
ci.ALLMed<-quantile(boot.xbar,c(0.025,0.975))
ci.ALLMed

############### 3c: AML: Bootstrap 95% CIs for median ############### 3c

data(golub, package = "multtest")
grep("Zyxin",golub.gnames[,2])
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

Zyxin<-golub[2124,gol.fac=="AML"]
n<-length(Zyxin)
nboot<-1000
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- Zyxin[sample(1:n,replace=TRUE)]
  boot.xbar[i]<-median(data.star)
}
ci.AMLMed<-quantile(boot.xbar,c(0.025,0.975))
ci.AMLMed


################################ 4a ########################### 4a

########### Formula -1 ###########

#Number of simulations
nsim <- 1000
lambda <-
  
#Simulated values and their mean
data <- matrix(rpois(50*nsim,10),nrow=nsim)
newlambda <- (apply(data,1,mean))

#4.1
tdist <- qt(.05,49) * sqrt(newlambda/50)

#90% CI
low = newlambda+tdist
High = newlambda-tdist

#Check coverage probabilities
sum(low<lambda & lambda<High)/1000

########### Formula -2 ###########

#Number of simulations
nsim <- 1000
lambda <-
  
#Simulated values and their mean
data <- matrix(rpois(50*nsim,10),nrow=nsim)
newlambda <- (apply(data,1,mean))

#90% CI
low = 49 *(newlambda)/qchisq(.95,49)
High = 49 *(newlambda)/qchisq(.05,49)

#Check coverage probabilities
sum(low<lambda & lambda<High)/1000


################################ 4b lambda = 0.1 ########################### 4b

########### Formula -1 # ##########

#Number of simulations
nsim<-1000
lambda<-0.1
  
#Simulated values and their mean
data<-matrix(rpois(50*nsim,lambda),nrow=nsim)
newlambda<-(apply(data,1,mean))

#4.1
tdist<- qt(.05,49) * sqrt(newlambda/50)

#90% CI
low=newlambda+tdist
High=newlambda-tdist

#Check coverage probabilities
sum(low<lambda & lambda<High)/1000

########### Formula -2 ###########

#Number of simulations
nsim <- 1000
lambda <- 0.1
  
#Simulated values and their mean
data <- matrix(rpois(50*nsim,lambda),nrow=nsim)
newlambda <- (apply(data,1,mean))

#90% CI
low = 49 *(newlambda)/qchisq(.95,49)
High = 49 *(newlambda)/qchisq(.05,49)

#Check coverage probabilities
sum(low<lambda & lambda<High)/1000

################################ 4b lambda = 1 ########################### 4b

########### Formula -1  ###########

#Number of simulations
nsim<-1000
lambda<-1
  
#Simulated values and their mean
data<-matrix(rpois(50*nsim,lambda),nrow=nsim)
newlambda<-(apply(data,1,mean))

#4.1
tdist<- qt(.05,49) * sqrt(newlambda/50)

#90% CI
low=newlambda+tdist
High=newlambda-tdist

#Check coverage probabilities
sum(low<lambda & lambda<High)/1000

########### Formula -2 ###########

#Number of simulations
nsim <- 1000
lambda <- 1
  
#Simulated values and their mean
data <- matrix(rpois(50*nsim,lambda),nrow=nsim)
newlambda <- (apply(data,1,mean))

#90% CI
low = 49 *(newlambda)/qchisq(.95,49)
High = 49 *(newlambda)/qchisq(.05,49)

#Check coverage probabilities
sum(low<lambda & lambda<High)/1000

################################ 4b lambda = 10 ########################### 4b

########### Formula -1 ###########

#Number of simulations
nsim<-1000
lambda<-10
  
#Simulated values and their mean
data<-matrix(rpois(50*nsim,lambda),nrow=nsim)
newlambda<-(apply(data,1,mean))

#4.1
tdist<- qt(.05,49) * sqrt(newlambda/50)

#90% CI
low=newlambda+tdist
High=newlambda-tdist

#Check coverage probabilities
sum(low<lambda & lambda<High)/1000

########### Formula -2 ###########

#Number of simulations
nsim <- 1000
lambda <- 10
  
#Simulated values and their mean
data <- matrix(rpois(50*nsim,lambda),nrow=nsim)
newlambda <- (apply(data,1,mean))

#90% CI
low = 49 *(newlambda)/qchisq(.95,49)
High = 49 *(newlambda)/qchisq(.05,49)

#Check coverage probabilities
sum(low<lambda & lambda<High)/1000