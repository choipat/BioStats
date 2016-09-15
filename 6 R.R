#####1a#####

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CD33",golub.gnames[,2])

#One-sided test:
t.test(golub[808,gol.fac=="ALL"],mu=-0.5, alternative = "greater")

#####1b#####

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CD33",golub.gnames[,2])

t.test(golub[808, gol.fac=="ALL"], golub[808, gol.fac=="AML"] )

#####1c#####

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CD33",golub.gnames[,2])
grep("CCND3 Cyclin D3",golub.gnames[,2])

t.test(golub[808, gol.fac=="ALL"], golub[1042, gol.fac=="ALL"],alternative = "less" )

#####1d#####

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CD33",golub.gnames[,2])
grep("CCND3 Cyclin D3",golub.gnames[,2])

ALLnum<-length(golub[808, gol.fac=="ALL"])
AMLnum<-length(golub[808, gol.fac=="AML"])

CD33.ALL<-golub[808, gol.fac=="ALL"]
CD33.ALL.gt<-sum(CD33.ALL>-.5)
prop.test(CD33.ALL.gt,ALLnum, alternative="greater")

#####1e#####if 1402#####

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CD33",golub.gnames[,2])
grep("CCND3 Cyclin D3",golub.gnames[,2])

ALLnum<-length(golub[808, gol.fac=="ALL"])
AMLnum<-length(golub[808, gol.fac=="AML"])

CD33.ALL<-golub[808, gol.fac=="ALL"]
CD33.ALL.gt<-sum(CD33.ALL>-.5)

CCND3.ALL<-golub[1402, gol.fac=="ALL"]
CCND3.ALL.gt<-sum(CCND3.ALL>-.5)

prop.test(x=c(CD33.ALL.gt,CCND3.ALL.gt), n=c(ALLnum,ALLnum), alternative="two.sided")

#####1e#####if 1042#####

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CD33",golub.gnames[,2])
grep("CCND3 Cyclin D3",golub.gnames[,2])

ALLnum<-length(golub[808, gol.fac=="ALL"])
AMLnum<-length(golub[808, gol.fac=="AML"])

CD33.ALL<-golub[808, gol.fac=="ALL"]
CD33.ALL.gt<-sum(CD33.ALL>-.5)

CCND3.ALL<-golub[1042, gol.fac=="ALL"]
CCND3.ALL.gt<-sum(CCND3.ALL>-.5)

prop.test(x=c(CD33.ALL.gt,CCND3.ALL.gt), n=c(ALLnum,ALLnum), alternative="two.sided")

#####1f#####

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CD33",golub.gnames[,2])
grep("CCND3 Cyclin D3",golub.gnames[,2])

ALLnum<-length(golub[808, gol.fac=="ALL"])
AMLnum<-length(golub[808, gol.fac=="AML"])

CD33.ALL<-golub[808, gol.fac=="ALL"]
CD33.ALL.gt<-sum(CD33.ALL>-.5)

CD33.AML<-golub[808, gol.fac=="AML"]
CD33.AML.gt<-sum(CD33.AML>-.5)

prop.test(x=c(CD33.ALL.gt,CD33.AML.gt), n=c(ALLnum,AMLnum), alternative="two.sided")

#####2b#####

pbinom(19,size=1000,prob=0.05)

#####3a#####

x.sim<-matrix(rnorm(10000*20, mean=3, sd=4), ncol=20)
tstat<-function(x) (mean(x)-3)/sd(x)*sqrt(length(x))
tstat.sim<-apply(x.sim,1,tstat) #t-test statistic for each data set
power.sim1<-mean(tstat.sim<qt(0.3,df=19))
power.sim2<-mean(tstat.sim<qt(0.4,df=19))
power.sim<-mean(power.sim2-power.sim1) #Calculate rejection rate
power.sim # Display rejection rate

#95%
power.sim+c(-1,0,1)*qnorm(0.975)*sqrt(power.sim*(1-power.sim)/10000)

#####4a#####

p.values <- apply(golub, 1, function(x) t.test(x ~ gol.fac)$p.value)
p.bon<-p.adjust(p=p.values,method="bonferroni")
p.fdr<-p.adjust(p=p.values,method="fdr")
sum(pt<0.05)
sum(p.bon<0.05)
sum(p.fdr<0.05)

#####4b#####

p.values <- apply(golub, 1, function(x) t.test(x ~ gol.fac)$p.value)
p.fdr<-p.adjust(p=p.values,method="fdr")
orderAML<-order(p.fdr, decreasing=FALSE)
golub.gnames[orderAML[1:3],2]

#####5a#####

#####Wald CI#####

wald.CI <- function(x,n,conf)
{
  k  <- qnorm(1-(1-conf)/2)
  p <- x/n
  sd <- sqrt((x+1)-(x)^2/(n+2))/(n+2)
  low <- p - k*sd
  up <- p + k*sd
  c(low, up)
}

#####Wilson CI#####

wilson.CI <-
  function(x, n, conf) {
    k <- qnorm((1/2)*(1 + conf))
    k1 <- k*sqrt((x*(n-x))/(n^2*n) + k^2/(4*n^2))
    low <- (n/(n+k^2))*(x/n + k^2/(2*n) - k1)
    up <- (n/(n+k^2))*(x/n + k^2/(2*n) + k1)
    c(low,up)
  }

#####Agresti-Coull CI#####

agresti.CI <- function(x,n,conf){
  k <- abs(qnorm((1-conf)/2))
  p <- (x+((k^2)/2))/(n+(k/2))
  stderr <- sqrt(p * (1-p)/(n+(k/2)))
  up <- p + k * stderr
  low <- p - k * stderr
  c(low,up)
}