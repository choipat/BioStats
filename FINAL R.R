getwd()
y <- as.numeric(t(read.table(file="DataPois.txt", header=TRUE)))

######################################## 1a ######################################## 1a

str(y) #120
mean(y) #1.816667

######################################## 1b ######################################## 1b

nlx <- function(x) -sum(log(dpois(y,x)))
optim(1, nlx) # 1.816797
log(1.816797)/log(2.71828)

######################################## 1c ######################################## 1c

n<-length(y)
nboot<-1000
boot.xbar <- rep(NA, nboot)
for (i in 1:nboot) {
  data.star <- y[sample(1:n,replace=TRUE)]
  boot.xbar[i]<-mean(data.star)
}
quantile(boot.xbar,c(0.025,0.975))

######################################## 2a ######################################## 2a

library(ISLR)
ncidata <- ncidatanew <- NCI60$data
ncilabs <- ncilabsnew <- NCI60$labs
filt <- NULL;k=1;
for (i in 1:64){
  if (sum(ncilabs==(ncilabs[i])) <= 3){
    filt[k] <- i;k=k+1;  
  }
}
ncidatanew <- ncidatanew[-filt,]
ncilabsnew <- ncilabsnew[-filt]

######################################## 2b ######################################## 2b

firstgene <- ncidatanew[,1]
anova(lm(firstgene ~ ncilabsnew))
summary(lm(firstgene ~ ncilabsnew))
pairwise.t.test(firstgene, ncilabsnew, p.adjust.method='fdr')

######################################## 2c ######################################## 2c

library(lmtest)
shapiro.test(residuals(lm(firstgene ~ ncilabsnew))) # It is normal
bptest(lm(firstgene ~ ncilabsnew), studentize = FALSE) # Accept homoscedasticity (Variance)

######################################## 2d ######################################## 2d

data2d = apply(ncidatanew, 2, function(x) anova(lm(x ~ ncilabsnew))$Pr[1])
p.fdr <- p.adjust(p=data2d, method="fdr")
sum(p.fdr<0.05)

######################################## 3a ######################################## 3a

data(state);head(state.x77)
state.data <- as.data.frame(state.x77[,1:8])
pairs(state.data) 
# Murder, HS Grad, Illeteracy, Frost 

######################################## 3b ######################################## 3b

state.data <- as.data.frame(state.x77[,c('Life Exp', 'Income', 'Illiteracy', 'Frost')]) 
names(state.data) <- c('Life.Exp', 'Income', 'Illiteracy', 'Frost')
lin.reg <- lm(Life.Exp~Income+Illiteracy+Frost, data=state.data)
summary(lin.reg)
# Life.Exp = 73 + 0.0002 Income - 1.6 Illiteracy -0.006 Frost
# Illiteracy

######################################## 3c ######################################## 3c

n <- dim(state.data)[1]
sqerr <- rep(NA, n)
for (i in 1:n) {
  cdata <- state.data[-i,]
  lin.reg <- lm(Life.Exp~Income+Illiteracy+Frost, data=cdata)
  coeffs <- lin.reg$coefficients
  pred <- coeffs[1] + sum(coeffs[-1] * state.data[i,-1])
  sqerr[i] <- (state.data[i,1] - pred)^2
}
err <- mean(sqerr)
err

######################################## 4a ######################################## 4a

library(ALL);data(ALL)
ALLB <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))]
data <- exprs(ALLB)
str(data)

######################################## 4b ######################################## 4b

meandata <- apply(data,1,mean)
sddata <- apply(data,1,sd)
data4b <- subset(sddata/meandata,sddata/meandata>0.2)
data.new <- data[names(data4b),]
sum(sddata/meandata>0.2)
row.names(data.new)

######################################## 4d ######################################## 4d

hc <- hclust(dist(t(data.new)),method="complete")
plot(hc, hang = -1, cex = .7, labels=ALL$BT[1:95])
rect.hclust(hc, k=5)
plot(hc, hang = -1, cex = .7, labels=ALL$mol.biol[1:95])
rect.hclust(hc, k=6)

######################################## 4e ######################################## 4e

library(gplots)
color.map <- function(T) {
  if (T=="B") "yellow"
  else if(T=="B1") "lightblue"
  else if(T=="B2") "darkolivegreen4"
  else if(T=="B3") "darkorange"
  else "purple"
}
clmy <- as.character(ALL$BT[1:95]) 
patientcolors <- unlist(lapply(clmy, color.map))
heatmap.2(data.new,
          col=greenred(75),
          scale="row",
          dendrogram="column",
          ColSideColors=patientcolors,
          margin=c(3, 12), cexRow=0.5,
          key=FALSE, trace="none", labRow=NA, labCol=NA)

color.map <- function(T) {
  if (T=="ALL1/AF4") "yellow"
  else if(T=="BCR/ABL") "lightblue"
  else if(T=="E2A/PBX1") "darkolivegreen4"
  else if(T=="NEG") "darkorange"
  else if(T=="NUP-98") "brown4"
  else "purple"
}
clmy <- as.character(ALL$mol.biol[1:95]) 
patientcolors <- unlist(lapply(clmy, color.map))
heatmap.2(data.new,
          col=greenred(75),
          scale="row",
          dendrogram="column",
          ColSideColors=patientcolors,
          margin=c(3, 12), cexRow=0.5,
          key=FALSE, trace="none", labRow=NA, labCol=NA)

######################################## 4f ######################################## 4f

levels <- ALL[,which(ALL$BT %in% c("B1","B2","B3","B4"))]$BT
levels(levels)[levels(levels)=="B3"] <- "B34";levels(levels)[levels(levels)=="B4"] <- "B34"
levels <- droplevels(levels)
library(limma)
design.ma <- model.matrix(~ 0 + factor(levels))
colnames(design.ma) <- c("B1","B2","B34") 
fit <- lmFit(ALL[,which(ALL$BT %in% c("B1","B2","B3","B4"))], design.ma)
fit <- eBayes(fit)
cont.ma <- makeContrasts(B1-B2,B2-B34,B1-B34, levels=levels)
fit1 <- contrasts.fit(fit, cont.ma)
fit1 <- eBayes(fit1)
dim(topTable(fit1, number=Inf, p.value=0.05, adjust.method="fdr"))[[1]]
a = topTable(fit1, number=Inf, p.value=0.05, adjust.method="fdr")
topTable(fit1, number=Inf, p.value=0.05, adjust.method="fdr")[1:10,]
names <- row.names(a)

######################################## 4g ######################################## 4g

library("hgu95av2.db");library(ALL);data(ALL)

######################################## 4h ######################################## 4h

sum(row.names(data.new) %in% names)
as.numeric(row.names(data.new) %in% names)
newgenes <- row.names(data.new)[row.names(data.new) %in% names]
newgenes

######################################## 5a ######################################## 5a

dat <- read.table(file="DataPoisReg.txt", header=TRUE)
nlx <- function(x) -sum(log(dpois(dat,x)))
optim(1, nlx)
