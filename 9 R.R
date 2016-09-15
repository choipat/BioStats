# 1 ######################################## 1 ######################################### 1

data(golub,package="multtest") 
grep("GRO2 GRO2",golub.gnames[,2])
grep("GRO3 GRO3",golub.gnames[,2])
GRO2data <- golub[2714,]
GRO3data <- golub[2715,]

# 1a ####################################### 1a ######################################## 1a

cor.results <- cor.test(GRO2data, GRO3data)
cor.results
# 0.7966283

# 1b ####################################### 1b ######################################## 1b

cor.results.conf <- cor.test(GRO2data, GRO3data, conf.level = 0.90)
cor.results.conf
# 0.6702984 0.8780861

# 1c ####################################### 1c ######################################## 1c

nboot <- 2000
boot.cor <- matrix(0, nrow=nboot, ncol = 1) 
data <- cbind(GRO2data, GRO3data)
for (i in 1:nboot){
  dat.star <- data[sample(1:nrow(data),replace=TRUE),] 
  boot.cor[i,] <- cor(dat.star[,1], dat.star[,2])
}
quantile(boot.cor[,1],c(0.05,0.95))

# 1d ####################################### 1d ######################################## 1d

nboot <- 2000
boot.cor <- matrix(0, nrow=nboot, ncol = 1) 
data <- cbind(GRO2data, GRO3data)
for (i in 1:nboot){
  dat.star <- data[sample(1:nrow(data),replace=TRUE),] 
  boot.cor[i,] <- cor(dat.star[,1], dat.star[,2])
}
quantile(boot.cor[,1],c(0.025,0.975))


# 2a ####################################### 2a ######################################## 2a

data(golub,package="multtest") 
grep("Zyxin",golub.gnames[,2])

Zyxin_data <- golub[2124,]
cor.data <- apply(golub, 1, function(x) cor.test(x, Zyxin_data)$estimate)

results.2a <- sum(cor.data < -0.5)
results.2a
# 85

# 2b ####################################### 2b ######################################## 2b

order.cor <- order(cor.data, decreasing=FALSE)
results.gnames <- golub.gnames[order.cor[1:5],2]
results.gnames

# 2c ####################################### 2c ######################################## 2c

cor.ttest <- apply(golub, 1, function(x) cor.test(x, Zyxin_data, alternative = "l")$p.value)
results.2c <- sum(cor.ttest < 0.05)
results.2c
# 572 p less than 0.05. That is, reject Null hyp that cor > 0. That is, accept cor < 0

cor.fdr <- p.adjust(p = cor.ttest, method = "fdr")
results.fdr.2c <- sum(cor.fdr < 0.05)
results.fdr.2c
#142


# 3a ####################################### 3a ######################################## 3a

data(golub,package="multtest") 

GRO2data <- golub[2714,]
GRO3data <- golub[2715,]

reg.fit <- lm(GRO3data ~ GRO2data)
reg.fit  
summary(reg.fit)
# as both p's are < 0.05, reject Null hypothesis (Beta = 0). Therefore, there is stat. sign. relationship b/w gro2,3
# Proportion = R squared = 0.6346

# 3b ####################################### 3b ######################################## 3b

confint(reg.fit, level = 0.95)
# Interval for slope parameter = 0.266 0.45; 0.3582 is within it.

# 3c ####################################### 3c ######################################## 3c

reg.fit <- lm(GRO3data ~ GRO2data)
newdata <- data.frame(GRO2data = 0)
predict(reg.fit, newdata, interval="prediction", level = 0.80)

# 3d ####################################### 3d ######################################## 3d

# Normal Assumption:

plot(reg.fit, which = 2)

# qq line seems to be normal.

# Confirmation:
shapiro.test(resid(reg.fit))
# Accept Null of normal dist.

plot(reg.fit,which = 1)
# Confirmation bptest(reg.fit)
# non-linear mean patterns
# Seems to violate Variance


# 4a ####################################### 4a ######################################## 4a

data(stackloss)
str(stackloss)

stackloss.data <- as.data.frame(stackloss[,c('Air.Flow', 'Water.Temp', 'Acid.Conc.', 'stack.loss')])
lin.reg <- lm(stack.loss~Air.Flow+Water.Temp+Acid.Conc., data=stackloss.data)
summary(lin.reg)
# Linear regression equation: -39.92 + 0.72 Air.flow + 1.3 Water.Temp -0.15 Acid Conc.

# 4b ####################################### 4b ######################################## 4b

# Air.FLow and Water.Temp are stat. sign. as p <0.05
# Proportion of var explained: 0.9136 or 91.36%

# 4c ####################################### 4c ######################################## 4c

given.data <- data.frame(Air.Flow = 60, Water.Temp = 20, Acid.Conc. = 90)
predict(lin.reg, given.data, interval="confidence", level = 0.90)
predict(lin.reg, given.data, interval="prediction", level = 0.90)

# END ###################################### END ####################################### END