# 1a ####################################### 1a ######################################## 1a

library("ALL"); data(ALL);data <- exprs(ALL)

ALL.fac <- factor(ALL$BT %in% c("B","B1","B2","B3","B4"), labels=c("1","2"))

# 1b ####################################### 1b ######################################## 1b

gene_1 <- data[1,]; gene_2 <- data[2,]; gene_3 <- data[3,]
par(mfrow=c(1,3))
hist(gene_1, main = "Gene 1"); hist(gene_2, main = "Gene 2"); hist(gene_3, main = "Gene 3")


# 1c ####################################### 1c ######################################## 1c

gene_4 <- data[4,]; gene_5 <- data[5,]
pairs(cbind(gene_1, gene_2, gene_3, gene_4, gene_5))

# 1d ####################################### 1d ######################################## 1d

install.packages("scatterplot3d")
require(scatterplot3d)
par(mfrow=c(1,1))
d <- rbind(exprs(ALL[c("39317_at","32649_at","481_at"),]))
scatterplot3d(t(d),color=ALL.fac)

# 1e ####################################### 1e ######################################## 1e

cl_1 <- kmeans(t(d),centers=2,nstart=10)
table(ALL.fac,cl_1$cluster)
cl_2 <- kmeans(t(d),centers=3,nstart=10)
table(ALL.fac,cl_2$cluster)

# 1f ####################################### 1f ######################################## 1f

pr.ALL <- prcomp(data, scale=TRUE)
summary(pr.ALL)
#93.6% and 0.95%

# 1g ####################################### 1g ######################################## 1g

par(mfrow=c(1,1))
biplot(pr.ALL, cex=0.5)
# PC1 is the average of patients

# 1h ####################################### 1h ######################################## 1h

o <- order(pr.ALL$x[,2], decreasing=T) 
dimnames(data)[[1]][[o[1]]]; dimnames(data)[[1]][[o[2]]]; dimnames(data)[[1]][[o[3]]]
dimnames(data)[[1]][[o[12623]]]; dimnames(data)[[1]][[o[12624]]]; dimnames(data)[[1]][[o[12625]]]

# 1i ####################################### 1i ######################################## 1i

annotation(ALL)
biocLite("hgu95av2.db"); biocLite("annotation")
library(help=hgu95av2.db)
ChrNrOfProbe <- as.list(hgu95av2CHR)
GNOfProbe <- as.list(hgu95av2GENENAME)
ChrNrOfProbe[o[1]]; GNOfProbe[o[1]] #biggest
ChrNrOfProbe[o[12625]]; GNOfProbe[o[12625]] #smallest


# 2a ####################################### 2a ######################################## 2a

iris_data <- iris[1:4]
mean <- mean(iris_data[,1])
sd <- sd(iris_data[,1])
Sepal.Length <- NULL
for (i in 1:150){Sepal.Length[i] <- (iris_data[i,1]-mean)/sd}

mean <- mean(iris_data[,2])
sd <- sd(iris_data[,2])
Sepal.Width <- NULL
for (i in 1:150){Sepal.Width[i] <- (iris_data[i,2]-mean)/sd}

mean <- mean(iris_data[,3])
sd <- sd(iris_data[,3])
Petal.Length <- NULL
for (i in 1:150){Petal.Length[i] <- (iris_data[i,3]-mean)/sd}

mean <- mean(iris_data[,4])
sd <- sd(iris_data[,4])
Petal.Width <- NULL
for (i in 1:150){Petal.Width[i] <- (iris_data[i,4]-mean)/sd}

scaled_data <- cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
xxx=data.frame(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
# 2b ####################################### 2b ######################################## 2b

cor_scaled <- cor(scaled_data)
cor_unscaled <- cor(iris_data)
cor_scaled;cor_unscaled
all.equal(cor_scaled, cor_unscaled)

# 2c ####################################### 2c ######################################## 2c

eucl_dist <- dist(scaled_data, method="euclidean")
eucl_dist^2
1-cor_scaled

# 2d ####################################### 2d ######################################## 2d

pca_unscaled <- prcomp(iris_data, scale=FALSE)
pca_unscaled
pca_scaled <- prcomp(scaled_data, scale=FALSE)
pca_scaled

# 2e ####################################### 2e ######################################## 2e

summary(pca_unscaled) # 97.77%
summary(pca_scaled) # 95.81%

# 2f ####################################### 2f ######################################## 2f

data2 <- iris_data; p <- ncol(data2); n <- nrow(data2) ; nboot<-1000
sdevs <- array(dim=c(nboot,p))
for (i in 1:nboot) {
  dat.star <- data2[sample(1:n,replace=TRUE),]
  sdevs[i,] <- prcomp(dat.star)$sdev
}
quantile(sdevs[,2], c(0.05,0.95))
