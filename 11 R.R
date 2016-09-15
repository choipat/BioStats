# 1 ######################################## 1 ######################################### 1

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("CCND3 Cyclin D3",golub.gnames[,2])

# 1a ####################################### 1a ######################################## 1a

clusdata <- golub[1042,]
# Conduct hierarchical clustering using single linkage and Ward linkage
# Plot cluster dendrogram
hclust.sing <- hclust(dist(clusdata, method="euclidian"), method="single") 
plot(hclust.sing,labels=gol.fac) 
hclust.ward <- hclust(dist(clusdata, method="euclidian"), method="ward.D2")
plot(hclust.ward,labels=gol.fac)

par(mfrow=c(1,2))
plot(hclust.sing,labels=gol.fac,cex=0.4) 
plot(hclust.ward,labels=gol.fac,cex=0.4)

# Get two clusters from each of the methods
clust.2sing <- cutree(hclust.sing,2)
clust.2ward <- cutree(hclust.ward,2)

# Use function table() to compare the clusters with the two patient groups ALL/AML
table(gol.fac, clust.2sing)
table(gol.fac, clust.2ward) # Better

# 1b ####################################### 1b ######################################## 1b

# Use k-means cluster analysis to get two clusters
clusters.2km <- kmeans(clusdata, centers=2)
table(gol.fac, clusters.2km$cluster)

# 1d ####################################### 1d ######################################## 1d

clusters.2km$centers
initial <- clusters.2km$centers
n <- length(clusters.2km$cluster)
nboot<-1000
boot.cl <- matrix(NA,nrow=nboot, ncol=2) 
for (i in 1:nboot){
  dat.star <- clusdata[sample(1:n,replace=TRUE)]
  cl <- kmeans(dat.star, centers=initial)
  boot.cl[i,] <- cl$centers
}
apply(boot.cl,2,mean)
quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))

# 1e ####################################### 1e ######################################## 1e

K <- (1:30)
sse.1e <- rep(NA,length(K)) 
for (k in K) {sse.1e[k] <- kmeans(clusdata, centers=k, nstart = 10)$tot.withinss}
plot(K, sse.1e, type='o', xaxt='n'); axis(1,at = K, las=2)


# 2 ######################################## 2 ######################################### 2

data(golub, package = "multtest")

# 2a ####################################### 2a ######################################## 2a

golub.onco <- agrep("^oncogene",golub.gnames[,2])
golub.angn <- agrep("^antigen",golub.gnames[,2])
golub.oncoangn <- unique(c(golub.onco, golub.angn))
gene.factor <- factor(c(rep("oncogene",length(golub.onco)),rep("antigen",length(golub.angn))))

# 2b ####################################### 2b ######################################## 2b

# Do the following four clustering analysis with K=2: 
# two clustering methods (K-means and K-Medoids) each using 
# two dissimilarity measures (Euclidean distance and 1-correlation)
data <- data.frame(golub[golub.oncoangn,])
clust.1 <- kmeans(dist(data, method='eucl'), centers=2) # k-means, Euclidean
clust.2 <- kmeans(as.dist(1-cor(t(data))),centers=2,nstart=10) # k-means, 1-correlation 
clust.3 <- pam(dist(data,method='eucl'),k=2) # k-medoids, Euclidean
clust.4 <- pam(as.dist(1-cor(t(data))),k=2) # k-medoids, 1-correlation

# Use table() to compare the resulting two clusters with the 
# two gene groups oncogenes and antigens for each of the four analysis
tab1 <- table(gene.factor, clust.1$cluster) # k-means, Euclidean
tab2 <- table(gene.factor, clust.2$cluster) # k-means, 1-correlation 
tab3 <- table(gene.factor, clust.3$cluster) # k-medoids, Euclidean
tab4 <- table(gene.factor, clust.4$cluster) # k-medoids, 1-correlation
tab1; tab2; tab3; tab4 # Print tables

# 2c ####################################### 2c ######################################## 2c

# Marginal Independence Test
chisq.test(tab1) # k-means, Euclidean
chisq.test(tab2) # k-means, 1-correlation 
chisq.test(tab3) # k-medoids, Euclidean
chisq.test(tab4) # k-medoids, 1-correlation

# 2d ####################################### 2d ######################################## 2d

hclust.sing.2 <- hclust(dist(data, method="euclidian"), method="single") 
plot(hclust.sing.2,labels=gene.factor, cex=0.62) # Eucledian, Single
hclust.comp.2 <- hclust(dist(data, method="euclidian"), method="complete")
plot(hclust.comp.2,labels=gene.factor, cex=0.62) # Eucledian, Complete 


# 3 ######################################## 3 ######################################### 3

install.packages('ISLR')
library(ISLR)
ncidata<-NCI60$data
ncilabs<-NCI60$labs

# 3a ####################################### 3a ######################################## 3a

#  Produce a plot of K versus SSE, for K=1,..., 30
K <- (1:30) 
sse <- rep(NA,length(K)) 
for (k in K) {sse[k]<-kmeans(ncidata, centers=k, nstart = 10)$tot.withinss}
plot(K, sse, type='o', xaxt='n'); axis(1,at = K, las=2)

# 3b ####################################### 3b ######################################## 3b

# Do K-means clustering (K=7) with 1-correlation as the dissimilarity measure on the data. 
# Compare the clusters with the cell lines.
clusters.7km <- kmeans(as.dist(1-cor(t(ncidata))),centers=7,nstart=10)
table(factor(ncilabs), clusters.7km$cluster)

# END ###################################### END ####################################### END