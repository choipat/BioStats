library("ALL"); data(ALL); install.packages(c('rpart','ROCR','VGAM')); library("hgu95av2.db"); library(ROCR) 

# 1a ####################################### 1a ######################################## 1a

IsB <- factor(ALL$BT %in% c("B","B1","B2","B3","B4"))

# 1b ####################################### 1b ######################################## 1b

names <- c("39317_at", "38018_g_at")
symb <- mget(names, env = hgu95av2SYMBOL)
ALLBTnames <- ALL[names, ]
probedat <- as.matrix(exprs(ALLBTnames))
row.names(probedat) <- unlist(symb)

require(rpart)
B.stage <- factor(IsB)
c.tr <- rpart(B.stage~., data = data.frame(t(probedat)))
plot(c.tr, branch=0,margin=0.1)
text(c.tr, digits=3, cex=0.7)
rpartpred <- predict(c.tr, type="class")
table(rpartpred, B.stage)

# ROC Curve

library(ROCR)
pred.prob <- predict(c.tr,type="prob")[,2]
pred <- prediction(pred.prob, IsB=="TRUE")
perf <- performance(pred,"tpr","fpr")
plot(perf, col="green")

# 1c ####################################### 1c ######################################## 1c

# misclassification = 
(11+2)/128
# fpr =
2/(31+2)
# fnr =
11/(84+11)
# Specificity = tnr =
31/(2+31)
# AUC = 0.922807
performance(pred,"auc")


# 1d ####################################### 1d ######################################## 1d

n <- length(rpartpred)
index <- 1:n
K <- 10
folds <- createFolds(index, k=K)
fnr.cv.raw <- rep(NA, K)
for (i in 1:K) {
  testID <- folds[[i]]
  c.tr <- rpart(IsB[-testID]~., data=data.frame(t(expr.data)[-testID,]))
  tr.pred <- predict(c.tr, newdata=data.frame(t(expr.data)[-testID,]), type="class")
  fnr.cv.raw[i] <- mean(tr.pred[testID] == 'TRUE' & IsB[testID] == 'TRUE')/mean(IsB[testID] == 'TRUE') 
}
fnr.cv <- mean(fnr.cv.raw)
fnr.cv


# 1e ####################################### 1e ######################################## 1e

prob.name <- c("39317_at", "38018_g_at")
expr.data <- exprs(ALL)[prob.name,]
data.lgr <- data.frame(IsB, t(expr.data))
fit.lgr <- glm(IsB~., family=binomial(link='logit'), data=data.lgr)
pred.prob <- predict(fit.lgr, data=data.lgr$expr.data, type="response")
pred.B1 <- factor(pred.prob> 0.5, levels=c(TRUE,FALSE), labels=c("B","not_B")) 
IsB_1 <- factor(IsB, levels=c(TRUE,FALSE), labels=c("B","not_B")) 
table(pred.B1, IsB_1) 

#CI

y <- as.numeric(IsB==TRUE) 
data.CI <- exprs(ALL)["39317_at", ]
data <- glm(data.frame(y, data.CI))
confint(data, level=0.8)

# 1f ####################################### 1f ######################################## 1f

install.packages('caret'); require(caret);
ALLB <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))]
data.lgr <- data.frame(IsB,t(expr.data))
n <- dim(data.lgr)[1]
index <- 1:n 
K <- 10
flds <- createFolds(index, k=K)
mcr.cv.raw <- rep(NA, K)
for (i in 1:K) {
  testID <- flds[[i]] 
  data.tr <- data.lgr[-testID,]
  data.test <- data.lgr[testID,] 
  fit.lgr <- glm(IsB~., family=binomial(link='logit'), data=data.tr)
  pred.prob <- predict(fit.lgr, newdata=data.test, type="response")
  pred <- (pred.prob> 0.5)
  mcr.cv.raw[i] <- sum(pred!=data.test$IsB)/length(pred)
}
mcr.cv <- mean(mcr.cv.raw) 
mcr.cv 

# 1g ####################################### 1g ######################################## 1g

pca.all <- prcomp(t(exprs(ALL)), scale=TRUE)
summary(pca.all)
PropVar <- summary(pca.all)$importance[2,] 
plot(1:length(PropVar), PropVar, xlab='number of principal components', ylab='proportion of variance explained',cex=0.3) 

# 1h ####################################### 1h ######################################## 1h

install.packages('e1071');library(e1071)
data.pca <- pca.all$x[,1:5]
svmest <- svm(IsB~data.pca,type="C-classification",kernal="linear")
svmpred <- predict(svmest,data.pca)
tpr.svm <- mean(svmpred==IsB&IsB==TRUE)/(mean(svmpred==IsB&IsB==TRUE)+mean(svmpred!=IsB&IsB==TRUE))
tpr.svm

# 1i ####################################### 1i ######################################## 1i

n <- length(IsB) 
mcr.cv.raw <-rep (NA,n) 
for (i in 1:n) {
  svmest <- svm(data.pca[-i,],IsB[-i],type="C-classification",kernal="linear") 
  svmpred <- predict(svmest,t(data.pca[i,])) 
  mcr.cv.raw[i] <-mean(svmpred!=IsB[i])
}
mcr.cv <- mean(mcr.cv.raw)
mcr.cv

# 2 ######################################## 2 ######################################### 2

# K=1
pca.iris <- prcomp(iris[,1:4], scale=TRUE)
Species <- iris$Species
data.pca <- pca.iris$x[,1,drop = F]
n <- length(Species)
iris2 <- data.frame(Species, data.pca)
# Logistic Regression
iris2.lgr <- vglm(Species~., family=multinomial, data=iris2) #Logistic Regression
pred.prob <- predict(iris2.lgr, iris2[,-1,drop=F], type="response") #get prediction
pred.lgr <- apply(pred.prob, 1, which.max) #Assign to the class with largest
pred.lgr <- factor(pred.lgr, levels=c("1","2","3"), labels=levels(iris2$Species)) #relabel 1,2,3 by species names
mcr.lgr <- mean(pred.lgr!=iris2$Species) #misclassification rate
### leave-one-out cross validation
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation
for (i in 1:n) {
  fit.lgr <- vglm(Species~., family=multinomial, data=iris2[-i,]) #fit logistic
  pred.prob <- predict(fit.lgr, iris2[i,-1, drop=F], type="response") #get prediction probability
  pred <- apply(pred.prob, 1, which.max) #Assign class
  pred <- factor(pred, levels=c("1","2","3"), labels=levels(iris2$Species))
  mcr.cv.raw[i] <- mean(pred!=Species[i]) #check misclassification
}
mcr.cv <- mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.lgr, mcr.cv)

# SVM
iris2.svm <- svm(data.pca, Species, type = "C-classification", kernel = "linear") #train SVM
svmpred <- predict(iris2.svm , data.pca) #get SVM prediction.
mcr.svm<- mean(svmpred!=Species) #misclassification rate
### leave-one-out cross validation
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation
for (i in 1:n) {
  svmest <- svm(data.pca[-i,], Species[-i], type = "C-classification", kernel =
                  "linear") #train SVM without i-th observation
  svmpred <- predict(svmest, t(data.pca[i,])) #predict i-th observation. Here
  mcr.cv.raw[i]<- mean(svmpred!=Species[i]) #misclassification rate
}
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.svm, mcr.cv)

# Classification Tree
fit <- rpart(Species ~ ., data = iris2, method = "class")
pred.tr<-predict(fit, iris2, type = "class")
mcr.tr <- mean(pred.tr!=Species) #misclassification rate
### leave-one-out cross validation
mcr.cv.raw <- rep(NA, n) #A vector to save mcr validation
for (i in 1:n) {
  fit.tr <- rpart(Species ~ ., data = iris2[-i,], method = "class") #train the tree without i-th observation
  pred <- predict(fit.tr, iris2[i,], type = "class")#use tree to predict i-th observation class
  mcr.cv.raw[i] <- mean(pred!=Species[i]) #check misclassifion
}
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
c(mcr.tr, mcr.cv)

# K =2,3,4
data(iris)
pca.iris <- prcomp(iris[,1:4], scale=TRUE)
Species <- iris$Species
Analysis <- function(k){
  print(c("When K=",k))
  data.pca <- pca.iris$x[,1:k]
  n <- length(Species)
  iris2 <- data.frame(Species, data.pca)
# Logistic Regression
iris2.lgr <- vglm(Species~., family=multinomial, data=iris2) #Logistic Regression
pred.prob <- predict(iris2.lgr, iris2[,-1], type="response") #get prediction
pred.lgr <- apply(pred.prob, 1, which.max) #Assign to the class with largest
pred.lgr <- factor(pred.lgr, levels=c("1","2","3"), labels=levels(iris2$Species)) #relabel 1,2,3 by species names
mcr.lgr <- mean(pred.lgr!=iris2$Species) #misclassification rate
### leave-one-out cross validation
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation
for (i in 1:n) {
  fit.lgr <- vglm(Species~., family=multinomial, data=iris2[-i,]) #fit logistic
  pred.prob <- predict(fit.lgr, iris2[i,-1], type="response") #get prediction probability
  pred <- apply(pred.prob, 1, which.max) #Assign class
  pred <- factor(pred, levels=c("1","2","3"), labels=levels(iris2$Species))
  #relabel 1,2,3 by species names
  mcr.cv.raw[i] <- mean(pred!=Species[i]) #check misclassification
}
mcr.cv <- mean(mcr.cv.raw) #average the mcr over all n rounds.
print(c("For Logistic regression: emperical mcr, leave one out mcr", mcr.lgr, mcr.cv))

# SVM
iris2.svm <- svm(data.pca, Species, type = "C-classification", kernel = "linear") #train SVM
svmpred <- predict(iris2.svm , data.pca) #get SVM prediction.
mcr.svm<- mean(svmpred!=Species) #misclassification rate
### leave-one-out cross validation
mcr.cv.raw<-rep(NA, n) #A vector to save mcr validation
for (i in 1:n) {
  svmest <- svm(data.pca[-i,], Species[-i], type = "C-classification", kernel =
                  "linear") #train SVM without i-th observation
  svmpred <- predict(svmest, t(data.pca[i,])) #predict i-th observation. Here
  mcr.cv.raw[i]<- mean(svmpred!=Species[i]) #misclassification rate
}
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
print(c("For SVM: emperical mcr, leave one out mcr", c(mcr.svm, mcr.cv)))

# Classification Tree
fit <- rpart(Species ~ ., data = iris2, method = "class")
pred.tr<-predict(fit, iris2, type = "class")
mcr.tr <- mean(pred.tr!=Species) #misclassification rate
### leave-one-out cross validation
mcr.cv.raw <- rep(NA, n) #A vector to save mcr validation
for (i in 1:n) {
  fit.tr <- rpart(Species ~ ., data = iris2[-i,], method = "class") #train the tree without i-th observation
  pred <- predict(fit.tr, iris2[i,], type = "class")#use tree to predict i-th observation class
  mcr.cv.raw[i] <- mean(pred!=Species[i]) #check misclassifion
  }
mcr.cv<-mean(mcr.cv.raw) #average the mcr over all n rounds.
print(c("For Classification Tree: emperical mcr, leave one out mcr", c(mcr.tr, mcr.cv)))
}

Analysis(2);Analysis(3);Analysis(4)