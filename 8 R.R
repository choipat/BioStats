install.packages("lmtest");library(ALL);data(ALL);library(lmtest);library(Biobase);
# Includes all the required libraries and packages for this module

######################################## 1a ######################################## 1a
# Conduct the one-way ANOVA. Do the disease stages affect the mean geneexpression value?

ALLB1234 <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))]
data <- exprs(ALLB1234)["109_at",]
anova(lm(data ~ ALLB1234$BT))

# p = 0.01082; Reject Null hypothesis as p-value is small and 
# conclude that 109_t affects the mean gene expression values

######################################## 1b ######################################## 1b
# Conduct the one-way ANOVA. Do the disease stages affect the mean gene expression value?

summary(lm(data ~ ALLB1234$BT-1))

######################################## 1c ######################################## 1c
# Which group’s mean gene expression value is different from that of group B?

pairwise.t.test(data, ALLB1234$BT)
# p > 0.05; Can't reject Null; So B1,2,3,4 are not different 
# from B

######################################## 1d ######################################## 1d
# Use the pairwise comparisons at FDR=0.05 to find which group means are different. 
# What is your conclusion?

pairwise.t.test(data,ALLB1234$BT,p.adjust.method='fdr')

######################################## 1e ######################################## 1e
# Check the ANOVA model assumptions with diagnostic tests? Do we need to
# apply robust ANOVA tests here? If yes, apply the appropriate tests and state
# your conclusion.

shapiro.test(residuals(lm(data ~ ALLB1234$BT)))
# Accept Null that it is Normal dist. Robust not required

bptest(lm(data ~ ALLB1234$BT), studentize = FALSE)
# Accept Null of equal variances. Robust not required

######################################## 2a ######################################## 2a
# Use FDR adjustments at 0.05 level. How many genes are expressed different in 
# some of the groups?

ALLB1234 <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4"))]
dataALL <- exprs(ALLB1234)[,]
ALLData <- apply(dataALL,1,function(x) kruskal.test(x ~ ALLB1234$BT)$p.value)
ALLData.fdr <- p.adjust(p=ALLData,method="fdr")
sum(ALLData.fdr<0.05)
# Reject null-hypothesis of equal distribution of gene expression for 423

######################################## 2b ######################################## 2b
# Find the gene names for the top five genes with smallest p-values. 

orderALLfdr <- order(ALLData.fdr, decreasing=FALSE)
k=1
names = NULL
for (i in orderALLfdr[1:5]){names[k] <- names(ALLData.fdr[i])
                            k=k+1}
print('Top five genes with smallest p-values ='); names

######################################## 3a ######################################## 3a
# Conduct the appropriate ANOVA analysis. Does any of the two factors affects the gene 
# expression values? Are there interaction between the two factors?

ALLdatasex <- ALL[,which(ALL$BT %in% c("B1","B2","B3","B4") & ALL$sex %in% c("M", "F"))]
ALLdataSex <- exprs(ALLdatasex)["38555_at",]
Bcell<-ALLdatasex$BT
sexg<-ALLdatasex$sex
anova(lm(ALLdataSex~ Bcell*sexg))

# BCell affects, sexg does not affect; No interaction as p = 0.9095

######################################## 3b ######################################## 3b
# Check the ANOVA model assumption with diagnostic tests? Are any of the 
# assumptions violated?

shapiro.test(residuals(lm(ALLdataSex ~ Bcell+sexg)))
#Reject Null hypothesis of normal distribution. Other robust tests required. Violated.

bptest(lm(ALLdataSex ~ Bcell+sexg), studentize = FALSE)
# Accept Null of equal variances. Not violated.

######################################## 4a ######################################## 4a
# Program this permutation test in R.

stages <- #Example: c("B1","B2")
  gene <-  #Example: "109_at"
  
  Statistics <- function(stages,gene){
    
    ALLB <- ALL[,which(ALL$BT %in% stages)]
    data <- exprs(ALLB)[gene,]
    group <- ALLB$BT[,drop=T]
    
    g <- length(stages)
    Means <- summary(lm(data ~ group-1))[["coefficients"]][1:g]
    TotalMean <- (1/g)*sum(Means)
    
    Uj_U <- NULL
    for (i in 1:g){
      Uj_U[i] <- (Means[i]-TotalMean)^2
    } 
    T.obs <- (1/(g-1))*sum(Uj_U) #Observed statistic
    
    n <- length(data)
    n.perm = 2000
    T.perm = NULL
    for(i in 1:n.perm) {
      data.perm = sample(data, n, replace=F)
      MeansPerm <- summary(lm(data.perm ~ ALLB$BT-1))[["coefficients"]][1:g]
      TotalMeanPerm <- (1/g)*sum(MeansPerm)
      Uj_U1 <- NULL
      for (k in 1:g){
        Uj_U1[k] <- (MeansPerm[k]-TotalMeanPerm)^2
      } 
      T.perm[i] = (1/(g-1))*sum(Uj_U1) #Permuted statistic
    }
    mean(T.perm>=T.obs) #p-value
  }

######################################## 4b ######################################## 4b
# Run this permutation test on the Ets2 repressor gene 1242_at on the patients
# in stage B1, B2, and B3 from the ALL data set.

stages <- c("B1","B2","B3")
gene <- "1242_at"

Statistics(stages,gene)  # Function "Statistics" defined in 4a
#or   Statistics(c("B1","B2","B3"), "1242_at")

# Accept Null: Equally distributed expression values

######################################## END ######################################## END