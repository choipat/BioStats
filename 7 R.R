################################### 1 ################################### 1
# Load data
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

# Wilcox 2-sample test:
data.wilcox= NULL
for(i in 1:3051){
  data.wilcox[i] = wilcox.test (golub[i,] ~gol.fac, paired=F, alternative="greater")$p.value
}

genes.wilcox <- data.wilcox<.05
sum(genes.wilcox)

wilcox.fdr <- p.adjust(p=data.wilcox, method="fdr")
sum(wilcox.fdr < .05)

#order non-fdr
orderAML1 <- order(data.wilcox, decreasing=FALSE)
golub.gnames[orderAML1[1:3],2]

#order fdr
orderAML <- order(wilcox.fdr, decreasing=FALSE)
golub.gnames[orderAML[1:3],2]

#order diff
meanALL = apply (golub[, gol.fac=="ALL"], 1, mean)
meanAML = apply (golub[, gol.fac=="AML"], 1, mean)
meanDIFF = meanALL - meanAML
orderDIFF <- order(meanDIFF, decreasing=TRUE)
golub.gnames[orderDIFF[1:3],2]

################################### 2 ################################### 2

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

data.p = apply (golub[, gol.fac=="AML"], 1, function(x) shapiro.test(x)$p.value)

sum(data.p>.05)
p.fdr <- p.adjust(p=data.p, method="fdr")

sum(p.fdr>0.05) #Pass
sum(p.fdr<0.05) #Fail

################################### 3 ################################### 3

data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
grep("HOXA9 Homeo box A9",golub.gnames[,2])
grep("CD33",golub.gnames[,2])

wilcox.test (x= golub[1391,], y= golub[808,], paired=T, alternative="two.sided")

################################### 4 ################################### 4

source("http://www.bioconductor.org/biocLite.R")
biocLite()
library(datasets); 

str(UCBAdmissions)
Dept <- c("Dept = A","Dept = B", "Dept = C", "Dept = D", "Dept = E", "Dept = F")

for (i in 1:6 ){
  print(Dept[i])
  DeptData <- matrix(c(UCBAdmissions[1,1,i], UCBAdmissions[2,1,i], UCBAdmissions[1,2,i], 
                       UCBAdmissions[2,2,i]), nrow=2, dimnames=list("Admit"=c("Admitted","Rejected"), 
                      "Gender"=c("Male","Female")))
  print(DeptData)
  print(chisq.test(DeptData))
  print(fisher.test(DeptData))
}

################################### 5 ################################### 5

# Formula

install.packages('gtools')
library(gtools)
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
gene <- grep("gene_name",golub.gnames[,2])

data <- golub[gene,]
n <- length(data)

T.obs <- var(data[gol.fac=="ALL"]) / var(data[gol.fac=="AML"])

n.perm = 2000
T.perm = NULL
for(i in 1:n.perm) {
  data.perm = sample(data, n, replace=F) 
  T.perm[i] = var(data.perm[gol.fac=="ALL"]) / var(data.perm[gol.fac=="AML"])
}

mean(T.perm <= T.obs)

# For CD33:

install.packages('gtools')
library(gtools)
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
gene <- grep("CD33",golub.gnames[,2])

data <- golub[gene,]
n <- length(data)

T.obs <- var(data[gol.fac=="ALL"]) / var(data[gol.fac=="AML"])

n.perm = 2000
T.perm = NULL
for(i in 1:n.perm) {
  data.perm = sample(data, n, replace=F) 
  T.perm[i] = var(data.perm[gol.fac=="ALL"]) / var(data.perm[gol.fac=="AML"])
}

mean(T.perm <= T.obs)

################################## END ################################## END