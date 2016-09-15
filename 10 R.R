# 1 ######################################## 1 ######################################### 1

biocLite("ArrayExpress");library(ArrayExpress)
yeast.raw <- ArrayExpress('E-MEXP-1551')

# 1a ####################################### 1a ######################################## 1a

eset.1a <- expresso(normData.1a,bgcorrect.method="mas",
                    normalize.method="quantiles",
                    pmcorrect.method="pmonly",
                    summary.method="medianpolish")
exprs.yeast.1a <- exprs(eset.1a)

# 1b ####################################### 1b ######################################## 1b

exprs.yeast.1a[1:5,]

# 1c ####################################### 1c ######################################## 1c

str(exprs.yeast.1a)


# 2 ######################################## 2 ######################################### 2

biocLite("ArrayExpress"); biocLite("annotate"); library(ArrayExpress)
yeast.raw <- ArrayExpress('E-MEXP-1551')

# 2a ####################################### 2a ######################################## 2a

annotation(yeast.raw)
biocLite("yeast2.db")

# 2b ####################################### 2b ######################################## 2b

library(yeast2.db)
go1769308_at <- get("1769308_at", env = yeast2GO)
library(annotate)
mf.go1769308_at <- getOntology(go1769308_at,"MF")
length(mf.go1769308_at) # 7

# 2c ####################################### 2c ######################################## 2c

biocLite("GO.db")
library(GO.db)

x <- get("1769308_at", env = yeast2GO)
gonr <- getOntology(go1769308_at, "MF")
gP <- getGOParents(gonr) #parents for GO numbers
pa <- sapply(gP,function(x) x$Parents) #extract GO numbers for parents
pa
length(pa)

# 2d ####################################### 2d ######################################## 2d

gC <- getGOChildren(gonr) #children for GO numbers
ch <- sapply(gC,function(x) x$Children) #extract GO numbers for children
ch
length(unlist(ch))


# 3 ######################################## 3 ######################################### 3

biocLite("genefilter");biocLite("genefilter");library("genefilter");library("ALL"); data(ALL);library("limma")

# 3a ####################################### 3a ######################################## 3a

patientB <- exprs(ALL)[,(ALL$BT %in% c("B2","B3"))]
factor <- droplevels(ALL$BT[ALL$BT %in% c("B2","B3")])
f1 <- function(x) (wilcox.test(x ~ factor, exact = F)$p.value < 0.001)
f2 <- function(x) (t.test(x ~ factor)$p.value < 0.001) 
wilcox <- genefilter(patientB, filterfun(f1))
t.test <- genefilter(patientB, filterfun(f2))

# 3b ####################################### 3b ######################################## 3b

x <- apply(cbind(wilcox,t.test), 2, as.integer)
vc <- vennCounts(x, include="both")
vennDiagram(vc)

# 3c ####################################### 3c ######################################## 3c

# From Venn diagram or we can use sum(wilcox); sum(selected)
# Wilcox = 48+297  345
# Both = 297

# 3d ####################################### 3d ######################################## 3d

annotation(ALL)
biocLite("hgu95av2.db")
library(hgu95av2.db)

GO.term <- function(term) {
  GTL <- eapply(GOTERM, function(x) {grep(term, x@Term, value=TRUE)}) 
  Gl <- sapply(GTL, length)         
  names(GTL[Gl>0])                  
}
oncogene.id<- GO.term("oncogene")
oncogene.id

# 3e ####################################### 3e ######################################## 3e

selected <- wilcox & t.test
sel<- patientB[selec,]
probes <- hgu95av2GO2ALLPROBES$"GO:0090402"
print(sum(probes %in% rownames(sel)))

# 4 ######################################## 4 ######################################### 4

library("limma");library("ALL");data(ALL);library("genefilter")

# 4a ####################################### 4a ######################################## 4a

allB <- ALL[,which(ALL$BT %in% c("B1","B2","B3"))]

# 4b ####################################### 4b ######################################## 4b

facB123 <- factor(allB$BT)
design.ma <- model.matrix(~ 0 + facB123)
colnames(design.ma) <- c("B1","B2","B3")
fit <- lmFit(allB, design.ma); fit <- eBayes(fit)
print(topTable(fit, number=5,adjust.method="fdr"), digits=4) # p is small. reject null that means = 0; they are different
sum(topTable(fit, number=Inf,adjust.method="fdr")$adj.P.Val<0.05) # differently expressed

# For B3 group:
print( topTable(fit, coef=3, number=5, adjust.method="fdr"), digits=4) 

# 4c ####################################### 4c ######################################## 4c

cont.ma <- makeContrasts(B1-B2,B2-B3, levels=factor(allB$BT))
fit1 <- contrasts.fit(fit, cont.ma);fit1 <- eBayes(fit1)
fdr.p.data <- topTable(fit1,number=Inf,adjust.method="fdr")$adj.P.Val
sum(fdr.p.data<0.01) # differently expressed

print(topTable(fit1, number=5,adjust.method="fdr"), digits=4) 

# END ###################################### END ####################################### END