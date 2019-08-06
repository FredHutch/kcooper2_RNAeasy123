#RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR
#https://www.bioconductor.org/help/workflows/RNAseq123/

#Packages needed
#library(readr)
#library(limma)
#library(Glimma)
#library(edgeR)
#library(Mus.musculus)
#library(RColorBrewer)
#library(Glimma)
#library(gplots)

#Data packaging
#Read in featureCounts count data, select DCs for analysis
#For Shiny- allow user to pick dataset to analyze (for this dataset, can analyze DP or DC populations)
library(readr)
counts <- read_delim("2019.07.12.counts.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
View(counts)
head(counts)
rownames(counts) <- counts$GeneSymbol
library(dplyr)
print(paste("counts class is", class(counts)))
DC <- select(counts, contains("DC")) #Dendritic cells
counts.m <- as.matrix(DC)
class(counts.m)
head(counts.m)
detach("package:dplyr", unload=TRUE) #Detatch because dplyr masks some required packages for limma

#Set up
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

data <- DGEList(counts.m)
class(data)
dim(data)
#[1] 19051     7

#Sample info for DC data
samplenames <- substring(colnames(data), 1, nchar(colnames(data)))
samplenames

colnames(data) <- samplenames
batch <- as.factor(c("09182018", "09182018", 
                     "09182018", "09182018", "09182018", 
                     "10312018", "10312018"))
data$samples$batch <- batch
group <- as.factor(rep(c("D0", "D4", "D7"), c(2,3,2)))
data$samples$group <- group
data$samples

#Organising gene annotations
geneid <- rownames(DC)
genes <- select(Mus.musculus, keys=geneid, columns=c("ENTREZID", "TXCHROM"), 
                keytype="SYMBOL")
head(genes)
genes <- genes[!duplicated(genes$SYMBOL),]
data$genes <- genes
data


#Transformations from the raw-scale
cpm <- cpm(data)
lcpm <- cpm(data, log=TRUE)
L <- mean(data$samples$lib.size) * 1e-6
M <- median(data$samples$lib.size) * 1e-6
c(L, M)
#[1] 35.51705 35.65764
summary(lcpm)

#edgeR pipeline
keep <- filterByExpr(data)
edgeR.data <- data[keep, , keep.lib.sizes=FALSE]
dim(edgeR.data)
#[1] 14172     7

#Removing genes that are lowly expressed
table(rowSums(data$counts==0)==7)
#FALSE  TRUE 
#18129   922 
keep.exprs <- filterByExpr(data, group=group)
data <- data[keep.exprs,, keep.lib.sizes=FALSE]
dim(data)
##[1] 14172     7

#Normalising gene expression distributions
data <- calcNormFactors(data, method = "TMM")
data$samples$norm.factors
#[1] 1.0742536 1.0634345 0.9199327 0.9146074 0.9555229 0.9862703 1.1039632

#Unsupervised clustering of samples, generate MDS plot
#For Shiny- generate plot for user
library(RColorBrewer)
cpm <- cpm(data)
lcpm <- cpm(data, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.batch <- batch
levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Set2")
col.batch <- as.character(col.batch)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. SL-TBI")
plotMDS(lcpm, labels=batch, col=col.batch, dim=c(3,4))
title(main="B. Batch")


#Creating a design matrix and contrasts
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  DC.D0vsD4 = D0 - D4, 
  DC.D0vsD7 = D0 - D7,
  DC.D4vsD7 = D4 - D7,
  levels = colnames(design))
contr.matrix

#Removing heteroscedascity from count data, plot voom results
par(mfrow=c(1,2))
v <- voom(data, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
tfit <- treat(vfit, lfc=log2(1.5))
plotSA(tfit, main="Final model: Mean-variance trend")

#Examining the number of DE genes
#For Shiny- print out summaries for user, allow user to decide test
summary(decideTests(efit))
#       DC.D0vsD4 DC.D0vsD7 DC.D4vsD7
#Down        1557       393      1084
#NotSig     10966     13413     11803
#Up          1649       366      1285
summary(decideTests(tfit))
#       DC.D0vsD4 DC.D0vsD7 DC.D4vsD7
#Down         649       151       286
#NotSig     13062     13893     13396
#Up           461       128       490
dt <- decideTests(tfit)
de <- decideTests(efit)
de.common <- which(de[,1]!=0 & de[,2]!=0)
length(de.common)
#[1] 566
head(efit$genes$SYMBOL[de.common], n=20)
