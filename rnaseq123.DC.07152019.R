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

#Generate venn diagrams for shared up-regulated and down-regulated genes
#For Shiny- generate venn diagrams for user, choice in color?
vennDiagram(de[,1:3], include = "up",
            circle.col=c("red", "blue", "green"), main= "Up-regulated genes in dendritic cells, eBays")
vennDiagram(dt[,1:3], include = "up",
            circle.col=c("red", "blue", "green"), main= "Up-regulated genes in dendritic cells, LFC")

vennDiagram(de[,1:3], include = "down",
            circle.col=c("red", "blue", "green"), main= "Down-regulated genes in dendritic cells, eBays")
vennDiagram(dt[,1:3], include = "down",
            circle.col=c("red", "blue", "green"), main= "Down-regulated genes in dendritic cells, LFC")

write.fit(tfit, dt, file="DEresults.lfc.txt")

#Examining individual DE genes from top to bottom
DC.d0vsd4 <- topTreat(tfit, coef=1, n=Inf)
DC.d0vsd7 <- topTreat(tfit, coef=2, n=Inf)
DC.d4vsd7 <- topTreat(tfit, coef=3, n=Inf)
head(DC.d0vsd4)

#MD plots of differential expression results
#For Shiny- generate MD plots
plotMD(efit, column=1, status=de[,1], main=colnames(efit)[1], 
       xlim=c(-5,15))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], col= c("blue", "red"),
       xlim=c(-5,15))
plotMD(tfit, column=2, status=dt[,2], main=colnames(tfit)[2], col= c("blue", "red"),
       xlim=c(-5,15))
plotMD(tfit, column=3, status=dt[,3], main=colnames(tfit)[3], col= c("blue", "red"),
       xlim=c(-5,15))

#Glimma MD plots
#For Shiny- can Glimma be integrated? Does CSS play well with Shiny?
library(Glimma)
glMDPlot(tfit, coef=3, status=dt, main=colnames(tfit)[3],
         side.main="SYMBOL", counts=lcpm, groups=group, launch=FALSE)

#Generate heatmap
#For Shiny- generate heatmap for user, allow user to determine how many genes plotted, what comparison to plot
library(gplots)
DC.d0vsd4.topgenes <- DC.d0vsd4$SYMBOL[1:100]
i <- which(v$genes$SYMBOL %in% DC.d0vsd4.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=rownames(v$targets), 
          col=mycol, trace="none", density.info="none", 
          margin=c(7,5), lhei=c(2,10), dendrogram="column",
          keysize = 1)

DC.d0vsd7.topgenes <- DC.d0vsd7$SYMBOL[1:100]
d7 <- which(v$genes$SYMBOL %in% DC.d0vsd7.topgenes)
#mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[d7,], scale="row",
          labRow=v$genes$SYMBOL[d7], labCol=rownames(v$targets), 
          col=mycol, trace="none", density.info="none", 
          margin=c(7,6), lhei=c(2,10), dendrogram="column",
          keysize = 1)

DC.d4vsd7.topgenes <- DC.d4vsd7$SYMBOL[1:100]
both <- which(v$genes$SYMBOL %in% DC.d4vsd7.topgenes)
#mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[both,], scale="row",
          labRow=v$genes$SYMBOL[both], labCol=rownames(v$targets), 
          col=mycol, trace="none", density.info="none", 
          margin=c(7,6), lhei=c(2,10), dendrogram="column",
          keysize = 1)

#Gene set testing with camera
#install.packages("msigdbr")
library(msigdbr)
#For Shiny- user selects gene set collection, view results of head() and choose direction
#Retrieve mouse curated collection (C2) gene sets
load("MSigDB/mouse_c2_v5p2.rdata")
idx <- ids2indices(Mm.c2,id=v$genes$ENTREZID)
cam.DC.d0vsd4 <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.DC.d0vsd4,5)
head(cam.DC.d0vsd4[cam.DC.d0vsd4$Direction == "Up", ], 5)
write.csv(cam.DC.d0vsd4, "camera.DC.d0vsd4.csv")

cam.DC.d0vsd7 <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.DC.d0vsd7,10)
write.csv(cam.DC.d0vsd7, "camera.DC.d0vsd7.csv")

cam.DC.d4vsd7 <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.DC.d4vsd7,5)
write.csv(cam.DC.d4vsd7, "camera.DC.d4vsd7.csv")

#Barcode example
#For Shiny- generate barcode, choose gene set
dev.off()
barcodeplot(tfit$t[,1], index=idx$NEMETH_INFLAMMATORY_RESPONSE_LPS_UP, 
            index2=idx$REACTOME_GENERIC_TRANSCRIPTION_PATHWAY, main="DC D0vsD4")

sessionInfo()