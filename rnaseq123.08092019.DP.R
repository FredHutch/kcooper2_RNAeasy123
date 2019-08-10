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
##Read in featureCounts count data
library(readr)
counts <- read_delim("2019.07.12.counts.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
View(counts)
head(counts)
rownames(counts) <- counts$GeneSymbol
library(dplyr)
DP <- select(counts, contains("DP")) #Double Positive Thymocytes
counts.m <- as.matrix(DP)
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
#[1] 19051     9


#Sample info for DP data
samplenames <- substring(colnames(data), 1, nchar(colnames(data)))
samplenames

colnames(data) <- samplenames
batch <- as.factor(c("08102018", "08102018", "08072018", 
                     "09182019", "09182019", "09182019", 
                     "08102019", "08102019", "08102019"))
data$samples$batch <- batch
group <- as.factor(rep(c("D0", "D4", "D7"), c(3,3,3)))
data$samples$group <- group
data$samples

#Organising gene annotations
geneid <- rownames(DP)
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
summary(lcpm)

#Removing genes that are lowly expressed
table(rowSums(data$counts==0)==9)
#FALSE  TRUE 
#17657  1394
keep.exprs <- filterByExpr(data, group=group)
data <- data[keep.exprs,, keep.lib.sizes=FALSE]
dim(data)
#[1] 12840     9

#Normalising gene expression distributions
data <- calcNormFactors(data, method = "TMM")
data$samples$norm.factors
#[1] 0.9444228 1.0200899 1.0185705 1.0308349 0.9885040 0.9969189 0.9897116 1.0059864
#[9] 1.0075715

#Unsupervised clustering of samples
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
  DP.D0vsD4 = D0 - D4, 
  DP.D0vsD7 = D0 - D7,
  DP.D4vsD7 = D4 - D7,
  levels = colnames(design))
contr.matrix

#Removing heteroscedascity from count data
par(mfrow=c(1,2))
v <- voom(data, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
tfit <- treat(vfit, lfc=log2(1.5))
plotSA(efit, main="Final model: Mean-variance trend")

#Examining the number of DE genes
summary(decideTests(efit))
#       DP.D0vsD4 DP.D0vsD7 DP.D4vsD7
#Down        1475         0      1244
#NotSig     10397     12840      9605
#Up           968         0      1991
summary(decideTests(tfit))
#       DP.D0vsD4 DP.D0vsD7 DP.D4vsD7
#Down         451         0       179
#NotSig     12294     12840     11909
#Up            95         0       752
dt <- decideTests(tfit)
de <- decideTests(efit)
de.common <- which(de[,1]!=0 & de[,3]!=0)
length(de.common)
#[1] 2027
head(efit$genes$SYMBOL[de.common], n=20)

#Generate venn diagrams for shared up-regulated and down-regulated genes
#For Shiny- generate venn diagrams for user, choice in color?
vennDiagram(de[,1:3], include = "up",
            circle.col=c("red", "blue", "green"), main= "Up-regulated genes in double-positive thymocytes, eBays")
vennDiagram(dt[,1:3], include = "up",
            circle.col=c("red", "blue", "green"), main= "Up-regulated genes in double-positive thymocytes, LFC")

vennDiagram(de[,1:3], include = "down",
            circle.col=c("red", "blue", "green"), main= "Down-regulated genes in double-positive thymocytes, eBays")
vennDiagram(dt[,1:3], include = "down",
            circle.col=c("red", "blue", "green"), main= "Down-regulated genes in double-positive thymocytes, LFC")

write.fit(tfit, dt, file="DEresults.lfc.txt")

#Examining individual DE genes from top to bottom
DP.d0vsd4 <- topTreat(tfit, coef=1, n=Inf)
DP.d0vsd7 <- topTreat(tfit, coef=2, n=Inf)
DP.d4vsd7 <- topTreat(tfit, coef=3, n=Inf)
head(DP.d0vsd4)

#MD plots of differential expression results
#For Shiny- generate MD plots
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
DP.d0vsd4.topgenes <- DP.d0vsd4$SYMBOL[1:100]
i <- which(v$genes$SYMBOL %in% DP.d0vsd4.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=rownames(v$targets), 
          col=mycol, trace="none", density.info="none", 
          margin=c(7,5), lhei=c(2,10), dendrogram="column",
          keysize = 1)

DP.d0vsd7.topgenes <- DP.d0vsd7$SYMBOL[1:100]
d7 <- which(v$genes$SYMBOL %in% DP.d0vsd7.topgenes)
#mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[d7,], scale="row",
          labRow=v$genes$SYMBOL[d7], labCol=rownames(v$targets), 
          col=mycol, trace="none", density.info="none", 
          margin=c(7,6), lhei=c(2,10), dendrogram="column",
          keysize = 1)

DP.d4vsd7.topgenes <- DP.d4vsd7$SYMBOL[1:100]
d4vsd7 <- which(v$genes$SYMBOL %in% DP.d4vsd7.topgenes)
#mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[d4vsd7,], scale="row",
          labRow=v$genes$SYMBOL[d4vsd7], labCol=rownames(v$targets), 
          col=mycol, trace="none", density.info="none", 
          margin=c(7,6), lhei=c(2,10), dendrogram="column",
          keysize = 1)

#Gene set testing with camera
#For Shiny- user selects gene set collection, view results of head() and choose direction
#Retrieve mouse curated collection (C2) gene sets
load("~/kcooper2_RNAeasy123/MSigDB/mouse_c2_v5p2.rdata")
idx <- ids2indices(Mm.c2,id=v$genes$ENTREZID)
cam.DP.d0vsd4 <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.DP.d0vsd4,5)
head(cam.DP.d0vsd4[cam.DP.d0vsd4$Direction == "Up", ], 5)
write.csv(cam.DP.d0vsd4, "camera.DP.d0vsd4.csv")

cam.DP.d0vsd7 <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.DP.d0vsd7,10)
write.csv(cam.DP.d0vsd7, "camera.DP.d0vsd7.csv")

cam.DP.d4vsd7 <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.DP.d4vsd7,5)
write.csv(cam.DP.d4vsd7, "camera.DP.d4vsd7.csv")

#Barcode example
#For Shiny- generate barcode, choose gene set
dev.off()
barcodeplot(tfit$t[,1], index=idx$REACTOME_RNA_POL_I_TRANSCRIPTION, 
            index2=idx$ICHIBA_GRAFT_VERSUS_HOST_DISEASE_D7_UP, main="DP D0vsD4")

sessionInfo()
