#This script has been tested to analyze the "GSE96592" dataset
#Guozhen Liu @Oct 29,  2017



#The analysis outlined in this workflow assumes that reads obtained from an RNA-seq 
#experiment have been aligned to an appropriate reference genome and summarised into 
#counts associated with gene-specific regions. In this instance, reads were aligned 
#to the mouse reference genome (mm10) using the R based pipeline available in the 
#Rsubread package (specifically the align function6 followed by featureCounts7 for 
#gene-level summarisation based on the in-built mm10 RefSeq-based annotation).

# Setup --------------------------------------------------------------
source("http://bioconductor.org/workflows.R")
workflowInstall("RNAseq123")

#run at the 1st time
#install.packages(c("gplots", "RColorBrewer", "R.utils"))

#------- Load the packages we will use----------------
library("limma")   # Linear models for differential expression
library("Glimma")  # Interactive plots for exploration
library("edgeR")   # Process count data from NGS experiments
library("Mus.musculus") # Gene annotations for the Mus musculus genome


library(RColorBrewer)
library("gplots")

#-------Set up working environment----------------
remove(list = ls())
GEOid <-  'GSE96592'
DESTINATIONfile <- paste(GEOid,  "_RAW.tar",  sep="") 
getwd()

#set a working folder for one analysis/dataset:
workingFolder <- paste('C:\\BigData\\RNAseq\\', GEOid, sep="")
setwd(workingFolder)
getwd()

# Download data ------------------------------------------------------

url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GEOid&format=file"
utils::download.file(url, destfile = DESTINATIONfile, mode = "wb")
utils::untar(DESTINATIONfile, exdir = ".")
  
 
files_gz <- Sys.glob("GSM*txt.gz")
files_gz <- Sys.glob("GSM*.gz")

for(f in files_gz)
   R.utils::gunzip(f, overwrite = TRUE)


#--- Import the data-------------------------
#Each of these text files contains the raw gene-level counts for a given sample.
files <- Sys.glob("GSM*.tsv")

#read the 1st file, only the top 5 rows. Each files has only 3 columns 
read.delim(files[1], skip = 3, sep = "\t", nrow = 5) 
x <- readDGE(files, skip = 3, columns = c(1, 3))
class(x)
dim(x)
names(x)
str(x)

 
x$samples
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames

group <- as.factor(c("crushed.control",
                       "crushed.dlki",
                       "uncrushed.dlki",
                       "crushed.control",
                       "crushed.dlki",
                       "uncrushed.control",
                       "crushed.dlki",
                       "crushed.dlki",
                       "crushed.control",
                       "uncrushed.control",
                       "uncrushed.dlki",
                       "crushed.dlki",
                       "uncrushed.dlki",
                       "uncrushed.dlki",
                       "uncrushed.dlki",
                       "uncrushed.control"))
#order the group
group <- factor(group,
                levels = c("crushed.dlki",
                           "uncrushed.dlki",
                           "crushed.control",
                           "uncrushed.control"))

x$samples$group <- group


# Annotate the genes.

head(x$counts)
dim(x$counts)


#Organising gene annotations
geneid <- rownames(x)
genes <- select(Mus.musculus, keys = geneid,
                columns = c("SYMBOL", "TXCHROM"),
                keytype = "ENTREZID")
head(genes)

genes <- genes[!duplicated(genes$ENTREZID), ]

x$genes <- genes
dim(x)

#optional, Entrez Ids that no longer have official gene symbols are dropped from the 
#analysis. The whole DGEList object, including annotation as well as counts, can be 
#subsetted by rows as if it was a matrix:

sum(is.na(x$genes$SYMOL))
dim(x)


# Data pre-processing -------------------------------------------------------
#-----Transformations from the raw-scale ------------------------------------
# Calculate (log) counts-per-million (cpm)

cpm <- cpm(x)
lcpm <- cpm(x, log = TRUE)


# Removing genes that are lowly expressed---------------------------------
#All datasets will include a mix of genes that are expressed and those that are not 
#expressed. Whilst it is of interest to examine genes that are expressed in one 
#condition but not in another, some genes are unexpressed throughout all samples. 

# Detect genes with zero counts across all 16 samples

table(rowSums(x$counts == 0) == 18)

# Visualize distribution of gene expression levels

plotDensities(lcpm, legend = FALSE, main = "Before filtering")
abline(v = 0, lty = 3)

#Only keep genes which have cpm greater than 1 in at least 4 samples.
keep.exprs <- rowSums(cpm > 1) > 2
#keep <- rowSums(x$counts) > 50   #another filtering method

x <- x[keep.exprs,, keep.lib.sizes=FALSE]

dim(x)


## ----filterplot1, fig.height=4, fig.width=8, fig.cap="The density of log-CPM values for raw pre-filtered data (A) and post-filtered data (B) are shown for each sample. Dotted vertical lines mark the log-CPM of zero threshold (equivalent to a CPM value of 1) used in the filtering step."----
library(RColorBrewer)
pdf(file="EffectofFilteringLowExpGenes.pdf", paper = "a4")

nsamples <- ncol(x)
nsamples <- 10 # reassign because brewer can only hand less than 12 samples in pair


col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples ){ # "nsamples" is too large, 12 is the maximum
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")


#realize the filtering effect in lcpm:
lcpm <- cpm(x, log=TRUE)

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){   # "nsamples" is too large, 12 is the maximum
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

dev.off()


# Normalization ------------------------------------------------------

# The default normalization provided by edgeR is TMM (trimmed mean of
# M-values), which prevents differences in highly expressed genes from
# biasing the entire distribution of gene expression. It often has a
# modest effect, as observed here.

pdf(file="Normalization_Effect.pdf", paper = "a4")

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors


x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[, 1] * 0.05)
x2$counts[,2] <- x2$counts[, 2] * 5


par(mfrow=c(1,2))
lcpm2 <- cpm(x2, log = TRUE)
boxplot(lcpm2, las = 2, col = col, main = "")
title(main="A. Example: Unnormalised data", ylab="Log-cpm")


x2 <- calcNormFactors(x2)
x2$samples$norm.factors


lcpm2 <- cpm(x2, log = TRUE)
boxplot(lcpm2, las=2, col=col, main = "")
title(main="B. Example: Normalised data",ylab="Log-cpm")

dev.off()


# Data Exploration --------------------------------------------------------
#---Unsupervised clustering of samples-----------------------
 
glMDSPlot(lcpm, labels=x$samples[,1], groups=group,   
          html = paste("unSupervisedSampleClustering", GEOid, sep="_"),
          launch=FALSE)


#- Differential expression analysis----------------------------
#-----Creating a design matrix and contrasts-------------------
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

design

contr.matrix <- makeContrasts(
  uncrushed.dlki_vs_uncrushed.control  = uncrushed.dlki - uncrushed.control,
  crushed.control_vs_uncrushed.control = crushed.control - uncrushed.control,
  crushed.dlki_vs_uncrushed.control    = crushed.dlki - uncrushed.control,
  crushed.dlki_vs_uncrushed.dlki       = crushed.dlki - uncrushed.dlki,
  levels = colnames(design))
contr.matrix

#-Removing heteroscedascity from count data-----------------------------

# Convert counts to be used in linear model --------------------------
pdf( file= "Mean-Variance-Trend.pdf",  paper = "a4")

v <- voom(x, design, plot = TRUE)
dim(v)


# Test for differential expression (DE) ------------------------------

# --Fit a linear model per gene.---------

vfit <- lmFit(v, design)

# Calculate the statistics for our specific contrasts of interest.

vfit <- contrasts.fit(vfit, contrasts = contr.matrix)

# Calculate levels of significance. Uses an empirical Bayes algorithm
# to shrink the gene-specific variances towards the average variance
# across all genes.

efit <- eBayes(vfit)

# As seen in this diagnostic plot of the residual variation versus the
# mean gene expression level, `voom` successfully removed the
# relationship between the mean and the variance.

plotSA(efit, main = "Final model: Mean-variance trend")

dev.off()

# Explore the results ------------------------------------------------

tfit <- decideTests(efit, p.value = 0.01)
# Tabulate the results
summary(decideTests(efit, p.value = 0.01))
summary(tfit)

#-Examining individual DE genes from top to bottom----------------
# Identify top DE genes.
comp5 <- topTable(efit, coef = NULL, number = 3000)
write.table(comp5, file = paste("MasterResult_", GEOid, ".txt", sep=""), row.names = FALSE, sep = "\t")

comp1 <- topTable(efit, coef = 1, number = 3000)
comp2 <- topTable(efit, coef = 2, number = 3000)
comp3 <- topTable(efit, coef = 3, number = 3000)
comp4 <- topTable(efit, coef = 4, number = 3000)
 
#comp3 <- topTreat(efit, coef = 3, n = Inf)
#comp4 <- topTreat(efit, coef = 4, n = Inf)
write.table(comp1, file = paste("topDEgenes", GEOid, colnames(efit)[1], ".txt", sep="_"), 
                  row.names = FALSE, sep = "\t")
write.table(comp2, file = paste("topDEgenes", GEOid, colnames(efit)[2], ".txt", sep="_"), 
            row.names = FALSE, sep = "\t")
write.table(comp3, file = paste("topDEgenes", GEOid, colnames(efit)[3], ".txt", sep="_"), 
            row.names = FALSE, sep = "\t")
write.table(comp4, file = paste("topDEgenes", GEOid, colnames(efit)[4], ".txt", sep="_"), 
            row.names = FALSE, sep = "\t")

#Create a joined table:
comp1s <- subset(comp1[ , c(-6, -9)])
colnames(comp1s)
comp1s <- comp1s[c(1,2,3,5,4,6,7)]
colnames(comp1s)
#change column names                        )
names(comp1s)[5] <- paste("logFC", colnames(efit)[1], sep="-" )
names(comp1s)[6] <- paste("P.Value", colnames(efit)[1], sep="-" )
names(comp1s)[7] <- paste("adj.P.Val", colnames(efit)[1], sep="-" )
colnames(comp1s)


colnames(comp2)
comp2s <- subset(comp2[, c(1,2,4,7,8)])
colnames(comp2s)
#names(comp2s) <- c("logFC"     "P.Value"   "adj.P.Val")
#names(comp2s) <- paste(names(comp2s), colnames(efit)[2], sep="-" )
names(comp2s)[3] <- paste("logFC", colnames(efit)[2], sep="-" )
names(comp2s)[4] <- paste("P.Value", colnames(efit)[2], sep="-" )
names(comp2s)[5] <- paste("adj.P.Val", colnames(efit)[2], sep="-" )


colnames(comp2s)

comp3s <- subset(comp3[, c(1,2,4,7,8)])
colnames(comp3s)
names(comp3s)[3] <- paste("logFC", colnames(efit)[3], sep="-" )
names(comp3s)[4] <- paste("P.Value", colnames(efit)[3], sep="-" )
names(comp3s)[5] <- paste("adj.P.Val", colnames(efit)[3], sep="-" )

#names(comp3s) <- paste(names(comp3s), colnames(efit)[3], sep="-" )
colnames(comp3s) 


comp4s <- subset(comp4[, c(1,2,4,7,8)])
colnames(comp4s)
names(comp4s)[3] <- paste("logFC", colnames(efit)[4], sep="-" )
names(comp4s)[4] <- paste("P.Value", colnames(efit)[4], sep="-" )
names(comp4s)[5] <- paste("adj.P.Val", colnames(efit)[4], sep="-" )



#names(comp4s) <- paste(names(comp4s), colnames(efit)[4], sep="-" )
colnames(comp4s) 

total12 <- merge(comp1s, comp2s,  by = c("ENTREZID","SYMBOL"), all = TRUE)
colnames(total12) 

head(total12)

total34 <- merge(comp3s, comp4s, by = c("ENTREZID","SYMBOL"), all = TRUE)
colnames(total34) 

total1234 <- merge(total12, total34, by = c("ENTREZID","SYMBOL"), all = TRUE)
colnames(total1234) 
#remove("total1234", "total_M")


write.table(total1234, file = paste("allDEgenes", GEOid, ".txt", sep="_"), 
            row.names = FALSE, sep = "\t")

dim(total1234)

comp1 <- topTable(efit, coef = 1, number = 500)
comp2 <- topTable(efit, coef = 2, number = 500)
comp3 <- topTable(efit, coef = 3, number = 500)
comp4 <- topTable(efit, coef = 4, number = 500)


head(comp1)
head(comp2)
head(comp3)

head(comp4)

#- Visualize DE genes.--------------------------------------

#To summarise results for all genes visually, mean-difference plots, which display 
#log-FCs from the linear model fit against the average log-CPM values can be generated
#using the plotMD function, with the differentially expressed genes highlighted.

pdf( file= "DEgenes-in-different_comparison.pdf",  paper = "a4")

plotMD(efit, column = 1, status = tfit[, 1], main = colnames(efit)[1],
       values=c(1,-1), col=c("red","blue"),legend="topright",
       xlim = c(-8, 13))

plotMD(efit, column = 2, status = tfit[, 2], main = colnames(efit)[2],
       values=c(1,-1), col=c("red","blue"),legend="topright",
       xlim = c(-8, 13))

plotMD(efit, column = 3, status = tfit[, 3], main = colnames(efit)[3],
       values=c(1,-1), col=c("red","blue"),legend="topright",
       xlim = c(-8, 13))

plotMD(efit, column = 4, status = tfit[, 4], main = colnames(efit)[4],
       values=c(1,-1), col=c("red","blue"),legend="topright",
       xlim = c(-8, 13))

dev.off()


glMDPlot(efit, coef = 1, status = tfit, main = colnames(efit)[1],
         side.main = "ENTREZID", counts = x$counts, groups = group,
#         id.column = "ENTREZID", counts = x$counts, groups = group,
#         html = "antiNGFDLKivsNGFdmso_TopDEgene",
         html = paste(colnames(efit)[1], "TopDEgene", sep = "_"),
         launch = FALSE)

glMDPlot(efit, coef = 2, status = tfit, main = colnames(efit)[2],
         side.main = "ENTREZID", counts = x$counts, groups = group,
         html = paste(colnames(efit)[2], "TopDEgene", sep = "_"),
         launch = FALSE)

glMDPlot(efit, coef = 3, status = tfit, main = colnames(efit)[3],
         side.main = "ENTREZID", counts = x$counts, groups = group,
         html = paste(colnames(efit)[3], "TopDEgene", sep = "_"),
         launch = FALSE)


glMDPlot(efit, coef = 4, status = tfit, main = colnames(efit)[4],
         side.main = "ENTREZID", counts = x$counts, groups = group,
         html = paste(colnames(efit)[4], "TopDEgene", sep = "_"),
         launch = FALSE)

#---Heat map clustering-------------------------
# View heatmap of top 30 DE genes according to the treat between the cmmparing 2 set cells.
library("gplots")

pdf( file= paste(colnames(efit)[1], "Heatmap_top40genes.pdf", sep = "_"),  paper = "a4")

topgene_list <- comp1$ENTREZID[1:40]
i <- which(v$genes$ENTREZID %in% topgene_list )
#mycol <- colorpanel(1000, "blue", "white", "red")
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(v$E[i, ], scale = "row",
           labRow = v$genes$SYMBOL[i], labCol = group,
           col = mycol, trace = "none", density.info = "none",
           margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()

#######################################

pdf( file= paste(colnames(efit)[2], "Heatmap_top40genes.pdf", sep = "_"),  paper = "a4")

topgene_list <- comp2$ENTREZID[1:40]
i <- which(v$genes$ENTREZID %in% topgene_list )
#mycol <- colorpanel(1000, "blue", "white", "red")
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(v$E[i, ], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group,
          col = mycol, trace = "none", density.info = "none",
          margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()

#######################################

pdf( file= paste(colnames(efit)[3], "Heatmap_top40genes.pdf", sep = "_"),  paper = "a4")

topgene_list <- comp3$ENTREZID[1:40]
i <- which(v$genes$ENTREZID %in% topgene_list )
#mycol <- colorpanel(1000, "blue", "white", "red")
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(v$E[i, ], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group,
          col = mycol, trace = "none", density.info = "none",
          margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()

#######################################

pdf( file= paste(colnames(efit)[4], "Heatmap_top40genes.pdf", sep = "_"),  paper = "a4")

topgene_list <- comp4$ENTREZID[1:40]
i <- which(v$genes$ENTREZID %in% topgene_list )
#mycol <- colorpanel(1000, "blue", "white", "red")
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(v$E[i, ], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group,
          col = mycol, trace = "none", density.info = "none",
          margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()


########################################------------------
#generate boxplot for special interest genes:

Apoptotic_Genes <- c("Ecel1", "Hrk", "Bbc3", "Puma", "Bcl2l11", "Bim", "Ddit3", "Chop") #upregulated by ONC
v_idx <- which(v$genes$SYMBOL %in% Apoptotic_Genes)

pdf( file= "BoxPlot-Apoptotic_Genes.pdf",  paper = "a4")
#g_list <- which(v$genes$ENTREZID %in% comp1$ENTREZID[1:10]) 

for (i in 1:length(v_idx) ) {
    boxplot(v$E[v_idx[i], ]~group,col="lightblue",
            main= paste(v$genes$SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
            par(cex.axis=0.9), ylab = "logCPM", xlab = "Sample Groups"
      )
}

dev.off()


Neuronal_and_Differentiation_Genes <- c(        #downregulated
                          "Calb2",
                          "Nrn1",
                          "Brn3b",
                          "Pou4f2",
                          "Kcnd2",
                          "Scn4b",
                          "Isl2",
                          "Vsnl1",
                          "Pvalb",
                          "Nefh" )

v_idx <- which(v$genes$SYMBOL %in% Neuronal_and_Differentiation_Genes)

pdf( file= "BoxPlot-Neuronal_and_Differentiation_Genes.pdf",  paper = "a4")
#g_list <- which(v$genes$ENTREZID %in% comp1$ENTREZID[1:10]) 
for (i in 1:length(v_idx)  ) {
  boxplot(v$E[v_idx[i], ]~group,col="lightblue",
          main= paste(v$genes$SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.9), ylab = "logCPM", xlab = "Sample Groups"
  )
}

dev.off()


#upregulated by ONC
Regeneration_Genes <- c("Sprr1a", "Sox11", "Klf6", "Atf3")
v_idx <- which(v$genes$SYMBOL %in% Regeneration_Genes)

pdf( file= "BoxPlot-Regeneration_Genes.pdf",  paper = "a4")
for (i in 1:length(v_idx)  ) {
  boxplot(v$E[v_idx[i], ]~group,col="lightblue",
          main= paste(v$genes$SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.9), ylab = "logCPM", xlab = "Sample Groups"
  )
}

dev.off()


#upregulated
Inflammation_Genes <- c("Aif1", "Gfap", "C1qc", "C4b", "Mmp12", "Cxcl10")
v_idx <- which(v$genes$SYMBOL %in% Inflammation_Genes)

pdf( file= "BoxPlot-Inflammation_Genes.pdf",  paper = "a4")
for (i in 1:length(v_idx)  ) {
  boxplot(v$E[v_idx[i], ]~group,col="lightblue",
          main= paste(v$genes$SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.9), ylab = "logCPM", xlab = "Sample Groups"
  )
}

dev.off()

DLK_pathway_Genes <- c("Dlk", "Map3k12", "Mlk1", "Map3k9", "Mlk2", "Map3k10", 
                       "Mlk3", "Map3k11"," Mkk4", "Map2k4", "Mkk7", "Map2k7", "Jun")
v_idx <- which(v$genes$SYMBOL %in% DLK_pathway_Genes)

pdf( file= "BoxPlot-DLK_pathway_Genes.pdf",  paper = "a4")
for (i in 1:length(v_idx)  ) {
  boxplot(v$E[v_idx[i], ]~group,col="lightblue",
          main= paste(v$genes$SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.9), ylab = "logCPM", xlab = "Sample Groups"
  )
}

dev.off()


# Test for enrichment of gene sets in Gene Ontology terms---------------------------

go <- goana(efit, coef =1, FDR = 0.05, species="Mm")

gos1 <- topGO(go, n=500, truncate.term = 40 )
gos2 <- gos1[ which (gos1$P.Up < 0.01 | gos1$P.Down < 0.01), ]


write.table(gos2, file= paste("GeneGOEnrich_", colnames(contr.matrix)[1], ".txt", 
      sep =""), col.names=TRUE, row.names=TRUE, sep="\t")
 
go <- goana(efit, coef =2, FDR = 0.05, species="Mm")
gos1 <- topGO(go, n=500,  truncate.term = 40 )
gos2 <- gos1[ which (gos1$P.Up < 0.01 | gos1$P.Down < 0.01), ]

write.table(gos2, file= paste("GeneGOEnrich_", colnames(contr.matrix)[2], ".txt", 
                                         sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

go <- goana(efit, coef =3, FDR = 0.05, species="Mm")
#topGO(go, n=15)
gos1 <- topGO(go, n=500,  truncate.term = 40 )
gos2 <- gos1[ which (gos1$P.Up < 0.01 | gos1$P.Down < 0.01), ]


write.table(gos2, file= paste("GeneGOEnrich_", colnames(contr.matrix)[3], ".txt", 
                                         sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

go <- goana(efit, coef =4, FDR = 0.05, species="Mm")
gos1 <- topGO(go, n=500,  truncate.term = 40 )
gos2 <- gos1[ which (gos1$P.Up < 0.01 | gos1$P.Down < 0.01), ]

write.table(gos2, file= paste("GeneGOEnrich5_", colnames(contr.matrix)[4], ".txt", 
                        sep =""), col.names=TRUE, row.names=TRUE, sep="\t")


keg <- kegga(efit, coef =1, FDR = 0.05, species="Mm")
#topKEGG(keg, n=15, truncate=34)
kegs1 <- topKEGG(keg, n=100,  truncate.path = 40 )
kegs2 <- kegs1[ which (kegs1$P.Up < 0.01 | kegs1$P.Down < 0.01), ]


write.table(kegs2, file= paste("KEGGenrich_", colnames(contr.matrix)[1], ".txt", 
                                    sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

keg <- kegga(efit, coef =2, FDR = 0.05, species="Mm")
#topKEGG(keg, n=15, truncate=34)
kegs1 <- topKEGG(keg, n=100,  truncate.path = 40 )
kegs2 <- kegs1[ which (kegs1$P.Up < 0.01 | kegs1$P.Down < 0.01), ]


write.table(kegs2, file= paste("KEGGenrich_", colnames(contr.matrix)[2], ".txt", 
                                    sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

keg <- kegga(efit, coef =3, FDR = 0.05, species="Mm")
#topKEGG(keg, n=15, truncate=34)
kegs1 <- topKEGG(keg, n=100,  truncate.path = 40 )
kegs2 <- kegs1[ which (kegs1$P.Up < 0.01 | kegs1$P.Down < 0.01), ]

write.table(kegs2, file= paste("KEGGenrich_", colnames(contr.matrix)[3], ".txt", 
                                   sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

keg <- kegga(efit, coef =4, FDR = 0.05,  species="Mm")
kegs1 <- topKEGG(keg, n=100,  truncate.path = 40 )
kegs2 <- kegs1[ which (kegs1$P.Up < 0.01 | kegs1$P.Down < 0.01), ]

write.table(kegs2, file= paste("KEGGenrich_", colnames(contr.matrix)[4], ".txt", 
                                  sep =""), col.names=TRUE, row.names=TRUE, sep="\t")



# Report session information -----------------------------------------
sessionInfo()

