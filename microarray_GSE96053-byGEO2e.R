# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Jun 7 12:11:07 EDT 2018

#this version change the label of the last sample 
#Gordon Liu @June 22, 2018
################################################################

getwd()
remove(list = ls())


#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

##############################

destGEOdir <- 'C:\\BigData\\TempGEO'
GEOid <- 'GSE96053'
destdir <- destGEOdir
#set a working folder for one analysis/dataset:
workingFolder <- paste('C:\\BigData\\ElikC\\GSE96053CorrectedLabel')
setwd(workingFolder)
getwd()
destdir

# load series and platform data from GEO
if (! exists("gset") ){
    gset <- getGEO("GSE96053", GSEMatrix =TRUE, AnnotGPL=FALSE,  destdir=destdir)
}

if (length(gset) > 1) idx <- grep("GPL11202", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))


Sample_list <- c("GSM2528016", "GSM2528017", "GSM2528018", "GSM2528019", "GSM2528020", 
                  "GSM2528021", "GSM2528022", "GSM2528023", "GSM2528024", "GSM2528025", 
                  "GSM2528026", "GSM2528027", "GSM2528028", "GSM2528029", "GSM2528030", 
                  "GSM2528031", "GSM2528032", "GSM2528033")


# str(gset)
# 
# 
# dim(gset)
# head(gset)
# dim(exprs(gset))
# class(gset)
# featureNames(gset)[1:5]
# colnames(fData(gset))
# dim(fData(gset))

HasSymbol <- !is.na(fData(gset)$GENE_SYMBOL)

gset <- gset[HasSymbol,] 

HasSymbol <- !is.na(fData(gset)$GENE)

gset <- gset[HasSymbol,] 

#gset <- gset[ , -c(11,18)] 

#colnames(fData(gset))
#fData(gset) <- fData(gset)[ , -c(1,2,3,5,8,9,10,11,12,14,16,17)]

dim(gset)

# y <- neqc(exprs(gset))
# 
# dim(y)

# library(genefilter)
# nsFilter(gset, require.entrez=TRUE,
#          require.GOBP=FALSE, require.GOCC=FALSE,
#          require.GOMF=FALSE, require.CytoBand=FALSE,
#          remove.dupEntrez=TRUE, var.func=IQR,
#          var.cutoff=0.5, var.filter=TRUE,
#          filterByQuantile=TRUE)  # feature.exclude="^AFFX", ...)


#Normalize between the arrays. Select a method appropriate for your data - see the user guide for details. 

ex <- exprs(gset)
y <- normalizeBetweenArrays(ex, method="quantile")


# #Use the avereps function to average replicate spots. 
# y.ave <- avereps(y, ID=y$genes$ProbeName)

# 
# vv <- gset[1:5, 1:3]
# dim(vv)
# gset <- gset[ , -c(7,11)] 
#   
# dim(gset)
# dim(fData(gset))
# sampleNames(gset)
# 
# length(sampleNames(gset))


# group names for all samples
# gsms <- "021321303020122011"
# sml <- c()
# for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
# 
# sml <- paste("G", sml, sep="")    # set group names
group <- as.factor(c(
  "Uncrushed.DlkNeg",
  "Crushed.DlkNeg",                     
  "Uncrushed.DlkEx",                       #
  "Crushed.DlkEx",
  "Crushed.DlkNeg",
  "Uncrushed.DlkEx",
  "Crushed.DlkEx",
  "Uncrushed.DlkNeg",
  "Crushed.DlkEx",
  "Uncrushed.DlkNeg",
  "Crushed.DlkNeg",
  "Uncrushed.DlkNeg",
  "Uncrushed.DlkEx",
  "Crushed.DlkNeg",
  "Crushed.DlkNeg",
  "Uncrushed.DlkNeg",
  "Uncrushed.DlkEx",
#  "Uncrushed.DlkEx" ))
  "Crushed.DlkEx"))

#order the group
group <- factor(group,
                levels = c("Crushed.DlkNeg",
                           "Uncrushed.DlkNeg",
                           "Crushed.DlkEx",
                           "Uncrushed.DlkEx"))



################################################################
#   Boxplot for selected GEO samples
# order samples by group
#ex <- exprs(gset)[ , order(group)]
exprs(gset) <- y[ , order(group)]

 
 

#######################################
#fl <- as.factor(sml)
gset$description <- group     #fl
design <- model.matrix(~ description + 0, gset)

colnames(design) <- levels(group)
###############################################
library(arrayQualityMetrics)

#colnames(pData(gset))
arrayQualityMetrics(expressionset = gset,
                    outdir = "QAreport_for_GSE96053",
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = "description"
                    )


######----------------------------
#get ready for linear modeling
###############################

fit <- lmFit(gset, design)

contr.matrix <- makeContrasts(
          Uncrushed.DlkNeg_vs_Uncrushed.DlkEx  = Uncrushed.DlkNeg - Uncrushed.DlkEx,
          Crushed.DlkEx_vs_Uncrushed.DlkEx   = Crushed.DlkEx - Uncrushed.DlkEx,
          Crushed.DlkNeg_vs_Uncrushed.DlkEx    = Crushed.DlkNeg - Uncrushed.DlkEx,
          Crushed.DlkNeg_vs_Uncrushed.DlkNeg     = Crushed.DlkNeg - Uncrushed.DlkNeg,
          levels=design)
contr.matrix
##############################################################
cfit <- contrasts.fit(fit, contrasts = contr.matrix)

summary(cfit)
#proportion	----the numeric value between 0 and 1, assumed proportion of genes which are differentially expressed
efit <- eBayes(cfit, proportion = 0.2)

summary(efit)

tfit <- decideTests(efit, p.value = 0.01)
# Tabulate the results
summary(tfit)




comp5 <- topTable(efit, coef = NULL, number = 1000)
write.table(comp5, file = paste("MasterResult2_", GEOid, ".txt", sep=""), row.names = FALSE, sep = "\t")

comp1 <- topTable(efit, coef = 1, number = 500)
comp2 <- topTable(efit, coef = 2, number = 500)
comp3 <- topTable(efit, coef = 3, number = 500)
comp4 <- topTable(efit, coef = 4, number = 500)

comp1s <- subset(comp1[ , c(6,7,13, 19, 18, 21, 22)])
colnames(comp1s)
comp2s <- subset(comp2[ , c(6,7,13, 19, 18, 21, 22)])
comp3s <- subset(comp3[ , c(6,7,13, 19, 18, 21, 22)])
comp4s <- subset(comp4[ , c(6,7,13, 19, 18, 21, 22)])

names(comp1s)[2] <- "SYMBOL"
names(comp1s)[1] <- "ENTREZID"
colnames(comp1s)
names(comp2s)[2] <- "SYMBOL"
names(comp2s)[1] <- "ENTREZID"

names(comp3s)[2] <- "SYMBOL"
names(comp3s)[1] <- "ENTREZID"

names(comp4s)[2] <- "SYMBOL"
names(comp4s)[1] <- "ENTREZID"

write.table(comp1s, file = paste("topDEgenes", GEOid, colnames(efit)[1], ".txt", sep="_"), 
            row.names = FALSE, sep = "\t")
write.table(comp2s, file = paste("topDEgenes", GEOid, colnames(efit)[2], ".txt", sep="_"), 
            row.names = FALSE, sep = "\t")
write.table(comp3s, file = paste("topDEgenes", GEOid, colnames(efit)[3], ".txt", sep="_"), 
            row.names = FALSE, sep = "\t")
write.table(comp4s, file = paste("topDEgenes", GEOid, colnames(efit)[4], ".txt", sep="_"), 
            row.names = FALSE, sep = "\t")


# colnames(comp1)
# [1] "ID"                   "SPOT_ID"              "CONTROL_TYPE"         "REFSEQ"              
# [5] "GB_ACC"               "GENE"                 "GENE_SYMBOL"          "GENE_NAME"           
# [9] "UNIGENE_ID"           "ENSEMBL_ID"           "TIGR_ID"              "ACCESSION_STRING"    
# [13] "CHROMOSOMAL_LOCATION" "CYTOBAND"             "DESCRIPTION"          "GO_ID"               
# [17] "SEQUENCE"             "logFC"                "AveExpr"              "t"                   
# [21] "P.Value"              "adj.P.Val"            "B"          

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


library("Glimma")  # Interactive plots for exploration
colnames(fData(gset))
fData(gset) <- fData(gset)[ , -c(1,2,3,5,8,9,10,11,12,14,16,17)]

glMDPlot(efit, coef = 1, status = tfit, main = colnames(efit)[1],
      #   side.main = "ENTREZID", 
         counts = exprs(gset), groups = group, labels=gset$GeneID, 
         #         id.column = "ENTREZID", counts = x$counts, groups = group,
         html = paste(colnames(efit)[1], "TopDEgene", sep = "_"),
         launch = FALSE)

 
glMDPlot(efit, coef = 2, status = tfit, main = colnames(efit)[2],
         #   side.main = "ENTREZID", 
         counts = exprs(gset), groups = group, labels=gset$GeneID, 
         #         id.column = "ENTREZID", counts = x$counts, groups = group,
         html = paste(colnames(efit)[2], "TopDEgene", sep = "_"),
         launch = FALSE)
 
glMDPlot(efit, coef = 3, status = tfit, main = colnames(efit)[3],
         #   side.main = "ENTREZID", 
         counts = exprs(gset), groups = group, labels=gset$GeneID, 
         #         id.column = "ENTREZID", counts = x$counts, groups = group,
         html = paste(colnames(efit)[3], "TopDEgene", sep = "_"),
         launch = FALSE)

glMDPlot(efit, coef = 4, status = tfit, main = colnames(efit)[4],
         #   side.main = "ENTREZID", 
         counts = exprs(gset), groups = group, labels=gset$GeneID, 
         #         id.column = "ENTREZID", counts = x$counts, groups = group,
         html = paste(colnames(efit)[4], "TopDEgene", sep = "_"),
         launch = FALSE)


#---Heat map clustering-------------------------
#######################################
fData.sorted = fData(gset)[order(rownames(fData(gset))),]

# Sort original data frame by the order of the new data frame
exprs(gset)[match(rownames(fData.sorted), rownames(exprs(gset))), ]


##########################################

# View heatmap of top 30 DE genes according to the treat between the cmmparing 2 set cells.
library("gplots")
# colnames(fData(gset))
# colnames(pData(gset))
# colnames(exprs(gset))

pdf( file= paste(colnames(efit)[1], "Heatmap_top20genes.pdf", sep = "_"),  paper = "a4")
topgene_list <- comp1$GENE[1:20]
topgene_list <- topgene_list[!is.na(topgene_list)]
topgene_list <- unique(topgene_list[topgene_list != ""]) 
topgene_list 

i <- which( fData(gset)$GENE %in% topgene_list )
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(exprs(gset)[i, ], scale = "row",
          labRow = fData(gset)$GENE_SYMBOL[i], labCol = gset$description,
          col = mycol, trace = "none", density.info = "none",
          margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()

#######################################
pdf( file= paste(colnames(efit)[2], "Heatmap_top340genes.pdf", sep = "_"),  paper = "a4")
topgene_list <- comp2$GENE[1:340]
topgene_list <- topgene_list[!is.na(topgene_list)]
topgene_list <- unique(topgene_list[topgene_list != ""]) 
topgene_list 

i <- which( fData(gset)$GENE %in% topgene_list )
#mycol <- colorpanel(1000, "blue", "white", "red")
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(exprs(gset)[i, ], scale = "row",
          labRow = fData(gset)$GENE_SYMBOL[i], labCol = gset$description,
          col = mycol, trace = "none", density.info = "none",
          margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()


#######################################
pdf( file= paste(colnames(efit)[2], "Heatmap_top25genes.pdf", sep = "_"),  paper = "a4")
topgene_list <- comp2$GENE[1:25]
topgene_list <- topgene_list[!is.na(topgene_list)]
topgene_list <- unique(topgene_list[topgene_list != ""]) 
topgene_list 

i <- which( fData(gset)$GENE %in% topgene_list )
#mycol <- colorpanel(1000, "blue", "white", "red")
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(exprs(gset)[i, ], scale = "row",
          labRow = fData(gset)$GENE_SYMBOL[i], labCol = gset$description,
          col = mycol, trace = "none", density.info = "none",
          margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()


#######################################
pdf( file= paste(colnames(efit)[3], "Heatmap_top20genes.pdf", sep = "_"),  paper = "a4")
topgene_list <- comp3$GENE[1:20]
topgene_list <- topgene_list[!is.na(topgene_list)]
topgene_list 

i <- which( fData(gset)$GENE %in% topgene_list )
#mycol <- colorpanel(1000, "blue", "white", "red")
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(exprs(gset)[i, ], scale = "row",
          labRow = fData(gset)$GENE_SYMBOL[i], labCol = gset$description,
          col = mycol, trace = "none", density.info = "none",
          margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()

########################

pdf( file= paste(colnames(efit)[4], "Heatmap_top20genes.pdf", sep = "_"),  paper = "a4")
topgene_list <- comp4$GENE_SYMBOL[1:20]
topgene_list <- topgene_list[!is.na(topgene_list)]
topgene_list 

# comp4$GENE[1:20]
# 
# dim(fData(gset)$GENE)
# head(fData(gset)$GENE)
# 
# head(fData(gset))
# 
# dim(exprs(gset))
# head(exprs(gset))
# colnames(exprs(gset))

i <- which( fData(gset)$GENE_SYMBOL %in% topgene_list )
#rownames(fData(gset))

#mycol <- colorpanel(1000, "blue", "white", "red")
mycol <- colorpanel(100, "blue", "white", "red")
heatmap.2(exprs(gset)[i, ], scale = "row",
          labRow = fData(gset)$GENE_SYMBOL[i], labCol = gset$description,
          col = mycol, trace = "none", density.info = "none",
          margin=c(12,10), lhei=c(2,10), dendrogram="both",cexRow=0.8, cexCol=1.2)
dev.off()
################################################################
########################################------------------
#generate boxplot for special interest genes:

Apoptotic_Genes <- c("Ecel1", "Hrk", "Bbc3", "Puma", "Bcl2l11", "Bim", "Ddit3", "Chop") #upregulated by ONC
v_idx <- which(fData(gset)$GENE_SYMBOL %in% Apoptotic_Genes)

pdf( file= "BoxPlot-Apoptotic_Genes.pdf",  paper = "a4")

for (i in 1:length(v_idx) ) {
  boxplot(exprs(gset)[v_idx[i], ]~group,col="lightblue",
          main= paste(fData(gset)$GENE_SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.7), ylab = "logEXP", xlab = "Sample Groups"
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

v_idx <- which(fData(gset)$GENE_SYMBOL %in% Neuronal_and_Differentiation_Genes)

pdf( file= "BoxPlot-Neuronal_and_Differentiation_Genes.pdf",  paper = "a4")
#g_list <- which(v$genes$ENTREZID %in% comp1$ENTREZID[1:10]) 
for (i in 1:length(v_idx)  ) {
  boxplot(exprs(gset)[v_idx[i], ]~group,col="lightblue",
          main= paste(fData(gset)$GENE_SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.7), ylab = "logEXP", xlab = "Sample Groups"
  )
}

dev.off()


#upregulated by ONC
Regeneration_Genes <- c("Sprr1a", "Sox11", "Klf6", "Atf3")
v_idx <- which(fData(gset)$GENE_SYMBOL %in% Regeneration_Genes)

pdf( file= "BoxPlot-Regeneration_Genes.pdf",  paper = "a4")
for (i in 1:length(v_idx)  ) {
  boxplot(exprs(gset)[v_idx[i], ]~group,col="lightblue",
          main= paste(fData(gset)$GENE_SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.7), ylab = "logEXP", xlab = "Sample Groups"
  )
}

dev.off()


#upregulated
Inflammation_Genes <- c("Aif1", "Gfap", "C1qc", "C4b", "Mmp12", "Cxcl10")
v_idx <- which(fData(gset)$GENE_SYMBOL %in% Inflammation_Genes)

pdf( file= "BoxPlot-Inflammation_Genes.pdf",  paper = "a4")
for (i in 1:length(v_idx)  ) {
  boxplot(exprs(gset)[v_idx[i], ]~group,col="lightblue",
          main= paste(fData(gset)$GENE_SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.7), ylab = "logEXP", xlab = "Sample Groups"
  )
}

dev.off()

DLK_pathway_Genes <- c("Dlk", "Map3k12", "Mlk1", "Map3k9", "Mlk2", "Map3k10", 
                       "Mlk3", "Map3k11"," Mkk4", "Map2k4", "Mkk7", "Map2k7", "Jun")
v_idx <- which(fData(gset)$GENE_SYMBOL %in% DLK_pathway_Genes)

pdf( file= "BoxPlot-DLK_pathway_Genes.pdf",  paper = "a4")
for (i in 1:length(v_idx)  ) {
  boxplot(exprs(gset)[v_idx[i], ]~group,col="lightblue",
          main= paste(fData(gset)$GENE_SYMBOL[v_idx[i]], "Gene Expression", sep=" "),
          par(cex.axis=0.7), ylab = "logEXP", xlab = "Sample Groups"
  )
}

dev.off()


# Test for enrichment of gene sets in Gene Ontology terms---------------------------
library("Mus.musculus") # Gene annotations for the Mus musculus genome

go <- goana(efit, coef =1, geneid = "GENE", FDR = 0.05, species="Mm")

gos1 <- topGO(go, n=500, truncate.term = 40 )
gos2 <- gos1[ which (gos1$P.Up < 0.01 | gos1$P.Down < 0.01), ]


write.table(gos2, file= paste("GeneGOEnrich_", colnames(contr.matrix)[1], ".txt", 
                              sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

go <- goana(efit, coef =2, geneid = "GENE", FDR = 0.05, species="Mm")
gos1 <- topGO(go, n=500,  truncate.term = 40 )
gos2 <- gos1[ which (gos1$P.Up < 0.01 | gos1$P.Down < 0.01), ]

write.table(gos2, file= paste("GeneGOEnrich_", colnames(contr.matrix)[2], ".txt", 
                              sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

go <- goana(efit, coef =3, geneid = "GENE", FDR = 0.05, species="Mm")
#topGO(go, n=15)
gos1 <- topGO(go, n=500,  truncate.term = 40 )
gos2 <- gos1[ which (gos1$P.Up < 0.01 | gos1$P.Down < 0.01), ]


write.table(gos2, file= paste("GeneGOEnrich_", colnames(contr.matrix)[3], ".txt", 
                              sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

go <- goana(efit, coef =4, geneid = "GENE", FDR = 0.05, species="Mm")
gos1 <- topGO(go, n=500,  truncate.term = 40 )
gos2 <- gos1[ which (gos1$P.Up < 0.01 | gos1$P.Down < 0.01), ]

write.table(gos2, file= paste("GeneGOEnrich5_", colnames(contr.matrix)[4], ".txt", 
                              sep =""), col.names=TRUE, row.names=TRUE, sep="\t")


#write.table(topGO(go, n=25, truncate.term = 34), file= paste("GeneGOEnrich_", colnames(contr.matrix)[4], ".txt", 
#                                         sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

#KEGG pathway enrichment analysis
#keg <- kegga(efit, coef =1, species="Mm")

keg <- kegga(efit, coef =1, geneid = "GENE", FDR = 0.05, species="Mm")
#topKEGG(keg, n=15, truncate=34)
kegs1 <- topKEGG(keg, n=100,  truncate.path = 40 )
kegs2 <- kegs1[ which (kegs1$P.Up < 0.01 | kegs1$P.Down < 0.01), ]


write.table(kegs2, file= paste("KEGGenrich_", colnames(contr.matrix)[1], ".txt", 
                               sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

keg <- kegga(efit, coef =2, geneid = "GENE", FDR = 0.05, species="Mm")
#topKEGG(keg, n=15, truncate=34)
kegs1 <- topKEGG(keg, n=100,  truncate.path = 40 )
kegs2 <- kegs1[ which (kegs1$P.Up < 0.01 | kegs1$P.Down < 0.01), ]


write.table(kegs2, file= paste("KEGGenrich_", colnames(contr.matrix)[2], ".txt", 
                               sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

keg <- kegga(efit, coef =3, geneid = "GENE", FDR = 0.05, species="Mm")
#topKEGG(keg, n=15, truncate=34)
kegs1 <- topKEGG(keg, n=100,  truncate.path = 40 )
kegs2 <- kegs1[ which (kegs1$P.Up < 0.01 | kegs1$P.Down < 0.01), ]

write.table(kegs2, file= paste("KEGGenrich_", colnames(contr.matrix)[3], ".txt", 
                               sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

keg <- kegga(efit, coef =4, geneid = "GENE", FDR = 0.05,  species="Mm")
kegs1 <- topKEGG(keg, n=100,  truncate.path = 40 )
kegs2 <- kegs1[ which (kegs1$P.Up < 0.01 | kegs1$P.Down < 0.01), ]

write.table(kegs2, file= paste("KEGGenrich_", colnames(contr.matrix)[4], ".txt", 
                               sep =""), col.names=TRUE, row.names=TRUE, sep="\t")

# Report session information -----------------------------------------
sessionInfo()

 
 
