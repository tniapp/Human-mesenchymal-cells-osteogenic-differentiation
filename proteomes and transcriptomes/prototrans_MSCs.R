####This it the code to reproduce the statistical data analysis from
####"Human mesenchymal cells retain the specificity of their embryonal origin both before and after osteogenic differentiation in vitro"
####code author: Lobov A.A.
####code includes the comparative analysis of proteomes and transcriptomes between cells with different embryological origin

#####Shotgun proteomics data----------------------------------------------------
###Opening the data
## Protein expression data (output from FragPipe)
dat_prot <- data.frame(read.table("combined_protein_mscs.tsv", sep = '\t', header = TRUE))
dat_prot <- dat_prot[-c(1:5),] #removing contaminants
dat_prot <- dat_prot[dat_prot$Combined.Total.Peptides > 1,] # removing the proteins identified by 1 peptide
rownames(dat_prot) <- make.names(dat_prot[,4], unique=TRUE)

#We will use MaxLFQ intensity for the analysis
dat_prot <- dat_prot[,c(215:264)]

#sample information file
library(readxl)
fact <- read_xlsx("prot_samples.xlsx")
fact <- as.data.frame(fact)

rownames(fact) <- fact$Sample
fact$Cell_type <- as.factor(fact$Cell_type)
fact$Cell_type

fact$Differentiation <- as.factor(fact$Differentiation)
fact$Differentiation

fact$Cluster <- as.factor(fact$Cluster)
fact$Cluster

fact$Cluster_2 <- as.factor(fact$Cluster_2)
fact$Cluster_2

fact$Dif_cells <- as.factor(fact$Dif_cells)
fact$Dif_cells

colnames(dat_prot) <- rownames(fact)

### Data preparation and multivariate analysis. We analyse control and differentiated cells separately
## Control
dat_prot_contr <- dat_prot[,fact$Differentiation == "Control"]
fact_prot_contr <- fact[fact$Differentiation == "Control",]
dat_prot_contr[dat_prot_contr == 0] <- NA

mean(complete.cases(dat_prot_contr))
dat_prot_contr1 <- dat_prot_contr[which(rowMeans(!is.na(dat_prot_contr[,fact_prot_contr$Cell_type == "ADSCs"])) >= 0.5), ]
dat_prot_contr1 <- dat_prot_contr1[which(rowMeans(!is.na(dat_prot_contr1[,fact_prot_contr$Cell_type == "DPSCs"])) >= 0.5), ]
dat_prot_contr1 <- dat_prot_contr1[which(rowMeans(!is.na(dat_prot_contr1[,fact_prot_contr$Cell_type == "GFs"])) >= 0.5), ]
dat_prot_contr1 <- dat_prot_contr1[which(rowMeans(!is.na(dat_prot_contr1[,fact_prot_contr$Cell_type == "OBs"])) >= 0.5), ]
dat_prot_contr1 <- dat_prot_contr1[which(rowMeans(!is.na(dat_prot_contr1[,fact_prot_contr$Cell_type == "PDLSCs"])) >= 0.5), ]
dat_prot_contr1 <- dat_prot_contr1[which(rowMeans(!is.na(dat_prot_contr1[,fact_prot_contr$Cell_type == "WJ-MSCs"])) >= 0.5), ]
mean(complete.cases(dat_prot_contr1))

#knn imputation of missng values
library(impute)
tdat <- t(dat_prot_contr1)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))

library(RColorBrewer)
pal <- brewer.pal(n = 6, name = "Set1")
cols <- pal[fact_prot_contr$Cell_type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact_prot_contr$Cell_type), fill = pal, bty = "n", xpd = T)

#log transformation 
dat_log <- log2(dat_knn+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact_prot_contr$Cell_type), fill = pal, bty = "n", xpd = T)

##Data normalization
library(vsn)
library(RColorBrewer)
meanSdPlot(as.matrix(dat_log))

library(limma)
#quantile normalization
dat_prot_contr <- normalizeQuantiles(dat_log)

boxplot(dat_prot_contr, outline = TRUE, col = cols, names=colnames(dat_prot_contr))
meanSdPlot(as.matrix(dat_prot_contr))

## multivariate analysis
#nMDS
library(vegan)
set.seed(125)
nmds <- metaMDS(t(dat_prot_contr), distance = "euclidean")
nmds

data.scores <- as.data.frame(scores(nmds)$sites)

data.scores$group <- fact_prot_contr$Cluster
head(data.scores)

library(ggplot2)
#tiff('dda_nMDS.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = group))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "group", y = "NMDS2") + scale_colour_manual(values = pal)
#dev.off()

#PCA
library(mixOmics)

PCA_contr <- pca(t(dat_prot_contr), ncomp = 8)
plot(PCA_contr)

PCA_contr <- pca(t(dat_prot_contr), ncomp = 3)
#tiff('dda_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(PCA_contr, comp = c(1,2), legend = TRUE,
          group = fact_prot_contr$Cluster,
          ellipse = TRUE,
          title = 'DDA PCA comp 1 - 2')
#dev.off()

#sPLS-DA
ordination.optimum.splsda <- splsda(t(dat_prot_contr), fact_prot_contr$Cluster, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('prot_contr_splsda_comp.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
#dev.off()
layout(1,1)

#tiff('prot_contr_splsda_graph.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, pch = fact_prot_contr$Cell_type, ind.names = F, ellipse = T, abline = TRUE, title = "PLS-DA ordination", legend=TRUE)
#dev.off()

## Differentiated cells
dat_prot_dif <- dat_prot[,fact$Differentiation == "Differentiation"]
fact_prot_dif <- fact[fact$Differentiation == "Differentiation",]
dat_prot_dif[dat_prot_dif == 0] <- NA

mean(complete.cases(dat_prot_dif))
dat_prot_dif1 <- dat_prot_dif[which(rowMeans(!is.na(dat_prot_dif[,fact_prot_dif$Cell_type == "ADSCs"])) >= 0.5), ]
dat_prot_dif1 <- dat_prot_dif1[which(rowMeans(!is.na(dat_prot_dif1[,fact_prot_dif$Cell_type == "DPSCs"])) >= 0.5), ]
dat_prot_dif1 <- dat_prot_dif1[which(rowMeans(!is.na(dat_prot_dif1[,fact_prot_dif$Cell_type == "GFs"])) >= 0.5), ]
dat_prot_dif1 <- dat_prot_dif1[which(rowMeans(!is.na(dat_prot_dif1[,fact_prot_dif$Cell_type == "OBs"])) >= 0.5), ]
dat_prot_dif1 <- dat_prot_dif1[which(rowMeans(!is.na(dat_prot_dif1[,fact_prot_dif$Cell_type == "PDLSCs"])) >= 0.5), ]
dat_prot_dif1 <- dat_prot_dif1[which(rowMeans(!is.na(dat_prot_dif1[,fact_prot_dif$Cell_type == "WJ-MSCs"])) >= 0.5), ]
mean(complete.cases(dat_prot_dif1))

#knn imputation of missng values
library(impute)
tdat <- t(dat_prot_dif1)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))

cols <- pal[fact_prot_dif$Cell_type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact_prot_dif$Cell_type), fill = pal, bty = "n", xpd = T)

#log transformation 
dat_log <- log2(dat_knn+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact_prot_dif$Cell_type), fill = pal, bty = "n", xpd = T)

##Data normalization
meanSdPlot(as.matrix(dat_log))

#quantile normalization
dat_prot_dif <- normalizeQuantiles(dat_log)

boxplot(dat_prot_dif, outline = TRUE, col = cols, names=colnames(dat_prot_dif))
meanSdPlot(as.matrix(dat_prot_dif))

## multivariate analysis
#nMDS
library(vegan)
set.seed(125)
nmds <- metaMDS(t(dat_prot_dif), distance = "euclidean")
nmds

data.scores <- as.data.frame(scores(nmds)$sites)

data.scores$group <- fact_prot_dif$Cluster
head(data.scores)

library(ggplot2)
#tiff('prot_dif_nMDS.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = group))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "group", y = "NMDS2") + scale_colour_manual(values = pal)
#dev.off()

#PCA
library(mixOmics)

PCA_dif <- pca(t(dat_prot_dif), ncomp = 8)
plot(PCA_dif)

PCA_dif <- pca(t(dat_prot_dif), ncomp = 3)
#tiff('dda_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(PCA_dif, comp = c(1,2), legend = TRUE,
          group = fact_prot_dif$Cluster,
          ellipse = TRUE,
          title = 'DDA PCA comp 1 - 2')
#dev.off()

#sPLS-DA 
ordination.optimum.splsda <- splsda(t(dat_prot_dif), fact_prot_dif$Cluster, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('PROT_dif_splsda1_comp.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", difib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
#dev.off()
layout(1,1)

#tiff('PROT_dif_splsda1_graph.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, pch = fact_prot_dif$Cell_type, ellipse = T, abline = TRUE, title = "PLS-DA ordination", legend=TRUE)
#dev.off()


library(VennDiagram)
venn.diagram(
  x = list(rownames(dat_prot_contr), rownames(dat_prot_dif)),
  category.names = c("Contr" , "Dif"),
  filename = '#unique_datasets_prot.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)



###Analysis of differential expression
##control
design = model.matrix(~0+fact_prot_contr$Cluster)
colnames(design) = c('Fetus', 'Mesenchyme', 'Neural_crest')

#To make all pair-wise comparisons between all four groups one could proceed

fit_contr <- lmFit(dat_prot_contr, design)
contrasts_groups = c('(Mesenchyme+Neural_crest)/2-Fetus','(Mesenchyme+Fetus)/2-Neural_crest','(Neural_crest+Fetus)/2-Mesenchyme')
contrast.matrix <- makeContrasts(contrasts = contrasts_groups,levels=design)

fit2_contr <- contrasts.fit(fit_contr, contrast.matrix)
fit2_contr <- eBayes(fit2_contr)

results_contr <- decideTests(fit2_contr)
#tiff('dif_expr_overlap_prot_contr.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
vennDiagram(results_contr)
#dev.off()

Anova_contr_prot <- topTable(fit2_contr, number=as.numeric(length(dat_prot_contr[,1])))
#write.csv(Anova_contr_prot, "Anova_contr_prot.csv")

Mesenchyme_vs_rest_contr <- topTable(fit2_contr, coef=3, adjust="BH", number = as.numeric(length(dat_prot_contr[,1])))
#write.csv(Mesenchyme_vs_rest_contr, "Mesenchyme_vs_rest_contr.csv")
Neural_crest_vs_resr_contr <- topTable(fit2_contr, coef=2, adjust="BH", number = as.numeric(length(dat_prot_contr[,1])))
#write.csv(Neural_crest_vs_resr_contr, "Neural_crest_vs_resr_contr.csv")
Fetus_vs_rest_contr <- topTable(fit2_contr, coef=1, adjust="BH", number = as.numeric(length(dat_prot_contr[,1])))
#write.csv(Fetus_vs_rest_contr, "Fetus_vs_rest_contr.csv")

Anova_contr_prot_top <- topTable(fit2_contr, number=40)

proteins_dif_prot_cont <- rownames(Anova_contr_prot_top)

proteins_dif_prot_cont <- dat_prot_contr[proteins_dif_prot_cont,]

library(pheatmap)
fact_prot_contr1 <- fact_prot_contr[,c(2,5)]
#tiff('heatmap_contr_full_prot.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
pheatmap::pheatmap(proteins_dif_prot_cont,scale="row", annotation_col = fact_prot_contr1,
                   color=colorRampPalette(c("navy", "white", "red"))(50))
#dev.off()

Anova_contr_prot_rel_full <- Anova_contr_prot[Anova_contr_prot$adj.P.Val <= 0.05,] #Will use for pathway analysis in the end

##differentiated cells
design = model.matrix(~0+fact_prot_dif$Cluster)
colnames(design) = c('Fetus', 'Mesenchyme', 'Neural_crest')

fit_dif <- lmFit(dat_prot_dif, design)
contrasts_groups = c('(Mesenchyme+Neural_crest)/2-Fetus','(Mesenchyme+Fetus)/2-Neural_crest','(Neural_crest+Fetus)/2-Mesenchyme')
contrast.matrix <- makeContrasts(contrasts = contrasts_groups,levels=design)

fit2_dif <- contrasts.fit(fit_dif, contrast.matrix)
fit2_dif <- eBayes(fit2_dif)

results_dif <- decideTests(fit2_dif)
#tiff('dif_expr_overlap_prot_dif.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
vennDiagram(results_dif)
#dev.off()

Anova_dif_prot <- topTable(fit2_dif, number=as.numeric(length(dat_prot_dif[,1])))
#write.csv(Anova_dif_prot, "Anova_dif_prot.csv")

Mesenchyme_vs_rest_dif <- topTable(fit2_dif, coef=3, adjust="BH", number = as.numeric(length(dat_prot_dif[,1])))
#write.csv(Mesenchyme_vs_rest_dif, "Mesenchyme_vs_rest_dif.csv")
Neural_crest_vs_resr_dif <- topTable(fit2_dif, coef=2, adjust="BH", number = as.numeric(length(dat_prot_dif[,1])))
#write.csv(Neural_crest_vs_resr_dif, "Neural_crest_vs_resr_dif.csv")
Fetus_vs_rest_dif <- topTable(fit2_dif, coef=1, adjust="BH", number = as.numeric(length(dat_prot_dif[,1])))
#write.csv(Fetus_vs_rest_dif, "Fetus_vs_rest_dif.csv")

Anova_dif_prot_top <- topTable(fit2_dif, number=40)

proteins_dif_prot_dif <- rownames(Anova_dif_prot_top)
proteins_dif_prot_dif <- dat_prot_dif[proteins_dif_prot_dif,]

library(pheatmap)
fact_prot_dif1 <- fact_prot_dif[,c(2,5)]
#tiff('heatmap_diff_full_prot.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
pheatmap::pheatmap(proteins_dif_prot_dif,scale="row", annotation_col = fact_prot_dif1,
                   color=colorRampPalette(c("navy", "white", "red"))(50))
#dev.off()

Anova_dif_prot_rel_full <- Anova_dif_prot[Anova_dif_prot$adj.P.Val <= 0.05,] #Will use for pathway analysis in the end

# Venn diagramm of comparison of dif. expr. proteins
VennDiagram::venn.diagram(
  x = list(rownames(Anova_dif_prot_rel_full), rownames(Anova_contr_prot_rel_full)),
  category.names = c("dif" , "contr"),
  filename = '#dif_prot_contr_and_dif.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)


####RNA-seq data----------------------------------------------------------------
## Opening the data
#upload initial table with counts
dat <- read.table("RNA-seq_counts.csv", sep = ',', header = TRUE)
dat <- dat[,-1]

#Convert ENSG to gene name (optional)
genes <- dat$GENE
library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'), 
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)
counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}
rownames(dat) <- make.names(genes_symbol$symbols, unique=TRUE)
dat <- dat[,-1]
head(dat)

#filtering null-value genes
mean(complete.cases(dat))
dat <- dat[rowMeans(dat) > 0,]
dat_f <- dat[rowSums(dat) > 50,] #filtering low-level genes

#upload a sample information
library(readxl)
fact <- as.data.frame(read_xlsx("RNA-seq_samples.xlsx", sheet = 1))

rownames(fact) <- fact$sample
fact$Cell_type <- as.factor(fact$Cell_type)
fact$Cell_type
fact$Cluster <- as.factor(fact$Cluster)
fact$Cluster
fact$Differentiation <- as.factor(fact$Differentiation)
fact$Differentiation
fact$condition <- as.factor(fact$condition)
fact$replicate <- as.factor(fact$replicate)
fact$Donor <- as.factor(fact$Donor)

fact$Mesenchyme_rest <- as.factor(fact$Mesenchyme_rest)
fact$Mesenchyme_rest
fact$Crest_rest <- as.factor(fact$Crest_rest)
fact$Crest_rest
fact$Fetus_rest <- as.factor(fact$Fetus_rest)
fact$Fetus_rest
#View(data.frame(colnames(dat), rownames(fact)))  # check names for silimarities
colnames(dat_f) <- rownames(fact)

#Data prep and split to three datasets
library(DESeq2)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = dat_f,
  colData = fact,
  design = ~ Cluster)
colData(ddsFullCountTable)[1:5,1:4]

#collapsing technical replicates 
ddsCountTable <- collapseReplicates(ddsFullCountTable,
                                    groupby = ddsFullCountTable$Sample_1,
                                    run = ddsFullCountTable$replicate)

#split dataset to three parts: control, dif 48 hours and dif 10 days
dds_contr <- ddsCountTable[ , ddsCountTable$Differentiation == "Control" ]
dds_48h <- ddsCountTable[ , ddsCountTable$Differentiation == "Dif_48h" ]
dds_10d <- ddsCountTable[ , ddsCountTable$Differentiation == "Dif_10d" ]

###Control
## multivariate analysis
rld_control <- rlog(dds_contr) # regularized log transformation

#PCA
library(ggplot2)
library(mixOmics)

pca_contr <- pca(t(assay(rld_control)), ncomp = 8)
plot(pca_contr)

pca_contr <- pca(t(assay(rld_control)), ncomp = 3)
#tiff('RNA-seq_contr_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(pca_contr, comp = c(1,2), legend = TRUE,
          group = rld_control$Cluster,
          pch = rld_control$Cell_type,
          ellipse = TRUE,
          title = 'DDA PCA comp 1 - 2')
#dev.off()

#PLS-DA 
ordination.optimum.splsda <- splsda(t(assay(rld_control)), rld_control$Cluster, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('RNAseq_contr_splsda1_comp.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", difib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
#dev.off()
layout(1,1)

#tiff('RNAseq_contr_splsda1_graph.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, pch = rld_control$Cell_type, ind.names = F, ellipse = T, abline = TRUE, 
          title = "PLS-DA ordination", legend=TRUE)
#dev.off()

##DEseq2
deseq_contr <- DESeq(dds_contr, test = "LRT", reduced = ~ 1) #analogous to ANOVA
res_contr <- results(deseq_contr)
res_contr
res_contr$padj[is.na(res_contr$padj)] <- 1 #Removing NA in padj

res_contr_rel_ANOVA <- res_contr[res_contr$padj <= 0.01,]
write.csv( as.data.frame(res_contr), file="RNA-seq_res_contr_ANOVA.csv" )

#Diagnostic plots
plotDispEsts(deseq_contr )

hist(res_contr$pvalue, breaks=20, col="grey" )

qs <- c( 0, quantile(res_contr$baseMean[res_contr$baseMean > 0], 0:7/7 ) )
bins <- cut(res_contr$baseMean, qs )
ratios <- tapply(res_contr$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

#Heatmap
res_contr_rel_ANOVA <- res_contr_rel_ANOVA[ order( res_contr_rel_ANOVA$padj ), ]
res_contr_rel_ANOVA_top <- res_contr_rel_ANOVA[1:40,]
res_contr_rel_ANOVA_top <- rownames(res_contr_rel_ANOVA_top)

count_contr <- data.frame(assay(dds_contr))
res_contr_rel_ANOVA_top <- count_contr[res_contr_rel_ANOVA_top,]
library(pheatmap)
fact_contr1 <- data.frame(dds_contr$Cell_type, dds_contr$Cluster)
rownames(fact_contr1) <- dds_contr$sample
colnames(fact_contr1) <- c("Cell_type", "Cluster")
#manually replace ENSG names of top-genes based on data from Ensembl, e.g. ENSG00000255438 is annotated as Novel Transcript, Antisense To SULF2 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=ENSG00000255438&keywords=ENSG00000255438; accessed 14.02.2023)
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000286502"] <- "COL5A1.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000255438"] <- "SULF2.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000249996"] <- "PPIC.AS1"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000227279"] <- "CDH2.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000266283"] <- "GATA6.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000270823"] <- "MEST.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000237286"] <- "CARD11.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000232295"] <- "OGFRL1.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000267774"] <- "CCBE1.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000287724"] <- "DHCR24.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000256013"] <- "EMP2.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000237380"] <- "HOXD.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000271774"] <- "TSHZ2.AS"
rownames(res_contr_rel_ANOVA_top)[rownames(res_contr_rel_ANOVA_top) == "ENSG00000285935"] <- "WASF3.AS"

#tiff('contr_RNA-seq_heatmap.tiff', units="in", width=7, height=9, res=300, compression = 'lzw')
pheatmap::pheatmap(res_contr_rel_ANOVA_top,scale="row", annotation_col = fact_contr1,
                   color=colorRampPalette(c("navy", "white", "red"))(50))
#dev.off()

# Identifying genes specific for some type against the others
res_contr_fetus_spec <- results(deseq_contr, contrast = c(1,-0.5,-0.5) )
res_contr_mes_spec <- results(deseq_contr, contrast = c(-0.5,1,-0.5) )
res_contr_crest_spec <- results(deseq_contr, contrast = c(-0.5,-0.5,1) )

res_contr_fetus_spec$padj[is.na(res_contr_fetus_spec$padj)] <- 1 #Removing NA in padj
res_contr_mes_spec$padj[is.na(res_contr_mes_spec$padj)] <- 1 #Removing NA in padj
res_contr_crest_spec$padj[is.na(res_contr_crest_spec$padj)] <- 1 #Removing NA in padj

res_contr_fetus_spec <- res_contr_fetus_spec[res_contr_fetus_spec$padj <= 0.05,]
res_contr_mes_spec <- res_contr_mes_spec[res_contr_mes_spec$padj <= 0.05,]
res_contr_crest_spec <- res_contr_crest_spec[res_contr_crest_spec$padj <= 0.05,]

res_contr_fetus_spec_up <- res_contr_fetus_spec[res_contr_fetus_spec$log2FoldChange >= 1 ,]
res_contr_fetus_spec_dn <- res_contr_fetus_spec[res_contr_fetus_spec$log2FoldChange <= -1 ,]

res_contr_mes_spec_up <- res_contr_mes_spec[res_contr_mes_spec$log2FoldChange >= 1 ,]
res_contr_mes_spec_dn <- res_contr_mes_spec[res_contr_mes_spec$log2FoldChange <= -1 ,]

res_contr_crest_spec_up <- res_contr_crest_spec[res_contr_crest_spec$log2FoldChange >= 1 ,]
res_contr_crest_spec_dn <- res_contr_crest_spec[res_contr_crest_spec$log2FoldChange <= -1 ,]

###Dif_48
## multivariate analysis
rld_d48 <- rlog(dds_48h)

#PCA
library(ggplot2)
library(mixOmics)

pca_48h <- pca(t(assay(rld_d48)), ncomp = 8)
plot(pca_48h)

pca_48h <- pca(t(assay(rld_d48)), ncomp = 3)
#tiff('RNA-seq_pca_48h.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(pca_48h, comp = c(1,2), legend = TRUE,
          group = rld_d48$Cluster,
          pch = rld_d48$Cell_type,
          ellipse = TRUE,
          title = 'pca_48h',
          ind.names = F)
#dev.off()


#PLS-DA 
ordination.optimum.splsda <- splsda(t(assay(rld_d48)), rld_d48$Cluster, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('RNAseq_48h_splsda1_comp.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", difib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
#dev.off()
layout(1,1)

#tiff('RNAseq_48h_splsda1_graph.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, pch = rld_d48$Cell_type, ind.names = F, ellipse = T, abline = TRUE, 
          title = "PLS-DA ordination", legend=TRUE)
#dev.off()

##DEseq2
deseq_48h <- DESeq(dds_48h, test = "LRT", reduced = ~ 1)
res_48h <- results(deseq_48h)
res_48h
res_48h$padj[is.na(res_48h$padj)] <- 1

res_48h_rel_ANOVA <- res_48h[res_48h$padj <= 0.01,]
#write.csv( as.data.frame(res_48h), file="RNA-seq_48h_ANOVA.csv" )

#Diagnostic plots
plotDispEsts(deseq_48h)

hist(res_48h$pvalue, breaks=20, col="grey" )

qs <- c( 0, quantile(res_48h$baseMean[res_48h$baseMean > 0], 0:7/7 ) )
bins <- cut(res_48h$baseMean, qs )
ratios <- tapply(res_48h$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

#Heatmap
res_48h_rel_ANOVA <- res_48h_rel_ANOVA[ order( res_48h_rel_ANOVA$padj ), ]
res_48h_rel_ANOVA_top <- res_48h_rel_ANOVA[1:50,]
res_48h_rel_ANOVA_top <- rownames(res_48h_rel_ANOVA_top)

count_48h <- data.frame(assay(dds_48h))
res_48h_rel_ANOVA_top <- count_48h[res_48h_rel_ANOVA_top,]
library(pheatmap)
fact_48h1 <- data.frame(dds_48h$Cell_type, dds_48h$Cluster)
rownames(fact_48h1) <- dds_48h$sample
colnames(fact_48h1) <- c("Cell_type", "Cluster")
#manually replace ENSG names of top-genes based on data from Ensembl, e.g. ENSG00000255438 is annotated as Novel Transcript, Antisense To SULF2 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=ENSG00000255438&keywords=ENSG00000255438; accessed 14.02.2023)
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000267787"] <- "ATP8B1.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000267774"] <- "CCBE1.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000256013"] <- "EMP2.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000261054"] <- "SYNM.AS2"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000232295"] <- "OGFRL1.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000245330"] <- "DEPTOR.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000287528"] <- "FZD8.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000285935"] <- "WASF3.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000232415"] <- "ELN.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000237380"] <- "HOXD.AS"
rownames(res_48h_rel_ANOVA_top)[rownames(res_48h_rel_ANOVA_top) == "ENSG00000227279"] <- "CDH2.AS"

#tiff('48h_RNA-seq_heatmap1.tiff', units="in", width=7, height=9, res=300, compression = 'lzw')
pheatmap::pheatmap(res_48h_rel_ANOVA_top,scale="row", annotation_col = fact_48h1,
                   color=colorRampPalette(c("navy", "white", "red"))(50))
#dev.off()

# Identifying genes specific for some type against the others
res_48h_fetus_spec <- results(deseq_48h, contrast = c(1,-0.5,-0.5) )
res_48h_mes_spec <- results(deseq_48h, contrast = c(-0.5,1,-0.5) )
res_48h_crest_spec <- results(deseq_48h, contrast = c(-0.5,-0.5,1) )

res_48h_fetus_spec$padj[is.na(res_48h_fetus_spec$padj)] <- 1 #Removing NA in padj
res_48h_mes_spec$padj[is.na(res_48h_mes_spec$padj)] <- 1 #Removing NA in padj
res_48h_crest_spec$padj[is.na(res_48h_crest_spec$padj)] <- 1 #Removing NA in padj

res_48h_fetus_spec <- res_48h_fetus_spec[res_48h_fetus_spec$padj <= 0.05,]
res_48h_mes_spec <- res_48h_mes_spec[res_48h_mes_spec$padj <= 0.05,]
res_48h_crest_spec <- res_48h_crest_spec[res_contr_crest_spec$padj <= 0.05,]

res_48h_fetus_spec_up <- res_48h_fetus_spec[res_48h_fetus_spec$log2FoldChange >= 1 ,]
res_48h_fetus_spec_dn <- res_48h_fetus_spec[res_48h_fetus_spec$log2FoldChange <= -1 ,]

res_48h_mes_spec_up <- res_48h_mes_spec[res_48h_mes_spec$log2FoldChange >= 1 ,]
res_48h_mes_spec_dn <- res_48h_mes_spec[res_48h_mes_spec$log2FoldChange <= -1 ,]

res_48h_crest_spec_up <- res_48h_crest_spec[res_48h_crest_spec$log2FoldChange >= 1 ,]
res_48h_crest_spec_dn <- res_48h_crest_spec[res_48h_crest_spec$log2FoldChange <= -1 ,]

###Dif_10d
## multivariate analysis
rld_10d <- rlog(dds_10d)

#PCA
library(ggplot2)
library(mixOmics)

pca_10d <- pca(t(assay(rld_10d)), ncomp = 8)
plot(pca_10d)

pca_10d <- pca(t(assay(rld_10d)), ncomp = 3)
#tiff('RNA-seq_pca_10d.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(pca_10d, comp = c(1,2), legend = TRUE,
          group = rld_10d$Cluster,
          pch = rld_10d$Cell_type,
          ellipse = TRUE,
          title = 'pca_10d',
          ind.names = F)
# dev.off()

#PLS-DA 
ordination.optimum.splsda <- splsda(t(assay(rld_10d)), rld_10d$Cluster, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('RNAseq_10d_splsda1_comp.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", difib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
#dev.off()
layout(1,1)

#tiff('RNAseq_10d_splsda1_graph.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, pch = rld_10d$Cell_type, ind.names = F, ellipse = T, abline = TRUE, 
          title = "PLS-DA ordination", legend=TRUE)
#dev.off()

##DEseq2
deseq_10d <- DESeq(dds_10d, test = "LRT", reduced = ~ 1)
res_10d <- results(deseq_10d)
res_10d
res_10d$padj[is.na(res_10d$padj)] <- 1

res_10d_rel_ANOVA <- res_10d[res_10d$padj <= 0.01,]
#write.csv( as.data.frame(res_10d), file="RNA-seq_res_10d_ANOVA.csv" )

#Diagnostic plots
plotDispEsts(deseq_10d)

hist(res_10d$pvalue, breaks=20, col="grey" )

qs <- c( 0, quantile(res_10d$baseMean[res_10d$baseMean > 0], 0:7/7 ) )
bins <- cut(res_10d$baseMean, qs )
ratios <- tapply(res_10d$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

#Heatmap
res_10d_rel_ANOVA <- res_10d_rel_ANOVA[ order( res_10d_rel_ANOVA$padj ), ]
res_10d_rel_ANOVA_top <- res_10d_rel_ANOVA[1:50,]
res_10d_rel_ANOVA_top <- rownames(res_10d_rel_ANOVA_top)

count_10d <- data.frame(assay(dds_10d))
res_10d_rel_ANOVA_top <- count_10d[res_10d_rel_ANOVA_top,]
library(pheatmap)
fact_10d1 <- data.frame(dds_10d$Cell_type, dds_10d$Cluster)
rownames(fact_10d1) <- dds_10d$sample
colnames(fact_10d1) <- c("Cell_type", "Cluster")
#manually replace ENSG names of top-genes based on data from Ensembl, e.g. ENSG00000255438 is annotated as Novel Transcript, Antisense To SULF2 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=ENSG00000255438&keywords=ENSG00000255438; accessed 14.02.2023)
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000232295"] <- "OGFRL1.AS"
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000227279"] <- "CDH2.AS"
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000267774"] <- "CCBE1.AS"
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000237380"] <- "HOXD.AS"
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000228925"] <- "MCFD2.AS"
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000287037"] <- "CXCL5.AS"
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000270823"] <- "MEST.AS2"
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000287387"] <- "ZFP36L2.AS"
rownames(res_10d_rel_ANOVA_top)[rownames(res_10d_rel_ANOVA_top) == "ENSG00000225670"] <- "CADM3.AS"

#tiff('10d_dif_RNA-seq_heatmap.tiff', units="in", width=7, height=9, res=300, compression = 'lzw')
pheatmap::pheatmap(res_10d_rel_ANOVA_top,scale="row", annotation_col = fact_10d1,
                   color=colorRampPalette(c("navy", "white", "red"))(50))
#dev.off()

# Identifying genes specific for some type against the others
res_10d_fetus_spec <- results(deseq_10d, contrast = c(1,-0.5,-0.5) )
res_10d_mes_spec <- results(deseq_10d, contrast = c(-0.5,1,-0.5) )
res_10d_crest_spec <- results(deseq_10d, contrast = c(-0.5,-0.5,1) )

res_10d_fetus_spec$padj[is.na(res_10d_fetus_spec$padj)] <- 1 #Removing NA in padj
res_10d_mes_spec$padj[is.na(res_10d_mes_spec$padj)] <- 1 #Removing NA in padj
res_10d_crest_spec$padj[is.na(res_10d_crest_spec$padj)] <- 1 #Removing NA in padj

res_10d_fetus_spec <- res_10d_fetus_spec[res_10d_fetus_spec$padj <= 0.05,]
res_10d_mes_spec <- res_10d_mes_spec[res_10d_mes_spec$padj <= 0.05,]
res_10d_crest_spec <- res_10d_crest_spec[res_10d_crest_spec$padj <= 0.05,]

res_10d_fetus_spec_up <- res_10d_fetus_spec[res_10d_fetus_spec$log2FoldChange >= 1 ,]
res_10d_fetus_spec_dn <- res_10d_fetus_spec[res_10d_fetus_spec$log2FoldChange <= -1 ,]

res_10d_mes_spec_up <- res_10d_mes_spec[res_10d_mes_spec$log2FoldChange >= 1 ,]
res_10d_mes_spec_dn <- res_10d_mes_spec[res_10d_mes_spec$log2FoldChange <= -1 ,]

res_10d_crest_spec_up <- res_10d_crest_spec[res_10d_crest_spec$log2FoldChange >= 1 ,]
res_10d_crest_spec_dn <- res_10d_crest_spec[res_10d_crest_spec$log2FoldChange <= -1 ,]

# Venn diagramm of comparison of dif. expr. proteins
VennDiagram::venn.diagram(
  x = list(rownames(res_contr_rel_ANOVA), rownames(res_10d_rel_ANOVA)),
  category.names = c("Contr" , "dif"),
  filename = '#dif_RNAseq_contr_and_dif10.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)


####Comparison of proteotranscriptomics data------------------------------------

## Venn diagramm of comparison of dif. expr. proteins
VennDiagram::venn.diagram(
  x = list(rownames(res_contr_rel_ANOVA), rownames(Anova_contr_prot_rel_full)),
  category.names = c("RNA-seq" , "prot"),
  filename = '#dif_RNAseq_vs_prot_contr.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)

VennDiagram::venn.diagram(
  x = list(rownames( res_10d_rel_ANOVA), rownames(Anova_dif_prot_rel_full)),
  category.names = c("RNA" , "Prot"),
  filename = '#dif_RNAseq_vs_prot_dif.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)

## Enrichment analysis
##Control
# Mesenchyme
Mesenchyme_vs_rest_contr_prot_up <- Mesenchyme_vs_rest_contr[Mesenchyme_vs_rest_contr$logFC >= 0.5, ]
Mesenchyme_vs_rest_contr_prot_dn <- Mesenchyme_vs_rest_contr[Mesenchyme_vs_rest_contr$logFC <= -0.5, ]

library("org.Hs.eg.db")
library("clusterProfiler")
mes_up_contr_prot <- bitr(rownames(Mesenchyme_vs_rest_contr_prot_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
mes_dn_contr_prot <- bitr(rownames(Mesenchyme_vs_rest_contr_prot_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

mes_up_contr_rna <- bitr(rownames(res_contr_mes_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
mes_dn_contr_rna <- bitr(rownames(res_contr_mes_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

mes_up_prot_egobp <- enrichGO(
  gene     = mes_up_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_up_rna_egobp <- enrichGO(
  gene     = mes_up_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_prot_egobp <- enrichGO(
  gene     = mes_dn_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_rna_egobp <- enrichGO(
  gene     = mes_dn_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('mes_up_prot_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_up_prot_egobp, showCategory = 15)
dev.off()
tiff('mes_up_rna_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_up_rna_egobp, showCategory = 15)
dev.off()
tiff('mes_dn_prot_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_dn_prot_egobp, showCategory = 15)
dev.off()
tiff('mes_dn_rna_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_dn_rna_egobp, showCategory = 15)
dev.off()

mes_up_prot_egocc <- enrichGO(
  gene     = mes_up_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_up_rna_egocc <- enrichGO(
  gene     = mes_up_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_prot_egocc <- enrichGO(
  gene     = mes_dn_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_rna_egocc <- enrichGO(
  gene     = mes_dn_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('mes_up_prot_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_up_prot_egocc, showCategory = 15)
dev.off()
tiff('mes_up_rna_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_up_rna_egocc, showCategory = 15)
dev.off()
tiff('mes_dn_prot_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_dn_prot_egocc, showCategory = 15)
dev.off()

# Neural_crest
Neural_crest_vs_resr_contr_prot_up <- Neural_crest_vs_resr_contr[Neural_crest_vs_resr_contr$logFC >= 0.5, ]
Neural_crest_vs_resr_contr_prot_dn <- Neural_crest_vs_resr_contr[Neural_crest_vs_resr_contr$logFC <= -0.5, ]

library("org.Hs.eg.db")
library("clusterProfiler")
NC_up_contr_prot <- bitr(rownames(Neural_crest_vs_resr_contr_prot_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
NC_dn_contr_prot <- bitr(rownames(Neural_crest_vs_resr_contr_prot_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

NC_up_contr_rna <- bitr(rownames(res_contr_crest_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
NC_dn_contr_rna <- bitr(rownames(res_contr_crest_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

NC_up_prot_egobp <- enrichGO(
  gene     = NC_up_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_up_rna_egobp <- enrichGO(
  gene     = NC_up_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_prot_egobp <- enrichGO(
  gene     = NC_dn_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_rna_egobp <- enrichGO(
  gene     = NC_dn_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('NC_up_prot_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_prot_egobp, showCategory = 15)
dev.off()
tiff('NC_up_rna_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_rna_egobp, showCategory = 15)
dev.off()
tiff('NC_dn_rna_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_dn_rna_egobp, showCategory = 15)
dev.off()

NC_up_prot_egocc <- enrichGO(
  gene     = NC_up_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_up_rna_egocc <- enrichGO(
  gene     = NC_up_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_prot_egocc <- enrichGO(
  gene     = NC_dn_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_rna_egocc <- enrichGO(
  gene     = NC_dn_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('NC_up_prot_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_prot_egocc, showCategory = 15)
dev.off()
tiff('NC_up_rna_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_rna_egocc, showCategory = 15)
dev.off()
tiff('NC_dn_prot_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_dn_prot_egocc, showCategory = 15)
dev.off()
tiff('NC_dn_rna_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_dn_rna_egocc, showCategory = 15)
dev.off()

# Fetus
Fetus_vs_rest_contr_prot_up <- Fetus_vs_rest_contr[Fetus_vs_rest_contr$logFC >= 0.5, ]
Fetus_vs_rest_contr_prot_dn <- Fetus_vs_rest_contr[Fetus_vs_rest_contr$logFC <= -0.5, ]

library("org.Hs.eg.db")
library("clusterProfiler")
Fet_up_contr_prot <- bitr(rownames(Fetus_vs_rest_contr_prot_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Fet_dn_contr_prot <- bitr(rownames(Fetus_vs_rest_contr_prot_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

Fet_up_contr_rna <- bitr(rownames(res_contr_fetus_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Fet_dn_contr_rna <- bitr(rownames(res_contr_fetus_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

Fet_up_prot_egobp <- enrichGO(
  gene     = Fet_up_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_up_rna_egobp <- enrichGO(
  gene     = Fet_up_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_prot_egobp <- enrichGO(
  gene     = Fet_dn_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_rna_egobp <- enrichGO(
  gene     = Fet_dn_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('Fet_up_prot_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_up_prot_egobp, showCategory = 15)
dev.off()
tiff('Fet_up_rna_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_up_rna_egobp, showCategory = 15)
dev.off()
tiff('Fet_dn_prot_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_dn_prot_egobp, showCategory = 15)
dev.off()
tiff('Fet_dn_rna_egobp.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_dn_rna_egobp, showCategory = 15)
dev.off()

Fet_up_prot_egocc <- enrichGO(
  gene     = Fet_up_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_up_rna_egocc <- enrichGO(
  gene     = Fet_up_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_prot_egocc <- enrichGO(
  gene     = Fet_dn_contr_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_rna_egocc <- enrichGO(
  gene     = Fet_dn_contr_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('Fet_up_prot_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_up_prot_egocc, showCategory = 15)
dev.off()
tiff('Fet_up_rna_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_up_rna_egocc, showCategory = 15)
dev.off()
tiff('Fet_dn_prot_egocc.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_dn_prot_egocc, showCategory = 15)
dev.off()

## dif_48h (RNA-seq only)
# Mesenchyme
library("org.Hs.eg.db")
library("clusterProfiler")
mes_up_48h_rna_up <- bitr(rownames(res_48h_mes_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
mes_dn_48h_rna_dn <- bitr(rownames(res_48h_mes_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

mes_up_rna_egobp_48h <- enrichGO(
  gene     = mes_up_48h_rna_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_rna_egobp_48h <- enrichGO(
  gene     = mes_dn_48h_rna_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_up_rna_egocc_48h <- enrichGO(
  gene     = mes_up_48h_rna_up[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_rna_egocc_48h <- enrichGO(
  gene     = mes_dn_48h_rna_dn[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('mes_dn_rna_egocc_48h.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_dn_rna_egocc_48h, showCategory = 15)
dev.off()
tiff('mes_dn_rna_egobp_48h.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_dn_rna_egobp_48h, showCategory = 15)
dev.off()

# Neural_crest
NC_up_48h_rna <- bitr(rownames(res_48h_crest_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
NC_dn_48h_rna <- bitr(rownames(res_48h_crest_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

NC_up_rna_egobp_48h <- enrichGO(
  gene     = NC_up_48h_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_rna_egobp_48h <- enrichGO(
  gene     = NC_dn_48h_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('NC_up_rna_egobp_48h.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_rna_egobp_48h, showCategory = 15)
dev.off()
tiff('NC_dn_rna_egobp_48h.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_dn_rna_egobp_48h, showCategory = 15)
dev.off()

NC_up_rna_egocc_48h <- enrichGO(
  gene     = NC_up_48h_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_rna_egocc_48h <- enrichGO(
  gene     = NC_dn_48h_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('NC_up_rna_egocc_48h.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_rna_egocc_48h, showCategory = 15)
dev.off()
tiff('NC_dn_rna_egocc_48h.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_dn_rna_egocc_48h, showCategory = 15)
dev.off()

# Fetus
Fet_up_48h_rna <- bitr(rownames(res_48h_fetus_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Fet_dn_48h_rna <- bitr(rownames(res_48h_fetus_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

Fet_up_rna_egobp_48h <- enrichGO(
  gene     = Fet_up_48h_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_rna_egobp_48h <- enrichGO(
  gene     = Fet_dn_48h_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_up_rna_egocc_48h <- enrichGO(
  gene     = Fet_up_48h_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_rna_egocc_48h <- enrichGO(
  gene     = Fet_dn_48h_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('Fet_up_rna_egocc_48h.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_up_rna_egocc_48h, showCategory = 15)
dev.off()

##Differentiation
# Mesenchyme
Mesenchyme_vs_rest_dif_prot_up <- Mesenchyme_vs_rest_dif[Mesenchyme_vs_rest_dif$logFC >= 0.5, ]
Mesenchyme_vs_rest_dif_prot_dn <- Mesenchyme_vs_rest_dif[Mesenchyme_vs_rest_dif$logFC <= -0.5, ]

mes_up_dif_prot <- bitr(rownames(Mesenchyme_vs_rest_dif_prot_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
mes_dn_dif_prot <- bitr(rownames(Mesenchyme_vs_rest_dif_prot_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

mes_up_dif_rna <- bitr(rownames(res_10d_mes_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
mes_dn_dif_rna <- bitr(rownames(res_10d_mes_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

mes_up_prot_egobp_10d <- enrichGO(
  gene     = mes_up_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_up_rna_egobp_10d <- enrichGO(
  gene     = mes_up_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_prot_egobp_10d <- enrichGO(
  gene     = mes_dn_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_rna_egobp_10d <- enrichGO(
  gene     = mes_dn_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('mes_up_prot_egobp_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_up_prot_egobp_10d, showCategory = 15)
dev.off()
tiff('mes_dn_prot_egobp_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_dn_prot_egobp_10d, showCategory = 15)
dev.off()


mes_up_prot_egocc_10d <- enrichGO(
  gene     = mes_up_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_up_rna_egocc_10d <- enrichGO(
  gene     = mes_up_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_prot_egocc_10d <- enrichGO(
  gene     = mes_dn_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

mes_dn_rna_egocc_10d <- enrichGO(
  gene     = mes_dn_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('mes_up_prot_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_up_prot_egocc_10d, showCategory = 15)
dev.off()
tiff('mes_dn_prot_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_dn_prot_egocc_10d, showCategory = 15)
dev.off()
tiff('mes_dn_rna_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(mes_dn_rna_egocc_10d, showCategory = 15)
dev.off()

# Neural_crest
Neural_crest_vs_resr_dif_prot_up <- Neural_crest_vs_resr_dif[Neural_crest_vs_resr_dif$logFC >= 0.5, ]
Neural_crest_vs_resr_dif_prot_dn <- Neural_crest_vs_resr_dif[Neural_crest_vs_resr_dif$logFC <= -0.5, ]

NC_up_dif_prot <- bitr(rownames(Neural_crest_vs_resr_dif_prot_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
NC_dn_dif_prot <- bitr(rownames(Neural_crest_vs_resr_dif_prot_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

NC_up_dif_rna <- bitr(rownames(res_10d_crest_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
NC_dn_dif_rna <- bitr(rownames(res_10d_crest_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

NC_up_prot_egobp_10d <- enrichGO(
  gene     = NC_up_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_up_rna_egobp_10d <- enrichGO(
  gene     = NC_up_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_prot_egobp_10d <- enrichGO(
  gene     = NC_dn_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_rna_egobp_10d <- enrichGO(
  gene     = NC_dn_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('NC_up_prot_egobp_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_prot_egobp_10d, showCategory = 15)
dev.off()


NC_up_prot_egocc_10d <- enrichGO(
  gene     = NC_up_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_up_rna_egocc_10d <- enrichGO(
  gene     = NC_up_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_prot_egocc_10d <- enrichGO(
  gene     = NC_dn_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

NC_dn_rna_egocc_10d <- enrichGO(
  gene     = NC_dn_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('NC_up_prot_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_prot_egocc_10d, showCategory = 15)
dev.off()
tiff('NC_dn_prot_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_dn_prot_egocc_10d, showCategory = 15)
dev.off()
tiff('NC_up_rna_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(NC_up_rna_egocc_10d, showCategory = 15)
dev.off()


# Fetus
Fetus_vs_rest_dif_prot_up <- Fetus_vs_rest_dif[Fetus_vs_rest_dif$logFC >= 0.5, ]
Fetus_vs_rest_dif_prot_dn <- Fetus_vs_rest_dif[Fetus_vs_rest_dif$logFC <= -0.5, ]

Fet_up_dif_prot <- bitr(rownames(Fetus_vs_rest_dif_prot_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Fet_dn_dif_prot <- bitr(rownames(Fetus_vs_rest_dif_prot_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

Fet_up_dif_rna <- bitr(rownames(res_10d_fetus_spec_up),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Fet_dn_dif_rna <- bitr(rownames(res_10d_fetus_spec_dn),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

Fet_up_prot_egobp_10d <- enrichGO(
  gene     = Fet_up_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_up_rna_egobp_10d <- enrichGO(
  gene     = Fet_up_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_prot_egobp_10d <- enrichGO(
  gene     = Fet_dn_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_rna_egobp_10d <- enrichGO(
  gene     = Fet_dn_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE)


tiff('Fet_dn_prot_egobp_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_dn_prot_egobp_10d, showCategory = 15)
dev.off()


Fet_up_prot_egocc_10d <- enrichGO(
  gene     = Fet_up_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_up_rna_egocc_10d <- enrichGO(
  gene     = Fet_up_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_prot_egocc_10d <- enrichGO(
  gene     = Fet_dn_dif_prot[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

Fet_dn_rna_egocc_10d <- enrichGO(
  gene     = Fet_dn_dif_rna[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  pvalueCutoff = 0.05,
  readable = TRUE)

tiff('Fet_up_prot_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_up_prot_egocc_10d, showCategory = 15)
dev.off()
tiff('Fet_dn_rna_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_dn_rna_egocc_10d, showCategory = 15)
dev.off()
tiff('Fet_dn_prot_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_dn_prot_egocc_10d, showCategory = 15)
dev.off()
tiff('Fet_up_rna_egocc_10d.tiff', units="in", width=6, height=7, res=300, compression = 'lzw')
dotplot(Fet_up_rna_egocc_10d, showCategory = 15)
dev.off()
