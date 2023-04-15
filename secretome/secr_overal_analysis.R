####This it the code to reproduce the statistical data analysis from
####"Human mesenchymal cells retain the specificity of their embryonal origin both before and after osteogenic differentiation in vitro"
####code author: Lobov A.A.
####code includes the comparative analysis of secretomes between cells with different embryological origin

###Opening the data
## Protein expression data (output from FragPipe)
dat <- data.frame(read.table("combined_protein_secr_fragpipe.tsv", sep = '\t', header = TRUE))
dat <- dat[-c(1:10),] # removing contaminants
dat <- dat[dat$Combined.Total.Peptides > 1,] # removing the proteins identified by 1 peptide
rownames(dat) <- make.names(dat[,4], unique=TRUE)

#We will use MaxLFQ intensity for the analysis
dat <- dat[,c(215:264)]

## sample information table
library(readxl)
fact <- read_xlsx("fact_secr.xlsx")
fact <- as.data.frame(fact)

rownames(fact) <- fact$Sample
fact$Cell_type <- as.factor(fact$Cell_type)
fact$Cell_type

fact$Differentiation <- as.factor(fact$Differentiation)
fact$Differentiation

fact$Donor <- as.factor(fact$Donor)
fact$Donor

fact$Cluster <- as.factor(fact$Cluster)
fact$Cluster

colnames(dat) <- rownames(fact)

### Data preparation and multimeric analysis. We analyse control and differentiated cells separately
## Control
dat_contr <- dat[,fact$Differentiation == "control"]
fact_contr <- fact[fact$Differentiation == "control",]
#replace 0 with NA
dat_contr[dat_contr == 0] <- NA

mean(complete.cases(dat_contr)) #we have a lot of NA
dat_contr1 <- dat_contr[which(rowMeans(!is.na(dat_contr[,fact_contr$Cell_type == "ADSCs"])) >= 0.5), ]
dat_contr1 <- dat_contr1[which(rowMeans(!is.na(dat_contr1[,fact_contr$Cell_type == "DPSCs"])) >= 0.5), ]
dat_contr1 <- dat_contr1[which(rowMeans(!is.na(dat_contr1[,fact_contr$Cell_type == "GFs"])) >= 0.6666667), ]
dat_contr1 <- dat_contr1[which(rowMeans(!is.na(dat_contr1[,fact_contr$Cell_type == "OBs"])) >= 0.5), ]
dat_contr1 <- dat_contr1[which(rowMeans(!is.na(dat_contr1[,fact_contr$Cell_type == "PDLSCs"])) >= 0.5), ]
dat_contr1 <- dat_contr1[which(rowMeans(!is.na(dat_contr1[,fact_contr$Cell_type == "WJ-MSCs"])) >= 0.5), ]
mean(complete.cases(dat_contr1))

#knn imputation of missng values
library(impute)
tdat <- t(dat_contr1)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))

library(ggplot2)
library(RColorBrewer)
pal <- brewer.pal(n = 6, name = "Set1")
cols <- pal[fact_contr$Cell_type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact_contr$Cell_type), fill = pal, bty = "n", xpd = T)

#log transformation 
dat_log <- log2(dat_knn+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact_contr$Cell_type), fill = pal, bty = "n", xpd = T)

#Data normalization
library(vsn)
library(RColorBrewer)
meanSdPlot(as.matrix(dat_log))

library(limma)
#quantile normalization
dat_contr <- normalizeQuantiles(dat_log)

boxplot(dat_contr, outline = TRUE, col = cols, names=colnames(dat_contr))
meanSdPlot(as.matrix(dat_contr))

## multivariate analysis
#nMDS
library(vegan)
set.seed(125)
nmds <- metaMDS(t(dat_contr), distance = "euclidean")
nmds
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores2 <- as.data.frame(scores(nmds)$sites)
data.scores$group <- fact_contr$Cell_type
data.scores2$group <- fact_contr$Cluster
head(data.scores)

library(ggplot2)
#tiff('contr_celltype_nMDS.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
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

#tiff('contr_embryo_nMDS.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
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

dda_pca <- pca(t(dat_contr), ncomp = 8)
plot(dda_pca)

contr_pca <- pca(t(dat_contr), ncomp = 3)
#tiff('contr_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(contr_pca, comp = c(1,2), legend = TRUE,
          group = fact_contr$Cell_type,
          ellipse = TRUE,
          title = 'contr_pca comp 1 - 2')
#dev.off()

#tiff('contr_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(contr_pca, comp = c(1,2), legend = TRUE,
          group = fact_contr$Cluster,
          pch = fact_contr$Cell_type,
          ellipse = TRUE,
          title = 'contr_pca comp 1 - 2')
#dev.off()


#PLS-DA 
ordination.optimum.splsda <- splsda(t(dat_contr), fact_contr$Cluster, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('secr_contr_splsda1_comp.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
#dev.off()

#tiff('secr_contr_splsda1_graph.tiff', units="in", width=7, height=6, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, pch = fact_contr$Cell_type, ind.names = F, ellipse = T, abline = TRUE,  title = "PLS-DA ordination", legend=TRUE)
#dev.off()

layout(1,1)


## Differentiation
dat_dif <- dat[,fact$Differentiation == "dif"]
fact_dif <- fact[fact$Differentiation == "dif",]
dat_dif[dat_dif == 0] <- NA
dat_dif <- dat_dif[,-9] # removing sample GF1 which referred to has technical problems during LC-MS
fact_dif <- fact_dif[-9,] # removing sample GF1 which referred to has technical problems during LC-MS

mean(complete.cases(dat_dif))
dat_dif1 <- dat_dif[which(rowMeans(!is.na(dat_dif[,fact_dif$Cell_type == "ADSCs"])) >= 0.5), ]
dat_dif1 <- dat_dif1[which(rowMeans(!is.na(dat_dif1[,fact_dif$Cell_type == "DPSCs"])) >= 0.5), ]
dat_dif1 <- dat_dif1[which(rowMeans(!is.na(dat_dif1[,fact_dif$Cell_type == "GFs"])) >= 0.5), ]
dat_dif1 <- dat_dif1[which(rowMeans(!is.na(dat_dif1[,fact_dif$Cell_type == "OBs"])) >= 0.5), ]
dat_dif1 <- dat_dif1[which(rowMeans(!is.na(dat_dif1[,fact_dif$Cell_type == "PDLSCs"])) >= 0.5), ]
dat_dif1 <- dat_dif1[which(rowMeans(!is.na(dat_dif1[,fact_dif$Cell_type == "WJ-MSCs"])) >= 0.5), ]
mean(complete.cases(dat_dif1))

#knn imputation of missng values
library(impute)
tdat <- t(dat_dif1)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))

cols <- pal[fact_dif$Cell_type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact_dif$Cell_type), fill = pal, bty = "n", xpd = T)

#log transformation 
dat_log <- log2(dat_knn+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact_dif$Cell_type), fill = pal, bty = "n", xpd = T)

##Data normalization
meanSdPlot(as.matrix(dat_log))

#quantile normalization
dat_dif <- normalizeQuantiles(dat_log)

boxplot(dat_dif, outline = TRUE, col = cols, names=colnames(dat_dif))
meanSdPlot(as.matrix(dat_dif))

## multivariate analysis
#nMDS
library(vegan)
set.seed(125)
nmds <- metaMDS(t(dat_dif), distance = "euclidean")
nmds
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$group <- fact_dif$Cell_type
data.scores2 <- as.data.frame(scores(nmds)$sites)
data.scores2$group <- fact_dif$Cluster
head(data.scores)

library(ggplot2)
#tiff('secr_dif_celltype_nMDS.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
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

#tiff('secr_dif_embryo_nMDS.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
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

dda_pca <- pca(t(dat_dif), ncomp = 8)
plot(dda_pca)

dda_pca <- pca(t(dat_dif), ncomp = 3)
#tiff('dif_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(dda_pca, comp = c(1,2), legend = TRUE,
          group = fact_dif$Cell_type,
          ellipse = TRUE,
          title = 'DDA PCA comp 1 - 2')
#dev.off()

#tiff('dif_pca.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
plotIndiv(dda_pca, comp = c(1,2), legend = TRUE,
          group = fact_dif$Cluster,
          pch = fact_dif$Cell_type,
          ellipse = TRUE,
          title = 'DDA PCA comp 1 - 2')
#dev.off()

#PLS-DA 
ordination.optimum.splsda <- splsda(t(dat_dif), fact_dif$Cluster, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('secr_dif_splsda1_comp.tiff', units="in", width=10, height=6, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", difib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
#dev.off()

#tiff('secr_dif_splsda1_graph.tiff', units="in", width=10, height=6, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, pch = fact_dif$Cell_type, ind.names = F, ellipse = T, abline = TRUE, 
          title = "PLS-DA ordination", legend=TRUE)
#dev.off()
layout(1,1)

### Qualititative comparison of control and differentiated secretomes
library(VennDiagram)
venn.diagram(
  x = list(rownames(dat_contr), rownames(dat_dif)),
  category.names = c("Contr" , "Dif"),
  filename = '#unique_datasets.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)

unique_secr <- VennDiagram::get.venn.partitions(list(contr=rownames(dat_contr), dif=rownames(dat_dif)))
Osteo_unique_secr <- unique_secr[2,4]
Osteo_unique_secr <- as.data.frame(Osteo_unique_secr)
Osteo_unique_secr <- Osteo_unique_secr$X2


library("org.Hs.eg.db")
library(clusterProfiler)
library(KEGG.db)
Osteo_unique_secr <- clusterProfiler::bitr(Osteo_unique_secr,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

kegg_osteo <- clusterProfiler::enrichKEGG(
  gene     = Osteo_unique_secr[[2]],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  use_internal_data = F)

#tiff('kegg_secr_ost_specific.tiff', units="in", width=7, height=5, res=300, compression = 'lzw')
clusterProfiler::dotplot(kegg_osteo, showCategory = 20)
#dev.off()

common_secr <- unique_secr[1,4]
common_secr <- as.data.frame(common_secr)
common_secr <- common_secr$X1
common_secr[89:91] <- c("CALD", "X", "NME")
common_secr <- clusterProfiler::bitr(common_secr[1:88],fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

kegg_common <- clusterProfiler::enrichKEGG(
  gene     = common_secr[[2]],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)


#tiff('secr_kegg_secr_common.tiff', units="in", width=8.5, height=6, res=300, compression = 'lzw')
clusterProfiler::dotplot(kegg_common, showCategory = 20)
#dev.off()

common_secr1 <- common_secr[,1]
Osteo_unique_secr1 <- Osteo_unique_secr[,1]

nnn <- max(length(common_secr1), length(Osteo_unique_secr1))
length(common_secr1) <- nnn                      
length(Osteo_unique_secr1) <- nnn

secr_over_unique_list <- cbind(common_secr1, Osteo_unique_secr1)
head(secr_over_unique_list)
#write.csv(secr_over_unique_list, "secr_over_unique_list.csv")

###Dif expression with limma
##control
design = model.matrix(~0+fact_contr$Cluster)
colnames(design) = c('Fetus', 'Mesenchyme','Neural_crest')
fit_contr <- lmFit(dat_contr, design)
contrasts_groups = c('Mesenchyme-Fetus','Mesenchyme-Neural_crest','Neural_crest-Fetus')
contrast.matrix <- makeContrasts(contrasts = contrasts_groups,levels=design)

fit2_contr <- contrasts.fit(fit_contr, contrast.matrix)
fit2_contr <- eBayes(fit2_contr)

results_contr <- decideTests(fit2_contr)
vennDiagram(results_contr)

Anova_contr_secr <- topTable(fit2_contr, number=95)
#write.csv(Anova_contr_secr, "Anova_contr_secr.csv")

Mesenchyme_Fetus_contr <- topTable(fit2_contr, coef=1, adjust="BH", number = 95)
#write.csv(Mesenchyme_Fetus_contr, "Mesenchyme_Fetus_contr.csv")
Mesenchyme_Neural_crest <- topTable(fit2_contr, coef=2, adjust="BH", number = 95)
#write.csv(Mesenchyme_Neural_crest, "Mesenchyme_Neural_crest_contr.csv")
Neural_crest_Fetus <- topTable(fit2_contr, coef=3, adjust="BH", number = 95)
#write.csv(Neural_crest_Fetus, "Neural_crest_Fetus_contr.csv")

Anova_contr_secr_rel <- topTable(fit2_contr, number=25)

proteins_dif_secr_cont <- rownames(Anova_contr_secr_rel)
proteins_dif_secr_cont <- proteins_dif_secr_cont[-c(3, 5, 8:12, 15, 21)] # removing keratins and collagens

proteins_dif_secr_cont <- dat_contr[proteins_dif_secr_cont,]
library(pheatmap)
fact_contr1 <- fact_contr[,c(2,5)]
#tiff('heatmap_contr_full.tiff', units="in", width=8, height=6, res=300, compression = 'lzw')
pheatmap::pheatmap(proteins_dif_secr_cont,scale="row", annotation_col = fact_contr1,
                   color=colorRampPalette(c("navy", "white", "red"))(50))
#dev.off()

##Differentiatied
design = model.matrix(~0+fact_dif$Cluster)
colnames(design) = c('Fetus', 'Mesenchyme','Neural_crest')

fit_dif <- lmFit(dat_dif, design)
contrasts_groups = c('Mesenchyme-Fetus','Mesenchyme-Neural_crest','Neural_crest-Fetus')
contrast.matrix <- makeContrasts(contrasts = contrasts_groups,levels=design)

fit2_dif <- contrasts.fit(fit_dif, contrast.matrix)
fit2_dif <- eBayes(fit2_dif)

Mesenchyme_Fetus_contr <- topTable(fit2_dif, coef=1, adjust="BH", number = 146)
#write.csv(Mesenchyme_Fetus_contr, "Mesenchyme_Fetus_dif.csv")
Mesenchyme_Neural_crest <- topTable(fit2_dif, coef=2, adjust="BH", number = 146)
#write.csv(Mesenchyme_Neural_crest, "Mesenchyme_Neural_crest_dif.csv")
Neural_crest_Fetus <- topTable(fit2_dif, coef=3, adjust="BH", number = 146)
#write.csv(Neural_crest_Fetus, "Neural_crest_Fetus_dif.csv")

Anova_dif_secr_rel <- topTable(fit2_dif, number=146)
#write.csv(Anova_dif_secr_rel, "Anova_dif_secr_rel.csv")

Anova_dif_secr_rel <- topTable(fit2_dif, number=20)

proteins_dif_secr_dif <- rownames(Anova_dif_secr_rel)

proteins_dif_secr_dif <- dat_dif[proteins_dif_secr_dif,]
fact_dif1 <- fact_dif[,c(2,5)]
tiff('heatmap_dif_secr.tiff', units="in", width=8, height=6, res=300, compression = 'lzw')
pheatmap::pheatmap(proteins_dif_secr_dif,scale="row", annotation_col = fact_dif1,
                   color=colorRampPalette(c("navy", "white", "red"))(50))
dev.off()


### Embryonic origin-specific proteins
## Control
dat_contr <- dat[,fact$Differentiation == "control"]
dat_contr[dat_contr == 0] <- NA

mean(complete.cases(dat_contr))
contr_ADSCs <- dat_contr[which(rowMeans(!is.na(dat_contr[fact_contr$Cell_type == "ADSCs"])) >= 0.5), ]
contr_DPSCs <- dat_contr[which(rowMeans(!is.na(dat_contr[fact_contr$Cell_type == "DPSCs"])) >= 0.5), ]
contr_GFs <- dat_contr[which(rowMeans(!is.na(dat_contr[fact_contr$Cell_type == "GFs"])) >= 0.6666667), ]
contr_ost <- dat_contr[which(rowMeans(!is.na(dat_contr[fact_contr$Cell_type == "OBs"])) >= 0.5), ]
contr_PDLSCs <- dat_contr[which(rowMeans(!is.na(dat_contr[fact_contr$Cell_type == "PDLSCs"])) >= 0.5), ]
contr_Fet <- dat_contr[which(rowMeans(!is.na(dat_contr[fact_contr$Cell_type == "WJ-MSCs"])) >= 0.5), ]

contr_Fet <- rownames(contr_Fet)

unique_crest <- VennDiagram::get.venn.partitions(list(dpscs=rownames(contr_DPSCs), gf=rownames(contr_GFs), pdlscs=rownames(contr_PDLSCs)))
unique_mes <- VennDiagram::get.venn.partitions(list(adscs=rownames(contr_ADSCs), ost=rownames(contr_ost)))

contr_crest <- unique_crest[1,5]
contr_crest <- as.data.frame(contr_crest)
contr_crest <- contr_crest$X1

contr_mes <- unique_mes[1,4]
contr_mes <- as.data.frame(contr_mes)
contr_mes <- contr_mes$X1

contr_mes1 <- contr_mes
contr_crest1 <- contr_crest
contr_Fet1 <- contr_Fet

nnn <- max(length(contr_mes1), length(contr_crest1), length(contr_Fet1))
length(contr_mes1) <- nnn                      
length(contr_crest1) <- nnn
length(contr_Fet1) <- nnn

control_unique_table <- cbind(contr_mes1, contr_crest1, contr_Fet1)
#write.csv(control_unique_table, "control_unique_table.csv")


VennDiagram::venn.diagram(
  x = list(contr_crest, contr_mes, contr_Fet),
  category.names = c("Crest" , "Mes", "Fet"),
  filename = '#contr_cluster_secr.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)


venn_contr_cluster <- VennDiagram::get.venn.partitions(list(crest=contr_crest, mes=contr_mes, fet=contr_Fet))

contr_Fet_spec <- venn_contr_cluster[4,5]
contr_Fet_spec <- as.data.frame(contr_Fet_spec)
contr_Fet_spec <- contr_Fet_spec$X4

contr_crest_spec <- venn_contr_cluster[7,5]
contr_crest_spec <- as.data.frame(contr_crest_spec)
contr_crest_spec <- contr_crest_spec$X7
contr_crest2_spec <- venn_contr_cluster[3,5]
contr_crest2_spec <- as.data.frame(contr_crest2_spec)
contr_crest2_spec <- contr_crest2_spec$X3
contr_crest_spec_full <- c(contr_crest_spec, contr_crest2_spec)

contr_crest_spec_full <- gsub("\\..*","",contr_crest_spec_full)
contr_Fet_spec <- gsub("\\..*","",contr_Fet_spec)

#enrichment

contr_crest_spec_full <- clusterProfiler::bitr(contr_crest_spec_full,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

kegg_crest <- clusterProfiler::enrichKEGG(
  gene     = contr_crest_spec_full[[2]],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)

library(clusterProfiler)
#tiff('kegg_secr_crest.tiff', units="in", width=8.5, height=6, res=300, compression = 'lzw')
clusterProfiler::dotplot(kegg_crest, showCategory = 20)
#dev.off()

egobp_crest <- clusterProfiler::enrichGO(
  gene     = contr_crest_spec_full[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

#tiff('egobp_secr_crest.tiff', units="in", width=8, height=7.5, res=300, compression = 'lzw')
clusterProfiler::dotplot(egobp_crest, showCategory = 20)
#dev.off()

contr_Fet_spec <- clusterProfiler::bitr(contr_Fet_spec,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

kegg_fet <- clusterProfiler::enrichKEGG(
  gene     = contr_Fet_spec[[2]],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)


#tiff('kegg_secr_fet.tiff', units="in", width=8.5, height=6, res=300, compression = 'lzw')
clusterProfiler::dotplot(kegg_fet, showCategory = 20)
#dev.off()

egobp_fet <- clusterProfiler::enrichGO(
  gene     = contr_Fet_spec[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

#tiff('egobp_secr_fet.tiff', units="in", width=8, height=11, res=300, compression = 'lzw')
clusterProfiler::dotplot(egobp_fet, showCategory = 20)
#dev.off()


## Differentiated
dat_dif <- dat[,fact$Differentiation == "dif"]
fact_dif <- fact[fact$Differentiation == "dif",]
dat_dif[dat_dif == 0] <- NA
dat_dif <- dat_dif[,-9]
fact_dif <- fact_dif[-9,]

mean(complete.cases(dat_dif))
dif_ADSCs <- dat_dif[which(rowMeans(!is.na(dat_dif[fact_dif$Cell_type == "ADSCs"])) >= 0.5), ]
dif_DPSCs <- dat_dif[which(rowMeans(!is.na(dat_dif[fact_dif$Cell_type == "DPSCs"])) >= 0.5), ]
dif_GFs <- dat_dif[which(rowMeans(!is.na(dat_dif[fact_dif$Cell_type == "GFs"])) >= 0.6666667), ]
dif_ost <- dat_dif[which(rowMeans(!is.na(dat_dif[fact_dif$Cell_type == "OBs"])) >= 0.5), ]
dif_PDLSCs <- dat_dif[which(rowMeans(!is.na(dat_dif[fact_dif$Cell_type == "PDLSCs"])) >= 0.5), ]
dif_Fet <- dat_dif[which(rowMeans(!is.na(dat_dif[fact_dif$Cell_type == "WJ-MSCs"])) >= 0.5), ]

dif_Fet <- rownames(dif_Fet)

unique_crest <- VennDiagram::get.venn.partitions(list(dpscs=rownames(dif_DPSCs), gf=rownames(dif_GFs), pdlscs=rownames(dif_PDLSCs)))
unique_mes <- VennDiagram::get.venn.partitions(list(adscs=rownames(dif_ADSCs), ost=rownames(dif_ost)))

dif_crest <- unique_crest[1,5]
dif_crest <- as.data.frame(dif_crest)
dif_crest <- dif_crest$X1

dif_mes <- unique_mes[1,4]
dif_mes <- as.data.frame(dif_mes)
dif_mes <- dif_mes$X1


dif_mes1 <- dif_mes
dif_crest1 <- dif_crest
dif_Fet1 <- dif_Fet

nnn <- max(length(dif_mes1), length(dif_crest1), length(dif_Fet1))
length(dif_mes1) <- nnn                      
length(dif_crest1) <- nnn
length(dif_Fet1) <- nnn

dif_unique_list <- cbind(dif_mes1, dif_crest1, dif_Fet1)
#write.csv(dif_unique_list, "dif_unique_list.csv")

VennDiagram::venn.diagram(
  x = list(dif_crest, dif_mes, dif_Fet),
  category.names = c("Crest" , "Mes", "Fet"),
  filename = '#dif_cluster_secr.png',
  resolution = 300,
  print.mode=c("raw","percent"),
  output=F
)


venn_dif_cluster <- VennDiagram::get.venn.partitions(list(crest=dif_crest, mes=dif_mes, fet=dif_Fet))

dif_Fet_spec <- venn_dif_cluster[4,5]
dif_Fet_spec <- as.data.frame(dif_Fet_spec)
dif_Fet_spec <- dif_Fet_spec$X4

dif_crest_spec <- venn_dif_cluster[7,5]
dif_crest_spec <- as.data.frame(dif_crest_spec)
dif_crest_spec <- dif_crest_spec$X7
dif_crest2_spec <- venn_dif_cluster[3,5]
dif_crest2_spec <- as.data.frame(dif_crest2_spec)
dif_crest2_spec <- dif_crest2_spec$X3
dif_crest_spec_full <- c(dif_crest_spec, dif_crest2_spec)

dif_crest_spec_full <- gsub("\\..*","",dif_crest_spec_full)
dif_Fet_spec <- gsub("\\..*","",dif_Fet_spec)

#enrichment
dif_crest_spec_full <- clusterProfiler::bitr(dif_crest_spec_full,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

kegg_crest <- clusterProfiler::enrichKEGG(
  gene     = dif_crest_spec_full[[2]],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)


#tiff('kegg_secr_crest.tiff', units="in", width=8.5, height=6, res=300, compression = 'lzw')
clusterProfiler::dotplot(kegg_crest, showCategory = 20)
#dev.off()

egobp_crest <- clusterProfiler::enrichGO(
  gene     = dif_crest_spec_full[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

#tiff('egobp_secr_crest.tiff', units="in", width=8, height=7.5, res=300, compression = 'lzw')
clusterProfiler::dotplot(egobp_crest, showCategory = 20)
#dev.off()

dif_Fet_spec <- clusterProfiler::bitr(dif_Fet_spec,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

kegg_fet <- clusterProfiler::enrichKEGG(
  gene     = dif_Fet_spec[[2]],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)


#tiff('kegg_secr_fet.tiff', units="in", width=8.5, height=6, res=300, compression = 'lzw')
clusterProfiler::dotplot(kegg_fet, showCategory = 20)
#dev.off()

egobp_fet <- clusterProfiler::enrichGO(
  gene     = dif_Fet_spec[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

#tiff('egobp_secr_fet.tiff', units="in", width=8, height=11, res=300, compression = 'lzw')
clusterProfiler::dotplot(egobp_fet, showCategory = 20)
#dev.off()

   