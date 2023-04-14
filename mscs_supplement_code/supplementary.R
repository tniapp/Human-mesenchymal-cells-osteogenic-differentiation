library(DEP)
library(readxl)
library(dplyr)
library(EnhancedVolcano)

### Proteome data analysis
data <- read.csv("./prot_input_files/ads_proteom_frag.csv", sep = ";")
data[data == 0] = NA

#filtering
data$na.count <- apply(data, 1, function(x) sum(is.na(x)))
data$na.count
nrow(data)

data <- subset(data, data$na.count < 3)  # < 3 for MSCs except GF (for GF < 2)
nrow(data)

#checking for duplicated proteins
data$Protein.ID %>% duplicated() %>% any()

#make unique names using the annotation
data_unique <- make_unique(data, "Protein.ID", "Gene", delim = ";")

#generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("MaxLFQ.Intensity", colnames(data_unique))
experimental_design <- read_xlsx("./prot_input_files/groups_proteom_frag.xlsx", sheet = 2)  # select list number in accordance with the MSCs type
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

#plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

#checking for missing values
data_filt <- filter_missval(data_se, thr = 1)
plot_frequency(data_filt)

#protein number
plot_numbers(data_filt)
dev.off()

#normalize the data (if necessary)
data_norm <- normalize_vsn(data_filt)

#visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)  # no differences

#plot a heatmap of proteins with missing values
plot_missval(data_filt)

#impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "knn")

#plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)


## Differential enrichment analysis  based on linear models and empherical Bayes statistics
#test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "ADSCs_c", design_formula = formula(~ 0 + condition + replicate))  # specify the MSCs type for the "control" option

#denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05)

#save dep results as a file
dep_results <- get_results(dep)
colnames(dep_results)
dep_results %>% filter(significant) %>% nrow()

write.csv(x = dep_results, file = "ADSCs_prot_dep.csv")  # change file name 

## PCA
#tiff(file="./adscs_pca_prot.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4, label = F) + scale_color_brewer(palette = "Set2")
dev.off()

#remove batch effect
library(SummarizedExperiment)
assay(dep) <- limma::removeBatchEffect(assay(dep),
                                       batch = colData(dep)[,'replicate'])


## PCA
#tiff(file="./adscs_pca_prot_no_batch.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4, label = F) + scale_color_brewer(palette = "Set2")
dev.off()

## Pearson correlation matrix
#tiff(file="./adscs_corr_plot_prot.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 4,
#     res=300)
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Blues")
dev.off()


## Heatmap of all significant proteins
#tiff(file="./adscs_heatmap.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 5,
#     res=300)
plot_heatmap(dep, type = "centered",
             kmeans = TRUE,
             k = 2,
             col_limit = 4, 
             show_row_names = TRUE,
             indicate = c("condition"),
             )
dev.off()

## Volcano plot for the contrast
#tiff(file="./adscs_volcano.tiff",  # change file name
#     units = "in",
#     width = 7,
#     height = 5,
#     res=300)
data <- dep_results
data_p <- data[, c(2, 7, 3, 4)]
names(data_p)[1] <- "GeneName"
names(data_p)[2] <- "logFC"
names(data_p)[3] <- "P.Value"
names(data_p)[4] <- "adj.P.Val"

EnhancedVolcano(data_p,
                lab = data_p$GeneName,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                title ="ADSCs_c vs ADSCs_d",  #change title
                labSize = 3.0,
                boxedLabels = F,
                col=c('black', 
                      '#CBD5E8', 
                      '#B3E2CD', 
                      '#FDCDAC'),
                colAlpha = 1)

dev.off()


### Secretome data analysis
data <- read.csv("./secr_input_files/ads_secr_fragpipe.csv", sep = ";")
data[data == 0] = NA

#filtering
data$na.count <- apply(data, 1, function(x) sum(is.na(x)))
data$na.count
nrow(data)

data <- subset(data, data$na.count < 3)  # < 3 for MSCs except GF (for GF < 2)
nrow(data)

#checking for duplicated proteins
data$Protein.ID %>% duplicated() %>% any()

#make unique names using the annotation
data_unique <- make_unique(data, "Protein.ID", "Gene", delim = ";")

#generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("MaxLFQ.Intensity", colnames(data_unique))
experimental_design <- read_xlsx("./secr_input_files/groups_secr.xlsx", sheet = 2)  # select list number in accordance with the MSCs type
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

#plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

#checking for missing values
data_filt <- filter_missval(data_se, thr = 1)
plot_frequency(data_filt)

#protein number
plot_numbers(data_filt)
dev.off()

#normalize the data (if necessary)
data_norm <- normalize_vsn(data_filt)

#visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm) 

#plot a heatmap of proteins with missing values
plot_missval(data_filt)

#impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "knn")

#plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

## Differential enrichment analysis  based on linear models and empherical Bayes statistics
#test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Ad_c", design_formula = formula(~ 0 + condition + replicate))  # specify the MSCs type for the "control" option

#denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05)

#save dep results as a file
dep_results <- get_results(dep)
colnames(dep_results)
dep_results %>% filter(significant) %>% nrow()

write.csv(x = dep_results, file = "ADSCs_secr_dep.csv")  # change file name 

## PCA
#tiff(file="./adscs_pca_prot.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
plot_pca(dep, x = 1, y = 2, n = 200, point_size = 4, label = F) + scale_color_brewer(palette = "Set2")
dev.off()

#remove batch effect
library(SummarizedExperiment)
assay(dep) <- limma::removeBatchEffect(assay(dep),
                                       batch = colData(dep)[,'replicate'])
## PCA
#tiff(file="./adscs_pca_prot_no_batch.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
plot_pca(dep, x = 1, y = 2, n = 200, point_size = 4, label = F) + scale_color_brewer(palette = "Set2")
dev.off()

## Pearson correlation matrix
#tiff(file="./adscs_corr_plot_prot.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 4,
#     res=300)
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Blues")
dev.off()

## Heatmap of all significant proteins
#tiff(file="./adscs_heatmap.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 5,
#     res=300)
plot_heatmap(dep, type = "centered",
             kmeans = TRUE,
             k = 1,
             col_limit = 4, 
             show_row_names = TRUE,
             indicate = c("condition"),
)
dev.off()

## Volcano plot for the contrast
#tiff(file="./adscs_volcano.tiff",  # change file name
#     units = "in",
#     width = 7,
#     height = 5,
#     res=300)
data <- dep_results
data_p <- data[, c(2, 7, 3, 4)]
names(data_p)[1] <- "GeneName"
names(data_p)[2] <- "logFC"
names(data_p)[3] <- "P.Value"
names(data_p)[4] <- "adj.P.Val"

EnhancedVolcano(data_p,
                lab = data_p$GeneName,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                title ="ADSCs_c vs ADSCs_d",  #change title
                labSize = 3.0,
                boxedLabels = F,
                col=c('black', 
                             '#CBD5E8', 
                             '#B3E2CD', 
                             '#FDCDAC'),
                             colAlpha = 1)

dev.off()


### Transcriptome data analysis
library(org.Hs.eg.db)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(edgeR)

count_matrix <- as.matrix(read.csv("./trans_input_files/fat_transcriptomes_total.csv", sep = ",", row.names = 1))  # path to the MSCs counts
head(count_matrix, 2)

#transcriptome metadata
snames <- colnames(count_matrix)
status <- substr(snames, 15, nchar(snames) - 2)  # substitute for correct column names (48h.undif, 48h.dif, and 10d.dif)
status <- as.factor(status)

coldata <- data.frame(
        sample = c(colnames(count_matrix)),
        replicates = substr(snames, 1, nchar(snames) - 2),
        status = status,
        batch = substr(snames, 1, 13),  # substitute manually for each MSCs
        row.names = "sample")

#check the order of samples in metadata and count_matrix
all(rownames(coldata) %in% colnames(count_matrix))


##Start analysis
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~batch + status)

dds <- collapseReplicates(dds, dds$replicates)

#set the reference condition (undif)
dds$status <- relevel(dds$status, ref = "48h.undif")

#median of ratios normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

#filtration of low-counted genes
dds <- dds[rowSums(counts(dds)) >= 10,]

#add gene symbols to gene properties
genes <- rownames(dds)
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

#change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds) 

## DESeq
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons 

#separate DESeq results 
res_dif_10d <- results(dds, contrast=c("status","10d.dif","48h.undif"))
res_dif_48h <- results(dds, contrast=c("status","48h.dif","48h.undif"))

#export DESeq results in the file
write.csv(as.data.frame(res_dif_10d[order(res_dif_10d$padj),] ), file="./fat_merged_replicates_dif10d_vs_undif.csv")  # change file name
write.csv(as.data.frame(res_dif_48h[order(res_dif_48h$padj),] ), file="./fat_merged_replicates_dif48h_vs_undif.csv")  # change file name

#summary of dif. expressed genes
summary(res_dif_10d)
summary(res_dif_48h)

#plot mean of normalized count
head(res_dif_10d[order(res_dif_10d$pvalue), ])
DESeq2::plotMA(res_dif_10d)
DESeq2::plotMA(res_dif_48h)
dev.off()

#plot counts
plotCounts(dds, gene=which.min(res_dif_10d$padj), intgroup="status")
dev.off()

## PCA
rlt <- rlog(dds)  #rlog Transformation

# remove batch effect
assay(rlt) <- limma::removeBatchEffect(assay(rlt),
                                       batch = colData(dds)[,'batch'])

library(ggplot2)
#tiff(file="./trans_adscs_pca.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
pcaData <- plotPCA(rlt, intgroup=c("status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
pcaData[,5] <- substr(names, 1, 13)  # change manually for each MSCs type
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")
ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        theme_bw() +
        scale_color_brewer(palette = "Set2")
dev.off()

## Heatmap by https://github.com/ACSoupir/Bioinformatics_YouTube/blob/master/Visualizing%20Counts/README.md
detectGroups <- function (x){  # x are col names
        tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
        #tem = gsub("_Rep|_rep|_REP","",tem)
        tem <- gsub("_$","",tem); # remove "_" from end
        tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
        tem <- gsub("_rep$","",tem); # remove "_rep" from end
        tem <- gsub("_REP$","",tem)  # remove "_REP" from end
        return( tem )
}

dist2 <- function(x, ...)   # distance function = 1-SCC (Spearman's correlation coefficient)
        as.dist(1-cor(t(x), method="spearman"))

hclust2 <- function(x, method="average", ...){  # average linkage in hierarchical clustering
        hclust(x, method=method, ...)
}

n = 30 # number of top genes by standard deviation

x = assay(rlt)
if(n>dim(x)[1]) n = dim(x)[1] # max as data

x = x[order(apply(x,1,sd),decreasing=TRUE),]  # sort genes by standard deviation

x = x[1:n,]   # only keep the n genes

# this will cutoff very large values, which could skew the color 
x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
cutoff = median(unlist(x)) + 4*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 4*sd (unlist(x)) 
x[x< cutoff] <- cutoff

groups = rlt@colData@listData[["status"]]
groups.colors = rainbow(length(unique(groups) ) )


lmat = rbind(c(5,4),c(0,1),c(3,2))
lwid = c(1.5,4)
lhei = c(1,.2,4)

#tiff(file="./trans_adscs_heatmap.tiff",  # change file name
#     units = "in",
#     width = 8,
#     height = 6,
#     res=300)
heatmap.2(x, distfun = dist2,hclustfun=hclust2,
          col=brewer.pal(11,"RdBu"), density.info="none", 
          trace="none", scale="none", keysize=.5
          ,key=T, symkey=F
          ,ColSideColors=groups.colors[ as.factor(groups)]
          ,margins=c(8,12) + 0.5
          ,cexRow=1
          ,srtCol=45
          ,cexCol=1.  # size of font for sample names
          ,lmat = lmat, lwid = lwid, lhei = lhei
)
dev.off()

## Volcano plot
library(EnhancedVolcano)

#tiff(file="./trans_adscs_48hdif_vs_undif_volcano.tiff", # change file name
#     units = "in",
#     width = 8,
#     height = 6,
#     res=300)
EnhancedVolcano(res_dif_48h,
                lab = rownames(res_dif_48h),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title ="MSC ADSCs (dif.48h vs undif.48h)", # change graph title
                labSize = 2,
                boxedLabels = F,
                col=c('black',
                             '#CBD5E8', 
                             '#B3E2CD', 
                             '#FDCDAC'),
                             colAlpha = 1)

dev.off()


### Total analysis of proteome, secretome and transcriptome data

## Enrichment analysis
# combine all upregulated genes from trans, prot and secr
library(org.Hs.eg.db)

#differentially expressed genes from transcriptome
trans_deg <- read.csv("./fat_merged_replicates_dif48h_vs_undif.csv")  # change file name
trans.up <- trans_deg[trans_deg$log2FoldChange > 1 & trans_deg$padj < 0.05 & !is.na(trans_deg$padj), 1] 
trans.dn <- trans_deg[trans_deg$log2FoldChange < -1 & trans_deg$padj < 0.05 & !is.na(trans_deg$padj), 1]
trans.up.entrez <- clusterProfiler::bitr(trans.up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
trans.dn.entrez <- clusterProfiler::bitr(trans.dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

trans_deg10 <- read.csv("./fat_merged_replicates_dif10d_vs_undif.csv")  # change file name
trans.up10 <- trans_deg10[trans_deg10$log2FoldChange > 1 & trans_deg10$padj < 0.05 & !is.na(trans_deg10$padj), 1] 
trans.dn10 <- trans_deg10[trans_deg10$log2FoldChange < -1 & trans_deg10$padj < 0.05 & !is.na(trans_deg10$padj), 1]
trans.up.entrez10 <- clusterProfiler::bitr(trans.up10,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
trans.dn.entrez10 <- clusterProfiler::bitr(trans.dn10,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#from proteome
prot_dep <- read.csv("./ADSCs_prot_dep.csv", row.names = 1)  # change file name
prot.up <- prot_dep[prot_dep[,7] > 1 & prot_dep[,4] < 0.05 & !is.na(prot_dep[,4]), 2] 
prot.dn <- prot_dep[prot_dep[,7] < -1 & prot_dep[,4] < 0.05 & !is.na(prot_dep[,4]), 2]
prot.up.entrez <- clusterProfiler::bitr(prot.up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
prot.dn.entrez <- clusterProfiler::bitr(prot.dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#from secretome
secr_dep <- read.csv("./ADSCs_secr_dep.csv", row.names = 1)  # change file name
secr.up <- secr_dep[secr_dep[,7] > 1 & secr_dep[,4] < 0.05 & !is.na(secr_dep[,4]), 2] 
secr.dn <- secr_dep[secr_dep[,7] < -1 & secr_dep[,4] < 0.05 & !is.na(secr_dep[,4]), 2]
secr.up.entrez <- clusterProfiler::bitr(secr.up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
secr.dn.entrez <- clusterProfiler::bitr(secr.dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

clusters.up <- list(trans_up_48h = trans.up.entrez$ENTREZID,
                    trans_up_10d = trans.up.entrez10$ENTREZID,
                    prot_up = prot.up.entrez$ENTREZID,
                    secr_up = secr.up.entrez$ENTREZID)


clusters.dn <- list(trans_down_48h = trans.dn.entrez$ENTREZID,
                    trans_down_10d = trans.dn.entrez10$ENTREZID,
                    prot_down = prot.dn.entrez$ENTREZID,
                    secr_down = secr.dn.entrez$ENTREZID)

library(ReactomePA)
library(clusterProfiler)

#enrichGO
ck1_1 <- compareCluster(geneCluster = clusters.up, fun = enrichGO, OrgDb=org.Hs.eg.db, ont = "MF")
ck1_2 <- compareCluster(geneCluster = clusters.dn, fun = enrichGO, OrgDb=org.Hs.eg.db, ont = "MF")

#enrichPathways
ck2_1 <- compareCluster(geneCluster = clusters.up, fun = enrichPathway)
ck2_2 <- compareCluster(geneCluster = clusters.dn, fun = enrichPathway)

#for enrichKEGG (if errors)
#library(R.utils)
#R.utils::setOption("clusterProfiler.download.method","auto")

library(ggplot2)
p1 <- dotplot(ck1_1, 
              title="MSC ADSCs upregulated, enrichGO (ONT = MF)", #change title
              includeAll=TRUE)+ scale_colour_distiller(palette="Set2")
p2 <- dotplot(ck1_2, 
              title="MSC ADSCs downregulated, enrichGO (ONT = MF)", #change title
              includeAll=TRUE)+ scale_colour_distiller(palette="Set2")
p3 <- dotplot(ck2_1, 
              title="MSC ADSCs upregulated, enrichPathway", #change title
              includeAll=TRUE)+ scale_colour_distiller(palette="Set2")
p4 <- dotplot(ck2_2, 
              title="MSC ADSCs downregulated, enrichPathway", #change title
              includeAll=TRUE)+ scale_colour_distiller(palette="Set2")

library(plotly)
#tiff(file="./adscs_ehrichGO_enrichPathway.tiff",  # change file name
#     units = "in",
#     width = 15,
#     height = 20,
#     res=300)
gridExtra::grid.arrange(p1,p2,p3,p4, ncol = 2)
dev.off()

#tiff(file="./adscs_ehrichGO_enrichPathway_up.tiff",  # change file name
#     units = "in",
#     width = 15,
#     height = 10,
#     res=300)
gridExtra::grid.arrange(p1,p3, ncol = 2)
dev.off()

#tiff(file="./adscs_ehrichGO_enrichPathway_dn.tiff", # change file name
#     units = "in",
#     width = 15,
#     height = 10,
#     res=300)
gridExtra::grid.arrange(p2,p4, ncol = 2)
dev.off()


## Check for overlaps between differentially expressed genes (prot, secr, trans)
library(VennDiagram)
library(RColorBrewer)

#extract genes which encode differentially expressed proteins (secretome)
secretom <- read.csv("./ADSCs_secr_dep.csv", row.names = 1)  # change file name
secretom_significant <- subset(secretom, secretom[,4] < 0.05)
secretom_genes_all <- secretom_significant$ID
secr_upreg <- subset(secretom_significant, secretom_significant[,7] > 1)
secr_upreg_genes <- secr_upreg$ID
secr_downreg <- subset(secretom_significant, secretom_significant[,7] < -1)
secr_downreg_genes <- secr_downreg$ID

#extract genes which encode differentially expressed proteins (secretome)
proteom <- read.csv("./ADSCs_prot_dep.csv", row.names = 1)  # change file name
proteom_significant <- subset(proteom, proteom[,4] < 0.05)
proteom_genes_all <- proteom_significant$ID
prot_upreg <- subset(proteom_significant, proteom_significant[,7] > 1)
prot_upreg_genes <- prot_upreg$ID
prot_downreg <- subset(proteom_significant, proteom_significant[,7] < -1)
prot_downreg_genes <- prot_downreg$ID

#extract differentially expressed genes (transcriptome)
transcriptome <- read.csv("./fat_merged_replicates_dif10d_vs_undif.csv")  # change file name
transcriptome_significant <- subset(transcriptome, transcriptome$padj < 0.05)
transcriptome_genes_all <- transcriptome$X
trans_upreg <- subset(transcriptome_significant, transcriptome_significant$log2FoldChange > 1)
trans_upreg_genes <- trans_upreg$X
trans_downreg <- subset(transcriptome_significant, transcriptome_significant$log2FoldChange < -1)
trans_downreg_genes <- trans_downreg$X

#visualize overlaps by Venn diagram
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
        x = list(secretom_significant$ID,
                 proteom_significant$ID,
                 transcriptome_significant$X),
        category.names = c("Secretome", 
                           "Proteome",
                           "Transcriptome (10d)"),
        filename = './adscs_padj005_10d_venn.png',  # change file name
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 700 , 
        width = 700 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .5,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(40, -10, 40),
        cat.dist = c(0.020, 0.040, 0.010),
        cat.fontfamily = "sans",
        rotation = 1
)

Reduce(intersect, list(proteom_significant$ID,
                       transcriptome_significant$X,
                       secretom_significant$ID))


## Compare LogFCs between proteome and secretome data
library(ggplot2)

proteom <- read.csv("./ADSCs_prot_dep.csv", row.names = 1)  # change file name
proteom.logFC <- proteom[,c(2,7)]
colnames(proteom.logFC)[2] <- "logFC.prot"

secretom <- read.csv("./ADSCs_secr_dep.csv", row.names = 1)  # change file name
secretom.logFC <- secretom[,c(2,7)]
colnames(secretom.logFC)[2] <- "logFC.secr"

logFC.dataframe <- merge(proteom.logFC,
                         secretom.logFC,
                         by = "ID")

#tiff(file="./adscs_prot_vs_secr.tiff",  # change file name
#     units = "in",
#     width = 7,
#     height = 7,
#     res=300)
ggplot(logFC.dataframe, 
       aes(x = logFC.prot,
           y = logFC.secr)) +
        geom_point(shape = 16) +
        geom_vline(xintercept = 1, linetype = 2) +
        geom_vline(xintercept = -1, linetype = 2) +
        geom_hline(yintercept = 1, linetype = 2) +
        geom_hline(yintercept = -1, linetype = 2) +
        labs(title = "Comparison of gene fold change between control and differentiated cells\n based on proteome and secretome data") + labs(x = "Proteome logFC") + labs(y = "Secretome logFC") +
        theme(legend.position = "top") +
        geom_text(aes(label=ifelse(logFC.prot>1 & logFC.secr>1,as.character(ID),'')),hjust=0,vjust=-0.5) + 
        geom_text(aes(label=ifelse(logFC.prot<(-1) & logFC.secr<(-1),as.character(ID),'')),hjust=0,vjust=-0.5)
dev.off()
