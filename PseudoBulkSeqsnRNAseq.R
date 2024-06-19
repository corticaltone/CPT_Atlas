# Load libraries
#library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
#library(topGO)

choroid.renamed <- readRDS("/omics/odcf/analysis/OE0519_projects/chptumor/FinalAnalysis/Manuscript/QC/choroid.renamed.RDS")

seurat <- choroid.renamed
#combine proliferative and epithelial
seurat <- RenameIdents(seurat, 'Proliferative' = "Epithelial")
seurat <- subset(seurat, subset = type %in% c("ChP", "CPP", "CPC"))
seurat$sample <- droplevels(seurat$sample)
seurat$type <- droplevels(seurat$type)
seurat@assays$rnded <- round(seurat@assays$RNA@counts)
seurat$sample_id <- seurat$sample
seurat$group_id <- seurat$type
seurat$met <- seurat$met.prof
counts <- seurat@assays$rnded
metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id", "age", "met", "status")]

# Explore the raw counts for the dataset
dim(colData(sce))

head(colData(sce))
## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))

## Explore the cellular metadata for the dataset
dim(colData(sce))

head(colData(sce))

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))

# Total number of samples 
ns <- length(sids)
ns


# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select(-"cluster_id")
ei




# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id", "met", "age")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)

pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)

# Explore the different components of list
str(pb)


# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$sample_id)


# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id", "met.prof", "age")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id, met.prof, age) 

metadata    

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

metadata  


# Generate vector of cluster IDs
clusters <- levels(as.factor(metadata$cluster_id))
clusters

clusters[2]

# Subset the metadata to only the Epithelial cells
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[2]), ]
head(cluster_metadata)
# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)

# Subset the counts to only the Epithelial cells
counts <- pb[[clusters[2]]]

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))   


#meth profile
# Run DESeq2 differential expression analysis
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ met.prof)
dds
dds <- DESeq(dds)
dds$met.prof <- factor(dds$met.prof, levels = c("ChP", "adult", "pedB"))
dds$met.prof <- relevel(dds$met.prof, ref = "ChP")
dds <- DESeq(dds)
res <- results(dds)
res


# Plot dispersion estimates
plotDispEsts(dds)


contrast <- c("met.prof", "pedB", "ChP")


res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, coef="met.prof_pedB_vs_ChP", type="apeglm")


# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Write all results to file
write.csv(res_tbl,
          paste0("results/", clusters[2], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
write.csv(res_tbl, "RES_pedBvChPEpithelial.csv", quote = FALSE, row.names = FALSE)
# Set thresholds

padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

topDegs <- filter(sig_res, log2FoldChange > 8.8)
topDegs <- topDegs[sort(topDegs$log2FoldChange, index.return = TRUE, decreasing = TRUE)$ix,]
topDegs <- dplyr::select(topDegs, c(gene, log2FoldChange, padj))
bottomDegs <- filter(sig_res, log2FoldChange > 8.8)
bottomDegs <- topDegs[sort(bottomDegs$log2FoldChange, index.return = TRUE, decreasing = TRUE)$ix,]
bottomDegs <- dplyr::select(bottomDegs, c(gene, log2FoldChange, padj))
# Write significant results to file
write.csv(sig_res,
          paste0("results", clusters[1], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
write.csv(sig_res, "RES_pedBvChPEpithelialSIG.csv", quote = FALSE, row.names = FALSE)

## ggplot of top genes
normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values and log2FoldChange
top40_sig_genes <- sig_res %>%
  filter(baseMean > 99) %>%
  dplyr::arrange(padj) %>%
  head(n=40) %>%
  dplyr::pull(gene) 

top20_sig_genes <- top40_sig_genes[c(1:11,13, 15:22)]

# order genes and remove non-coding genes (that are absent in Thomas et. al bulk seq data)
SigOrder <- c("PRLR", "GPR155", "CAB39L", "DIRC3", "APOD", "TMEM64", "SPECC1", "CLDN2", "LMCD1", "SPTB", "LRAT", "NABP1", "GNB3", "SDK2", "SLC1A7", "SAMD3","KIAA1211", "WNT5B", "MAML3", "IL1RAPL2", "TMTC2", "MID1")




top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)

gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top20_sig <- inner_join(ei[, c("sample_id", "met.prof", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))

## plot using ggplot2
gathered_top20_sig <- mutate(gathered_top20_sig, normalized_counts = normalized_counts+1)


df <- gathered_top20_sig %>% filter(met.prof %in% c("pedB", "ChP"))
df$gene <- factor(df$gene, levels = SigOrder)
ggplot(df) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = met.prof), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes snPEDB v snChP") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(text = element_text(size = 18)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#56B4E9", "#E69F00", "#33FF00"))
p+scale_fill_manual(
  ...,
  values,
  aesthetics = "fill",
  breaks = waiver(),
  na.value = "grey50"
)

# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)
# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata[, c("group_id", "met.prof")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 

pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata[, c("group_id", "age")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 


## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table_thres <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

#Venn
res1 <- read.table("RES_pedBvChPEpithelial.csv", sep = ",", header = TRUE)
CPCsigs <- filter(res1, padj < 0.05)

res2 <- read.table("RES_adultvChPEpithelial.csv", sep = ",", header = TRUE)
CPPsigs <- filter(res2, padj < 0.05)

CPCup <- filter(CPCsigs, log2FoldChange > 1)[,1]
CPCdn <- filter(CPCsigs, log2FoldChange < -1)[,1]
CPPup <- filter(CPPsigs, log2FoldChange > 1)[,1]
CPPdn <- filter(CPPsigs, log2FoldChange < -1)[,1]

x <- list(
  'pedB up'=CPCup,
  'pedB dn'=CPCdn,
  'adult up'=CPPup,
  'adult dn'=CPPdn
)

ggvenn(x, c('pedB up', 'adult up'), fill_color = c("red", "blue"), show_percentage = FALSE)
ggvenn(x, c('pedB dn', 'adult dn'), fill_color = c("red", "blue"), show_percentage = FALSE)


CPPdegsUP <- filter(CPPdegsUP, log2FoldChange > 1)
CPPbackground <- filter(res, baseMean > 2)
CPCdegsUP <- filter(res, padj > 0.05)
CPCdegsDN <- filter(CPCdegsUP, log2FoldChange < -1)
CPCdegsUP <- filter(CPCdegsUP, log2FoldChange > 1)
CPCbackground <- filter(res, baseMean > 1.5)

write.table(CPCdegsDN, file = "CPCdegsDN.csv", sep = ",")
write.table(CPCdegsDN$gene, file = "CPCdeggenesDN.csv", sep = ",")
write.table(CPCdegsUP, file = "CPCdegsUP.csv", sep = ",")
write.table(CPCdegsUP$gene, file = "CPCdeggenesUP.csv", sep = ",")






#Volcano
res <- read.table("RES_CPCvChPEpithelial.csv", sep = ",", header = TRUE)
CPCsigs <- filter(res, padj < 0.1)

res2 <- read.table("RES_CPPvChPEpithelial.csv", sep = ",", header = TRUE)
CPPsigs <- filter(res, padj < 0.1)

library(EnhancedVolcano)
EnhancedVolcano(res_tbl,
                lab = res_tbl$gene,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'CPP v ChP',
                pCutoff = 0.001,
                FCcutoff = 1.5,
                pointSize = 1.0,
                labSize = 3.0,
                col=c('grey', 'black', 'blue', 'red3'),
                colAlpha = 1, 
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'))




# Enrichment analysis
library(EnsDb.Hsapiens.v79)

#load genes
CPCcomp <- read.table("RES_CPCvChPEpithelial.csv", sep = ",", header = TRUE)
CPPcomp <- read.table("RES_CPCvChPEpithelial.csv", sep = ",", header = TRUE)
#filter poorly expressed genes
CPCcomp <- filter(CPCcomp, baseMean > 1.5)
CPPcomp <- filter(CPPcomp, baseMean > 1.5)
#extract gene names
CPCbaseGenes <- CPCcomp$gene
CPPbaseGenes <- CPPcomp$gene
CPCupGenes <- CPCcomp %>% filter(padj < 0.1) %>% filter(log2FoldChange > 1) %>% dplyr::select(gene)
CPPupGenes <- CPPcomp %>% filter(padj < 0.1) %>% filter(log2FoldChange > 1) %>% select(gene)
CPCupGenes5 <- CPCcomp %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1) %>% select(gene)
CPPupGenes5 <- CPPcomp %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1) %>% select(gene)
CPCdnGenes <- CPCcomp %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1) %>% select(gene)
CPPdnGenes <- CPPcomp %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1) %>% select(gene)
CPCdnGenes5 <- CPCcomp %>% filter(padj < 0.05) %>% filter(log2FoldChange< -1) %>% select(gene)
CPPdnGenes5 <- CPPcomp %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1) %>% select(gene)
# 2. Convert from gene.symbol to ensembl.gene
geneSymbols <-  c('DDX26B','CCDC83',  'MAST3', 'RPL11', 'ZDHHC20',  'LUC7L3',  'SNORD49A',  'CTSH', 'ACOT8')

CPCupGenes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPCupGenes[[1]], keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPCupGenes5 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPCupGenes5[[1]], keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPPupGenes5 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPPupGenes5[[1]], keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPPupGenes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPPupGenes[[1]], keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPCdnGenes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPCdnGenes[[1]], keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPCdnGenes5 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPCdnGenes5[[1]], keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPPdnGenes5 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPPdnGenes[[1]], keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPPdnGenes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPPdnGenes5[[1]], keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPCbaseGenes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPCbaseGenes, keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
CPPbaseGenes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CPPbaseGenes, keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))






CPCupGenes5 <- CPCcomp %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1)
CPPupGenes5 <- CPPcomp %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1) 
CPPupGenes <- CPPcomp %>% filter(padj < 0.10) %>% filter(log2FoldChange > 1) 
CPCdnGenes <- CPCcomp %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1) 
CPPdnGenes <- CPPcomp %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1)
CPCdnGenes5 <- CPCcomp %>% filter(padj < 0.05) %>% filter(log2FoldChange< -1) 
CPPdnGenes5 <- CPPcomp %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1) 
CPCupGenes <- CPCupGenes %>% select(gene)

setwd("/omics/odcf/analysis/OE0519_projects/chptumor/FinalAnalysis/CPT/BulkSeq/DegEnrich")

write.table(CPCupGenes, file = "CPCupGens.csv", sep = ",")
write.table(CPCupGenes5, file = "CPCupGens5.csv", sep = ",")
write.table(CPCdnGenes, file = "CPCdnGens.csv", sep = ",")
write.table(CPCdnGenes5, file = "CPCdnGens5.csv", sep = ",")
write.table(CPPupGenes, file = "CPPupGens.csv", sep = ",")
write.table(CPPupGenes5, file = "CPPupGens5.csv", sep = ",")
write.table(CPPdnGenes, file = "CPPdnGens.csv", sep = ",")
write.table(CPPdnGenes5, file = "CPPdnGens5.csv", sep = ",")
write.table(CPCbaseGenes, file = "CPCbaseGenes.csv", sep = ",")
write.table(CPPbaseGenes, file = "CPPbaseGenes.csv", sep = ",")




write.table(as.matrix(GetAssayData(object = choroid.renamed, slot = "counts")), 
            'choroid.renamedCounts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)


#adult

contrast <- c("met.prof", "adult", "ChP")


res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, coef="met.prof_adult_vs_ChP", type="apeglm")


# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Write all results to file
write.csv(res_tbl,
          paste0("results/", clusters[2], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
write.csv(res_tbl, "RES_adultvChPEpithelial.csv", quote = FALSE, row.names = FALSE)
# Set thresholds

padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

topDegs <- filter(sig_res, log2FoldChange > 8.8)
topDegs <- topDegs[sort(topDegs$log2FoldChange, index.return = TRUE, decreasing = TRUE)$ix,]
topDegs <- dplyr::select(topDegs, c(gene, log2FoldChange, padj))
bottomDegs <- filter(sig_res, log2FoldChange > 8.8)
bottomDegs <- topDegs[sort(bottomDegs$log2FoldChange, index.return = TRUE, decreasing = TRUE)$ix,]
bottomDegs <- dplyr::select(bottomDegs, c(gene, log2FoldChange, padj))
# Write significant results to file
write.csv(sig_res,
          paste0("results", clusters[1], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
write.csv(sig_res, "RES_adultvChPEpithelialSIG.csv", quote = FALSE, row.names = FALSE)

# Plot PCA
rld <- rlog(dds, blind=TRUE)
DESeq2::plotPCA(rld, intgroup = "met.prof")


pcaData <- plotPCA(rld, intgroup=c("met.prof", "group_id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group_id, shape=met)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("ChP" = "green", "CPC" = "red", "CPP" = "blue"))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id", "met.prof"), drop=F], legend = FALSE) 



## ggplot of top genes from snRNAseq

  
  
normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top40_sig_genes <- sig_res %>%
  filter(baseMean > 100) %>%
  dplyr::arrange(padj) %>%
  head(n=40) %>%
 dplyr::pull(gene) 

top20_sig_genes <- top40_sig_genes[1:25]
top20_sig_genes <- top20_sig_genes[c(1,3:11,14:15, 17, 19:25)]


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)



gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top20_sig <- inner_join(ei[, c("sample_id", "group_id", "met.prof" )], gathered_top20_sig, by = c("sample_id" = "samplename"))
## plot using ggplot2
gathered_top20_sig <- mutate(gathered_top20_sig, normalized_counts = normalized_counts+1)

gathered_top20_sig <- subset(gathered_top20_sig, subset = met.prof %in% c("ChP", "adult"))
gathered_top20_sig$gene <- factor(gathered_top20_sig$gene, levels = c("PHACTR2", "GPR155", "CRYBG1", "CA12", "SLC16A12",  "LINC00598", "SLC15A4", "TGFBR2", "SLC19A1", "LRAT", "CACNA2D4", "SAMD3", "AFF2", "TMPRSS7", "TMEM163", "WNT5B", "MAML3", "HOTAIRM1", "MAMDC2", "ADCYAP1R1", "PRR5L", "LY75", "RCAN3AS", "NR2E3", "MAOA", "CXXC4", "MEGF11", "IQGAP2", "NPNT", "MID1"))
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = met.prof), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes snADULT v snChP") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))



