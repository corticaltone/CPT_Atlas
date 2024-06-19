library(dplyr)
library(DESeq2)
library(org.Hs.eg.db)
library(scales)

data_dir <- "/omics/odcf/analysis/OE0519_projects/chptumor/BulkSeq/"
load(paste0(data_dir, "scripts/bsDDS-METH.Robj"))
gene_table <- readRDS(paste0(data_dir, "scripts/GeneTable.RDS"))
#make sample table
md_table <- read.table(paste0(data_dir, "Samples_Paper.txt"), sep = "\t", header = TRUE)
coldata <- md_table[1:53,] %>% dplyr::select("KRYO_ID", "Paper_ID", "DIAGNOSIS", "AGE_YEARS", "SEX", "Family.Member", "Total.CNVs", "recurrence", "JAK2", "PRKCA", "TERT", "TP53", "P53.IH")
colnames(coldata) <- c("tumor", "sample", "type", "age", "sex", "Meth", "CNVs", "recurrence", "JAK2", "PRKCA", "TERT", "TP53", "P53.IH")
coldata$sample <- gsub("-", ".", coldata$sample)
#add controls
coldata$sample[c(51:53)] <- c("ChP.1", "ChP.2", "ChP.3")
coldata$age[c(51:53)] <- c(85, 47, 78)
coldata$sex[c(51:53)] <- c("male", "female", "male")
coldata$type[c(51:53)] <- c("ChP", "ChP", "ChP")
coldata$Meth[c(51:53)] <- rep("ChP", 3)
coldata$status <- c(rep("CPT", 50), rep("ChP", 3))
samples <- coldata$sample

ct_table <- read.table(paste0(data_dir, "cpt.counts.txt"), sep = "\t", header = TRUE)
ct_table <- ct_table %>% filter(!gene %in% c("__too_low_aQual", "__not_aligned", "__alignment_not_unique"))
gene <- ct_table$gene
cts <- ct_table[,samples]
rownames(cts) <- gene

#get gene symbols




#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~Meth)

dds$Meth <- factor(dds$Meth, levels = c("PED_B", "PED_A","AD", "ChP"))
dds$Meth <- relevel(dds$Meth, ref = "ChP")
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients
resUS <- results(dds, name="Meth_PED_B_vs_ChP")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="Meth_PED_B_vs_ChP", type="apeglm")
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

resOrdered <- res[order(res$pvalue),]

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)


write.csv(as.data.frame(resOrdered), 
          file="condition_pedBvChP_results.csv")
resSig <- subset(resOrdered, padj < 0.05)
resSig
write.csv(as.data.frame(resSig), 
          file="condition_pedBvChP_resultsSIG.csv")


top_Degs <- data.frame(resSig) %>% filter(padj < 0.5) %>% filter(abs(lfcSE) > 1)
data <- rownames(top_Degs)
data <- substr(data,1,15)

genes <- select(org.Hs.eg.db, keys=data, 
                 columns="SYMBOL", keytype="ENSEMBL")
genes <- genes[!duplicated(genes[,1]),]
rownames(genes) <- genes[,1]
genes$revision <- 
symbol <- genes[data,2]
top_Degs$symbol <- symbol
CPCvChPdegs <- na.omit(top_Degs)$symbol
write.table(top_Degs, file = "PEDBvChPtopTable.txt", sep = "/t")
save(CPCvChPdegs, file = "PEDBvChPdegGenes.Robj")

#Repeat for CPP
#create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~type)

dds$type <- factor(dds$type, levels = c("CPC","aCPP", "CPP", "ChP"))
dds$type <- relevel(dds$type, ref = "ChP")
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients
resUS <- results(dds, name="type_CPP_vs_ChP")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="type_CPP_vs_ChP", type="apeglm")
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

resOrdered <- res[order(res$pvalue),]

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)


write.csv(as.data.frame(resOrdered), 
          file="condition_CPPvChP_results.csv")
resSig <- subset(resOrdered, padj < 0.1)
resSig
write.csv(as.data.frame(resSig), 
          file="condition_CPPvChP_resultsSIG.csv")

