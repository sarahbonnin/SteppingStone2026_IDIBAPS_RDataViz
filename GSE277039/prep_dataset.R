# Prepare GSE277039 dataset for IDIBAPS course
# April 2026

BiocManager::install("GEOquery")

library(GEOquery)
library(dplyr)
library(stringr)

gse <- "GSE277039"

# get series matrix
gse <- getGEO(gse)

head(as(gse[[1]]@phenoData, "data.frame"))

gse_annot <- as(gse[[1]]@phenoData, "data.frame") 

# simplify, for course, and rename selected columns
gse_annot_simple <- gse_annot %>% 
  rename_with(~str_remove(., ':ch1|_ch1')) %>% 
  rename(cell_type=`cell type`) %>%
  dplyr::select(geo_accession, title, source_name, genotype, treatment, cell_type) %>%
  mutate(group=gsub("[0-9]*", "", title))

readr::write_csv(gse_annot_simple, "GSE277039/metadata.csv") 

# Read in raw counts
raw <- read.table("GSE277039/GSE277039_raw_counts.txt.gz", header = T, row.names = 1)

# format metadata properly so it can be imported in DESeq
meta <- gse_annot_simple
rownames(meta) <- meta$title


# Create DESeq object
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = raw, colData = meta, design = ~ group)

# filter
smallestGroupSize <- min(table(meta$group))
keep <- rowSums(counts(dds) >= 25) >= smallestGroupSize
dds <- dds[keep,]


## Run analysis
dds <- DESeq(dds)

## Check the size factors
sizeFactors(dds)
# WT1_bas   WT2_bas   WT3_bas   WT4_bas   WT5_bas  WT1_stim  WT2_stim  WT3_stim  WT4_stim  WT5_stim   KO1_bas 
# 1.0669821 0.9217579 1.1425440 1.1654818 0.7724094 1.3006658 1.3029119 0.9226698 0.8458719 1.0212337 1.3649310 
# KO2_bas   KO3_bas   KO4_bas   KO5_bas   KO6_bas   KO7_bas   KO8_bas  KO1_stim  KO2_stim  KO3_stim  KO4_stim 
# 0.9676636 1.0753153 1.2682063 1.0462741 0.8471065 0.9444999 0.9475567 0.8664291 0.9028415 0.8477821 1.1732781 
# KO5_stim  KO6_stim  KO7_stim  KO8_stim 
# 0.8513128 0.7892406 1.4241475 0.8711420 

## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)


#res1 <- results(dds, contrast=c("group", "KO_stim", "KO_bas"))
#summary(res1, alpha=0.05)

res2 <- results(dds, contrast=c("group", "KO_bas", "WT_bas"))
summary(res2, alpha=0.05)

# add and select columns
res2_sel <- as.data.frame(res2) %>% 
  filter(complete.cases(.)) %>%
  dplyr::select(log2FoldChange, padj) %>%
  mutate(GeneSymbol=rownames(.), 
         DE=ifelse(padj < 0.05 & log2FoldChange > 0, "UP", ifelse(padj < 0.05 & log2FoldChange < 0, "DOWN", "NO"))) %>%
  dplyr::relocate(GeneSymbol)

# write to file
readr::write_csv(res2_sel, "GSE277039/KO_vs_WT_stats.csv")

norm_counts <- log2(counts(dds, normalized=TRUE)+1)
# 19205 genes

# select only genes that were kept in res2_sel
norm_counts <- norm_counts[rownames(norm_counts) %in% rownames(res2_sel),]
# 14363 genes
# select only KO_bas and WT_base samples, and rename columns to simplify dataset
norm_counts <- as.data.frame(norm_counts) %>% 
  dplyr::select(ends_with("bas")) %>%
  rename_with(~str_remove(., '_bas')) %>%
  tibble::rownames_to_column("GeneSymbol")
  
# save counts
readr::write_csv(norm_counts, "GSE277039/counts.csv")

# Put together counts and DEA into a single file
counts_de <- merge(res2_sel, norm_counts, by="GeneSymbol")

readr::write_csv(counts_de, "GSE277039/DEG_counts.csv")

# also write to Excel for the read Excel example
readr::write_excel_csv(counts_de, "GSE277039/DEG_counts.xlsx")

# Write smaller random file for the first plots
set.seed(11)

sample_de <- counts_de[sample(1:nrow(counts_de))[1:50],]

readr::write_csv(sample_de, "GSE277039/DEG_counts_sample.csv")
readr::write_excel_csv(sample_de, "GSE277039/DEG_counts-sample.xlsx")


