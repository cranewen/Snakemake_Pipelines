library(dplyr)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(data.table)
library(RColorBrewer)
library(EnhancedVolcano)

# setwd(utils::getSrcDirectory())

# clean up the bioMart cache
#biomartCacheClear()

# Set arguments for Rscript commend, args length are flexiable depending on how many you entered in the command line
args <- commandArgs(trailingOnly = TRUE)

# args[1] is the counts dir where consolidated_counts.csv and col_data.csv stored
# args[2] is the sample_pool
# args[3] is the treatment, cellline etc. combination. Format: "treatment:cellline:..."
# args[4] is the results() function parameters
# args[5] is the diffexp_genes dir
# args[6] is the rlog_expression dir
# args[7] is the PCA dir
# args[8] is the heatmaps dir
# args[9] is the volcano plot dir
# args[10] is ensemble genome name

mm10_dir <- "/mnt/storage/dept/pedonc/Reference/Mouse_gencode_M24_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans"
hg19_dir <- "/mnt/storage/dept/pedonc/Reference/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans"
hg38_dir <- "/mnt/storage/dept/pedonc/Reference/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans"

countData_df <- NULL
if (tail(strsplit(args[1], split = "/", fixed = TRUE)[[1]], n=1) == "sorted"){
  countData_df <- read.csv(paste0(args[1], "consolidated_counts_sorted.csv"), header = TRUE, check.names = FALSE)
}
if (tail(strsplit(args[1], split = "/", fixed = TRUE)[[1]], n=1) == "sorted_rmdup"){
  countData_df <- read.csv(paste0(args[1], "consolidated_counts_sorted_rmdup.csv"), header = TRUE, check.names = FALSE)
}
countData_matrix <- as.matrix(countData_df)
colnames(countData_matrix) <- colnames(countData_df)
rownames(countData_matrix) <- countData_df[, 1]
countData_matrix <- countData_matrix[, -1]
countData_matrix <- countData_matrix[, sort(colnames(countData_matrix))]
mode(countData_matrix) <- "numeric"


# sample_pool is a str, format: "mm160x|mm160z"
# args[2]
sample_pool <- args[2]
sample_pool_list <- as.list(strsplit(sample_pool, split = "|", fixed = TRUE))

countData_row_names <- rownames(countData_matrix)
countData_column_names <- colnames(countData_matrix)
# get the indexes for cellline, eg. mm160X is in column range from 1-6
countData_column_range_index <- which(countData_column_names %in% sample_pool_list[[1]])
# get individual cellline data  
countData <- countData_matrix[, countData_column_range_index]

# colData 
colData <- NULL

if (tail(strsplit(args[1], split = "/", fixed = TRUE)[[1]], n=1) == "sorted"){
  colData <- read.csv(paste0(args[1], "col_data_sorted.csv"), header = TRUE, check.names = FALSE)
}
if (tail(strsplit(args[1], split = "/", fixed = TRUE)[[1]], n=1) == "sorted_rmdup"){
  colData <- read.csv(paste0(args[1], "col_data_sorted_rmdup.csv"), header = TRUE, check.names = FALSE)
}

colData <- colData[match(colnames(countData), colData[,1]),]

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

# using formula to make "design = ~x" dynamically with a variable
design_param <- ''
if (grepl(":", args[3], fixed = TRUE)) {
  design_param <- strsplit(args[3], split = ":", fixed = TRUE)[[1]][1]
} else {
  design_param <- args[3]
}


dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = formula(paste("~", design_param)))
# input args[3], format "treatment:cellline"
# ddsGroupFunc function to detect if there are multi-input 
ddsGroupFunc <- function(cmdArgs) {
  if (grepl(":", cmdArgs, fixed = TRUE)) {
    group_list <- strsplit(cmdArgs, split = ":", fixed = TRUE)[[1]]
    if (length(group_list) == 2) {
      return(factor(paste0(dds[[group_list[1]]], ".", dds[[group_list[2]]])))
    }
    if (length(group_list) == 3) {
      return(factor(paste0(dds[[group_list[1]]], ".", dds[[group_list[2]]], ".", dds[[group_list[3]]])))
    }
    if (length(group_list) == 4) {
      return(factor(paste0(dds[[group_list[1]]], ".", dds[[group_list[2]]], ".", dds[[group_list[3]]], ".", dds[[group_list[4]]])))
    }
  } else {
      return(factor(paste0(dds[[cmdArgs]])))
  }
}

# e.g. dds$group <- factor(paste0(dds[["treatment"]], ".", dds[["cellline"]]))
dds$group <- ddsGroupFunc(args[3])
design(dds) <- ~ group
dds <- DESeq(dds, betaPrior = TRUE)

# Don't write csv file
vsd <- vst(dds, blind=TRUE)
vsdSampleDists <- dist(t(assay(vsd)))
vsdSampleDistMatrix <- as.matrix(vsdSampleDists)
rownames(vsdSampleDistMatrix) <- dds$group
colnames(vsdSampleDistMatrix) <- NULL

png(filename = paste0(args[7], ".png"),  width = 1600, height = 1600)
plotPCA(vsd, intgroup=strsplit(args[3], split = ":", fixed = TRUE)[[1]]) + geom_text(aes(label=name),vjust=2)
dev.off()

png(filename = paste0(args[8], ".png"))
pheatmap(vsdSampleDistMatrix, clustering_distance_rows=vsdSampleDists, clustering_distance_cols=vsdSampleDists, col=colors)
dev.off()

ensembl_genome_local_file <- ""
ensembl_genome_id_name <- ""
if (grepl("mm", args[10], fixed = TRUE)) {
  ensembl_genome_local_file <- mm10_dir
  ensembl_genome_id_name <- read.table(ensembl_genome_local_file, sep = '\t', header = FALSE)
}
if (grepl("hg38", args[10], fixed = TRUE)) {
  ensembl_genome_local_file <- hg38_dir
  ensembl_genome_id_name <- read.table(ensembl_genome_local_file, sep = '\t', header = FALSE)
}
if (grepl("hg19", args[10], fixed = TRUE)) {
  ensembl_genome_local_file <- hg19_dir
  ensembl_genome_id_name <- read.table(ensembl_genome_local_file, sep = '\t', header = FALSE)
}


# Use rlog for count analysis (data range/format easier to work with downstream)
rldAnalysis <- rlog(dds, blind=FALSE)
rldAnalysisSampleDists <- dist(t(assay(rldAnalysis)))
rldAnalysisSampleDistMatrix <- as.matrix(rldAnalysisSampleDists)
rownames(rldAnalysisSampleDistMatrix) <- dds$group
colnames(rldAnalysisSampleDistMatrix) <- NULL
rldAnalysisSampleDist_df <- as.data.frame(assay(rldAnalysis))

#write.csv(as.data.frame(assay(rldAnalysis)), file=paste0(args[6], ".csv"))

# Calculate results using all samples (remember that dds2 is treatment-only and dds is treatment and cellline)

resd <- results(dds, contrast = strsplit(args[4], split = ":", fixed = TRUE)[[1]])
resdsig <- resd[which(resd$padj < 0.1), ]

ensembl_genome_name <- ""
if (grepl("mm", args[10], fixed = TRUE)) {
  ensembl_genome_name = "mmusculus_gene_ensembl"
}
if (grepl("hg", args[10], fixed = TRUE)) {
  ensembl_genome_name = "hsapiens_gene_ensembl"
}

#ensembl <- useEnsembl(biomart = "ensembl", dataset = ensembl_genome_name, mirror = "useast")
resd$ensembl <- sapply( strsplit( rownames(resd), split="\\+" ), "[", 1 )
resd$ensembl <- gsub("\\..*","",resd$ensembl )
#genemap <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = resd$ensembl, mart = ensembl, useCache = FALSE)

idx <- match(rownames(resd), ensembl_genome_id_name$V1)
##################################
rldAnalysisSampleDist_df$gene <- ensembl_genome_id_name$V6[idx]
rldAnalysisSampleDist_df$biotype <- ensembl_genome_id_name$V7[idx]
write.csv(rldAnalysisSampleDist_df, file=paste0(args[6], ".csv"))
##################################
resd$gene <- ensembl_genome_id_name$V6[idx]
resd$biotype <- ensembl_genome_id_name$V7[idx]
# file name placeholder
write.csv(as.data.frame(resd), file=paste0(args[5], ".csv"))

png(filename = paste0(args[9], ".png"))
EnhancedVolcano(resd, lab = resd$gene, x = "log2FoldChange", y = "padj", FCcutoff = 1, pCutoff = 0.05)
dev.off()


