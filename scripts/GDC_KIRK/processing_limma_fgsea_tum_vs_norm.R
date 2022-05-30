# rm(list = ls(all.names = TRUE))
# gc()
# .rs.restartR() ##Only if Rstudio

library(readr)
library(vsn)
library(limma)
library(biomaRt)
library(fgsea)
library(clusterProfiler)
library(pheatmap)
library(GSEABase)
library(ggplot2)
library(ggrepel)

source("scripts/support_functions.R")

df <- as.data.frame(read_csv("data/GDC_KIRK/GDC_KIRK_raw_20191024.csv")) ##Assembled from download https://portal.gdc.cancer.gov/projects/TCGA-KIRC // clinical can be downloaded after adding files to cart
row.names(df) <- df$gene
df <- df[,-which(names(df) == "gene")]

df <- log2(df)

plot(density(as.matrix(df)))

df[df < 7.5] <- NA

plot(colSums(is.na(df)))
df <- df[,colSums(is.na(df)) <= 49000]

df_filtered <- df[rowSums(is.na(df)) < 350,]

plot(density(as.matrix(df_filtered), na.rm = T))

##########

targets <- as.data.frame(
  read_delim("support/gdc_sample_sheet.2019-10-24.tsv", 
             "\t", escape_double = FALSE, trim_ws = TRUE))
targets$`File Name` <- tolower(gsub("[.]gz","",targets$`File Name`))

targets <- targets[targets$`File Name` %in% names(df_filtered),]
targets <- targets[targets$`Sample Type` != "Additional - New Primary",]

file_to_id <- targets[,1]
names(file_to_id) <- targets[,2]

for(i in 1:length(df_filtered[1,]))
{
  if(names(df_filtered[i]) %in% names(file_to_id))
  {
    names(df_filtered[i]) <- file_to_id[names(df_filtered[i])]
  }
}

clinical <- as.data.frame(read_delim("support/clinical.tsv", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE))

clinical <- unique(clinical)

clinical_reduced <- clinical[,c("case_submitter_id","tumor_stage")]
names(clinical_reduced)[1] <- "Case ID"

targets <- merge(targets,clinical_reduced, by = "Case ID")
targets <- unique(targets)


df_filtered <- 2^df_filtered

fit <- vsnMatrix(as.matrix(df_filtered))
meanSdPlot(fit)
df_filtered_vsn <- as.data.frame(vsn::predict(fit,as.matrix(df_filtered)))

# pca_plots <- nicePCA(df_filtered_vsn, targets_simple, no_label = T)
# plot(pca_plots)
targets_simple <- targets[,c(3,8)]
names(targets_simple) <- c("sample","condition")

df_filtered_vsn <- df_filtered_vsn[,targets_simple$sample]

limmaRes <- runLimma(df_filtered_vsn, targets_simple, comparisons = list(c(1,-2)))
ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 13418, adjust.method = "fdr"))
ttop$ID <- gsub("[.].*","",ttop$ID)

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_to_symbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ttop$ID, mart = mart)  

ensembl_to_symbol_vec <- ensembl_to_symbol[,2]
names(ensembl_to_symbol_vec) <- ensembl_to_symbol[,1]

for(i in 1:length(ttop[,1]))
{
  ttop[i,1] <- ensembl_to_symbol_vec[ttop[i,1]]
}

pathways <- import_gmt("support/c2.cp.v7.0.symbols.gmt")
pathways_list <- list()
i <- 1
for(pathway in unique(pathways$term))
{
  pathways_list[[i]] <- pathways[pathways$term == pathway,"gene"]
  i <- i+1
}
names(pathways_list) <- unique(pathways$term)
pathways_list <- pathways_list[grepl("KEGG",names(pathways_list))]


stats_vec <- ttop$t 
names(stats_vec) <- ttop$ID
fgsea_res <- fgsea(pathways = pathways_list, stats = stats_vec, nperm = 1000, nproc = 3)
fgsea_res <- fgsea_res[,-8]
fgsea_res$pathway <- tolower(gsub("_"," ",fgsea_res$pathway))
fgsea_res$pathway <- tolower(gsub("kegg ","",fgsea_res$pathway))

fgsea_res <- fgsea_res[order(fgsea_res$padj, decreasing = F),]
fgsea_res_top <- fgsea_res[1:30,]
fgsea_res_top <- fgsea_res_top[order(fgsea_res_top$NES),]
fgsea_res_top$pathway <- factor(fgsea_res_top$pathway, levels = fgsea_res_top$pathway)

fgsea_res_top$pval <- -log10(as.numeric(fgsea_res_top$pval))

ggplot(fgsea_res_top, aes(x = NES, y = pathway, color = NES)) + 
  geom_point(aes(size = abs(pval))) + 
  scale_color_gradient2(low="blue", high="red", midpoint = 0, mid = "white") +
  theme_minimal()
