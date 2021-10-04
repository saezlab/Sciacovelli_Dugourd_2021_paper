rm(list = ls(all.names = TRUE))
gc()
.rs.restartR() ##Only if Rstudio

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

targets$tumor_stage <- ifelse(targets$`Sample Type` == "Solid Tissue Normal", "stage 0", targets$tumor_stage)

# targets <- targets[targets$`Sample Type` != "Solid Tissue Normal",]

df_filtered <- df_filtered[,targets$`File Name`]
names(df_filtered) <- targets$`File ID`

targets_simple <- targets[,c(2,9)]
names(targets_simple) <- c("sample","condition")

df_filtered <- 2^df_filtered

fit <- vsnMatrix(as.matrix(df_filtered))
meanSdPlot(fit)
df_filtered_vsn <- as.data.frame(vsn::predict(fit,as.matrix(df_filtered)))

targets_simple <- targets_simple[targets_simple$condition != "not reported",]

# targets_simple$condition <- ifelse(targets_simple$condition == "stage i" | targets_simple$condition == "stage ii", "stage 1 and 2", "stage 3 and 4")

df_filtered_vsn <- df_filtered_vsn[,names(df_filtered_vsn) %in% targets_simple$sample]
# pca_plots <- nicePCA(df_filtered_vsn, targets_simple, no_label = T)
# plot(pca_plots)

unique(targets_simple$condition)
comparisons <- list("stage i"  = c(1,-4),
                    "stage ii"  = c(3,-4),
                    "stage iii"  = c(2,-4),
                    "stage iv"  = c(5,-4))

limmaRes <- runLimma(df_filtered_vsn, targets_simple, comparisons = comparisons)

t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(df_filtered_vsn[,1]),
                         adjust.method = "fdr"))

t_table$ID <- gsub("[.].*","",t_table$ID)

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_to_symbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = t_table$ID, mart = mart)  

ensembl_to_symbol_vec <- ensembl_to_symbol[,2]
names(ensembl_to_symbol_vec) <- ensembl_to_symbol[,1]

t_table_symbol <- t_table

for(i in 1:length(t_table_symbol[,1]))
{
  t_table_symbol[i,1] <- ensembl_to_symbol_vec[t_table_symbol[i,1]]
}
