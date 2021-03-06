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

targets <- targets[targets$`Sample Type` != "Solid Tissue Normal",]

df_filtered <- df_filtered[,targets$`File Name`]
names(df_filtered) <- targets$`File ID`

targets_simple <- targets[,c(2,9)]
names(targets_simple) <- c("sample","condition")

df_filtered <- 2^df_filtered

fit <- vsnMatrix(as.matrix(df_filtered))
meanSdPlot(fit)
df_filtered_vsn <- as.data.frame(vsn::predict(fit,as.matrix(df_filtered)))

targets_simple <- targets_simple[targets_simple$condition != "not reported",]

targets_simple$condition <- ifelse(targets_simple$condition == "stage i" | targets_simple$condition == "stage ii", "stage 1 and 2", "stage 3 and 4")

df_filtered_vsn <- df_filtered_vsn[,names(df_filtered_vsn) %in% targets_simple$sample]

limmaRes <- runLimma(df_filtered_vsn, targets_simple, comparisons = list(c(2,-1)))
ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 13418, adjust.method = "fdr"))
ttop$ID <- gsub("[.].*","",ttop$ID)

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_to_symbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ttop$ID, mart = mart)  

ensembl_to_symbol_vec <- ensembl_to_symbol[,2]
names(ensembl_to_symbol_vec) <- ensembl_to_symbol[,1]

ttop_symbol <- ttop

for(i in 1:length(ttop_symbol[,1]))
{
  ttop_symbol[i,1] <- ensembl_to_symbol_vec[ttop_symbol[i,1]]
}

df_filtered_vsn_late <- df_filtered_vsn[,names(df_filtered_vsn) %in% targets_simple[targets_simple$condition == "stage 3 and 4","sample"]]

ASS1_exp <- as.numeric(df_filtered_vsn_late[grepl("ENSG00000130707",row.names(df_filtered_vsn_late)),])
ASS1_exp <- ASS1_exp[complete.cases(ASS1_exp)]
plot(hist(ASS1_exp,breaks = 1000))

length(ASS1_exp[ASS1_exp < 10])

ASS1_exp_df <- as.data.frame(t(data.frame(df_filtered_vsn_late[grepl("ENSG00000130707",row.names(df_filtered_vsn_late)),])))
row.names(ASS1_exp_df) <- names(df_filtered_vsn_late[grepl("ENSG00000130707",row.names(df_filtered_vsn_late)),])
ASS1_exp_df <- ASS1_exp_df[complete.cases(ASS1_exp_df),,drop = F]
ASS1_exp_df <- ASS1_exp_df[ASS1_exp_df$ENSG00000130707.16 > 10,, drop = F]

library(mclust)
plot(hist(ASS1_exp_df$ENSG00000130707.16, breaks = 100))

BIC <- mclustBIC(ASS1_exp_df)

GMM_class <- Mclust(ASS1_exp_df, G = 2, x = BIC, modelNames = "E")
GMM_groups <- GMM_class$classification
GMM_uncertainty <- GMM_class$uncertainty
GMM_z <- as.data.frame(GMM_class$z)

ASS1_exp_df_G1 <- ASS1_exp_df[GMM_z[,1] > 0.5,, drop = F]
ASS1_exp_df_G2 <- ASS1_exp_df[GMM_z[,2] > 0.5,, drop = F]

plot(hist(ASS1_exp_df_G1$ENSG00000130707.16, breaks = 100))
plot(hist(ASS1_exp_df_G2$ENSG00000130707.16, breaks = 100))

patients_G2 <- targets[targets$`File ID` %in% row.names(ASS1_exp_df_G2),]

write(as.character(row.names(ASS1_exp_df_G1)),"results/GDC_KIRK/G1_members.txt")
write(as.character(row.names(ASS1_exp_df_G2)),"results/GDC_KIRK/G2_members.txt")
write_csv(patients_G2, 'results/GDC_KIRK/patients_G2.csv')

targets_simple$condition_2 <- ifelse(targets_simple$sample %in% row.names(ASS1_exp_df_G1),"late_G1",
                                     ifelse(targets_simple$sample %in% row.names(ASS1_exp_df_G2),"late_G2",
                                            ifelse(targets_simple$sample %in% targets[targets$tumor_stage == "stage i",2],"early","less_early")))

unique(targets_simple$condition_2)
limmaRes <- runLimma(df_filtered_vsn, targets_simple[,c(1,3)], comparisons = list(c(2,-1),c(4,-1)))
ttop_G1 <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 13418, adjust.method = "fdr"))
ttop_G1$ID <- gsub("[.].*","",ttop_G1$ID)

ttop_G2 <- ttopFormatter(topTable(limmaRes[[1]], coef = 2, number = 13418, adjust.method = "fdr"))
ttop_G2$ID <- gsub("[.].*","",ttop_G2$ID)

write_csv(ttop_G1, "results/GDC_KIRK/ttop_lateVsEarly_noASS1.csv")
write_csv(ttop_G2, "results/GDC_KIRK/ttop_lateVsEarly_ASS1_reExp.csv")

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_to_symbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ttop$ID, mart = mart)

ensembl_to_symbol_vec <- ensembl_to_symbol[,2]
names(ensembl_to_symbol_vec) <- ensembl_to_symbol[,1]

for(i in 1:length(ttop[,1]))
{
  ttop[i,1] <- ensembl_to_symbol_vec[ttop[i,1]]
}

for(i in 1:length(ttop_G1[,1]))
{
  ttop_G1[i,1] <- ensembl_to_symbol_vec[ttop_G1[i,1]]
}

for(i in 1:length(ttop_G2[,1]))
{
  ttop_G2[i,1] <- ensembl_to_symbol_vec[ttop_G2[i,1]]
}

write_csv(ttop_G1, "results/GDC_KIRK/ttop_lateVsEarly_noASS1.csv")
write_csv(ttop_G2, "results/GDC_KIRK/ttop_lateVsEarly_ASS1_reExp.csv")

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

plotEnrichment(pathway = pathways_list[["KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION"]],
               stats_vec, ticksSize = 0.4)

plotGseaTable(pathways = pathways_list, stats = stats_vec, fgseaRes = fgsea_res)

fgsea_res_top$pval <- -log10(as.numeric(fgsea_res_top$pval))

ggplot(fgsea_res_top, aes(x = NES, y = pathway, color = NES)) + 
  geom_point(aes(size = abs(pval))) + 
  scale_color_gradient2(low="blue", high="red", midpoint = 0, mid = "white") +
  theme_minimal()


#######

metabolic_enzymes <- as.data.frame(read_csv("support/metabolic_enzymes.csv"))

ttop_metab <- ttop[ttop$ID %in% metabolic_enzymes$gene,]

volcano_nice(ttop_metab, FCIndex = 2, pValIndex = 5, IDIndex = 1,label = T, nlabels = 50) + ggtitle("stage III/IV vs stage II/I")

ttop_metab <- ttop_G1[ttop_G1$ID %in% metabolic_enzymes$gene,]

volcano_nice(ttop_metab, FCIndex = 2, pValIndex = 5, IDIndex = 1,label = T, nlabels = 50) + ggtitle("stage III/IV noASS1 vs stage II/I")

ttop_metab <- ttop_G2[ttop_G2$ID %in% metabolic_enzymes$gene,]

volcano_nice(ttop_metab, FCIndex = 2, pValIndex = 5, IDIndex = 1,label = T, nlabels = 50) + ggtitle("stage III/IV ASS1 reExp vs stage II/I")

stats_vec <- ttop_G1$t 
names(stats_vec) <- ttop_G1$ID
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
  geom_point(aes(size = abs(NES))) + 
  scale_color_gradient2(low="blue", high="red", midpoint = 0, mid = "white") +
  theme_minimal() + ggtitle("ASS1 noreExp late Vs early")

stats_vec <- ttop_G2$t 
names(stats_vec) <- ttop_G2$ID
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
  geom_point(aes(size = abs(NES))) + 
  scale_color_gradient2(low="blue", high="red", midpoint = 0, mid = "white") +
  theme_minimal() + ggtitle("ASS1 reExp late Vs early")

