library(readr)
library(decoupleR)
library(GSEABase)
library(fgsea)

source("scripts/support_functions.R")

fc_table <- as.data.frame(read_csv("results/RNA/fc_table.csv"))

DA_tumorVsHealthy <- as.data.frame(read_csv("~/Dropbox/marco_kidney_singlecell/results/DA_tumorVsHealthy.csv"))
DA_ASS1vsNoASS1 <- as.data.frame(read_csv("~/Dropbox/marco_kidney_singlecell/results/DA_ASS1vsNoASS1.csv"))

DA_tumorVsHealthy <- DA_tumorVsHealthy[which(DA_tumorVsHealthy$contrast == "Epithelium_and_Vascular"),]
DA_ASS1vsNoASS1 <- DA_ASS1vsNoASS1[which(DA_ASS1vsNoASS1$contrast),]

DA_tumorVsHealthy <- DA_tumorVsHealthy[,c(3,4)]
DA_ASS1vsNoASS1 <- DA_ASS1vsNoASS1[,c(3,4)]

names(DA_tumorVsHealthy) <- c("ID","tumorVsHealthy")
names(DA_ASS1vsNoASS1) <- c("ID","tumor_ASS1vsTumor")

fc_table <- merge(fc_table,DA_tumorVsHealthy, by = "ID")
fc_table <- merge(fc_table,DA_ASS1vsNoASS1, by = "ID")

fc_table <- fc_table[which(fc_table$tumorVsHealthy != 0),]
fc_table <- fc_table[which(fc_table$tumor_ASS1vsTumor != 0),]


fc_table_ranks <- fc_table
fc_table_ranks[,-1] <- apply(fc_table_ranks[,-1], 2, rank)

pheatmap::pheatmap(fc_table_ranks[,-1])

pathways <- import_gmt("support/c2.cp.v7.0.symbols.gmt")

row.names(DA_tumorVsHealthy) <- DA_tumorVsHealthy$ID
row.names(DA_ASS1vsNoASS1) <- DA_ASS1vsNoASS1$ID
DA_tumorVsHealthy <- DA_tumorVsHealthy[,-1,drop = F]
DA_ASS1vsNoASS1 <- DA_ASS1vsNoASS1[,-1,drop = F]

pathways <- as.data.frame(
  decoupleR::intersect_regulons(mat = DA_tumorVsHealthy,network = pathways, .source = "term", .target = "gene", minsize = 20))

pathways_list <- list()
i <- 1
for(pathway in unique(pathways$term))
{
  pathways_list[[i]] <- pathways[pathways$term == pathway,"gene"]
  i <- i+1
}
names(pathways_list) <- unique(pathways$term)

stats_vec <- DA_tumorVsHealthy$tumorVsHealthy
names(stats_vec) <- row.names(DA_tumorVsHealthy)
fgsea_res_tumorVsHealthy <- fgsea(pathways = pathways_list, stats = stats_vec, nperm = 1000, nproc = 3)
fgsea_res_tumorVsHealthy <- fgsea_res_tumorVsHealthy[,-8]

stats_vec <- DA_ASS1vsNoASS1$tumor_ASS1vsTumor
names(stats_vec) <- row.names(DA_ASS1vsNoASS1)
fgsea_res_ASS1vsNoASS1 <- fgsea(pathways = pathways_list, stats = stats_vec, nperm = 1000, nproc = 3)
fgsea_res_ASS1vsNoASS1 <- fgsea_res_ASS1vsNoASS1[,-8]

### HALLMARKS
hallmarkss <- import_gmt("support/h.all.v7.5.1.symbols.gmt")

hallmarkss <- as.data.frame(
  decoupleR::intersect_regulons(mat = DA_tumorVsHealthy,network = hallmarkss, .source = "term", .target = "gene", minsize = 20))

hallmarkss_list <- list()
i <- 1
for(hallmarks in unique(hallmarkss$term))
{
  hallmarkss_list[[i]] <- hallmarkss[hallmarkss$term == hallmarks,"gene"]
  i <- i+1
}
names(hallmarkss_list) <- unique(hallmarkss$term)

stats_vec <- DA_tumorVsHealthy$tumorVsHealthy
names(stats_vec) <- row.names(DA_tumorVsHealthy)
fgsea_hallmarks_tumorVsHealthy <- fgsea(pathways = hallmarkss_list, stats = stats_vec, nperm = 1000, nproc = 3)
fgsea_hallmarks_tumorVsHealthy <- fgsea_hallmarks_tumorVsHealthy[,-8]

stats_vec <- DA_ASS1vsNoASS1$tumor_ASS1vsTumor
names(stats_vec) <- row.names(DA_ASS1vsNoASS1)
fgsea_hallmarks_ASS1vsNoASS1 <- fgsea(pathways = hallmarkss_list, stats = stats_vec, nperm = 1000, nproc = 3)
fgsea_hallmarks_ASS1vsNoASS1 <- fgsea_hallmarks_ASS1vsNoASS1[,-8]
