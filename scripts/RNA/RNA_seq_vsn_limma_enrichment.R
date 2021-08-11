rm(list = ls(all.names = TRUE))
gc()
.rs.restartR() ##Only if Rstudio

library(vsn)
library(limma)
library(readr)
library(biomaRt)
library(dorothea)
library(ggplot2)
library(reshape2)
library(viper)
library(GSEABase)

source("scripts/support_functions.R")

count_df_integrated <- as.data.frame(read_csv("data/RNA/counts.csv")) #Rsubread, default parameters, with Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa
count_df_integrated$geneID <- as.character(count_df_integrated$geneID)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_list=unique(count_df_integrated$geneID)
mapping_df=getBM(attributes = c("entrezgene_id", "hgnc_symbol"), filters = "entrezgene_id", values = gene_list, bmHeader = T, mart = mart)
colnames(mapping_df)=c("entrezgene_id", "hgnc_symbol")
mapping_df$entrezgene_id <- as.character(mapping_df$entrezgene_id)
mapping_vec <- mapping_df$hgnc_symbol
names(mapping_vec) <- mapping_df$entrezgene_id

count_df_integrated <- count_df_integrated[count_df_integrated$geneID %in% mapping_df$entrezgene_id,]

for(i in 1:length(count_df_integrated[,1]))
{
  if(count_df_integrated[i,1] %in% names(mapping_vec))
  {
    count_df_integrated[i,1] <- mapping_vec[count_df_integrated[i,1]]
  }
}

count_df_integrated <- count_df_integrated[count_df_integrated$geneID != "",]

doublons <- count_df_integrated[duplicated(count_df_integrated$geneID),1]
count_df_integrated <- count_df_integrated[!(count_df_integrated$geneID %in% doublons),]

row.names(count_df_integrated) <- count_df_integrated$geneID
count_df_integrated <- count_df_integrated[,-1]

count_df_integrated[is.na(count_df_integrated)] <- 0

mean_count_df <- sapply(c(1,4,7), function(x, count_df_integrated){
  rowMeans(count_df_integrated[,c(x,x+1,x+2)])
}, count_df_integrated = count_df_integrated)

to_long <- melt(log2(mean_count_df))
ggplot(to_long, aes(x = value, group = Var2, fill = Var2)) + geom_density(alpha = 0.3)

mean_count_df <- as.data.frame(apply(mean_count_df, 2, function(x){
  ifelse(log2(x) < 5, FALSE, TRUE)
}))

mean_count_df$hgnc_symbol <- row.names(mean_count_df)
mean_count_df <- merge(mean_count_df, mapping_df)
names(mean_count_df)[c(2,3,4)] <- c("M1A","O","HK2")
mean_count_df <- mean_count_df[,c(4,3,2,1,5)]

write_csv(mean_count_df, "data/metabolomic/expressed_genes_plasmax.csv")

count_df_integrated <- count_df_integrated[rowMeans(count_df_integrated) > 50,]

count_df_integrated[count_df_integrated == 0] <- 0.5

fit <- vsnMatrix(as.matrix(count_df_integrated))
meanSdPlot(fit)
count_df_integrated_vsn <- as.data.frame(vsn::predict(fit,as.matrix(count_df_integrated)))

targets <- as.data.frame(matrix(NA,length(count_df_integrated[1,]),2))
names(targets) <- c("sample","condition")
targets$sample <- names(count_df_integrated)
targets$condition <- gsub("_[0-9]$","",targets$sample)

to_write <- count_df_integrated_vsn
to_write$gene <- row.names(to_write)

write_csv(to_write,'data/RNA/counts_processed.csv')
write_csv(targets,'support/targets_RNA.csv')

limmaRes <- runLimma(count_df_integrated_vsn, targets, comparisons = list(c(1,-3),c(2,-3),c(1,-2)))

ttop_786vHK2 <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 12219, adjust.method = "fdr"))
ttop_OvHK2 <- ttopFormatter(topTable(limmaRes[[1]], coef = 2, number = 12219, adjust.method = "fdr"))
ttop_786vO <- ttopFormatter(topTable(limmaRes[[1]], coef = 3, number = 12219, adjust.method = "fdr"))


dorothea <- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c('A','B','C','D'),c(3,1,4)]) #ABC

names(dorothea) <- c("target_genesymbol", "source_genesymbol", "sign")

regulon_viper <- df_to_viper_regulon(dorothea)

eset <- merge(ttop_OvHK2[,c('ID','t')], ttop_786vHK2[,c('ID','t')], by = 'ID')
eset <- merge(eset, ttop_786vO[,c('ID','t')], by = 'ID')
row.names(eset) <- eset$ID

t_table <- eset
names(t_table) <- c("ID",'OvHK2','786vHK2','786vO')
write_csv(t_table, file = "results/RNA/t_table.csv")

eset <- eset[eset$ID != "HIF1A",] #This is a genetic deletion, shouldn't be used for TF activity estimation

eset <- eset[,-1]
names(eset) <- c('OvHK2','786vHK2','786vO')

TF_activity <- as.data.frame(
  viper(eset = eset, 
        regulon = regulon_viper, 
        nes = T, 
        minsize = 5, 
        eset.filter = F, 
        pleiotropy = T)) #minsize 20
TF_activity$ID <- row.names(TF_activity)

write_csv(TF_activity, file = "results/RNA/TF_activity.csv")


TF_activity <- TF_activity[abs(TF_activity$OvHK2) >= 2 | abs(TF_activity$`786vHK2`) >= 2,] #1.5

TF_activity_long <- melt(TF_activity)

TF_activity_long_HK2 <- TF_activity_long[grep('786vHK2',TF_activity_long$variable),]

temp <- TF_activity_long[grep('OvHK2',TF_activity_long$variable),]

TF_activity_long_HK2 <- as.data.frame(rbind(temp[order(temp$value, decreasing = F),], TF_activity_long_HK2[order(TF_activity_long_HK2$value, decreasing = F),]))
TF_activity_long_HK2$ID <- factor(TF_activity_long_HK2$ID, levels = unique(TF_activity_long_HK2$ID))
TF_activity_long_HK2 <- TF_activity_long_HK2[c(126,c(1:125),c(127:length(TF_activity_long_HK2$ID))),]
TF_activity_long_HK2$variable <- factor(TF_activity_long_HK2$variable, levels = unique(TF_activity_long_HK2$variable))

ggplot(TF_activity_long_HK2, aes(x = value, y = ID, group = variable, color = variable, size = abs(value))) + geom_point() +
  theme_minimal() + xlab('NES (activity score)') + ylab('TF')

TF_activity_long_M1AvO <- TF_activity_long[grep('786vO',TF_activity_long$variable),]
TF_activity_long_M1AvO <- TF_activity_long_M1AvO[order(temp$value, decreasing = F),]

TF_activity_long_M1AvO$ID <- factor(TF_activity_long_M1AvO$ID, levels = TF_activity_long_M1AvO$ID)

ggplot(TF_activity_long_M1AvO, aes(x = value, y = ID, group = variable, size = abs(value))) + geom_point(color = 'blue') +
  theme_minimal() + xlab('NES (activity score)') + ylab('TF')

##### pathways

library(fgsea)

pathways <- import_gmt("support/c2.cp.v7.0.symbols.gmt")
pathways_list <- list()
i <- 1
for(pathway in unique(pathways$term))
{
  pathways_list[[i]] <- pathways[pathways$term == pathway,"gene"]
  i <- i+1
}
names(pathways_list) <- unique(pathways$term)

stats_vec <- ttop_OvHK2$t 
names(stats_vec) <- ttop_OvHK2$ID
fgsea_res <- fgsea(pathways = pathways_list, stats = stats_vec, nperm = 1000, nproc = 3)
fgsea_res <- fgsea_res[,-8]

write_csv(fgsea_res, "results/RNA/fgsea_cp_OvHK2.csv")

stats_vec <- ttop_786vHK2$t 
names(stats_vec) <- ttop_786vHK2$ID
fgsea_res <- fgsea(pathways = pathways_list, stats = stats_vec, nperm = 1000, nproc = 3)
fgsea_res <- fgsea_res[,-8]

write_csv(fgsea_res, "results/RNA/fgsea_cp_M1AvHK2.csv")

stats_vec <- ttop_786vO$t 
names(stats_vec) <- ttop_786vO$ID
fgsea_res <- fgsea(pathways = pathways_list, stats = stats_vec, nperm = 1000, nproc = 3)
fgsea_res <- fgsea_res[,-8]

write_csv(fgsea_res, "results/RNA/fgsea_cp_M1AvO.csv")
