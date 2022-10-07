rm(list = ls(all.names = TRUE))
gc()
.rs.restartR() ##Only if Rstudio

library(readr)
library(dplyr)
library(vsn)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(limma)
library(fgsea)
library(GSEABase)
library(pheatmap)

runfgsea_ttop <- function(ttop, 
                          pathways,
                          nperm, 
                          minSize = 1, 
                          maxSize = Inf, 
                          nproc = 0,
                          gseaParam = 1, 
                          BPPARAM = NULL)
{
  stats_vec <- ttop$t 
  names(stats_vec) <- ttop$ID
  fgsea_res <- fgsea(pathways = pathways, stats = stats_vec, nperm = nperm, nproc = nproc, minSize = minSize,
                     maxSize = maxSize, gseaParam = gseaParam, BPPARAM = BPPARAM)
  
  return(fgsea_res)
}

source('scripts/support_functions.R')

proteomic <- as.data.frame(
  read_delim("data/proteomic/proteomic_VHL.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE)) #formerly named proteinGroups raw.txt

proteomic$unique_gene_name <- gsub(";.*","",proteomic$`Gene names`)
proteomic <- proteomic[,c(44,1:24)]

proteomic <- proteomic %>% group_by(unique_gene_name) %>% summarise_each(funs(sum(., na.rm = TRUE)))
proteomic <- as.data.frame(proteomic)
proteomic <- proteomic[complete.cases(proteomic),]

row.names(proteomic) <- proteomic[,1]
proteomic <- proteomic[,-1]
proteomic[proteomic == 0] <- NA

targets <- as.data.frame(matrix(NA,dim(proteomic)[2],2))
names(targets) <- c("sample","condition")
targets$sample <- names(proteomic)
targets$condition <- c(rep("HK2",5),rep("786_O",5),rep("786_O+VHL",5),rep("786_M1A",5),rep("786_M1A+VHL",4))

plots <- magicPlotMakerLight(log2(proteomic),targets)

plot(plots[[1]])
plot(plots[[2]])

#now we can normalise the cleaned dataframe using vsn
fit <- vsnMatrix(as.matrix(proteomic)) #train vsn parameters

#make sure the mean/sd trend is not going crazy
meanSdPlot(fit)

#if good, normalise data with the trained parameters of vsn
proteomic_vsn <- as.data.frame(vsn::predict(fit,as.matrix(proteomic)))

plots_vsn <- magicPlotMakerLight(proteomic_vsn,targets)

plot(plots_vsn[[1]])
plot(plots_vsn[[2]])

#first check the conditions order
unique(targets$condition)

#we want to compare the KO condition with the WT condition so we build a
#comparison list
comparisons <- list("786_OvsHK2" = c(2,-1),
                    "786_O+VHLvsHK2" = c(3,-1),
                    "786_M1AvsHK2" = c(4,-1),
                    "786_M1A+VHLvsHK2" = c(5,-1),
                    "786_O+VHLvs786_O" = c(3,-2),
                    "786_M1A+VHLvs786_M1A" = c(5,-4)) 

#now that the comparisons are defined, we can run limma
limmaRes <- runLimma(measurements = proteomic_vsn, 
                     targets = targets, 
                     comparisons = comparisons)

#once limma has run, we extract the statistics dataframe to summarise the
#differential analysis
t_table <- merge(ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(proteomic_vsn[,1]), adjust.method = "fdr"))[,c(1,4)],
                 ttopFormatter(topTable(limmaRes[[1]], coef = 2, number = length(proteomic_vsn[,1]), adjust.method = "fdr"))[,c(1,4)],
                 by = 'ID')

for( i in 3:6)
{
  t_table <- merge(t_table,
                   ttopFormatter(topTable(limmaRes[[1]], coef = i, number = length(proteomic_vsn[,1]), adjust.method = "fdr"))[,c(1,4)],
                   by = 'ID')
}

names(t_table) <- c("ID",names(comparisons))

ttop_list <- sapply(1:6,function(x,limmaRes,proteomic_vsn){
  return(ttopFormatter(topTable(limmaRes[[1]], coef = x, number = length(proteomic_vsn[,1]), adjust.method = "fdr")))
}, limmaRes = limmaRes, proteomic_vsn = proteomic_vsn, USE.NAMES = T, simplify = F)
names(ttop_list) <- names(comparisons)

####### FGSEA

pathways <- import_gmt("~/Dropbox/GDC_KIRK/c2.cp.v7.0.symbols.gmt")
pathways <- pathways[grep("KEGG",pathways$term),]
pathways_list <- list()
i <- 1
for(pathway in unique(pathways$term))
{
  pathways_list[[i]] <- pathways[pathways$term == pathway,"gene"]
  i <- i+1
}
names(pathways_list) <- unique(pathways$term)

fgsea_res_list <- sapply(ttop_list, function(x, pathways_list){
  return(as.data.frame(runfgsea_ttop(x, pathways_list, 10000)[,-8]))
}, pathways_list = pathways_list, USE.NAMES = T, simplify = F)

pathway_table <- merge(fgsea_res_list[[1]][,c(1,5)],fgsea_res_list[[2]][,c(1,5)], by = "pathway")

for( i in 3:6)
{
  pathway_table <- merge(pathway_table,
                         fgsea_res_list[[i]][,c(1,5)],
                         by = 'pathway')
}

names(pathway_table) <- c("pathway",names(comparisons))

pathway_table_top_dn <- pathway_table[order(pathway_table$`786_OvsHK2`),]
pathway_table_top_dn <- pathway_table_top_dn[1:25,]

pathway_table_top_up <- pathway_table[order(pathway_table$`786_OvsHK2`, decreasing = T),]
pathway_table_top_up <- pathway_table_top_up[1:25,]

pathway_table_top <- as.data.frame(rbind(pathway_table_top_dn, pathway_table_top_up))
row.names(pathway_table_top) <- pathway_table_top$pathway
pathway_table_top <- pathway_table_top[,-1]

pheatmap(pathway_table_top, cluster_cols = F)

#######
t_table_matrix <- t_table
row.names(t_table_matrix) <- t_table_matrix[,1]
t_table_matrix <- t_table_matrix[,-1]

to_gg <- as.data.frame(t(t_table_matrix['ASS1',]))
to_gg$condition <- row.names(to_gg)
names(to_gg)[1] <- "stat"
to_gg$condition <- factor(to_gg$condition, levels = unique(to_gg$condition))

ggplot(to_gg, aes(x = condition, y = stat, fill = condition)) + 
  geom_bar(stat = 'identity') +
  theme_minimal() + ggtitle('ASS1 expression change')

volcano_nice(ttop_list$`786_O+VHLvsHK2`, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = c("ASS1","VHL")) + ggtitle("786+VHLvsHK2") + ylab("-log10(p-value)")
volcano_nice(ttop_list$`786_O+VHLvs786_O`, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = c("ASS1","VHL")) + ggtitle("786+VHLvs786-O") + ylab("-log10(p-value)")
volcano_nice(ttop_list$`786_M1A+VHLvsHK2`, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = c("ASS1","VHL")) + ggtitle("M1A+VHLvsHK2") + ylab("-log10(p-value)")
volcano_nice(ttop_list$`786_M1A+VHLvs786_M1A`, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = c("ASS1","VHL")) + ggtitle("M1A+VHLvsM1A-O") + ylab("-log10(p-value)")

ttop_786VHLvs786O_BCCA <- ttop_list$`786_O+VHLvs786_O`
ttop_786VHLvs786O_BCCA <- ttop_786VHLvs786O_BCCA[ttop_786VHLvs786O_BCCA$ID %in% pathways[pathways$term == "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","gene"],]

ttop_M1AVHLvsM1AO_BCCA <- ttop_list$`786_M1A+VHLvs786_M1A`
ttop_M1AVHLvsM1AO_BCCA <- ttop_M1AVHLvsM1AO_BCCA[ttop_M1AVHLvsM1AO_BCCA$ID %in% pathways[pathways$term == "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","gene"],]

volcano_nice(ttop_786VHLvs786O_BCCA, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = "ASS1") + ggtitle("BCCA_786+VHLvs786-O") + ylab("-log10(p-value)")
volcano_nice(ttop_M1AVHLvsM1AO_BCCA, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = "ASS1") + ggtitle("BCCA_M1A+VHLvsM1A-O") + ylab("-log10(p-value)")

###############

ttop_786VHLvsHK2_BCCA <- ttop_list$`786_O+VHLvsHK2`
ttop_786VHLvsHK2_BCCA <- ttop_786VHLvsHK2_BCCA[ttop_786VHLvsHK2_BCCA$ID %in% pathways[pathways$term == "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","gene"],]

ttop_M1AVHLvsHK2_BCCA <- ttop_list$`786_M1A+VHLvsHK2`
ttop_M1AVHLvsHK2_BCCA <- ttop_M1AVHLvsHK2_BCCA[ttop_M1AVHLvsHK2_BCCA$ID %in% pathways[pathways$term == "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","gene"],]

volcano_nice(ttop_786VHLvsHK2_BCCA, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = "ASS1") + ggtitle("BCCA_786+VHLvsHK2") + ylab("-log10(p-value)")
volcano_nice(ttop_M1AVHLvsHK2_BCCA, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = "ASS1") + ggtitle("BCCA_M1A+VHLvsHK2") + ylab("-log10(p-value)")

###############

ttop_786VHLvsHK2_LIMONENE <- ttop_list$`786_O+VHLvsHK2`
ttop_786VHLvsHK2_LIMONENE <- ttop_786VHLvsHK2_LIMONENE[ttop_786VHLvsHK2_LIMONENE$ID %in% pathways[pathways$term == "KEGG_LIMONENE_AND_PINENE_DEGRADATION","gene"],]

ttop_M1AVHLvsHK2_LIMONENE <- ttop_list$`786_M1A+VHLvsHK2`
ttop_M1AVHLvsHK2_LIMONENE <- ttop_M1AVHLvsHK2_LIMONENE[ttop_M1AVHLvsHK2_LIMONENE$ID %in% pathways[pathways$term == "KEGG_LIMONENE_AND_PINENE_DEGRADATION","gene"],]

volcano_nice(ttop_786VHLvsHK2_LIMONENE, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = "ASS1") + ggtitle("LIMONENE_786+VHLvs786-O") + ylab("-log10(p-value)")
volcano_nice(ttop_M1AVHLvsHK2_LIMONENE, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 50, label = T, manual_labels = "ASS1") + ggtitle("LIMONENE_M1A+VHLvsM1A-O") + ylab("-log10(p-value)")
