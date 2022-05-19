rm(list = ls(all.names = TRUE))
gc()
.rs.restartR() ##Only if Rstudio

library(ocean)
library(readxl)
library(reshape2)
library(ggplot2)
library(readr)
library(limma)
library(visNetwork)
library(dplyr)
library(pheatmap)

library(igraph)

source("scripts/support_functions.R")

labelling  <- as.data.frame(read_excel("data/metabolomic/GLS_inhibitor.xlsx"))

names(labelling) <- tolower(names(labelling))
names(labelling) <- gsub(",",'-',names(labelling))

metabs <- names(labelling)[-1]
labelling$sample <- gsub(" ", "",labelling$sample)
samples <- labelling$sample

batches <- as.data.frame(t(labelling[,-1]))
names(batches) <- samples
names(batches) <- paste(names(batches), c(1:length(batches)), sep = "_")

targets <- as.data.frame(matrix(NA,length(batches[1,]),2))
names(targets) <- c("sample","condition")

targets$sample <- names(batches)
targets$condition <- gsub("_[0-9]+$","",names(batches))

non_zero_min <- sort(as.matrix(batches))[sort(as.matrix(batches)) != 0][1]

# batches[batches == 0] <- non_zero_min/2
batches[batches == 0] <- NA
batches <- log2(batches)

unique(targets$condition)
comparisons <- list("HK2_CB839_100nM" = c(2,-1), 
                    "HK2_CB839_500nM" = c(3,-1), 
                    "HK2_CB839_1μM" = c(4,-1),
                    "786O_CB839_100nM" = c(6,-5),
                    "786O_CB839_500nM" = c(7,-5),
                    "786O_CB839_1μM" = c(8,-5),
                    "786M1A_CB839_100nM" = c(10,-9),
                    "786M1A_CB839_500nM" = c(11,-9),
                    "786M1A_CB839_1μM" = c(12,-9),
                    "HK2_BMS_5μM" = c(13,-1),
                    "HK2_BMS_10μM" = c(14,-1),
                    "786O_BMS_5μM" = c(15,-5),
                    "786O_BMS_10μM" = c(16,-5),
                    "786M1A_BMS_5μM" = c(17,-9),
                    "786M1A_BMS_10μM" = c(18,-9))

limmaRes <- runLimma(batches, targets, comparisons = comparisons)

t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(batches[,1]),
                         adjust.method = "fdr"))

t_table_names <- t_table

GLS_inhib_metab_to_kegg <- as.data.frame(
  read_csv("support/GLS_inhib_metab_to_kegg.txt"))

t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = GLS_inhib_metab_to_kegg,
                                              affixes = c("c","l","x","m","e","n","r"))

##Prepare the metabolic enzyme sets
penalty_min <- 8 #minimum 1 and integer
penalty_max <- 8 #maximum 9 and integer

######## SUBNETWORKS

expressed_genes <- as.data.frame(read_csv("data/metabolomic/expressed_genes_plasmax.csv"))

# View(unique(recon2_redhuman$pathway))

all_pathways <- unique(recon2_redhuman$pathway)
sub_network <- model_to_pathway_sif(pathway_to_keep = all_pathways$X1)

sub_network <- translate_complexes(sub_network)

sub_network_nocofact <- remove_cofactors(sub_network)

non_expressed_genes <- expressed_genes[rowSums(expressed_genes[,c(1,2,3)]) == 0,c(4,5)] #HK2 and O
non_expressed_genes <- c(non_expressed_genes$hgnc_symbol,non_expressed_genes$entrezgene_id)

tokeep <- list()
for(i in 1:length(sub_network_nocofact$reaction_network[,1]))
{
  elements <- gsub("[><].*","",sub_network_nocofact$reaction_network[i,])
  elements <- elements[!grepl("cpd:",elements)]
  elements <- unlist(strsplit(elements, "_"))
  tokeep[[i]] <- sum(elements %in% non_expressed_genes) == 0
}
tokeep <- unlist(tokeep)
sub_network_nocofact$reaction_network <- sub_network_nocofact$reaction_network[tokeep,]

sub_network_nocofact <- compress_transporters(sub_network_nocofact = sub_network_nocofact)

sub_network_nocofact <- split_transaminases(sub_network_nocofact = sub_network_nocofact)

# plot_network(sub_network_nocofact$reaction_network)

enzymes <- unique(sub_network_nocofact$attributes$V1)
enzymes <- enzymes[!grepl("_[clxmenr]$",enzymes)]


sub_forest <- forestMaker(enzymes, sub_network_nocofact$reaction_network, branch_length = c(3,3))

###################

reaction_set_list <- prepare_metabolite_set(penalty_range = penalty_min:penalty_max,   
                                            # forest = tree_without_cofactors,
                                            forest = sub_forest,
                                            measured_metabolites = t_table$KEGG)

reaction_set_list_merged <- condense_metabolite_set(reaction_set_list = reaction_set_list)

penalty <- 8 #has be between penalty_min and penalty_max and integer

regulons_df <- prepare_regulon_df(reaction_set_list_merged, penalty, filter_imbalance = c(0,1))

t_table <- t_table[!grepl("cpd:C00438",t_table$KEGG),]
##Compute metabolic enzme enrichment score
metactivity_res <- metactivity(metabolomic_t_table = t_table, 
                               regulons_df = regulons_df, 
                               compartment_pattern = "_[a-z]$", 
                               k = 10000)

mean_ES_df <- metactivity_res$ES

mean_NES_df <- metactivity_res$NES


##translate the metabolic ids back to names
translated_results <- translate_results(regulons_df = regulons_df, t_table = t_table, mapping_table = mapping_table)

translated_regulons_df <- translated_results$regulons_df
translated_regulons_df$ID <- paste(translated_regulons_df$set, gsub("_[a-z]$","",translated_regulons_df$targets), sep = "___")
translated_regulons_df <- translated_regulons_df[,-c(1,2)]

translated_regulons_df <- translated_regulons_df %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
translated_regulons_df <- as.data.frame(translated_regulons_df)

translated_regulons_df$set <- gsub("___.*","",translated_regulons_df$ID)
translated_regulons_df$targets <- gsub(".*___","",translated_regulons_df$ID)
translated_regulons_df <- translated_regulons_df[,c(3,4,2)]

##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = 'GLS_glugln', 
                                    stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 9, stat_name = 't', scaling_factor = 1, nLabels =  20)

plot(plots$scatter)
plot(plots$cumsumPlot)

##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = 'AACS', 
                                    stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 9, stat_name = 't', scaling_factor = 1, nLabels =  20)

plot(plots$scatter)
plot(plots$cumsumPlot)

##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = 'ACLY', 
                                    stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 9, stat_name = 't', scaling_factor = 1, nLabels =  20)

plot(plots$scatter)
plot(plots$cumsumPlot)


##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = 'NAGS', 
                                    stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 9, stat_name = 't', scaling_factor = 1, nLabels =  20)

plot(plots$scatter)
plot(plots$cumsumPlot)
############

to_hm <- t_table_names[,c(4,7,10)]
row.names(to_hm) <- t_table_names[,1]

pheatmap(to_hm)
##Visualise results at pathway level pathways

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 6, vis.height = 2000) %>% 
  visSave(file = "results/networks/GLSi/O_GLSi_full_NES_p8.html")
plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 9, vis.height = 2000) %>% 
  visSave(file = "results/networks/GLSi/M1A_GLSi_full_NES_p8.html")

ignet <- graph_from_data_frame(sub_network_nocofact$reaction_network)
shortest_paths(ignet, from = "ACLY",to = "cpd:C00049_c")
