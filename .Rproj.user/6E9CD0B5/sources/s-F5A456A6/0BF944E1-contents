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

library(igraph)

source("scripts/support_functions.R")

labelling  <- as.data.frame(read_excel("data/metabolomic/Re-run HK2_ccRCC metabolomics +BCAT inhibitor.xlsx"))


labelling <- labelling[,-1]
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

batches <- as.data.frame(rbind(batches,batches["kmv+kic",]))
row.names(batches)[length(batches[,1])] <- "kmv"
batches <- as.data.frame(rbind(batches,batches["kmv+kic",]))
row.names(batches)[length(batches[,1])] <- "kic"
batches <- batches[-which(row.names(batches) == "kmv+kic"),]

# magicPlotMaker(batches, "~/Documents/Marco/Marco/results/metabolomic/labelling/l1/log2",targets)

unique(targets$condition)
comparisons <- list("HK2_BCATi" = c(2,-1), 
                    "O_BCATi" = c(4,-3), 
                    "M1A_BCATi" = c(6,-5),
                    "M2A_BCATi" = c(8,-7))

limmaRes <- runLimma(batches, targets, comparisons = comparisons)

t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(batches[,1]),
                         adjust.method = "fdr"))

t_table[which(t_table$ID == "succinyladenosine"),1] <-"adenylosuccinate"

t_table_names <- t_table

t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
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

##Compute metabolic enzme enrichment score
metactivity_res <- metactivity(metabolomic_t_table = t_table, 
                               regulons_df = regulons_df, 
                               compartment_pattern = "_[a-z]$", 
                               k = 10000)

mean_ES_df <- metactivity_res$ES

mean_NES_df <- metactivity_res$NES


##translate the metabolic ids back to names
translated_results <- translate_results(regulons_df = regulons_df, t_table = t_table, mapping_table = mapping_table)

##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = 'BCAT1>780', 
                                    stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_results$regulons_df, 
                                    contrast_index = 2, stat_name = 't', scaling_factor = 3, nLabels =  20)

plot(plots$scatter)
plot(plots$cumsumPlot)

##Visualise results at pathway level pathways

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 1, vis.height = 2000) %>% 
  visSave(file = "results/networks/BCATi/HK2_BCATi_full_NES_p8.html")
plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 2, vis.height = 2000) %>% 
  visSave(file = "results/networks/BCATi/O_BCATi_full_NES_p8.html")
plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 3, vis.height = 2000) %>% 
  visSave(file = "results/networks/BCATi/M1A_BCATi_full_NES_p8.html")
plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 4, vis.height = 2000) %>% 
  visSave(file = "results/networks/BCATi/M1A_BCATi_full_NES_p8.html")

ignet <- graph_from_data_frame(sub_network_nocofact$reaction_network)
shortest_paths(ignet, from = "BCAT1",to = "cpd:C00037_m",)