rm(list = ls(all.names = TRUE))
gc()
.rs.restartR() ##Only if Rstudio

library(ocean)
library(readxl)
library(reshape2)
library(ggplot2)
library(readr)
library(dplyr)
library(visNetwork)

library(igraph)

# source("script/support_pheatmap_colors.R")
source("scripts/support_functions.R")
source("scripts/support_labdiff_functions.R")

labelling  <- as.data.frame(read_excel("data/metabolomic/15N_leu_iso_combined.xlsx", 
                                       sheet = "Sheet1"))

labelling <- labelling[,-1]
names(labelling) <- tolower(names(labelling))

metabs <- names(labelling)[-1]
labelling$sample <- gsub(" ", "",labelling$sample)
samples <- labelling$sample
samples <- paste(samples,1:length(samples),sep = "_")

batches <- as.data.frame(t(labelling[,-1]))
names(batches) <- samples
# names(batches) <- paste(names(batches), c(1:length(batches)), sep = "_")

batches$metab <- row.names(batches)
batches_label <- batches[grepl("_n1",row.names(batches)),]
batches_label$metab <- gsub("_n1","",row.names(batches_label))

batches_label <- batches_label %>% group_by(metab) %>% summarise_each(funs(sum(., na.rm = TRUE)))
batches_label <- as.data.frame(batches_label)
row.names(batches_label) <- batches_label$metab
batches_label <- batches_label[,-1]

batches <- batches[-which(row.names(batches) == "n-acetylaspartate"),]
batches <- batches[!grepl("_n1",row.names(batches)),]
batches <- batches[,-length(batches[1,])]
batches <- batches[row.names(batches_label),]

batches <- batches_label / (batches + batches_label)
# batches <- batches - batches_label
batches[is.na(batches)] <- 0
batches <- batches[rowSums(batches) != 0,]
batches <- as.data.frame(batches)

targets <- as.data.frame(matrix(NA,length(batches[1,]),2))
names(targets) <- c("sample","condition")

targets$sample <- names(batches)
targets$condition <- gsub("_[0-9]+$","",names(batches))

# non_zero_min <- sort(as.matrix(batches))[sort(as.matrix(batches)) != 0][1]

# batches[batches == 0] <- non_zero_min/2
# batches[batches == 0] <- NA
# batches <- log2(batches)

row.names(batches) <- gsub(" ","",row.names(batches))

sum(row.names(batches) %in% mapping_table$metab)
row.names(batches)[!(row.names(batches) %in% mapping_table$metab)]

batches <- as.data.frame(rbind(batches,batches["leucine+isoleucine",]))
row.names(batches)[length(batches[,1])] <- "leucine"
batches <- as.data.frame(rbind(batches,batches["leucine+isoleucine",]))
row.names(batches)[length(batches[,1])] <- "isoleucine"
batches <- batches[-which(row.names(batches) == "leucine+isoleucine"),]

####
#

# batches[batches == 0] <- NA
# rowMins_div2 <- apply(batches,1,function(x){min(x, na.rm = T)}) / 2
# batches[is.na(batches)] <- 0
# batches <- log2(batches + rowMins_div2)

#
####

row.names(batches)[!(row.names(batches) %in% mapping_table$metab)]

unique(targets$condition)

avg_lab_prop <- as.data.frame(average_labelling_proportion(batches, targets))
avg_lab_prop[avg_lab_prop == 0] <- NA
rowMins_div2 <- apply(avg_lab_prop,1,function(x){min(x, na.rm = T)}) / 2
avg_lab_prop[is.na(avg_lab_prop)] <- 0
# avg_lab_prop <- log2(avg_lab_prop + rowMins_div2)

comparisons <- list("OvHK2" = c(2,-1), 
                    "M1AvHK2" = c(3,-1), 
                    "M2AvHK2" = c(4,-1),
                    "M1AvsO" = c(3,-2))

avg_lab_diff <- average_labelling_differences(avg_lab_prop, comparisons)

t_table <- as.data.frame(do.call(cbind,avg_lab_diff))
names(t_table) <- names(comparisons)

t_table$ID <- row.names(t_table)
t_table <- t_table[,c(length(t_table[1,]), 1:(length(t_table[1,])-1))]

#####
#

# limmaRes <- runLimma(batches, targets, comparisons = comparisons)
# 
# t_table <- ttop_list_to_t_table(
#   limma_res_to_ttop_list(limma_res = limmaRes,
#                          comp_names = names(comparisons),
#                          number = length(batches[,1]),
#                          adjust.method = "fdr"))

#
######

t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
                                              affixes = c("c","l","x","m","e","n","r"))

t_table <- t_table[complete.cases(t_table),]

# write_csv(t_table, file = "results/metabolomic/plasmax/t_table.csv")

######## SUBNETWORKS

expressed_genes <- as.data.frame(read_csv("support/expressed_genes.csv"))
# View(unique(recon2_redhuman$pathway))

all_pathways <- unique(recon2_redhuman$pathway)
sub_network <- model_to_pathway_sif(pathway_to_keep = all_pathways$X1)

sub_network <- translate_complexes(sub_network)

sub_network_nocofact <- remove_cofactors(sub_network)

non_expressed_genes <- expressed_genes[rowSums(expressed_genes[,c(1,2,3)]) == 0,c(4,5)] #HK2 and O
non_expressed_genes <- c(non_expressed_genes$hgnc_symbol,non_expressed_genes$entrezgene_id,"HMR3832") #Arginine to citruline should be catalised by nos, not expressed here

tokeep <- list()
for(i in 1:length(sub_network_nocofact$reaction_network[,1]))
{
  elements <- gsub("[><].*","",sub_network_nocofact$reaction_network[i,])
  elements <- elements[!grepl("cpd:",elements)]
  elements <- gsub("^HMR_","HMR",elements)
  elements <- unlist(strsplit(elements, "_"))
  tokeep[[i]] <- sum(elements %in% non_expressed_genes) == 0
}
tokeep <- unlist(tokeep)
sub_network_nocofact$reaction_network <- sub_network_nocofact$reaction_network[tokeep,]

sub_network_nocofact <- compress_transporters(sub_network_nocofact = sub_network_nocofact)

sub_network_nocofact <- nitrogen_tracking(sub_network_nocofact = sub_network_nocofact)

enzymes <- unique(sub_network_nocofact$attributes$V1)
enzymes <- enzymes[!grepl("_[clxmenr]$",enzymes)]

sub_forest <- forestMaker(enzymes, sub_network_nocofact$reaction_network, branch_length = c(1,1), remove_reverse = T)

###################
##Prepare the metabolic enzyme sets
penalty_min <- 6 #minimum 1 and integer
penalty_max <- 8 #maximum 9 and integer

reaction_set_list <- prepare_metabolite_set(penalty_range = penalty_min:penalty_max,   
                                            forest = sub_forest,
                                            measured_metabolites = t_table$KEGG)

reaction_set_list_merged <- condense_metabolite_set(reaction_set_list = reaction_set_list)

penalty <- 8 #has be between penalty_min and penalty_max and integer

regulons_df <- prepare_regulon_df(reaction_set_list_merged, penalty, filter_imbalance = c(0,1))

##Compute metabolic enzme enrichment score
metactivity_res <- metactivity(metabolomic_t_table = t_table, 
                               regulons_df = regulons_df, 
                               compartment_pattern = "_[a-z]$", 
                               k = 1000)

mean_ES_df <- metactivity_res$ES

mean_NES_df <- metactivity_res$NES

# write_csv(mean_NES_df, file = "results/metabolomic/plasmax/mean_nes_df_p8.csv")

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
plots <- plotMetaboliteContribution(enzyme = 'ASS1', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 1, nLabels =  30)

plot(plots$scatter)
# plot(plots$cumsumPlot)

##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = 'ACO1>1218', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 4, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 1, nLabels =  30)

plot(plots$scatter)

####################



plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 1, vis.height = 2000) %>%
  visSave(file = "results/networks/15N/786OvsHK2_full_NES_p8.html")

network <- generate_sifAtt(sub_network_nocofact, t_table, mean_NES_df, column_index = 1)
write_csv(network$sif, file = "results/networks/15N/786OvsHK2_full_NES_p8_sif.csv")
write_csv(network$att[,-7], file = "results/networks/15N/786OvsHK2_full_NES_p8_att.csv")

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 2, vis.height = 2000) %>%
  visSave(file = "results/networks/15N/786M1AvsHK2_full_NES_p8.html")

network <- generate_sifAtt(sub_network_nocofact, t_table, mean_NES_df, column_index = 2)
write_csv(network$sif, file = "results/networks/15N/786M1AvsHK2_full_NES_p8_sif.csv")
write_csv(network$att[,-7], file = "results/networks/15N/786M1AvsHK2_full_NES_p8_att.csv")

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 3, vis.height = 2000) %>%
  visSave(file = "results/networks/15N/786M2AvsHK2_full_NES_p8.html")

network <- generate_sifAtt(sub_network_nocofact, t_table, mean_NES_df, column_index = 3)
write_csv(network$sif, file = "results/networks/15N/786M2AvsHK2_full_NES_p8_sif.csv")
write_csv(network$att[,-7], file = "results/networks/15N/786M2AvsHK2_full_NES_p8_att.csv")

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 4, vis.height = 2000) %>%
  visSave(file = "results/networks/15N/786M1Avs786O_full_NES_p8.html")

network <- generate_sifAtt(sub_network_nocofact, t_table, mean_NES_df, column_index = 4)
write_csv(network$sif, file = "results/networks/15N/786M1Avs786O_full_NES_p8_sif.csv")
write_csv(network$att[,-7], file = "results/networks/15N/786M1Avs786O_full_NES_p8_att.csv")