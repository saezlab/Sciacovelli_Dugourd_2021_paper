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
library(decoupleR)
library(cowplot)
library(gridExtra)

library(igraph)




source("scripts/support_functions.R")


PLASMAX <- as.data.frame(
  read_excel("data/metabolomic/PLASMAX.xlsx"))

batches <- as.data.frame(t(PLASMAX[,-c(1,2)]))

names(batches) <- paste(PLASMAX$Sample, 1:length(names(batches)), sep = "_")

targets <- as.data.frame(cbind(names(batches),NA))
names(targets) <- c("sample","condition")
targets$condition <- gsub("_.*","",targets$sample)

batches <- batches[rowSums(batches) != 0,]
batches[batches == 0] <- NA 
batches <- as.data.frame(log2(batches))

batches <- as.data.frame(rbind(batches,batches["methylmalonylcarnitine C3-DC-M/succinylcarnitine C4-DC",]))
row.names(batches)[length(batches[,1])] <- "methylmalonylcarnitine C3-DC-M"
# batches <- as.data.frame(rbind(batches,batches["methylmalonylcarnitine C3-DC-M/succinylcarnitine C4-DC",]))
# row.names(batches)[length(batches[,1])] <- "succinylcarnitine C4-DC"
batches <- batches[-which(row.names(batches) == "methylmalonylcarnitine C3-DC-M/succinylcarnitine C4-DC"),]

sub_batches <- batches[row.names(batches) %in% c("leucine",
                                         "isoleucine",
                                         "KMV",
                                         "isovalerylcarnitine C5",
                                         "propionylcarnitine C3",
                                         "methylmalonylcarnitine C3-DC-M"),]

sub_batches <- 2^sub_batches

to_plot <- melt(batches)
to_plot$condition <- gsub("_.*","",to_plot$variable)

ggplot(to_plot, aes(x = variable, y = value, group = variable, fill = condition)) + geom_violin()

plot(nicePCA(batches, targets))

unique(targets$condition)

comparisons <- list("OvHK2" = c(2,-1), 
                    "M1AvHK2" = c(3,-1), 
                    "M2AvHK2" = c(4,-1),
                    "M1AvO" = c(3,-2))

limmaRes <- runLimma(batches, targets, comparisons = comparisons)

ttop_O <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(batches[,1]), adjust.method = "fdr"))
ttop_M1A <- ttopFormatter(topTable(limmaRes[[1]], coef = 2, number = length(batches[,1]), adjust.method = "fdr"))

t_table <- ttop_list_to_t_table(
  limma_res_to_ttop_list(limma_res = limmaRes,
                         comp_names = names(comparisons),
                         number = length(batches[,1]),
                         adjust.method = "fdr"))

sub_t_table <- t_table[t_table$ID %in% c("leucine",
                                         "isoleucine",
                                         "KMV",
                                         "isovalerylcarnitine C5",
                                         "propionylcarnitine C3",
                                         "methylmalonylcarnitine C3-DC-M"),]

to_dotplot <- t_table
to_dotplot <- to_dotplot[order(to_dotplot$OvHK2, decreasing = F),]
to_dotplot_long <- melt(to_dotplot[,c(1,2,3)])
to_dotplot_long <- to_dotplot_long[c(180:358,1:179),]

to_dotplot_long$ID <- factor(to_dotplot_long$ID, levels = unique(to_dotplot_long$ID))
to_dotplot_long$variable <- factor(to_dotplot_long$variable, levels = unique(to_dotplot_long$variable))
to_dotplot_long <- to_dotplot_long[abs(to_dotplot_long$value) > 11,]

ggplot(to_dotplot_long, aes(x = value, y = ID, size = abs(value))) +
  geom_point(aes(colour = variable)) +
  # scale_color_gradient2(low="blue", high="blue", midpoint = 0, mid = "blue") +
  theme_minimal() + ggtitle("786-M1A and 786-O versus HK2 metabolomic") #+
  # geom_vline(xintercept = -2) +
  # geom_vline(xintercept = 2)

to_dotplot <- to_dotplot[order(to_dotplot$M1AvO),]
to_dotplot$ID <- factor(to_dotplot$ID, levels = to_dotplot$ID)
to_dotplot <- to_dotplot[abs(to_dotplot$M1AvO) > 3.5,]

ggplot(to_dotplot, aes(x = M1AvO, y = ID, size = abs(M1AvO))) +
  geom_point(aes(colour = M1AvO)) +
  # scale_color_gradient2(low="blue", high="red", midpoint = 0, mid = "white") +
  scale_color_gradient2(low="blue", high="blue", midpoint = 0, mid = "blue") +
  theme_minimal() + ggtitle("786-M1A versus 786-O metabolomic") #+
  # geom_vline(xintercept = -2) +
  # geom_vline(xintercept = 2)

t_table$ID <- tolower(gsub("[, ]","",t_table$ID))
sum(t_table$ID %in% mapping_table$metab)
t_table[!(t_table$ID %in% mapping_table$metab),"ID"]

# mapping_table <- as.data.frame(rbind(mapping_table, c("methylmalonylcarnitinec3-dc-m","methylmalonylcarnitine-C3_m")))

t_table <- t_table_metactivity_input_formater(metabolomic_t_table = t_table,
                                              mapping_table = mapping_table,
                                              affixes = c("c","l","x","m","e","n","r"))

write_csv(t_table, file = "results/metabolomic/plasmax/t_table.csv")

######## SUBNETWORKS

expressed_genes <- as.data.frame(read_csv("data/metabolomic/expressed_genes_plasmax.csv"))
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

sub_network_nocofact <- split_transaminases(sub_network_nocofact = sub_network_nocofact)

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
plots <- plotMetaboliteContribution(enzyme = 'ASNS_glugln', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 1, nLabels =  30)

plot(plots$scatter)
# plot(plots$cumsumPlot)

##Visualise results for single enzmes
plots <- plotMetaboliteContribution(enzyme = 'GRHPR', stat_df = translated_results$t_table, 
                                    metabolite_sets = translated_regulons_df, 
                                    contrast_index = 1, stat_name = 'Abundance Down <==> Up (t-value)', 
                                    scaling_factor = 1, nLabels =  30)

plot(plots$scatter)

ig_net <- graph_from_data_frame(sub_network_nocofact$reaction_network)

shortest_paths(ig_net, from = "ADSL>2", to = "cpd:C00152_c", mode = "out")
shortest_paths(ig_net, from = "BCAT1>780", to = "cpd:methylmalonylcarnitine-C3_m", mode = "out")

shortest_paths(ig_net, from = "BCAT1>780_gluakg", to = "cpd:C03406_c", mode = "out")
shortest_paths(ig_net, from = "BCAT1>780", to = "cpd:C03406_c", mode = "out")

####################

#Visualise results at pathway level pathways
hm <- pathway_HM(mean_NES_df = mean_NES_df, pathway_name = 'KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION', pathways = kegg_pathways)
plot(hm)

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 1, vis.height = 2000) %>%
  visSave(file = "results/networks/PLASMAX/786OvsHK2_full_NES_p8.html")

network <- generate_sifAtt(sub_network_nocofact, t_table, mean_NES_df, column_index = 1)
write_csv(network$sif, file = "results/networks/PLASMAX/786OvsHK2_full_NES_p8_sif.csv")
write_csv(network$att[,-7], file = "results/networks/PLASMAX/786OvsHK2_full_NES_p8_att.csv")

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 2, vis.height = 2000) %>%
  visSave(file = "results/networks/PLASMAX/786M1AvsHK2_full_NES_p8.html")

network <- generate_sifAtt(sub_network_nocofact, t_table, mean_NES_df, column_index = 2)
write_csv(network$sif, file = "results/networks/PLASMAX/786M1AvsHK2_full_NES_p8_sif.csv")
write_csv(network$att[,-7], file = "results/networks/PLASMAX/786M1AvsHK2_full_NES_p8_att.csv")

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 3, vis.height = 2000) %>%
  visSave(file = "results/networks/PLASMAX/786M2AvsHK2_full_NES_p8.html")

network <- generate_sifAtt(sub_network_nocofact, t_table, mean_NES_df, column_index = 3)
write_csv(network$sif, file = "results/networks/PLASMAX/786M2AvsHK2_full_NES_p8_sif.csv")
write_csv(network$att[,-7], file = "results/networks/PLASMAX/786M2AvsHK2_full_NES_p8_att.csv")

plot_reaction_network(sub_network_nocofact, t_table, mean_NES_df, column_index = 4, vis.height = 2000) %>%
  visSave(file = "results/networks/PLASMAX/786M1Avs786O_full_NES_p8.html")

network <- generate_sifAtt(sub_network_nocofact, t_table, mean_NES_df, column_index = 4)
write_csv(network$sif, file = "results/networks/PLASMAX/786M1Avs786O_full_NES_p8_sif.csv")
write_csv(network$att[,-7], file = "results/networks/PLASMAX/786M1Avs786O_full_NES_p8_att.csv")