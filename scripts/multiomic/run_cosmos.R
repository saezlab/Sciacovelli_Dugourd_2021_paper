library(cosmosR)
library(readr)
library(dplyr)

data("meta_network")

#probably erroneous interaction
meta_network <- meta_network[-which(meta_network$source == "PRKCA" & meta_network$target == "SRC"),]

#Likely mistake in STICH database
meta_network <- meta_network[-which(meta_network$source == "Metab__HMDB0000517_c" & meta_network$target == "AKT1"),]

TF_activity <- as.data.frame(
  read_csv("results/RNA/TF_activity.csv"))

t_table <- as.data.frame(
  read_csv("results/RNA/t_table.csv"))

sig_input <- 1
names(sig_input) <- "Metab__HMDB0003378_c"
#HMDB0003378 NO
#HMDB0000052 argininosucc

metab_input <- TF_activity[,3]
names(metab_input) <- TF_activity$ID

RNA_input <- t_table[,4]
names(RNA_input) <- t_table$ID

##Filter significant inputs
metab_input <- metab_input[abs(metab_input) > 3]

#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options(solver = "cplex")

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
# my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable
my_options$solverPath <- "./cplex"
# my_options$solver <- "cplex" #or cbc
my_options$solver <- "cplex"
my_options$timelimit <- 3600*0.5
my_options$mipGAP <- 0.05
my_options$threads <- 6

metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network,
                                                      signaling_data = sig_input,
                                                      metabolic_data = metab_input,
                                                      diff_expression_data = RNA_input,
                                                      maximum_network_depth = 3,
                                                      remove_unexpressed_nodes = T,
                                                      filter_tf_gene_interaction_by_optimization = T,
                                                      CARNIVAL_options = my_options)

##https://pubs.acs.org/doi/abs/10.1021/jm00109a032
##Synthesis of NO downstream of arginine -> go angiogenesis and invasion
##https://www.frontiersin.org/articles/10.3389/fcell.2021.658861/full

my_options$timelimit <- 3600*5

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)

formatted_res <- format_COSMOS_res(test_result_for)

SIF <- formatted_res[[1]]
ATT <- formatted_res[[2]]

SIF <- SIF[which(SIF$Weight != 0),]

RNA_input_df <- data.frame(Nodes = names(RNA_input), t = RNA_input)
ATT <- merge(ATT, RNA_input_df, all.x = T)
ATT <- ATT[ATT$AvgAct != 0,]

write_csv(SIF, file = "results/cosmos/ASS_SIF.csv")
write_csv(ATT, file = "results/cosmos/ASS_ATT.csv")
