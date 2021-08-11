rm(list = ls(all.names = TRUE))
gc()
.rs.restartR() ##Only if Rstudio

library(readr)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(cowplot)
library(grid)
library(gridExtra)
library(limma)
library(GSEABase)

source('scripts/support_functions.R')

proteomic_plasmax <- as.data.frame(
  read_delim("data/proteomic/proteomic_plasmax.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE))

proteomic_plasmax <- proteomic_plasmax[,c(1:20,59)]
proteomic_plasmax <- proteomic_plasmax[!is.na(proteomic_plasmax$`Gene names`),]

row.names(proteomic_plasmax) <- proteomic_plasmax$`Gene names`
proteomic_plasmax <- proteomic_plasmax[,c(21,1:20)]
names(proteomic_plasmax)[1] <- "ID" 

targets <- as.data.frame(cbind(names(proteomic_plasmax),names(proteomic_plasmax)))
targets <- targets[-1,]
names(targets) <- c("sample","condition")
targets$condition <- gsub("LFQ intensity ","",targets$condition)
targets$condition <- gsub("_.*","",targets$condition)

plots <- magicPlotMakerLight(proteomic_plasmax, targets)

plot(plots[[2]])

#we want to compare the KO condition with the WT condition so we build a
#comparison list
comparisons <- list("786_OvsHK2" = c(2,-1),
                    "786_M1AvsHK2" = c(3,-1),
                    "786_M2AvsHK2" = c(4,-1),
                    "786_M1Avs786_O" = c(3,-2)) 

#now that the comparisons are defined, we can run limma
limmaRes <- runLimma(measurements = proteomic_plasmax[,-1], 
                     targets = targets, 
                     comparisons = comparisons)

ttop_786_OvsHK2 <- ttopFormatter(topTable(limmaRes[[1]], 
                                          coef = 1, number = length(proteomic_plasmax[,1]), 
                                          adjust.method = "fdr"))

ttop_786_M1AvsHK2 <- ttopFormatter(topTable(limmaRes[[1]], 
                                          coef = 2, number = length(proteomic_plasmax[,1]), 
                                          adjust.method = "fdr"))

ttop_786_M1Avs786_O <- ttopFormatter(topTable(limmaRes[[1]], 
                                          coef = 4, number = length(proteomic_plasmax[,1]), 
                                          adjust.method = "fdr"))

write_csv(ttop_786_OvsHK2, file ="results/proteomic/plasmax/ttop_786_OvsHK2.csv")
write_csv(ttop_786_M1AvsHK2, file ="results/proteomic/plasmax/ttop_786_M1AvsHK2.csv")
write_csv(ttop_786_M1Avs786_O, file ="results/proteomic/plasmax/ttop_786_M1Avs786_O.csv")

pathways <- import_gmt("support/c2.cp.v7.0.symbols.gmt")
pathways <- pathways[grep("KEGG",pathways$term),]

ttop_786_OvsHK2_BCCA <- ttop_786_OvsHK2[
  ttop_786_OvsHK2$ID %in% pathways[
    pathways$term == "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","gene"],]

volcano_nice(ttop_786_OvsHK2_BCCA, FCIndex = 2, IDIndex = 1, pValIndex = 5, nlabels = 30, hAss = 0.5, vAss = 0.1, label = T) + ylab("-log10(p-value)")
