rm(list = ls(all.names = TRUE))
gc()
.rs.restartR() ##Only if Rstudio

library(readr)
library(fgsea)

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

pathways <- import_gmt("support/c2.cp.v7.0.symbols.gmt")
pathways <- pathways[grep("KEGG",pathways$term),]
pathways_list <- list()
i <- 1
for(pathway in unique(pathways$term))
{
  pathways_list[[i]] <- pathways[pathways$term == pathway,"gene"]
  i <- i+1
}
names(pathways_list) <- unique(pathways$term)


ttop_786vHK2 <- as.data.frame(read_csv("results/proteomic/plasmax/ttop_786_OvsHK2.csv"))
ttop_M1AvHK2 <- as.data.frame(read_csv("results/proteomic/plasmax/ttop_786_M1AvsHK2.csv"))
ttop_M1Av786 <- as.data.frame(read_csv("results/proteomic/plasmax/ttop_786_M1Avs786_O.csv"))

fgsea_res_786vHK2 <- as.data.frame(runfgsea_ttop(ttop_786vHK2, pathways_list, 10000)[,-8])
fgsea_res_M1AvHK2 <- as.data.frame(runfgsea_ttop(ttop_M1AvHK2, pathways_list, 10000)[,-8])
fgsea_res_M1Av786 <- as.data.frame(runfgsea_ttop(ttop_M1Av786, pathways_list, 10000)[,-8])

fgsea_res_786vHK2 <- fgsea_res_786vHK2[order(fgsea_res_786vHK2$NES, decreasing = F),]
fgsea_res_786vHK2$contrast <- "OvHK2"
# fgsea_res_M1AvHK2 <- fgsea_res_M1AvHK2[order(fgsea_res_786vHK2$NES, decreasing = F),]
fgsea_res_M1AvHK2$contrast <- "M1AvHK2"
O_order <- as.vector(sapply(fgsea_res_M1Av786$pathway,function(x)
{
  which(fgsea_res_786vHK2$pathway == x)
}))
fgsea_res_M1Av786 <- fgsea_res_M1Av786[order(O_order,decreasing = F),]


combined_fgseaRes <- as.data.frame(rbind(fgsea_res_786vHK2,fgsea_res_M1AvHK2))
# combined_fgseaRes <- combined_fgseaRes[combined_fgseaRes$padj < 0.05,]
combined_fgseaRes$pathway <- gsub("KEGG_","",combined_fgseaRes$pathway)
combined_fgseaRes$pathway <- tolower(gsub("_"," ",combined_fgseaRes$pathway))
combined_fgseaRes$pathway <- factor(combined_fgseaRes$pathway,levels = unique(combined_fgseaRes$pathway))

library(ggplot2)

ggplot(combined_fgseaRes[combined_fgseaRes$padj < 0.05,], aes(x = NES, y = pathway)) +
  geom_point(aes(colour = contrast, size = -log10(padj))) +
  # scale_color_gradient2(low="blue", high="red", midpoint = 0, mid ="white") +
  theme_minimal() + ggtitle("Metastasis and tumor versus HK2 FGSEA")

fgsea_res_M1Av786$pathway <- gsub("KEGG_","",fgsea_res_M1Av786$pathway)
fgsea_res_M1Av786$pathway <- tolower(gsub("_"," ",fgsea_res_M1Av786$pathway))
fgsea_res_M1Av786 <- fgsea_res_M1Av786[fgsea_res_M1Av786$pathway %in% combined_fgseaRes[combined_fgseaRes$padj < 0.05,"pathway"],]

fgsea_res_M1Av786 <- fgsea_res_M1Av786[order(fgsea_res_M1Av786$NES, decreasing = F),]

fgsea_res_M1Av786$pathway <- factor(fgsea_res_M1Av786$pathway, levels = fgsea_res_M1Av786$pathway)

ggplot(fgsea_res_M1Av786, aes(x = NES, y = pathway, size = -log10(padj), colour = NES)) +
  geom_point() +
  scale_color_gradient2(low="blue", high="blue", midpoint = 0, mid = "blue") +
  theme_minimal() + ggtitle("Metastasis versus tumor FGSEA")

write_csv(fgsea_res_786vHK2, "results/proteomic/plasmax/fgsea_OvHK2.csv")
write_csv(fgsea_res_M1AvHK2, "results/proteomic/plasmax/fgsea_M1AvHK2.csv")
write_csv(fgsea_res_M1Av786, "results/proteomic/plasmax/fgsea_M1AvO.csv")

