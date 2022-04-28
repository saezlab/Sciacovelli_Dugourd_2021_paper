library(readxl)
library(pheatmap)
library(vsn)
library(readr)
library(limma)

Plasmax_metabolomics_vs_standard_medium <- as.data.frame(read_excel("data/metabolomic/Plasmax metabolomics vs standard medium.xlsx"))

targets <- Plasmax_metabolomics_vs_standard_medium[,c(1,2)]
names(targets) <- c("sample","condition")

batches <- as.data.frame(t(Plasmax_metabolomics_vs_standard_medium[,-c(1,2)]))
names(batches) <- Plasmax_metabolomics_vs_standard_medium[,1]
row.names(batches) <- gsub(" ", "_", row.names(batches))
row.names(batches) <- tolower(row.names(batches))

fit <- vsnMatrix(as.matrix(batches))
meanSdPlot(fit)
batches <- as.data.frame(vsn::predict(fit,as.matrix(batches)))

###################
#Import data
raw_metabolomic <- as.data.frame(read_excel("data/metabolomic/cosmos_paper/raw_metabolomic.xlsx", 
                                            na = "NF"))

metab_id_correspondance <- as.data.frame(read_delim("data/metabolomic/cosmos_paper/metab_id_correspondance", 
                                                    "\t", escape_double = FALSE, trim_ws = TRUE))

names(raw_metabolomic) <- c("ID", metab_id_correspondance$ID)
names(raw_metabolomic) <- gsub(" T","TU", names(raw_metabolomic))
names(raw_metabolomic) <- gsub(" H","KI", names(raw_metabolomic))

row.names(raw_metabolomic) <- raw_metabolomic[,1]
raw_metabolomic <- raw_metabolomic[,-1]

samples_to_keep <- c(paste(c("11", "15", "16", "24", "29", "32", "38", "40", "35"), "TU", sep = ""), paste(c("11", "15", "16", "24", "29", "32", "38", "40", "35"), "KI", sep = ""))
samples_to_keep <- samples_to_keep[samples_to_keep %in% names(raw_metabolomic)]

raw_metabolomic <- raw_metabolomic[,samples_to_keep]

row.names(raw_metabolomic) <- gsub(" ", "_", row.names(raw_metabolomic))
row.names(raw_metabolomic) <- tolower(row.names(raw_metabolomic))

targets_1 <- as.data.frame(matrix(NA,length(raw_metabolomic[1,]),2))
names(targets_1) <- c("sample","condition")
targets_1$sample <- names(raw_metabolomic)
targets_1$condition <- gsub("[0-9]*","",targets_1$sample)

fit <- vsnMatrix(as.matrix(raw_metabolomic))
meanSdPlot(fit)
raw_metabolomic <- as.data.frame(vsn::predict(fit,as.matrix(raw_metabolomic)))
##################

targets_cell_lines <- targets
targets_patients <- targets_1

cell_lines <- batches
names(cell_lines) <- paste(targets_cell_lines$condition, c(1:length(targets_cell_lines[,1])), sep = "")
targets_cell_lines$sample <- names(cell_lines) 

patient_samples <- raw_metabolomic

cell_lines$metabolites <- row.names(cell_lines)


patient_samples$metabolites <- row.names(patient_samples)

all_samples <- merge(cell_lines, patient_samples, by = "metabolites")

all_samples <- all_samples[complete.cases(all_samples),]
cols <- colorRampPalette(c("white","blue"))(100)
pheatmap(cor(all_samples[,-1], method = "spearman"), display_numbers = F, color = cols)
##### CONTROL_PLASMAX ########

controls_plasmax <- merge(cell_lines[,names(cell_lines) == "metabolites" | names(cell_lines) %in% targets_cell_lines[targets_cell_lines$condition == "HK2_Plasmax","sample"]],
                          patient_samples[,names(patient_samples) == "metabolites" | names(patient_samples) %in% targets_patients[targets_patients$condition == "KI", "sample"]],
                          by = "metabolites")
row.names(controls_plasmax) <- controls_plasmax$metabolites
controls_plasmax <- controls_plasmax[,-1]
controls_plasmax <- controls_plasmax[complete.cases(controls_plasmax),]
# controls_plasmax <- as.data.frame(t(controls_plasmax))

cor_controls_plasmax <- as.data.frame(cor(controls_plasmax, method = "spearman"))

mean(as.matrix(cor_controls_plasmax[c(7:14),c(1:6)]))

######## CONTROL RPMI ##########

controls_RPMI <- merge(cell_lines[,names(cell_lines) == "metabolites" | names(cell_lines) %in% targets_cell_lines[targets_cell_lines$condition == "HK2_RPMI","sample"]],
                       patient_samples[,names(patient_samples) == "metabolites" | names(patient_samples) %in% targets_patients[targets_patients$condition == "KI", "sample"]],
                       by = "metabolites")
row.names(controls_RPMI) <- controls_RPMI$metabolites
controls_RPMI <- controls_RPMI[,-1]
controls_RPMI <- controls_RPMI[complete.cases(controls_RPMI),]
# controls_RPMI <- as.data.frame(t(controls_RPMI))

cor_controls_RPMI <- as.data.frame(cor(controls_RPMI, method = "spearman"))

mean(as.matrix(cor_controls_RPMI[c(7:14),c(1:6)]))

t.test(as.vector(cor_controls_plasmax[c(7:14),c(1:6)]),as.vector(cor_controls_RPMI[c(7:14),c(1:6)]))

######### TUMOR PLASMAX ##########
tumor_plasmax <- merge(cell_lines[,names(cell_lines) == "metabolites" | names(cell_lines) %in% targets_cell_lines[targets_cell_lines$condition == "786-O_Plasmax","sample"]],
                       patient_samples[,names(patient_samples) == "metabolites" | names(patient_samples) %in% targets_patients[targets_patients$condition == "TU", "sample"]],
                       by = "metabolites")
row.names(tumor_plasmax) <- tumor_plasmax$metabolites
tumor_plasmax <- tumor_plasmax[,-1]
tumor_plasmax <- tumor_plasmax[complete.cases(tumor_plasmax),]

cor_tumor_plasmax <- as.data.frame(cor(tumor_plasmax, method = "spearman"))
mean(as.matrix(cor_tumor_plasmax[c(7:14),c(1:6)]))

########## TUMOR RPMI ############
tumor_RPMI <- merge(cell_lines[,names(cell_lines) == "metabolites" | names(cell_lines) %in% targets_cell_lines[targets_cell_lines$condition == "786-O_RPMI","sample"]],
                    patient_samples[,names(patient_samples) == "metabolites" | names(patient_samples) %in% targets_patients[targets_patients$condition == "TU", "sample"]],
                    by = "metabolites")
row.names(tumor_RPMI) <- tumor_RPMI$metabolites
tumor_RPMI <- tumor_RPMI[,-1]
tumor_RPMI <- tumor_RPMI[complete.cases(tumor_RPMI),]

cor_tumor_RPMI <- as.data.frame(cor(tumor_RPMI, method = "spearman"))
mean(as.matrix(cor_tumor_RPMI[c(7:14),c(1:6)]))

t.test(as.vector(cor_tumor_plasmax[c(7:14),c(1:6)]),as.vector(cor_tumor_RPMI[c(7:14),c(1:6)]))

########## CONTROL TUMOR PLASMAX ###########
control_tumor_plasmax <- merge(cell_lines[,names(cell_lines) == "metabolites" | names(cell_lines) %in% targets_cell_lines[targets_cell_lines$condition == "786-O_Plasmax","sample"]],
                               patient_samples[,names(patient_samples) == "metabolites" | names(patient_samples) %in% targets_patients[targets_patients$condition == "KI", "sample"]],
                               by = "metabolites")
row.names(control_tumor_plasmax) <- control_tumor_plasmax$metabolites
control_tumor_plasmax <- control_tumor_plasmax[,-1]
control_tumor_plasmax <- control_tumor_plasmax[complete.cases(control_tumor_plasmax),]

cor_control_tumor_plasmax <- as.data.frame(cor(control_tumor_plasmax, method = "spearman"))
mean(as.matrix(cor_control_tumor_plasmax[c(7:14),c(1:6)]))

############ TUMOR CONTROL PLASMAX ##############
tumor_control_plasmax <- merge(cell_lines[,names(cell_lines) == "metabolites" | names(cell_lines) %in% targets_cell_lines[targets_cell_lines$condition == "HK2_Plasmax","sample"]],
                               patient_samples[,names(patient_samples) == "metabolites" | names(patient_samples) %in% targets_patients[targets_patients$condition == "TU", "sample"]],
                               by = "metabolites")
row.names(tumor_control_plasmax) <- tumor_control_plasmax$metabolites
tumor_control_plasmax <- tumor_control_plasmax[,-1]
tumor_control_plasmax <- tumor_control_plasmax[complete.cases(tumor_control_plasmax),]

cor_tumor_control_plasmax <- as.data.frame(cor(tumor_control_plasmax, method = "spearman"))
mean(as.matrix(cor_tumor_control_plasmax[c(7:14),c(1:6)]))
