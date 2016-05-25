#Prova seleccio atributs

library("SummarizedExperiment")
library("edgeR")
library("geneplotter")
library("ggplot2")
library("sva")
old.par <- par(mar = c(0, 0, 0, 0)) # to get CLI 

prac.se <- readRDS("data/sePRAD.rds")
n_occur <- data.frame(table(colData(prac.se)$bcr_patient_barcode))
prac.se.unique <- prac.se[, colData(prac.se)$bcr_patient_barcode 
                          %in% n_occur$Var1[n_occur$Freq == 1] & 
                            !is.na(colData(prac.se)$bcr_patient_barcode) |
                            (colData(prac.se)$bcr_patient_barcode %in% 
                               n_occur$Var1[n_occur$Freq > 1] & 
                               colData(prac.se)$type == "normal" & 
                               !is.na(colData(prac.se)$bcr_patient_barcode))]

prac.dge <- DGEList(counts = assays(prac.se)$counts, genes = mcols(prac.se), group = prac.se$type)


assays(prac.se)$logCPM <- cpm(prac.dge, log=TRUE, prior.count=3.5)

prac.dge.unique <- DGEList(counts = assays(prac.se.unique)$counts, genes = mcols(prac.se.unique), group = prac.se.unique$type)

# library size, non-paired samples
prac.dge.unique.filtlib <- prac.dge.unique[,(prac.dge.unique$samples$lib.size/1e6) > 40 ] # 50e06 threshold value

# adding to prac.se the filtering by threshold value.
prac.se.sub <- prac.se.unique[,rownames(prac.dge.unique.filtlib$samples)]

# TSS table of the complete values
tss <- substr(colnames(prac.se.sub), 6, 7)
tss_df <- data.frame(TYPE=prac.se.sub$type, TSS=tss)

selection  <- c('CH', 'EJ', 'G9', 'HC')
prac.dge.unique.filtlib <- prac.dge.unique.filtlib[,tss_df$TSS %in% selection]
prac.se.sub <- prac.se.unique[,rownames(prac.dge.unique.filtlib$samples)]


GetID <- function(number, se){
  id <- grep(number,rownames(colData(se)))
  return(rownames(colData(prac.se.sub))[id])
}
number.id <- c("7737","7745","7740","7738","7211","7747","8258") #from MDS plot
pruning.vec <- sapply(number.id, GetID, prac.se.sub)

prac.dge.unique.pruned <- prac.dge.unique.filtlib[,
                                                  !rownames(prac.dge.unique.filtlib$samples) %in% pruning.vec]
prac.dge.unique.pruned$counts

#genes, samples,counts

###########################################################3