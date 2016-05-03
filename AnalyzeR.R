# david path setwd("/Users/davidmasp/GD_UPF-bioinformatica/IEO/IEO_PRAC_analysis/")

# PRAC Analysis

#### Authors: 
#### ADRIA AULADELL
#### JOAN MARTI
#### DAVID MAS


## 0.1 LOADING PACKAGES

library("SummarizedExperiment")
library("edgeR")
library("geneplotter")
library("ggplot2")

## 0.2 READING DATA

pracse <- readRDS("data/sePRAD.rds")
pracse

# check if it is the new data
metadata(pracse)$objectCreationDate #the result should be [1] "2016-04-25"

mcols(colData(pracse), use.names=TRUE)

### bcr_patient_barcode




### elements of the data

feature.info <- mcols(pracse)
sample.info <- colData(pracse)
sample.info.sum <- sapply(sample.info, summary)
head(sample.info)
TSS <- substr(rownames(sample.info),6,7) # Tissue Sort Site
center <- substr(rownames(sample.info),27,28)
summary(as.data.frame(center))
summary(as.data.frame(TSS))

# to filter for NA values (currently not working for NOT AVAILABLE)
filter.NAS <- function(i) {
  if (sum(is.na(i)) / length(i) < 0.20) {
    TRUE
  } else {
    FALSE
  }
}

# to show 
colData(pracse)$project_code

# to apply the filter.NAS function to the 
sample.info.mask <- sapply(colData(pracse), filter.NAS)

# to get a list of sample variants
colnames(sample.info[sample.info.mask])

# to show a summary of a sample variable
summary(sample.info[,"tissue_source_site"])

# to check the CDEID code for a sample variable
mcols(colData(pracse), use.names=TRUE)["tissue_source_site",]


### POSSIBLE BATCH EFFECT
# prospective_collection
# tissue_source_site
summary(as.data.frame(sample.info[,"tissue_source_site"])) # as all of them are andenocarcinomas we can check TSS as the place of processing.

# dim(sample.info)
# colnames(sample.info) # get the sample info fields

sample.info.df <- as.data.frame(sample.info)

## Ploting some descriptives from samples

qplot(race, data=sample.info.df, geom="bar", fill=ethnicity)
qplot(gleason_score, data=sample.info.df, geom="bar", fill=ethnicity)

### NORMALITZATION PIPELINE

dge <- DGEList(counts = assays(pracse)$counts, genes = mcols(pracse), group = pracse$type)
names(dge)
summary(pracse$type)
head(dge$samples$lib.size)

mean(dge$samples$lib.size)
dge.filtred <- dge[,dge$samples$lib.size > mean(dge$samples$lib.size) ]

dge.filtred 
mean(dge.filtred$samples$lib.size)

logCPM <- cpm(dge, log = TRUE, prior.count = 3.25)

multidensity(as.list(as.data.frame(logCPM)), xlab = "log2 CPM", legend = NULL, main = "Count density vs expression level")

plotSmear(dge, lowess = TRUE)
abline(h = 0, col = "blue", lwd = 2)

dgenorm <- calcNormFactors(dge)

dgenormfilt <- dgenorm[, !rownames(dgenorm$samples) %in% "NA19172"]


#########DAVID_BEGIN
class(duplicated(colData(pracse)$bcr_patient_barcode))
length(duplicated(colData(pracse)$bcr_patient_barcode))
class(unique(colData(pracse)$bcr_patient_barcode))
class(colData(pracse)$bcr_patient_barcode)

pracse.filt <- pracse[, duplicated(colData(pracse)$bcr_patient_barcode)]
tumor.list <- pracse.filt$colData(pracse)$type == "tumor" 

duplicated(colData(pracse.filt)$bcr_patient_barcode)
length(colData(pracse.filt)$bcr_patient_barcode)

dup.list <- colData(pracse)$bcr_patient_barcode[duplicated(colData(pracse)$bcr_patient_barcode)]

pracse.filt.paired <- pracse[, colData(pracse)$bcr_patient_barcode %in% dup.list]

duplicated(colData(pracse.filt.paired)$bcr_patient_barcode)

pracse.filt.paired

ALL_RECORDS <- df[df$ID==df$ID[duplicated(df$ID)],]
pracse.filt.cases.dup <- pracse[,duplicated(colData(pracse)$bcr_patient_barcode)]

duplicated(colData(pracse.filt.cases.dup)$bcr_patient_barcode)
unique(colData(pracse.filt.cases.dup)$bcr_patient_barcode)


pracse.filt.controls.dup <- pracse[,duplicated(colData(pracse)$bcr_patient_barcode) & colData(pracse)$type == "normal" ]

unique(colData(pracse.filt.controls.dup)$bcr_patient_barcode)


barcode <- colData(pracse)$bcr_patient_barcode

pracse.filt.unique <- pracse[, !duplicated(barcode) | (!duplicated(barcode) & colData(pracse)$type == "normal") ]
pracse.filt.unique

dup.list <- barcode[duplicated(barcode)]
length(dup.list)
pracse.filt.paired <- pracse[,barcode
                             %in% dup.list & !is.na(barcode)]
duplicated(colData(pracse.filt.dupl)$bcr_patient_barcode)

unique <- unique(colData(pracse)$bcr_patient_barcode)
tumor_unique <- pracse[ , colData(pracse)$bcr_patient_barcode %in%  unique]
tumor_unique


n_occur <- data.frame(table(colData(pracse)$bcr_patient_barcode))

pracse.filt.dupl <- pracse[,colData(pracse)$bcr_patient_barcode %in% n_occur$Var1[n_occur$Freq > 1] & !is.na(colData(pracse)$bcr_patient_barcode)]

pracse.filt.dupl

dup.list.tumor <- colData(pracse)$bcr_patient_barcode[]

pracse.filt.unique <- pracse[,colData(pracse)$bcr_patient_barcode %in% n_occur$Var1[n_occur$Freq == 1] & !is.na(colData(pracse)$bcr_patient_barcode) ]

pracse.filt.unique

table(colData(pracse.filt.unique)$type)


#########DAVID_END

plotMDS(dge, col = c("red", "blue")[as.integer(dgenorm$samples$group)], cex = 0.7)
legend("topleft", c("female", "male"), fill = c("red", "blue"), inset = 0.05, cex = 0.7)



# Batch Effect Detection

## Check BE in prospective collection

pc <- pracse$prospective_collection

table(data.frame(TYPE = pracse$type, PC = pc))

d <- as.dist(1 - cor(logCPM, method = "spearman")) 
sampleClustering <- hclust(d)

batch <- as.integer(pracse$prospective_collection)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- colnames(pracse)
outcome <- as.character(pracse$type)
names(outcome) <- colnames(pracse)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) {
  ## for every node in the dendrogram if it is a leaf node if (is.leaf(x)) {
  attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")])) ## c
  attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome }
  x
}, batch, outcome)  ## these are the second and third arguments in the function

plot(sampleDendrogram)
legend("topright", paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

## Result Negative


## Check for TSS

pc <- TSS

table(data.frame(TYPE = pracse$type, TSS = pc))

d <- as.dist(1 - cor(logCPM, method = "spearman")) 
sampleClustering <- hclust(d)

ConvertNamesToColor <- function(x){
  x <- as.factor(x)
  levels(x) <- 1:length(levels(x))
  x <- as.numeric(x)
  return(x)
}


batch <- ConvertNamesToColor(TSS)

sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- colnames(pracse)
outcome <- as.character(pracse$type)
names(outcome) <- colnames(pracse)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) {
  ## for every node in the dendrogram if it is a leaf node if (is.leaf(x)) {
  attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")])) ## c
  attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome }
  x
}, batch, outcome)  ## these are the second and third arguments in the function

plot(sampleDendrogram)
legend("topright", paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

## NOT SURE -> go to sva


library("sva")

mod <- model.matrix(~ type + tissue_source_site, data = colData(pracse)) 
head(mod)
mod0 <- model.matrix(~tissue_source_site, data = colData(pracse))

sv <- sva(logCPM, mod, mod0)
sv




ord <- order(dge$sample$lib.size/1e6)
barplot(dge$sample$lib.size[ord]/1e6, las=1, ylab="Millions of reads",
        xlab="Samples", col=c("blue", "red")[(pracse$type[ord] == "tumor") + 1])
legend("topleft", c("tumor", "normal"), fill=c("red", "blue"), inset=0.01)





##########JOAN_BEGIN
	
plot(density(dge$samples$lib.size/1e6))

dge.filtred <- dge[,(dge$samples$lib.size/1e6) > 50 ]

ord <- order(dge.filtred$samples$lib.size/1e6)
barplot(dge.filtred$sample$lib.size[ord]/1e6, las=1, ylab="Millions of reads",
        xlab="Samples", col=c("blue", "red")[(pracse$type[ord] == "tumor") + 1])
legend("topleft", c("tumor", "normal"), fill=c("red", "blue"), inset=0.01)
# normal 25 of 52, tumor 112 of 502
	
##########JOAN_END

##########ADRIA_BEGIN

pracse[, !duplicated(colData(pracse)$bcr_patient_barcode) & duplicated(colData(pracse)$type == 'normal')]
unique <- unique(colData(pracse)$bcr_patient_barcode)
tumor_unique <- pracse[ , colnames(pracse) %in%  unique]
# tumor_unique <- unique(pracse[,(colData(pracse)$type == 'tumor')])
normal
normal_dupl <- pracse[, colData(pracse)$type == 'normal']
normal_dupl
normal_dupl[,!duplicated(colData(normal_dupl)$bcr_patient_barcode)]


# Non duplicated values 

pracse.unique <- pracse[,!duplicated(colData(pracse)$bcr_patient_barcode) ] 

table(colData(pracse.unique)$type)

dup.list <- colData(pracse)$bcr_patient_barcode[duplicated(colData(pracse)$bcr_patient_barcode)]
length(dup.list)
pracse.filt.dupl <- pracse[,colData(pracse)$bcr_patient_barcode %in% dup.list & !is.na(colData(pracse)$bcr_patient_barcode)]
duplicated(colData(pracse.filt.dupl)$bcr_patient_barcode)

table(pracse.filt.dupl$type)
# IMPLEMENTACIÓ TAUL

n_occur <- data.frame(table(colData(pracse)$bcr_patient_barcode))

pracse.filt.dupl <- pracse[,colData(pracse)$bcr_patient_barcode %in% n_occur$Var1[n_occur$Freq > 1] & !is.na(colData(pracse)$bcr_patient_barcode)]

pracse.filt.unique <- pracse[,colData(pracse)$bcr_patient_barcode %in% n_occur$Var1[n_occur$Freq == 1] & !is.na(colData(pracse)$bcr_patient_barcode)]


#########ADRIA_END


#########DAVID_BEGIN


dge <- DGEList(counts = assays(pracse)$counts, genes = mcols(pracse), group = pracse$type)
head(dge$samples)

assays(pracse)$logCPM <- cpm(dge, log=TRUE, prior.count=3.5)
assays(pracse)$logCPM[1:5, 1:5]

n_occur <- data.frame(table(colData(pracse)$bcr_patient_barcode))

pracse.filt.dupl <- pracse[,colData(pracse)$bcr_patient_barcode %in% n_occur$Var1[n_occur$Freq > 1] & !is.na(colData(pracse)$bcr_patient_barcode)]

pracse.filt.unique <- pracse[,colData(pracse)$bcr_patient_barcode %in% n_occur$Var1[n_occur$Freq == 1] & !is.na(colData(pracse)$bcr_patient_barcode)]

dge
assays(pracse.filt.unique)$logCPM[1:5, 1:5]



dge.filtred <- assays(pracse)$logCPM[,(dge$samples$lib.size/1e6) > 50 ]
#########DAVID_END

