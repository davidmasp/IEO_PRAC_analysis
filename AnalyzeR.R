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
head(dge$samples)

logCPM <- cpm(dge, log = TRUE, prior.count = 3.25)

multidensity(as.list(as.data.frame(logCPM)), xlab = "log2 CPM", legend = NULL, main = "Count density vs expression level")

plotSmear(dge, lowess = TRUE)
abline(h = 0, col = "blue", lwd = 2)

dgenorm <- calcNormFactors(dge)

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


library(sva)

mod <- model.matrix(~ type + tissue_source_site, data = colData(pracse)) 
head(mod)
mod0 <- model.matrix(~tissue_source_site, data = colData(pracse))

sv <- sva(logCPM, mod, mod0)
sv






##########JOAN_BEGIN

##########JOAN_END


##########ADRIA_BEGIN

#########ADRIA_END


#########DAVID_BEGIN

#########DAVID_END
