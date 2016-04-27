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

filter.NAS <- function(i){
  if (sum(is.na(i)) / length(i) < 0.5) {
    TRUE
  } else {
    FALSE
  }
}

sample.info.mask <- sapply(colData(pracse), filter.NAS)


# dim(sample.info)
# colnames(sample.info) # get the sample info fields

sample.info.df <- as.data.frame(sample.info)

## Ploting some descriptives from samples

qplot(race, data=sample.info.df, geom="bar", fill=ethnicity)
qplot(gleason_score, data=sample.info.df, geom="bar", fill=ethnicity)

### NORMALITZATION PIPELINE

dge <- DGEList(counts = assays(pracse)$counts, genes = mcols(pracse), group = pracse$type)
names(dge)

head(dge$samples)

logCPM <- cpm(dge, log = TRUE, prior.count = 3.25)

multidensity(as.list(as.data.frame(logCPM)), xlab = "log2 CPM", legend = NULL, main = "Count density vs expression level")

plotSmear(dge, lowess = TRUE)
abline(h = 0, col = "blue", lwd = 2)

dgenorm <- calcNormFactors(dge)

plotMDS(dge, col = c("red", "blue")[as.integer(dgenorm$samples$group)], cex = 0.7)
legend("topleft", c("female", "male"), fill = c("red", "blue"), inset = 0.05, cex = 0.7)
