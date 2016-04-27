# david path setwd("/Users/davidmasp/GD_UPF-bioinformatica/IEO/IEO_PRAC_analysis/data/")

### Loading packages:
library("SummarizedExperiment")
library("edgeR", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
library("geneplotter")
library("ggplot2")

### Reading data

pracse <- readRDS("sePRAD.rds")
pracse

#check if it is the new data
metadata(pracse)$objectCreationDate #the result should be [1] "2016-04-25"

mcols(colData(pracse), use.names=TRUE)


### elements of the data

feature.info <- mcols(pracse)
sample.info <- colData(pracse)
dim(sample.info)
colnames(sample.info) # get the sample info fields
sample.info.df <- as.data.frame(sample.info)

## Ploting some descriptives from samples

qplot(race, data=sample.info.df, geom="bar", fill=ethnicity)

### NORMALITZATION PIPELINE

dge <- DGEList(counts = assays(pracse)$counts, genes = mcols(pracse), group = pracse$type)
names(dge)

head(dge$samples)

logCPM <- cpm(dge, log = TRUE, prior.count = 3.25)

multidensity(as.list(as.data.frame(logCPM)), xlab = "log2 CPM", legend = NULL, main = "Count density vs expression level")

plotSmear(dge, lowess = TRUE)
abline(h = 0, col = "blue", lwd = 2)

plotMDS(dge, col = c("red", "blue")[as.integer(dgenorm$samples$group)], cex = 0.7)
legend("topleft", c("female", "male"), fill = c("red", "blue"), inset = 0.05, cex = 0.7)
