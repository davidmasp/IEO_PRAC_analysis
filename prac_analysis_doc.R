## ----install, eval=FALSE-------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## library(BiocInstaller)
## biocLite(c("edgeR","SummarizedExperiment","geneplotter","sva"))
## install.packages(ggplot2)

## ----loading, message=FALSE----------------------------------------------
library("SummarizedExperiment")
library("edgeR")
library("geneplotter")
library("ggplot2")
library("sva")
old.par <- par(mar = c(0, 0, 0, 0)) # to get CLI 

## ----reading-------------------------------------------------------------
prac.se <- readRDS("data/sePRAD.rds")
prac.se

## ----dim-----------------------------------------------------------------
dim(colData(prac.se))

## ----data----------------------------------------------------------------
head(names(colData(prac.se)))
colData(prac.se)[1:5, 1:5]
mcols(colData(prac.se), use.names=TRUE)

## ----exploring some variables, fig.align='center', echo=FALSE, message=FALSE, fig.width=12----
sampleinfo.df <- as.data.frame(colData(prac.se))
# require(gridExtra)
plot1 <- ggplot(sampleinfo.df, aes(x=vital_status, fill =gleason_score)) + geom_bar()
plot2 <- ggplot(sampleinfo.df, aes(x=race, fill =gleason_score)) + geom_bar() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
# grid.arrange(plot1, plot2, ncol=2)

## ----table types---------------------------------------------------------
table(prac.se$type)


## ----row data------------------------------------------------------------
rowRanges(prac.se)

## ----gc content hist, fig.align='center', echo=FALSE, fig.width=12-------
geneinfo.df <- as.data.frame(mcols(prac.se))
plot1 <- ggplot(geneinfo.df, aes(x=txlen)) + geom_histogram(bins = 50)
plot2 <- ggplot(geneinfo.df, aes(x=txgc)) + geom_histogram(bins = 50)
# grid.arrange(plot1, plot2, ncol=2)

## ----metadata------------------------------------------------------------
metadata(prac.se)$objectCreationDate #the result should be [1] "2016-04-25"

## ----subset--------------------------------------------------------------
# Subset the paired samples
## we build a table from the bcr barcode occurences in our dataset
n_occur <- data.frame(table(colData(prac.se)$bcr_patient_barcode))

## we filter out the na patients (having a na in the barcode)
## and the ones that occur only 1 time (they are unique)
prac.se.dupl <- prac.se[,colData(prac.se)$bcr_patient_barcode
                 %in% n_occur$Var1[n_occur$Freq > 1] & 
                 !is.na(colData(prac.se)$bcr_patient_barcode)]
## See the results.
table(colData(prac.se.dupl)$type)

## ----nonpaired subset----------------------------------------------------
# Subset the non-paired 
prac.se.unique <- prac.se[, colData(prac.se)$bcr_patient_barcode 
                             %in% n_occur$Var1[n_occur$Freq == 1] & 
                                !is.na(colData(prac.se)$bcr_patient_barcode) |
                                (colData(prac.se)$bcr_patient_barcode %in% 
                                n_occur$Var1[n_occur$Freq > 1] & 
                                colData(prac.se)$type == "normal" & 
                                !is.na(colData(prac.se)$bcr_patient_barcode))]
table(colData(prac.se.unique)$type)


## ----dge, warning= FALSE-------------------------------------------------
prac.dge <- DGEList(counts = assays(prac.se)$counts, genes = mcols(prac.se), group = prac.se$type)
head(prac.dge$samples)

## ----log2----------------------------------------------------------------
assays(prac.se)$logCPM <- cpm(prac.dge, log=TRUE, prior.count=3.5)
assays(prac.se)$logCPM[1:5, 1:5]

## ----normalitzation of the filtred dataset, warning= FALSE---------------
prac.dge.unique <- DGEList(counts = assays(prac.se.unique)$counts, genes = mcols(prac.se.unique), group = prac.se.unique$type)

## ----library size all after filter unique, echo=FALSE, fig.align='center', fig.cap='Fig. 1: Library Size per sample in the  non-paired design dataset. '----

ord <- order(prac.dge.unique$samples$lib.size/1e6)
barplot((prac.dge.unique$sample$lib.size/1e06)[ord], las=1, ylab="Millions of reads",
        xlab="Samples", col=c("blue", "red")[(prac.se$type[ord] == "tumor") + 1])
legend("topleft", c("tumor", "normal"), fill=c("red", "blue"), inset=0.01)


## ----lib.unique----------------------------------------------------------
# library size, non-paired samples
prac.dge.unique.filtlib <- prac.dge.unique[,(prac.dge.unique$samples$lib.size/1e6) > 40 ] # 50e06 threshold value
table(prac.dge.unique.filtlib$samples$group)

## ----filter--------------------------------------------------------------
# not filtering by TSS
prac.se.sub <- prac.se.unique[,rownames(prac.dge.unique.filtlib$samples)]

# filter by TSS block
tss <- substr(colnames(prac.se.sub), 6, 7)
tss_df <- data.frame(TYPE=prac.se.sub$type, TSS=tss)
table(data.frame(TYPE=prac.se.sub$type, TSS=tss))

## ----selec---------------------------------------------------------------
selection  <- c('CH', 'EJ', 'G9', 'HC')
prac.dge.unique.filtlib <- prac.dge.unique.filtlib[,tss_df$TSS %in% selection]
prac.se.sub <- prac.se.unique[,rownames(prac.dge.unique.filtlib$samples)]
table(prac.dge.unique.filtlib$samples$group)

## ----figure library size, echo=FALSE, fig.align='center', fig.cap='Fig. 2: Density sample distribution using already filtred data by sequencing depth.'----

ord.unique <- order(prac.dge.unique.filtlib$samples$lib.size/1e6)

libsize <- data.frame(libsize = prac.dge.unique.filtlib$sample$lib.size/1e6, type = prac.se.sub$type, samples = rownames(prac.dge.unique.filtlib$samples) )

summary(libsize)


ggplot(libsize) +
  geom_density(aes(x=libsize, fill=type), alpha=0.5) +
  ylab("Density\n") + xlab("\nMillions of Reads") +
  theme_bw()


## ----barplot libsize samples filtred 2, echo=FALSE, fig.align='center',fig.cap="Fig. 3: Barplot of the final selected samples with its library size."----


ggplot(libsize, aes(x=reorder(samples, libsize), y=libsize, fill=type)) + geom_bar(stat = "identity") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())

#barplot((prac.dge.unique.filtlib$sample$lib.size/1e06)[ord.unique], las=1, ylab="Millions of reads",
#        xlab="Samples", col=c("blue", "red")[(prac.dge.unique.filtlib$sample$group[ord.unique] == "tumor") + 1])

#legend("topleft", c("tumor", "normal"), fill=c("red", "blue"), inset=0.01)

## ----Multydensity, fig.align='center', fig.cap='Fig 4: Multidensity plot representing the clusters of the non expressed genes at the left and the expressed genes at the right.', echo=FALSE, fig.width=12----
MDS_normal <- prac.dge.unique.filtlib[,prac.dge.unique.filtlib$samples$group == 'normal']
MDS_tumor <- prac.dge.unique.filtlib[,prac.dge.unique.filtlib$samples$group == 'tumor']

logCPM.MDS_normal <- cpm(MDS_normal, log = TRUE, prior.count = 3.25)
logCPM.MDS_tumor <- cpm(MDS_tumor, log = TRUE, prior.count = 3.25)

par(mfrow=c(1,2))
multidensity(as.list(as.data.frame(logCPM.MDS_normal)), xlab = "log2 CPM", legend = NULL, main = "Normal samples")
multidensity(as.list(as.data.frame(logCPM.MDS_tumor)), xlab = "log2 CPM", legend = NULL, main = "Tumor samples")


## ----avgexp gene, fig.align='center', fig.cap='Fig. 5: Histogram presenting gene frequency expression for different logCPM values'----
assays(prac.se.sub)$logCPM <- cpm(prac.dge.unique.filtlib, log=TRUE, prior.count=0.5)

avgexp <- rowMeans(assays(prac.se.sub)$logCPM)

hist(avgexp, xlab="log2 CPM", main="", las=1)

abline(v=1, col="red", lwd=2)


## ----gene filt-----------------------------------------------------------
mask <- avgexp > 1
dim(prac.se.sub)

prac.se.sub <- prac.se.sub[mask, ]
dim(prac.se.sub)

dim(prac.dge.unique.filtlib )


prac.dge.unique.filtlib <- prac.dge.unique.filtlib[mask, ]
dim(prac.dge.unique.filtlib)


## ----calcNOrm------------------------------------------------------------
prac.dge.unique.filtlib <- calcNormFactors(prac.dge.unique.filtlib)
assays(prac.se.sub)$logCPM <- cpm(prac.dge.unique.filtlib, log=TRUE, prior.count=0.5)

## ----maPlotsTumor, fig.align="center", fig.height= 35, fig.width=12, dpi=100, echo=FALSE, fig.cap="Fig. 6: MA-plots of the tumor samples."----

setmp <- prac.se.sub[, prac.se.sub$type == "tumor"]
dgetmp <- prac.dge.unique.filtlib[, prac.dge.unique.filtlib$type == "tumor"]

samples.number <- length(rownames(colData(setmp)))
if (samples.number %% 6 == 0) { 
  nr <- samples.number / 6
} else {nr <- samples.number %/% 6 + 1}
par(mfrow=c(nr, 6), mar=c(4, 5, 3, 1))

for (i in 1:ncol(setmp)) {
  A <- rowMeans(assays(setmp)$logCPM)
  M <- assays(setmp)$logCPM[, i] - A
  samplename <- substr(as.character(setmp$bcr_patient_barcode[i]), 1, 12)
  smoothScatter(A, M, main=samplename, las=1)
  abline(h=0, col="blue", lwd=2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col="red", lwd=2)
}

## ----maPlotsNormal, fig.align="center", fig.height= 12.5, fig.width=12, dpi=100, echo=FALSE, fig.cap="Fig. 7: MA-plots of the tumor samples."----

setmp <- prac.se.sub[, prac.se.sub$type == "normal"]
dgetmp <- prac.dge.unique.filtlib[, prac.dge.unique.filtlib$type == "normal"]

samples.number <- length(rownames(colData(setmp)))
if (samples.number %% 6 == 0) { 
  nr <- samples.number / 6
} else {nr <- samples.number %/% 6 + 1}
par(mfrow=c(nr, 6), mar=c(4, 5, 3, 1))

for (i in 1:ncol(setmp)) {
  A <- rowMeans(assays(setmp)$logCPM)
  M <- assays(setmp)$logCPM[, i] - A
  samplename <- substr(as.character(setmp$bcr_patient_barcode[i]), 1, 12)
  smoothScatter(A, M, main=samplename, las=1)
  abline(h=0, col="blue", lwd=2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col="red", lwd=2)
}

## ----getting the batch labels, collapse=TRUE-----------------------------
tss <- substr(colnames(prac.se.sub), 6, 7)
table(tss)
center <- substr(colnames(prac.se.sub), 27, 28)
table(center)
plate <- substr(colnames(prac.se.sub), 22, 25)
table(plate)
portionanalyte <- substr(colnames(prac.se.sub), 18, 20)
table(portionanalyte)
samplevial <- substr(colnames(prac.se.sub), 14, 16)
table(samplevial)

## ----table TSS/TYPE------------------------------------------------------
table(data.frame(TYPE=prac.se.sub$type, TSS=tss))

## ----hirechical clustering, fig.align='center', fig.width=14, dpi=100, fig.cap='Fig. 8: Hirechical clustering of the prac.se.sub dataset', fig.height=8----
  
par(old.par)
logCPM <- cpm(prac.dge.unique.filtlib, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(tss))
colors <- palette()
palette(c(colors,"darkorange1","darkmagenta","darkolivegreen","red3","seagreen","pink1","yellow3","sienna1"))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(prac.se.sub)
outcome <- paste(substr(colnames(prac.se.sub), 9, 12), as.character(prac.se.sub$type), sep="-")
names(outcome) <- colnames(prac.se.sub)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)

plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))), fill=sort(unique(batch)))


## ----mds, echo=FALSE, fig.align='center', fig.width=14, dpi=100, fig.cap='Fig. 9: MDS of the prac.se.sub dataset', fig.height=8----
par(mfrow=c(1, 1), mar=c(0,0,0,0))
plotMDS(prac.dge.unique.filtlib, labels=outcome, col=batch)
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))),
       fill=sort(unique(batch)), inset=0.05)

## ----mds2, echo=FALSE, fig.align='center', fig.width=14, dpi=100, fig.cap='Fig. 10: MDS of the prac.se.sub dataset', fig.height=8----
type <- as.integer(factor(prac.se.sub$type))
plotMDS(prac.dge.unique.filtlib, pch = 19, col=type)
legend("topright", paste("Type", sort(unique(type)), levels(factor(prac.se.sub$type))),
       fill=sort(unique(type)), inset=0.05)

## ----function to prune, echo= FALSE--------------------------------------
GetID <- function(number, se){
  id <- grep(number,rownames(colData(se)))
  return(rownames(colData(prac.se.sub))[id])
}


## ----pruning the data set------------------------------------------------
number.id <- c("7737","7745","7740","7738","7211","7747","8258") #from MDS plot
pruning.vec <- sapply(number.id, GetID, prac.se.sub)

prac.dge.unique.pruned <- prac.dge.unique.filtlib[,
                   !rownames(prac.dge.unique.filtlib$samples) %in% pruning.vec]
table(prac.dge.unique.pruned$samples$group)
prac.se.pruned <- prac.se.sub[,rownames(prac.dge.unique.pruned$samples)]
v <- voom(dge, design, plot=TRUE)
## ----MDS after prunning, echo=FALSE--------------------------------------
tss <- substr(colnames(prac.se.pruned), 6, 7)
batch <- as.integer(factor(tss))
outcome <- paste(substr(colnames(prac.se.pruned), 9, 12), as.character(prac.se.pruned$type), sep="-")
names(outcome) <- colnames(prac.se.pruned)

## ----mds3, echo=FALSE, fig.align='center', fig.width=14, dpi=100, fig.cap='Fig. 11: MDS of the prac.se.sub dataset', fig.height=8----
plotMDS(prac.dge.unique.pruned, labels=outcome, col=batch)
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))),
       fill=sort(unique(batch)), inset=0.05)



## ----mds4, echo=FALSE, fig.align='center', fig.width=14, dpi=100, fig.cap='Fig. 12: MDS of the prac.se.sub dataset', fig.height=8----

type <- as.integer(factor(prac.se.pruned$type))
plotMDS(prac.dge.unique.pruned, pch = 19, col=type)
legend("topright", paste("Type", sort(unique(type)), levels(factor(prac.se.sub$type))),
       fill=sort(unique(type)), inset=0.05)

## ----diferential gene expr1----------------------------------------------
mod <- model.matrix(~ prac.se.pruned$type, colData(prac.se.pruned))
mod0 <- model.matrix(~ 1, colData(prac.se.pruned))
pv <- f.pvalue(assays(prac.se.pruned)$logCPM, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.01)

## ----pdist, echo=FALSE, out.width="600px", fig.cap="Fig. 13: Distribution of raw p-values for an F-test on every gene between tumor and normal samples.", fig.align='center'----


hist(pv, main="", xlab= 'p-value', las=1)


## ----sva 1---------------------------------------------------------------
sv <- sva(assays(prac.se.pruned)$logCPM, mod, mod0)
sv$n

## ----sva2----------------------------------------------------------------
modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(prac.se.pruned)$logCPM, modsv, mod0sv)
sum(p.adjust(pvsv, method="fdr") < 0.01)

## ----psvdist, echo=FALSE, out.width="600px", fig.cap="Fig. 14: Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA.", fig.align='center'----
hist(pvsv, main="", xlab= 'p-value', las=1)

## ----sva tss-------------------------------------------------------------
# Introduce the extracted data in the se obj
colData(prac.se.pruned)$tissue_source_site <- tss
table(data.frame(TYPE=prac.se.pruned$type, TSS = tss))
mod <- model.matrix(~ type + as.factor(tissue_source_site), data = colData(prac.se.pruned))
mod0 <- model.matrix(~ as.factor(tissue_source_site), data = colData(prac.se.pruned))
pValues <- f.pvalue(assays(prac.se.pruned)$logCPM, mod, mod0) 
sum(p.adjust(pValues, method = "BH") < 0.01)

## ----pvales plot tss, echo=FALSE, fig.align='center', fig.cap= 'Fig. 15: p-values of the differential expression analysis using TSS as null model.'----
hist(pValues, xlab= 'p-value')

## ----sva plate-----------------------------------------------------------
# introduce the extracted data in the se obj
plate <- substr(colnames(prac.se.pruned), 22, 25)
colData(prac.se.pruned)$plate <- plate
mod <- model.matrix(~ type + as.factor(plate), data = colData(prac.se.pruned))
mod0 <- model.matrix(~ as.factor(plate), data = colData(prac.se.pruned))
pValues <- f.pvalue(assays(prac.se.pruned)$logCPM, mod, mod0) 
sum(p.adjust(pValues, method = "BH") < 0.01)

## ----pvales plot p, echo=FALSE, fig.align='center', fig.cap= 'p-values of the differential expression analysis using plate as null model.'----
hist(pValues, xlab= 'p-value')

## ------------------------------------------------------------------------
sessionInfo()

### DE analysis

### How to basic
assays(prac.se.sub)$logCPM <- cpm(prac.dge.unique.filtlib, log=TRUE, prior.count=30)

design <- model.matrix(~ type, data=colData(prac.se.sub)) 
fit <- lmFit(assays(prac.se.sub)$logCPM, design)
fit <- eBayes(fit) 
tt <- topTable(fit, coef=2, n=Inf)
par(mfrow=c(1,1))
qqt(fit$t[, 2], df=fit$df.prior+fit$df.residual, main="", pch=".", cex=3, ylim= c(-50,50)) 
abline(0,1)


#################


# Modifying NAs values adding gleason score 2 as 'normal'.

na.mask <- is.na((colData(prac.se.pruned)$gleason_score))
colData(prac.se.pruned)$gleason_score <- as.character(colData(prac.se.pruned)$gleason_score)
colData(prac.se.pruned)$gleason_score[na.mask] <- "2"
colData(prac.se.pruned)$gleason_score <- as.factor(colData(prac.se.pruned)$gleason_score)

gs.mask <-  (colData(prac.se.pruned)$gleason_score %in% c("6","9"))
             
#------LETS APPLY the masking -------

prac.se.pruned <- prac.se.pruned[,gs.mask]
colData(prac.se.pruned)$gleason_score <- droplevels(colData(prac.se.pruned)$gleason_score)
design <- model.matrix(~ gleason_score, data=colData(prac.se.pruned)) 

head(design)
fit <- lmFit(assays(prac.se.pruned)$logCPM, design)
names(fit)
fit <- eBayes(fit) 
class(fit)

names(fit)
dim(prac.dge.unique.pruned)

res <- decideTests(fit)
summary(res)

FDRcutoff <- 0.01
res <- decideTests(fit, p.value=FDRcutoff) 
summary(res)
tt <- topTable(fit, coef=2, n=Inf) 
class(tt)
head(tt)

library("SummarizedExperiment")
genesmd <- data.frame(chr=as.character(seqnames(rowRanges(prac.se.pruned))),
                      symbol=as.character(rowRanges(prac.se.pruned)[, 1]), stringsAsFactors=FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef=2, n=Inf) 
head(tt)
chr.distro.df <- data.frame(chr = tt$chr[tt$adj.P.Val < FDRcutoff])
chr.distro.df$chr <- sort(chr.distro.df$chr)
ggplot(chr.distro.df, aes(x=chr, fill=chr)) + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

DEgenes <- rownames(tt)[tt$adj.P.Val < FDRcutoff]
length(DEgenes)



# to diagnose the DE

par(mfrow=c(1,1), mar=c(4, 5, 2, 2))
hist(tt$P.Value, xlab="Raw P-values", main="", las=1)
qqt(fit$t[, 2], df=fit$df.prior+fit$df.residual, main="", pch=".", cex=3, ylim=c(-20,20)) 
abline(0, 1, lwd=2)
qqline(tt$P.Value, col = 2,lwd=2,lty=2)# ojo! Sembla que son diferents amb abline

# VOOM
gs.mask <-  (colData(prac.se.pruned)$gleason_score %in% c("6","9"))
prac.dge.unique.pruned <- prac.dge.unique.pruned[,gs.mask]

par(mfrow=c(1,1))
v <- voom(prac.dge.unique.pruned, design, plot=TRUE)##
dim(prac.dge.unique.pruned)
dim(design)

# redo the process with the voom weight

fit2 <- lmFit(v, design)
fit2 <- eBayes(fit2)
res2 <- decideTests(fit2, p.value=FDRcutoff)

fit2$genes <- genesmd
tt2 <- topTable(fit2, coef=2, n=Inf)
head(tt2)
chr.distro2.df <- data.frame(chr = tt2$chr[tt2$adj.P.Val < FDRcutoff])
chr.distro2.df$chr <- sort(gsub("chr","",chr.distro2.df$chr))
ggplot(chr.distro2.df, aes(x=chr, fill=chr)) + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

par(mfrow=c(1,2), mar=c(4, 5, 2, 2))
hist(tt2$P.Value, xlab="Raw P-values", main="", las=1)
qqt(fit2$t[, 2], df=fit2$df.prior+fit2$df.residual, main="", pch=".", cex=3)
abline(0, 1, lwd=2)
qqline(fit2$t[, 2], col = 2,lwd=2,lty=2)# ojo! Sembla que son diferents amb abline

par(mfrow=c(1,1))


head(colData(prac.se.pruned), n=8)



## The results are shit

## We should get them working with covariates


design <- model.matrix(~ factor(type) + gleason_score, data=colData(prac.se.pruned))
head(design)
dim(v)
dim(design)


fit3 <- lmFit(v, design)
fit3 <- eBayes(fit3)
tt3 <- topTable(fit3, coef=2, n=Inf)


par(mfrow=c(1,2), mar=c(4, 5, 2, 2))
hist(tt3$P.Value, xlab="Raw P-values", main="", las=1)
qqt(fit3$t[, 2], df=fit3$df.prior+fit3$df.residual, main="", pch=".", cex=3)
abline(0, 1, lwd=2)
qqline(fit3$t[, 2], col = 2,lwd=2,lty=2)# ojo! Sembla que son diferents amb abline

par(mfrow=c(1,1))


library(sva)
mod0 <- model.matrix(~ 1, colData(prac.se.pruned))
sv <- sva(v$E, mod=design, mod0=mod0)

sv$n

design <- cbind(design, sv$sv)
colnames(design) <- c(colnames(design)[1:2], paste0("SV", 1:sv$n))

fit4 <- lmFit(v, design)

fit4 <- eBayes(fit4)

res4 <- decideTests(fit4, p.value=FDRcutoff)

summary(res4)


fit4$genes <- genesmd
tt4 <- topTable(fit4, coef=2, n=Inf) 
head(tt4)

sort(table(tt4$chr[tt4$adj.P.Val < FDRcutoff]), decreasing=TRUE)

par(mfrow=c(1,2), mar=c(4, 5, 2, 2))
hist(tt4$P.Value, xlab="Raw P-values", main="", las=1)
qqt(fit4$t[, 2], df=fit4$df.prior+fit4$df.residual, main="after SV", pch=".", cex=3) 
qqline(fit4$t[, 2], col = 2,lwd=2,lty=2)# ojo! Sembla que son diferents amb abline
abline(0, 1, lwd=2)

par(mfrow=c(1,2), mar=c(4, 5, 2, 2))
qqt(fit4$t[, 2], df=fit4$df.prior+fit4$df.residual, main="after SV", pch=".", cex=3,ylim=c(-20,20)) 
abline(0, 1, lwd=2)
qqt(fit4$t[, 2], df=fit4$df.prior+fit4$df.residual, main="original", pch=".", cex=3,ylim=c(-20,20)) 
abline(0, 1, lwd=2)



library(SummarizedExperiment)
lclse <- readRDS("data/sePRAD.rds")

library(edgeR)
dge <- DGEList(counts=assays(lclse)$counts, group=lclse$type, genes=mcols(lclse))
dge <- calcNormFactors(dge)
assays(lclse)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
mask <- rowMeans(assays(lclse)$logCPM) > 1
sum(mask)
lclse <- lclse[mask, ]
dim(lclse)
dge <- dge[mask, ]
dim(dge)
design <- model.matrix(~ type, data=colData(lclse))
head(design)
fit <- lmFit(assays(lclse)$logCPM, design)
class(fit)
names(fit)
fit <- eBayes(fit)
class(fit)
names(fit)
res <- decideTests(fit)
summary(res) #always 3 rows
FDRcutoff <- 0.01
res <- decideTests(fit, p.value=FDRcutoff)
summary(res)
tt <- topTable(fit, coef=2, n=Inf)
genesmd <- data.frame(chr=as.character(seqnames(rowRanges(lclse))),
                      symbol=as.character(rowRanges(lclse)[, 1]),
                      stringsAsFactors=FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef=2, n=Inf)
head(tt)
sort(table(tt$chr[tt$adj.P.Val < FDRcutoff]), decreasing=TRUE)


par(mfrow=c(1,2), mar=c(4, 5, 2, 2))
hist(tt$P.Value, xlab="Raw P-values", main="", las=1)
qqt(fit$t[, 2], df=fit$df.prior+fit$df.residual, main="", pch=".", cex=3)
abline(0, 4, lwd=2)

