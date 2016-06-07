#Prova seleccio atributs

library("SummarizedExperiment")
library("edgeR")
library("geneplotter")
library("ggplot2")
library("sva")
library("limma")
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

prac.se.pruned <- prac.se.sub[,rownames(prac.dge.unique.pruned$samples)]

logCPM <- cpm(prac.dge.unique.filtlib, log=TRUE, prior.count=3)

##### ALL NA's at Gleason Score are NORMAL SAMPLES ########




# normal <- (colData(prac.se.sub)$type[colData(prac.se.sub)$type == 'normal'])

# tumor <- (colData(prac.se.sub)$type[colData(prac.se.sub)$type == 'tumor'])

# names(gleason.score.na) %in% names(normal)

# names(gleason.score.na) %in% names(tumor)

na.mask <- is.na((colData(prac.se.sub)$gleason_score))
colData(prac.se.sub)$gleason_score <- as.character(colData(prac.se.sub)$gleason_score)
colData(prac.se.sub)$gleason_score[na.mask] <- "2"
colData(prac.se.sub)$gleason_score <- as.numeric(colData(prac.se.sub)$gleason_score)

###########################################################

gs.mask <-  (colData(prac.se.pruned)$gleason_score %in% c("6","9"))

#------LETS APPLY the masking -------
dim(prac.se.pruned)
prac.se.pruned <- prac.se.pruned[,gs.mask]
dim(prac.se.pruned)
colData(prac.se.pruned)$gleason_score <- droplevels(colData(prac.se.pruned)$gleason_score)
design <- model.matrix(~ gleason_score, data=colData(prac.se.pruned)) 
dim(design)
fit <- lmFit(assays(prac.se.pruned)$logCPM, design)
names(fit)
fit <- eBayes(fit) 
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

# to diagnose the DE

par(mfrow=c(1,1), mar=c(4, 5, 2, 2))
hist(tt$P.Value, xlab="Raw P-values", main="", las=1)
qqt(fit$t[, 2], df=fit$df.prior+fit$df.residual, main="", pch=".", cex=3, ylim=c(-20,20)) 
abline(0, 1, lwd=2)
qqline(tt$P.Value, col = 2,lwd=2,lty=2)# ojo! Sembla que son diferents amb abline

# VOOM
gs.mask <-  (colData(prac.se.pruned)$gleason_score %in% c("6","9"))
dim(prac.dge.unique.pruned)
prac.dge.unique.pruned <- prac.dge.unique.pruned[,gs.mask]
dim(prac.dge.unique.pruned)

design <- model.matrix(~ gleason_score ,data=colData(prac.se.pruned)) 
dim(design)
design
par(mfrow=c(1,1))
v <- voom(prac.dge.unique.pruned, design, plot=TRUE)##
dim(design)
dim(prac.dge.unique.pruned)

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


design_mix <- model.matrix(~ type + gleason_score, data=colData(prac.se.pruned))
length(colData(prac.se.pruned))
design_mix
dim(v)
dim(design_mix)


fit3 <- lmFit(design_mix)

fit3 <- eBayes(fit3) #can't calculate eBayes, too much 0s
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

