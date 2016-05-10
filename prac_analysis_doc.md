# Analysis of a TCGA RNA-seq data set on Prostate Adenocarcinoma 
[Joan Martí](mailto:joan.marti02@estudiant.upf.edu), [David Mas](mailto:david.mas01@estudiant.upf.edu), [Adrià Auladell](mailto:adria.auladell01@estudiant.upf.edu)  
10th May 2016  





# Introduction

[Prostate cancer](https://en.wikipedia.org/wiki/Prostate_cancer)
is a disease of the prostate, a walnut-size gland in the male
reproductive system.  Nearly all prostate cancer is prostate adenocarcinoma.
 Prostate cancer is graded based on its Gleason score, which is how the cells
 look under a microscope and ranges from two to ten. A low Gleason score means
 that the cancer tissue is similar to normal cells and unlikely to spread. A
 high Gleason score means that the cancer cells are very different from normal
 cells and are likely to spread. For patients whose cancer has spread, their
 survival time is usually one to three years. It was estimated that for 2011,
 240,890 men would be diagnosed and 33,720 would die from prostate cancer.

The Cancer Genome Atlas already studied this dataset for a wide characterization of the prostate adenocarcinoma subtypes. TCGA established the following points:

* 74% of all tumors being assignable to one of seven molecular classes based on distinct oncogenic drivers: 
    - fusions involving (1) **ERG**, (2) **ETV1**, (3) **ETV4**, or (4) **FLI1** (46%, 8%, 4%, and 1%, respectively)
    - mutations in (5) **SPOP** or (6) **FOXA1**; or (7) **IDH1**. (11%, 3%, and 1%, respectively).
    
* 25% of the prostate cancers had a presumed actionable lesion in the PI3K or MAPK signaling pathways, and DNA repair genes were inactivated in 19%.

You can find the full article [here](http://dx.doi.org/10.1016/j.cell.2015.10.025).

# Data import

## Package Installation and Import
First, we install the packages using BiocInstaller:


```r
source("http://www.bioconductor.org/biocLite.R")
library(BiocInstaller)
biocLite(c("edgeR","SummarizedExperiment","geneplotter","sva"))
install.packages(ggplot2)
```





```r
library("SummarizedExperiment")
library("edgeR")
library("geneplotter")
library("ggplot2")
library("sva")
old.par <- par(mar = c(0, 0, 0, 0)) # to get CLI 
```

## Data Reading
After installing and loading the packages, we use them to import the [dataset](http://functionalgenomics.upf.edu/courses/IEO/projects/datasets/sePRAD.rds). This dataset is provided by our professor [Robert Castelo](mailto:robert.castelo@upf.edu). 

The data is contained in a rds object and when reading it by the function `readRDS` we extract the `Summarized Experiment`object that contains the data we are going to use. This class is from the `SummarizedExperiment` package loaded above. 


```r
prac.se <- readRDS("data/sePRAD.rds")
prac.se
```

```
class: RangedSummarizedExperiment 
dim: 20115 554 
metadata(5): experimentData annotation cancerTypeCode
  cancerTypeDescription objectCreationDate
assays(1): counts
rownames(20115): 1 2 ... 102724473 103091865
rowRanges metadata column names(3): symbol txlen txgc
colnames(554): TCGA.2A.A8VL.01A.21R.A37L.07
  TCGA.2A.A8VO.01A.11R.A37L.07 ... TCGA.HC.8262.11A.01R.2263.07
  TCGA.J4.A83J.11A.11R.A36G.07
colData names(549): type bcr_patient_uuid ...
  lymph_nodes_aortic_pos_by_ihc lymph_nodes_aortic_pos_total
```

## Nomenclature

In the this and following reports we are going to try naming variables and objects using a nomenclature. We believe this nomenclature follows the Google's R Style Guidelines. As recomended by this guidelines, different parts of the variable name are separated by points. 

* The first part of the variable defines what it contains in a clear as possible. 
* Second part defines the `class` of the object. This is missed when the object have a common class (such as character vectors) or it is a temporary variable. 
* Third pard defines the filters applied to the original object. 

Based on this nomenclature, below you can find variables like `<name>.<type_var>.<other>`:

```
prac.se 
# where prac is the item it represent (cancer type) and
# se is the class (summarized Experiment)
prac.se.unique 
# where a third labeled is applied in order to note that the data have 
# been filtred. 
```
We also try to mantain our code chunks below 80 chrs.

## Exploring the DataSet

As defined in the `SummarizedExperiment` documentation there are 3 data types in a SE. The **Column data** contains info about the samples used in the experiment and it is refered as column data because represent the columns of the original expression table. Then we have the **metadata** that brings us information about the assays that were performed. Finally we got the **row data** that contains the gene information. This is called row data because in the expression table genes are in the rows. 

### Data Size

```r
dim(colData(prac.se))
```

```
[1] 554 549
```

Our dataset presents 554 different samples with 549 clinical variables analyzed. 


### Column Data

In order to acces the column data from our `SE`object we can use its own implemented function, `colData()`.


```r
head(names(colData(prac.se)))
```

```
[1] "type"                     "bcr_patient_uuid"        
[3] "bcr_patient_barcode"      "form_completion_date"    
[5] "prospective_collection"   "retrospective_collection"
```

```r
colData(prac.se)[1:5, 1:5]
```

```
DataFrame with 5 rows and 5 columns
                                 type                     bcr_patient_uuid
                             <factor>                             <factor>
TCGA.2A.A8VL.01A.21R.A37L.07    tumor 49197847-CC83-4CE1-8397-D09CEA4C4928
TCGA.2A.A8VO.01A.11R.A37L.07    tumor 91C0D161-2B59-4B7A-8C19-6D26DEA31849
TCGA.2A.A8VT.01A.11R.A37L.07    tumor 931B549F-B9F2-4E8D-83ED-FF663671883C
TCGA.2A.A8VV.01A.11R.A37L.07    tumor 75A7AFB5-66D5-47E3-8A8A-3E3A1E749A96
TCGA.2A.A8VX.01A.11R.A37L.07    tumor 942F1788-D977-4AC0-A177-659F9D4CD077
                             bcr_patient_barcode form_completion_date
                                        <factor>             <factor>
TCGA.2A.A8VL.01A.21R.A37L.07        TCGA-2A-A8VL            2014-3-29
TCGA.2A.A8VO.01A.11R.A37L.07        TCGA-2A-A8VO            2014-3-30
TCGA.2A.A8VT.01A.11R.A37L.07        TCGA-2A-A8VT            2014-3-29
TCGA.2A.A8VV.01A.11R.A37L.07        TCGA-2A-A8VV            2014-3-29
TCGA.2A.A8VX.01A.11R.A37L.07        TCGA-2A-A8VX            2014-3-29
                             prospective_collection
                                           <factor>
TCGA.2A.A8VL.01A.21R.A37L.07                     NO
TCGA.2A.A8VO.01A.11R.A37L.07                     NO
TCGA.2A.A8VT.01A.11R.A37L.07                     NO
TCGA.2A.A8VV.01A.11R.A37L.07                     NO
TCGA.2A.A8VX.01A.11R.A37L.07                     NO
```

```r
mcols(colData(prac.se), use.names=TRUE)
```

```
DataFrame with 549 rows and 2 columns
                                                         labelDescription
                                                              <character>
type                                           sample type (tumor/normal)
bcr_patient_uuid                                         bcr patient uuid
bcr_patient_barcode                                   bcr patient barcode
form_completion_date                                 form completion date
prospective_collection            tissue prospective collection indicator
...                                                                   ...
lymph_nodes_pelvic_pos_total                               total pelv lnp
lymph_nodes_aortic_examined_count                           total aor lnr
lymph_nodes_aortic_pos_by_he                          aln pos light micro
lymph_nodes_aortic_pos_by_ihc                                 aln pos ihc
lymph_nodes_aortic_pos_total                                total aor-lnp
                                        CDEID
                                  <character>
type                                       NA
bcr_patient_uuid                           NA
bcr_patient_barcode                   2673794
form_completion_date                       NA
prospective_collection                3088492
...                                       ...
lymph_nodes_pelvic_pos_total          3151828
lymph_nodes_aortic_examined_count     3104460
lymph_nodes_aortic_pos_by_he          3151832
lymph_nodes_aortic_pos_by_ihc         3151831
lymph_nodes_aortic_pos_total          3151827
```

By observing the clinical variables, from the 549 present some have information collected and others present NA values. Our interest is in the informative variables.Some of them are useful to study the data, like bcr\_patient\_barcode, race, ethnicity, age\_at\_initial\_pathologic\_diagnosis...There are variables related directly with prostate cancer, like gleason score, the measurement of the development of the adenocarcinoma.

<img src="prac_analysis_doc_files/figure-html/exploring some variables-1.png" title="" alt="" style="display: block; margin: auto;" />

Looking at the metadata, we observe the common structure followed in the datasets of the TCGA project. The first column is the clinical variable abreviated, the second a explanation of the variable and the last one contains the  [CDEID](https://www.nlm.nih.gov/cde/glossary.html#cdedefinition ) code to obtain more information in the [CDE website](https://www.nlm.nih.gov/cde). 

As an example we can see in the figures, there are clinical variables and sample variables giving information about the sampels. 

### Tumor / Normal Samples

```r
table(prac.se$type)
```

```

normal  tumor 
    52    502 
```

The dataset presents 554 samples. There are 52 labeled as normal samples and 504 labeled as tumor samples. The first ones are not from healthy individuals but from affected ones. This mean that normal samples do not have tumor phenotype but the same individual background than the tumor samples. In fact, most of the normal samples (up to 50) have a paired sample (2 samples from the same individual) in the tumor side of the table. We can check this information using the BCR patient barcode that identifies individuals. 


### Row Data

We can get information about the genes using the following. 


```r
rowRanges(prac.se)
```

```
GRanges object with 20115 ranges and 3 metadata columns:
            seqnames               ranges strand   |      symbol     txlen
               <Rle>            <IRanges>  <Rle>   | <character> <integer>
          1    chr19 [58345178, 58362751]      -   |        A1BG      3322
          2    chr12 [ 9067664,  9116229]      -   |         A2M      4844
          9     chr8 [18170477, 18223689]      +   |        NAT1      2280
         10     chr8 [18391245, 18401218]      +   |        NAT2      1322
         12    chr14 [94592058, 94624646]      +   |    SERPINA3      3067
        ...      ...                  ...    ... ...         ...       ...
  100996331    chr15 [20835372, 21877298]      -   |       POTEB      1706
  101340251    chr17 [40027542, 40027645]      -   |    SNORD124       104
  101340252     chr9 [33934296, 33934376]      -   |   SNORD121B        81
  102724473     chrX [49303669, 49319844]      +   |      GAGE10       538
  103091865    chr21 [39313935, 39314962]      +   |   BRWD1-IT2      1028
                 txgc
            <numeric>
          1 0.5644190
          2 0.4882329
          9 0.3942982
         10 0.3895613
         12 0.5249429
        ...       ...
  100996331 0.4308324
  101340251 0.4903846
  101340252 0.4074074
  102724473 0.5055762
  103091865 0.5924125
  -------
  seqinfo: 455 sequences (1 circular) from hg38 genome
```

For each gene we have present the chromosome, the range of position, the strand, its symbol, the length of the transcript and the GC content. 
<img src="prac_analysis_doc_files/figure-html/gc content hist-1.png" title="" alt="" style="display: block; margin: auto;" />

### Metadata
From the metadata we can check that we have the current version of the data. 

```r
metadata(prac.se)$objectCreationDate #the result should be [1] "2016-04-25"
```

```
[1] "2016-04-25"
```

# Subsetting the Data

In order to reduce possible batch effects and reduce the computational costs of the analysis we have decided to subset our data. We will use several approaches in order to do perform this action. Finally, we will apply one or other subset if we want to avoid the disadvantages of a defined strategy.

## **Paired Subsetting** 

In the paired subsetting we try to get only the patients (based on the bcr barcode) that have a normal sample and a tumor one. This approach is suitable to distinguish exactly the gene expression changes that drive the tumour. However, it also drives some disadvantages as we have been taugh in class.

* As seen in the resulting numbers, we have filter out 454 patients that were not paired. We have same number of normal samples and tumor samples (50) because for each individual we have 2 samples, each one in a different group.* 


```r
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
```

```

normal  tumor 
    50     50 
```
 
## **Non-Paired Subsetting** 

With this strategy we want to discard all the paired data in order to obtain only unique individuals. Our aim is to only get the normal samples of the replicated individuals. As the initial filter strategy get rid of most of the normal samples and kept the tumor ones we had to adapt it and generate a more complex filter as you can see below. 

* The replicated samples are 2 samples extracted from the same individual in 2 different conditions (tumor and normal). We want to distinguish general variation that is associated to one or the other subset and introducing 2 samples for the same individual could include batch effect in our samples.

* In the resuting table we can see that all the normal samples remain meanwhile we have filter out 67 tumor samples that contained `NA`values in the bcr barcode field or that they were extracted from the same individuals that provided the normal samples.
    

```r
# Subset the non-paired 
prac.se.unique <- prac.se[, colData(prac.se)$bcr_patient_barcode 
                             %in% n_occur$Var1[n_occur$Freq == 1] & 
                                !is.na(colData(prac.se)$bcr_patient_barcode) |
                                (colData(prac.se)$bcr_patient_barcode %in% 
                                n_occur$Var1[n_occur$Freq > 1] & 
                                colData(prac.se)$type == "normal" & 
                                !is.na(colData(prac.se)$bcr_patient_barcode))]
table(colData(prac.se.unique)$type)
```

```

normal  tumor 
    52    435 
```
  

# Quality data control

Having understood completely our dataset, to perform efficiently quality assesment and normalization of the data, we have to convert our values to counts per million values (CPM). In order to do that, the package edgeR is used, creating a 'DGEList Object'. 
Digital Gene Expression data class (DGE) has been implemented for storing read counts and associated information from digital gene expression or sequencing technologies. The user specifies the counts, the samples, and optional components include the genes and the groups.


```r
prac.dge <- DGEList(counts = assays(prac.se)$counts, genes = mcols(prac.se), group = prac.se$type)
head(prac.dge$samples)
```

```
                             group lib.size norm.factors
TCGA.2A.A8VL.01A.21R.A37L.07 tumor 36714836            1
TCGA.2A.A8VO.01A.11R.A37L.07 tumor 26883656            1
TCGA.2A.A8VT.01A.11R.A37L.07 tumor 42636097            1
TCGA.2A.A8VV.01A.11R.A37L.07 tumor 39132257            1
TCGA.2A.A8VX.01A.11R.A37L.07 tumor 36293081            1
TCGA.2A.A8W1.01A.11R.A37L.07 tumor 36955738            1
```

Finally, from each CPM value we calculate the $\log_2$ measure and we include it in our dataset in order to use it in the following comparisons/normalizations. 


```r
assays(prac.se)$logCPM <- cpm(prac.dge, log=TRUE, prior.count=3.5)
assays(prac.se)$logCPM[1:5, 1:5]
```

```
   TCGA.2A.A8VL.01A.21R.A37L.07 TCGA.2A.A8VO.01A.11R.A37L.07
1                      2.030317                   0.01585209
2                      8.111182                  10.23848110
9                     -3.623878                  -3.62387802
10                    -3.623878                  -3.62387802
12                     5.322306                   5.50185474
   TCGA.2A.A8VT.01A.11R.A37L.07 TCGA.2A.A8VV.01A.11R.A37L.07
1                     0.3793317                   -0.1133961
2                     8.3965136                    8.9351599
9                    -3.6238780                   -3.6238780
10                   -3.6238780                   -3.6238780
12                    5.7006877                    5.6144791
   TCGA.2A.A8VX.01A.11R.A37L.07
1                      1.065240
2                      8.936584
9                     -3.623878
10                    -3.623878
12                     5.682520
```

## Filtering by library size
Now we subset our dataset by library size, which is a measure of sequencing deepness or how robust are the RNA-seqs by samples. We renormalize the data with the new `SE`filtred element. 


```r
prac.dge.unique <- DGEList(counts = assays(prac.se.unique)$counts, genes = mcols(prac.se.unique), group = prac.se.unique$type)
```

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/library size all after filter unique-1.png" alt="Fig. 1: Library Size per sample in the  non-paired design dataset. "  />
<p class="caption">Fig. 1: Library Size per sample in the  non-paired design dataset. </p>
</div><p class="caption">Fig. 1: Library Size per sample in the  non-paired design dataset. </p>

In the fig. 1 , we cannot distinguish the different samples due the large number of samples. We will filter in a high treshold of sequencing depth to recude the number of samples erasing the ones with less sequencing quality. 



```r
# library size, non-paired samples
prac.dge.unique.filtlib <- prac.dge.unique[,(prac.dge.unique$samples$lib.size/1e6) > 40 ] # 50e06 threshold value
table(prac.dge.unique.filtlib$samples$group)
```

```

normal  tumor 
    49    224 
```
<!--- Adri START --->


Since we still have a lot of samples, we look at the distribution of the Tissue Soucre Site (TSS) for selecting the batches with a good distribution of normal/tumor samples in order to obtain a robust dataset.



```r
# not filtering by TSS
prac.se.sub <- prac.se.unique[,rownames(prac.dge.unique.filtlib$samples)]

# filter by TSS block
tss <- substr(colnames(prac.se.sub), 6, 7)
tss_df <- data.frame(TYPE=prac.se.sub$type, TSS=tss)
table(data.frame(TYPE=prac.se.sub$type, TSS=tss))
```

```
        TSS
TYPE     2A 4L CH EJ FC G9 H9 HC HI J4 KC KK M7 QU V1 VN VP XJ XK YJ YL ZG
  normal  0  0  2 22  0 11  0 13  0  1  0  0  0  0  0  0  0  0  0  0  0  0
  tumor   4  1 18 46  2 26  2 32  3 15  5 19  3  1 11  3  6  1  5  1 11  9
```


We can see some TSS presenting only tumor samples with a low number. These type of samples can be a potential batch factor, so we choose a selection of balanced samples (for example CH, EJ, G9 and HC).


```r
selection  <- c('CH', 'EJ', 'G9', 'HC')
prac.dge.unique.filtlib <- prac.dge.unique.filtlib[,tss_df$TSS %in% selection]
prac.se.sub <- prac.se.unique[,rownames(prac.dge.unique.filtlib$samples)]
table(prac.dge.unique.filtlib$samples$group)
```

```

normal  tumor 
    48    122 
```

 <!--- Adri END --->


```
    libsize          type                             samples   
 Min.   :40.05   normal: 48   TCGA.CH.5737.01A.11R.1580.07:  1  
 1st Qu.:45.10   tumor :122   TCGA.CH.5738.01A.11R.1580.07:  1  
 Median :50.64                TCGA.CH.5739.01A.11R.1580.07:  1  
 Mean   :52.89                TCGA.CH.5743.01A.21R.1580.07:  1  
 3rd Qu.:56.08                TCGA.CH.5744.01A.11R.1580.07:  1  
 Max.   :96.68                TCGA.CH.5745.01A.11R.1580.07:  1  
                              (Other)                     :164  
```

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/figure library size-1.png" alt="Fig. 2: Density sample distribution using already filtred data by sequencing depth."  />
<p class="caption">Fig. 2: Density sample distribution using already filtred data by sequencing depth.</p>
</div><p class="caption">Fig. 2: Density sample distribution using already filtred data by sequencing depth.</p>

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/barplot libsize samples filtred 2-1.png" alt="Fig. 3: Barplot of the final selected samples with its library size."  />
<p class="caption">Fig. 3: Barplot of the final selected samples with its library size.</p>
</div><p class="caption">Fig. 3: Barplot of the final selected samples with its library size.</p>

After filtering, we obtain a dataset of 48 normal samples and 122 tumor (fig. 2 and 3). We introduce it again in the prac.se dataset to obtain the final filtered version.
Bear in mind that it seems there is an outlayer, a normal sample with an unusual or unexpected good quality. The converage is unusually high but since it is theorically a good attribute, we decided to keep it. 

## Sample expression distribution

Once we have our filtered dataset, we observe the density logCPM distribution of the tumor and normal samples separatedly. It will give us an impression of the possible problems of the sampling, since we expect to have a similar distribution of the samples between tumor / normal. 

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/Multydensity-1.png" alt="Fig 4: Multidensity plot representing the clusters of the non expressed genes at the left and the expressed genes at the right."  />
<p class="caption">Fig 4: Multidensity plot representing the clusters of the non expressed genes at the left and the expressed genes at the right.</p>
</div><p class="caption">Fig 4: Multidensity plot representing the clusters of the non expressed genes at the left and the expressed genes at the right.</p>

Analysing the density graphs from fig. 4, we cannot establish differences between the tumoral/normal logCPM distribution.

## Gene expression distribution

Ending the sample distribution analysis, we can observe the distribution of the log2CPM by gene. We will erase the genes presenting a logCPM below 1, since for these values the CPM the expression level is very low, we cannot establish real conclusions for that genes (spontaneus transcription, non-specific of tissue transcription...).



```r
assays(prac.se.sub)$logCPM <- cpm(prac.dge.unique.filtlib, log=TRUE, prior.count=0.5)

avgexp <- rowMeans(assays(prac.se.sub)$logCPM)

hist(avgexp, xlab="log2 CPM", main="", las=1)

abline(v=1, col="red", lwd=2)
```

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/avgexp gene-1.png" alt="Fig. 5: Histogram presenting gene frequency expression for different logCPM values"  />
<p class="caption">Fig. 5: Histogram presenting gene frequency expression for different logCPM values</p>
</div><p class="caption">Fig. 5: Histogram presenting gene frequency expression for different logCPM values</p>
As can be observed from fig. 5, there are values for logCPM inferior to 1 (which is considered the minumum standard to consider from). Those values inferior to 1 should not be taken into account when performing the analysis, as those values fall within the Grey Zone (difficulties for differencing data artifacts from the techniques and data obtained from the samples).

## Filtering by gene expression levels

In order to eliminate those values from our, a mask is creating, using as a factor of exclusiong avergae of expression higher than 1. 


```r
mask <- avgexp > 1
dim(prac.se.sub)
```

```
[1] 20115   170
```

```r
prac.se.sub <- prac.se.sub[mask, ]
dim(prac.se.sub)
```

```
[1] 11722   170
```

```r
dim(prac.dge.unique.filtlib )
```

```
[1] 20115   170
```

```r
prac.dge.unique.filtlib <- prac.dge.unique.filtlib[mask, ]
dim(prac.dge.unique.filtlib)
```

```
[1] 11722   170
```


<!--- JOAN START --->
# Normalization

First of all, we calculate the normalized factors between sample and we introduce it to prac.se.sub dataset. 


```r
prac.dge.unique.filtlib <- calcNormFactors(prac.dge.unique.filtlib)
assays(prac.se.sub)$logCPM <- cpm(prac.dge.unique.filtlib, log=TRUE, prior.count=0.5)
```

From this values, we create an MA plot of each sample to observe the distribution of expression. Unusual samples will be filtered of the dataset later on. 
## MA-plots

<center><h4>Tumor samples MA-plot</h4></center>

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/maPlotsTumor-1.png" alt="Fig. 6: MA-plots of the tumor samples."  />
<p class="caption">Fig. 6: MA-plots of the tumor samples.</p>
</div><p class="caption">Fig. 6: MA-plots of the tumor samples.</p>

<center><h4>Normal samples MA-plot</h4></center>

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/maPlotsNormal-1.png" alt="Fig. 7: MA-plots of the tumor samples."  />
<p class="caption">Fig. 7: MA-plots of the tumor samples.</p>
</div><p class="caption">Fig. 7: MA-plots of the tumor samples.</p>
As can be percived in fig. 6 and 7, the tendency line representing the mean of normalised values (red line) is reasonably close to the expected tendency line (blue line), which is centered at mean 0.
Despite there are some samples with line tails diverging from the extected values, it is not a big deal, usually.
<!--- JOAN END --->

<!--- DAVID START --->


# Batch Effect

The batch effect occurs when there is a non-biological variable that accounts for the variability between samples. Usually this is related to handling, storing and processing the samples. This source of sample variability might obscure the analysis of variation and comparison between samples.

First of all, the batch labels need to be obtained from the sample codes of the `SE` object as the coldata column refering to it contains some NAs. The Tissue Source Site (TSS) is the center in which the samples where processed and the tissue where they were extracted. The correlation between TSS codes and centers can be found [here](https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm).

Our selected samples were collected and analysed in the University of Pittsburgh, the Roswell Park, the International Genomics Consortium and the company Indivumed.

We also extract from the sample code the plate, its portion of analyte and sample_vial.


```r
tss <- substr(colnames(prac.se.sub), 6, 7)
table(tss)
tss
CH EJ G9 HC 
20 68 37 45 
center <- substr(colnames(prac.se.sub), 27, 28)
table(center)
center
 07 
170 
plate <- substr(colnames(prac.se.sub), 22, 25)
table(plate)
plate
1580 1789 1858 1965 2118 2263 2403 A29R A30B A311 A31N A32O A33R A352 A36G 
  48   15    1   17   44   16    3    2    7    1    1    3    4    2    2 
A41O 
   4 
portionanalyte <- substr(colnames(prac.se.sub), 18, 20)
table(portionanalyte)
portionanalyte
01R 02R 11R 12R 13R 21R 22R 31R 
 72   4  80   3   1   7   1   2 
samplevial <- substr(colnames(prac.se.sub), 14, 16)
table(samplevial)
samplevial
01A 01B 11A 11B 
119   3  47   1 
```
<!--- ADRI START --->

By looking at the table between type and tss, we will observe if the data chosen presents batch effect. We have corrected this possibility at the first selections by taking only the biggets batches. 


```r
table(data.frame(TYPE=prac.se.sub$type, TSS=tss))
```

```
        TSS
TYPE     CH EJ G9 HC
  normal  2 22 11 13
  tumor  18 46 26 32
```

<!--- ADRI END --->

A good way to see if there is batch effect on the samples is to study the data clustering.

In the following dendogram, samples are being clustered by Gene Expression levels and colored by TSS.


```r
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
```

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/hirechical clustering-1.png" alt="Fig. 8: Hirechical clustering of the prac.se.sub dataset"  />
<p class="caption">Fig. 8: Hirechical clustering of the prac.se.sub dataset</p>
</div><p class="caption">Fig. 8: Hirechical clustering of the prac.se.sub dataset</p>
When in a dendogram some samles cluster together, it means that they are closer to each other, than the rest of the sample. If we expect no batch effect, the batch labels should be more or less randomly clustered together.

As we can see from the dendogram in fig. 8, some red samples (batch EJ) and blue samples (batch HC) seem to cluster together. This can happen as a result of 2 causes, which are non-mutualy exclusive:

* There are more red and blue samples, thus, it is more likely that they cluster together.
* Red and blue samples are similiar between them, and different between the other samples, hence, clustering together. 

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/mds-1.png" alt="Fig. 9: MDS of the prac.se.sub dataset"  />
<p class="caption">Fig. 9: MDS of the prac.se.sub dataset</p>
</div><p class="caption">Fig. 9: MDS of the prac.se.sub dataset</p>

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/mds2-1.png" alt="Fig. 10: MDS of the prac.se.sub dataset"  />
<p class="caption">Fig. 10: MDS of the prac.se.sub dataset</p>
</div><p class="caption">Fig. 10: MDS of the prac.se.sub dataset</p>

MDS-plots help to see also how data cluster together projecting the variability of each sample, in a Multiple Dimension Scaling plots.

As we  can see in fig. 9 and 10, there is a clustering of normal samples, belonging to batch HC (blue at left-inferior corner). For the sake of keeping the dataset as homogenous as possible, it might be interessting to do not use those samples. Maybe, that mini-cluster of samples is partially responsable of the overall disposition of samples (which differentates poorly between normal and cancer tissues).


<!--- DAVID END --->



<!--- DAVID START ---> 

## Filtering after MDS ands MA plots

As some unexpected results have been found, regarding normal samples, one viable option for the analysis is to remove those samples.




```r
number.id <- c("7737","7745","7740","7738","7211","7747","8258") #from MDS plot
pruning.vec <- sapply(number.id, GetID, prac.se.sub)

prac.dge.unique.pruned <- prac.dge.unique.filtlib[,
                   !rownames(prac.dge.unique.filtlib$samples) %in% pruning.vec]
table(prac.dge.unique.pruned$samples$group)
```

```

normal  tumor 
    41    122 
```

```r
prac.se.pruned <- prac.se.sub[,rownames(prac.dge.unique.pruned$samples)]
```

After pruning, we observe again the distribution of the samples (Fig.11). 



<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/mds3-1.png" alt="Fig. 11: MDS of the prac.se.sub dataset"  />
<p class="caption">Fig. 11: MDS of the prac.se.sub dataset</p>
</div><p class="caption">Fig. 11: MDS of the prac.se.sub dataset</p>

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/mds4-1.png" alt="Fig. 12: MDS of the prac.se.sub dataset"  />
<p class="caption">Fig. 12: MDS of the prac.se.sub dataset</p>
</div><p class="caption">Fig. 12: MDS of the prac.se.sub dataset</p>



As we can see in fig. 11 and 12, after filter those samples out, the distribution gains centrality and uniformity generating a more reliable data set. 



# Differential Gene Expression

After having corrected by those disturbing samples, the differential gene expression can be analysed. The null model created in this case is not adjusted for any attribute (mod0 with intercept 1).



```r
mod <- model.matrix(~ prac.se.pruned$type, colData(prac.se.pruned))
mod0 <- model.matrix(~ 1, colData(prac.se.pruned))
pv <- f.pvalue(assays(prac.se.pruned)$logCPM, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.01)
```

```
[1] 5610
```

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/pdist-1.png" alt="Fig. 13: Distribution of raw p-values for an F-test on every gene between tumor and normal samples." width="600px" />
<p class="caption">Fig. 13: Distribution of raw p-values for an F-test on every gene between tumor and normal samples.</p>
</div><p class="caption">Fig. 13: Distribution of raw p-values for an F-test on every gene between tumor and normal samples.</p>


The histogram presents a peak in the 0.0-0.1 region, indicating a high presence of DE genes. The distribution of the other values is quite uniform and low in comparison with the 0.1-0.0 values. 


In order to focus on the DE analysis we use only the SV affecting our gene expression variability. Using this adjustment we can see how the number increases significantly. This may be related to a unkown source of variability that we haven't been able to detect or to the intrinsec variability of individuals. 


```r
sv <- sva(assays(prac.se.pruned)$logCPM, mod, mod0)
```

```
Number of significant surrogate variables is:  20 
Iteration (out of 5 ):1  2  3  4  5  
```

```r
sv$n
```

```
[1] 20
```



```r
modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(prac.se.pruned)$logCPM, modsv, mod0sv)
sum(p.adjust(pvsv, method="fdr") < 0.01)
```

```
[1] 7145
```

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/psvdist-1.png" alt="Fig. 14: Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA." width="600px" />
<p class="caption">Fig. 14: Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA.</p>
</div><p class="caption">Fig. 14: Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA.</p>

The total of DE genes has increased a 27% with respect to the initial number (from fig. 14). The distribution of p-values has remained equal, having a decrease in the frequency of each interval (from 400 ocurrences to 300). The uniformity of the distribution is mantained. 

## Checking for batch effect from TSS

As it have not been 100% clear if there is batch effect at the previous step we will use again the sva package to stimate if there exists actual batch effect by using the TSS as a null model when building the `mod.matrix`. We will also perform the same analysis with the plate attribute of the samples, also to recheck for some batch effect. 

### TSS

The analysis is quite similar from the one performed above, however, here we consider the TSS variation as a null model. If the variation in the gene expression was explained by the variation in the TSS the p values will be uniform, meaning that there is not a peak in the significant genes. 

In the following figure we can observe that this is not true and the p-value distribution is mantained. Moreover, the decrease of DE genes between the uniform null model (5610) and the TSS null model (5546) is very small.


```r
# Introduce the extracted data in the se obj
colData(prac.se.pruned)$tissue_source_site <- tss
table(data.frame(TYPE=prac.se.pruned$type, TSS = tss))
```

```
        TSS
TYPE     CH EJ G9 HC
  normal  2 22 11  6
  tumor  18 46 26 32
```

```r
mod <- model.matrix(~ type + as.factor(tissue_source_site), data = colData(prac.se.pruned))
mod0 <- model.matrix(~ as.factor(tissue_source_site), data = colData(prac.se.pruned))
pValues <- f.pvalue(assays(prac.se.pruned)$logCPM, mod, mod0) 
sum(p.adjust(pValues, method = "BH") < 0.01)
```

```
[1] 5546
```

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/pvales plot tss-1.png" alt="Fig. 15: p-values of the differential expression analysis using TSS as null model."  />
<p class="caption">Fig. 15: p-values of the differential expression analysis using TSS as null model.</p>
</div><p class="caption">Fig. 15: p-values of the differential expression analysis using TSS as null model.</p>
As can be seen in fig. 15, the histogram shows a distribution similar to fig. 14 (high proportion of genes with p-values < 0.01 with respect to all the others).


### Plate
Here we perform the same analysis and we see that the number aso decrease also by a small extend. This would impliy that there not exists significant batch effect by this variable. 


```r
# introduce the extracted data in the se obj
plate <- substr(colnames(prac.se.pruned), 22, 25)
colData(prac.se.pruned)$plate <- plate
mod <- model.matrix(~ type + as.factor(plate), data = colData(prac.se.pruned))
mod0 <- model.matrix(~ as.factor(plate), data = colData(prac.se.pruned))
pValues <- f.pvalue(assays(prac.se.pruned)$logCPM, mod, mod0) 
sum(p.adjust(pValues, method = "BH") < 0.01)
```

```
[1] 5518
```

<div class="figure" style="text-align: center">
<img src="prac_analysis_doc_files/figure-html/pvales plot p-1.png" alt="p-values of the differential expression analysis using plate as null model."  />
<p class="caption">p-values of the differential expression analysis using plate as null model.</p>
</div><p class="caption">p-values of the differential expression analysis using plate as null model.</p>

<!--- DAVID END --->


# Conclusions


<!--- DAVID START --->
- In our data summary we can see a lot of missing values for clinical traits that we should check and consider in the following analysis. However, gene information seems to be quite robust.  

- We have subset our data set using a non-paired analysis design to avoid assesing the extend of differential expression for non-independent samples. In order to do that we have used complementary criteria exluding the tumor samples that were duplicated as normal. 
- After that subset to define the analysis approach we have also include another filter layer on top by selecting the individuals with larger library size values.

- We have detected some problematic normal samples in the MDS plot from an unique batch. They have been filtered out to avoid further problems. 

- The main source of variation seems to be driven by the sample type (tumor or normal) as seen in the DE analysis. 

- When adjusting this DE analysis by the detected Surrogate Variants the number of DE genes increases from 5800 up to 7000. This may indicate there is still some sort of heterogeneity we cannot detect or may be just the individual variability of the gene expression.  

- We cannot detect any significant Batch Effect after checking for the **TSS** and the **plate** attributes. For the TSS we have checked using a dendogram plot, a MDS and a SVA considering the TSS variablity as the null model. For the Plate variables we have only performed an SVA that suggests there is no Batch Effect. 

- When using SVA to check Batch Effect of plate and tss attributes, we observe a decrease in the quantity of DE genes. Although this migth imply that indeed exists this batch effect, the change is not significant and it not affect the p-value distribution. Therefore, although the variability of some genes may be explained by plate or TSS, this genes are a few hundreds that won't really affect the conclusions of this analysis. 


<!--- DAVID END--->


<!--- ADRI START --->


- There are present 7145 genes with a Differental Expression with a FDR of 1% (having as null model with a matrix of intercept term ~ 1 and adjusting by the detected SVs).

 
<!--- ADRI END --->
- Wide projects like The Cancer Genome Atlas imply a really difficult process of filtering and assesing the data to obtain relevant and real results. Even with that treatment, we can always have a lot of variability making difficult the extraction of real or relevant conclusions. 

# Session information


```r
sessionInfo()
```

```
R version 3.2.4 (2016-03-10)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.11.4 (El Capitan)

locale:
[1] es_ES.UTF-8/es_ES.UTF-8/es_ES.UTF-8/C/es_ES.UTF-8/es_ES.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] gridExtra_2.2.1            sva_3.18.0                
 [3] genefilter_1.52.1          mgcv_1.8-12               
 [5] nlme_3.1-125               ggplot2_2.1.0             
 [7] geneplotter_1.48.0         annotate_1.48.0           
 [9] XML_3.98-1.4               AnnotationDbi_1.32.3      
[11] lattice_0.20-33            edgeR_3.12.1              
[13] limma_3.26.9               SummarizedExperiment_1.0.2
[15] Biobase_2.30.0             GenomicRanges_1.22.4      
[17] GenomeInfoDb_1.6.3         IRanges_2.4.8             
[19] S4Vectors_0.8.11           BiocGenerics_0.16.1       
[21] knitr_1.12.3              

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.4        formatR_1.3        RColorBrewer_1.1-2
 [4] plyr_1.8.3         XVector_0.10.0     tools_3.2.4       
 [7] zlibbioc_1.16.0    digest_0.6.9       evaluate_0.8.3    
[10] RSQLite_1.0.0      gtable_0.2.0       Matrix_1.2-4      
[13] DBI_0.3.1          yaml_2.1.13        stringr_1.0.0     
[16] grid_3.2.4         survival_2.38-3    rmarkdown_0.9.5   
[19] magrittr_1.5       codetools_0.2-14   splines_3.2.4     
[22] scales_0.4.0       htmltools_0.3.5    xtable_1.8-2      
[25] colorspace_1.2-6   labeling_0.3       KernSmooth_2.23-15
[28] stringi_1.0-1      munsell_0.4.3     
```
