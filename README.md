# IEO_PRAC_analysis

## Contributors

* Adrià Auladell ( `@adriaaula` )
* Joan Marti ( `@joanmarticarreras` )
* David Mas(`@davidmasp`)

## What is prostate cancer?

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

Complete Name: Prostate Adenocarcinoma

source: [Cancer Genome - NIH](http://cancergenome.nih.gov/cancersselected/prostatecancer)

##Tips & Tricks
Open R with terminal instead of Rstudio (Ubuntu users). It shuts down your other applications, a side from Rstudio itself.

##TO DO

- [x] Choose a RNA-seq data set among the ones offered in the data sets folder.
- [ ]Read introductory material about the corresponding type of tumor you are going to analyze. Try to understand the major tumorigenic mechanisms that participate in the growth and proliferation of such a tumor type. This essentially means to gather information about what are the most relevant genes to this cancertype. The goal of the project is not to reproduce previous findings about this cancer type, although you may choose to do so, but rather to find some simple question related to this cancer type that you can answer focusing on the contrast between tumor and normal samples, or on some other simple contrast of interest using the analysis techinques we have seen in class. In this respect, you can consider analyzing some of the available clinical variables such as the tumor stage encoded in 'ajcc_pathologic_tumor_stage'. Each clinical variable has a metadata with a so-called 'CDE' identifier which you can use to fetch futher information about at https://cdebrowser.nci.nih.gov.
- [ ]Try to figure out factors that generate variability unrelated to the outcome of interest. By means of the diagnostics we have seen the lecture on batch identification, try to ensure that you do not have a major confounding with the outcome of interest.
- [ ]To speed up some parts of the analysis, you may choose to work with a subset of the samples.
- [ ]Carry out quality assessment and normalization of the data.
- [ ]Search for differentially expressed (DE) genes using the simple F-test implemented in the package SVA to do a two-group comparison. For this first part of the project *do not attempt* to interpret the list of DE genes, just report how many of them do you find and how the distribution of p-values looks like. Consider estimating surrogate variables with SVA to see whether the number of expression changes increases or decreases. Bear in mind that the actual DE analysis you will do in the second part will be more sophisticated since you will have to take into account aspects such as variance heterogeneity of log CPM values or the fact that normal samples were derived from the same pool of individuals as a fraction of the tumor samples. For this first part of the project, there is no need to address these issues.
- [ ]The project template that is provided is not comprehensive and it just tries to help you in quickly learning how to work with R Markdown files, so you should try to make your supplementary material more readable and complete than the template provided.

