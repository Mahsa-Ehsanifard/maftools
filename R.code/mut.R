library(TCGAbiolinks)
library(dplyr)
library(DT)
library(tidyverse)
luad_maf <- GDCquery(project = "TCGA-LUAD", 
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     access = "open",
                     workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(luad_maf, method = "client")
luadmaf <- GDCprepare(luad_maf,summarizedExperiment = T)
load("C:/Users/Lenovo/Desktop/.RData")

library(maftools)
#
luad= read.maf(maf = luad_maf) # reads MAF files, summarizes it in various ways and stores it as an MAF object. 
                              # recommended to provide annotations associated with samples in MAF.
                              # One can also integrate copy number data if available.

sample<-getSampleSummary(luad) #sample summary.
gene <- getGeneSummary(luad) #genes summary
clin <- getClinicalData(luad) # clinical data associated with samples
field <- getFields(luad)
write.mafSummary(luad, basename = "luad")# save all summary functions to the current directory

#*********************************************************************************
#Analysis
Somaticinter <- somaticInteractions(maf = luad, top = 20, pvalue = c(0.05,0.01))
print(Somaticinter)

#Detecting cancer driver genes
luadsig = oncodrive(maf = luad, minMut = 5, pvalMethod = "zscore")
sig <- luadsig[1:20,]
plotOncodrive(res = sig, fdrCutOff = 0.1, useFraction = TRUE,
              labelSize = 0.5)

##########
#Adding and summarizing pfam domains
pfam = pfamDomains(maf = luad, top = 10)
prSum <- pfam$proteinSummary
dSum <- pfam$domainSummary

#Survival analysis
cohort <- tcgaLoad(study = "LUAD")
data <- cohort@data
clin <- getClinicalData(cohort)
clin$days_to_last_followup[clin$days_to_last_followup=="[Not Available]"] <- NA
clin$vital_status[clin$vital_status=="Alive"] <- 0
clin$vital_status[clin$vital_status=="Dead"] <- 1


luad= read.maf(maf = data, clinicalData = clin) 
mafSurvival(maf = luad, genes = "TP53", time = "days_to_last_followup", 
                   Status = "vital_status")


plotmafSummary(maf = LUAD, rmOutlier = T, addStat = "mean", dashboard = T, titvRaw = F)
mafbarplot(maf = LUAD,)

oncoplot(maf = LUAD, top = 10)

laml.titv = titv(maf = LUAD, plot = F, useSyn = TRUE)
plotTiTv(res = laml.titv)

lollipopPlot(
  maf = LUAD,
  gene = 'TP53',AACol = "HGVSp_Short",labelPos = 272,
  showMutationRate = TRUE)

plotProtein(gene = "TP53", refSeqID = "NM_001126117")

rainfallPlot(maf = LUAD, detectChangePoints = TRUE, pointSize = 0.4,#most mutated sample,
             tsb = "TCGA-05-4382-01A-01D-1931-08") 


par(mfrow=c(1,1))
par(pin=c(2,2))
luad.mutload = tcgaCompare(maf = LUAD, cohortName = 'Example-LUAD', logscale = TRUE, 
                           capture_size = 50)

plotVaf(maf = LUAD)

all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gistic <- system.file("extdata", "scores.gistic", package = "maftools")
luad.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, gisticScoresFile = scores.gistic, 
                         isTCGA = TRUE)

#getSampleSummary, getGeneSummary

gisticChromPlot(gistic = luad.gistic, markBands = "all")
gisticBubblePlot(gistic = luad.gistic)

gisticOncoPlot(gistic = luad.gistic, clinicalData = getClinicalData(x = luad),
               annotationFontSize = 0.9,fontSize = 0.8,
               clinicalFeatures = "CDR_ajcc_pathologic_tumor_stage", 
               top = 10, sortByAnnotation = T)

gisticOncoPlot(gistic = luad.gistic, clinicalData = clin,
               clinicalFeatures = 'CDR_ajcc_pathologic_tumor_stage', 
               showTumorSampleBarcodes = T,
               sortByAnnotation = T, top = 10)
