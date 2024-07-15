library(TCGAbiolinks)
library(DT)
library(dplyr)
library(tidyverse)
library(SummarizedExperiment)

luad_maf <- GDCquery(project = "TCGA-LUAD", 
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     access = "open")
GDCdownload(luad_maf, method = "client")
luadmaf <- GDCprepare(luad_maf, summarizedExperiment = T)
##***************************************************************************
BiocManager::install("maftools")
library(maftools)
#read.maf
# reads MAF files, summarizes in various ways and stores as an MAF object.
# recommended => providing annotations with samples in MAF
# can also integrate copy number data; if available

LUAD <- read.maf(maf = luadmaf)
sample <- getSampleSummary(LUAD) # sample frequency summary
gene <- getGeneSummary(LUAD) # Gene frequency summary
clin <- getClinicalData(LUAD)
field <- getFields(LUAD)
write.mafSummary(LUAD, basename = "LUAD")
#*************************************************************************
# Analysis
# Somatic Interactions
# ♣ Co-occurance and ♣ Mutually exclusive => pair genes
# Useful for identification of relationships between genes based on their mutation patterns 
### across a popularion of samples.
##
#♠ Mutually exlusive => One gene is mutated, the other gene is not
#♠ Co-accurring (co-occurance frequently) => both genes are often mutated together
# Top ranked paired genes
# highlights potential functional connections between the genes 
##
#Significant => Fisher's Exact Test. 
##costumization, P.Value < 0.05 or 0.01 

somaticinter <- somaticInteractions(maf = LUAD, top = 20, pvalue = 0.01)

#0,1
#0 => mutated , 1 => non-mutated
#*****************************************************************************
#Detecting cancer driver genes
#oncodrive => driver genes => specific loci 

luadSig <- oncodrive(maf = LUAD, minMut = 5, pvalMethod = "zscore")
plotOncodrive(res = luadsig, fdrCutOff = 0.05, useFraction = T, labelSize = 0.5)
sig <- luadsig[1:20,]
plotOncodrive(res = sig, fdrCutOff = 0.05, bubbleSize = 25,useFraction = T, labelSize = 0.5)

