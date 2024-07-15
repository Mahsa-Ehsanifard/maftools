library(TCGAbiolinks)
library(DT)
library(dplyr)
library(maftools)
#visualization
# MAF summary
#varient_clasification
plotmafSummary(maf = luad, rmOutlier = T, addStat = "median", dashboard = T, titvRaw = F)
mafbarplot(maf = luad, genes = "TP53")
#Oncoplot
oncoplot(maf = luad, genes = c("TP53", "TTN"))
#Transition and Transversion
luad.titv <- titv(maf = luad,plot = F, useSyn = T)
plotTiTv(luad.titv)
#Lollipop plot
#changes in amino acids
lollipopPlot(maf = luad, gene = "TP53", AACol = "HGVSp_Short",
             showMutationRate = T, showDomainLabel = T, cBioPortal = T)
#rainfall plot
rainfallPlot(maf = luad, detectChangePoints = T, pointSize = 0.4) #most mutated sample
rainfallPlot(maf = luad, detectChangePoints = T, pointSize = 0.7,
             tsb = "TCGA-05-4244-01A-01D-1105-08")
#Compare mutation against TCGA cohorts
par(mfrow=c(1,1))
par(pin=c(4,2))
luad.compare <- tcgaCompare(maf = luad, cohortName = "my_luad",
                            tcga_cohorts = c("LUAD","LUSC", "SKCM", "BRCA"),
                            logscale = T, capture_size = 50, axisFontSize = 0.5)
#TMB=Tumor mutation Burden
#total numbers of mutations per megabase (mb)