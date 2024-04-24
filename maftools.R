library(TCGAbiolinks)
library(maftools)
install.packages("maftools")
BiocManager::install("maftools")
library(maftools)
library(dplyr)
tcga_maf <- GDCquery(project = "TCGA-SKCM", 
                     data.category = "Simple Nucleotide Variation", # Simple nucleotide variation if legacy
                     data.type = "Masked Somatic Mutation",
                     access = "open")
GDCdownload(tcga_maf)
maf <- GDCprepare(tcga_maf)


maf_mutect2 <- GDCquery_Maf("SKCM", pipelines = "mutect") 
mutect_braf <- maf[maf$=="BRAF", ]
p.V600E <- mutect_braf[mutect_braf$HGVSp_Short=="p.V600E", ]
maf_mutect2 <- data.frame(maf_mutect2)
getGeneSummary(maf_mutect2)
plotmafSummary(maf = maf_mutect2, rmOutlier = T, addStat = "median", dashboard = T)

maf = read.maf(maf = maf_mutect2)
sample_summery <- getSampleSummary(maf)
gene_summary <- getGeneSummary(maf)
clinic_summary <- getClinicalData(maf)
getFields(maf)
write.mafSummary(maf = maf, basename = 'maf')

plotmafSummary(maf, rmOutlier = T, addStat = "median", dashboard = T)
oncoplot(maf = maf, top = 10)

PV <- read.maf(maf = p.V600E)
ss <- getSampleSummary(PV)
gs<- getGeneSummary(PV)
cs <- getClinicalData(PV)
plotmafSummary(PV, rmOutlier = T, dashboard = T, addStat = "median")

lollipopPlot(
  maf = maf,
  gene = "BRAF",
  AACol = "HGVSp_Short",
  showMutationRate = TRUE,labelPos = "all", labPosSize = 0.7,showDomainLabel = F,cBioPortal = T,
  roundedRect = T,repel = T,collapsePosLabel = T ,labPosAngle = 45,
  domainLabelSize = 0.6, colors = c("navy","red3","green","yellow"),domainAlpha = 1,domainBorderCol = NA,
  bgBorderCol = NA, labelOnlyUniqueDoamins = T,pointSize = 2)

rainfallPlot(maf = maf, detectChangePoints = TRUE, pointSize = 0.5, tsb = NULL)

skcm.mutload <- tcgaCompare(maf = maf, cohortName = "my_SKCM", logscale = T, capture_size = 50,
                            tcga_capture_size =  35.8,cohortFontSize = 0.6)
skcm.mutload = tcgaCompare(maf = maf, cohortName = "my_SKCM", logscale = T, capture_size = 50,
                           tcga_cohorts = "SKCM")

somaticInteractions(maf = maf, top = 25, pvalue = c(0.05,0.01))

# gistic
skcm.maf <- system.file("maf_mutect2", package = "maftools")
SKCM <- read.maf(skcm.maf)

all.lesions <- system.file("maf", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("maf", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("maf", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- system.file("maf", "scores.gistic", package = "maftools")
skcm.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)

PlotOncogenicPathways(maf = maf,fullPathway = F, pathways = "Cell_Cycle")

