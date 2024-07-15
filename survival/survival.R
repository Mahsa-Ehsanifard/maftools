library(TCGAbiolinks)
library(DT)
library(dplyr)
library(maftools)
#*********************
#Adding pfam domain changes
## which domain of gene or protien is affected
## pfam annotation
pfam <- pfamDomains(maf = LUAD, top = 10)
prSum <- pfam$proteinSummary
Dsum <- pfam$domainSummary
#****************************
#*
#Survival analysis => Kaplan_miere (KM)
# required clinical input => Tumor_sample_Barcode -> matches to those in MAF file
### binary event -> 1,0 => Status (vital_status): Alive = 0, Dead = 1
### time to event -> last followup => days 
# Clinical annotation 

cohort <- tcgaLoad(study = "LUAD")
clin <- getClinicalData(cohort)
clin$days_to_last_followup[clin$days_to_last_followup=="[Not Available]"] <- NA
clin$vital_status[clin$vital_status=="Alive"] <- 0
clin$vital_status[clin$vital_status=="Dead"] <- 1

data <- cohort@data

luad <- read.maf(data, clinicalData = clin)

mafSurvival(maf = luad, genes = "DNMT3A", 
            time = "days_to_last_followup", Status = "vital_status")

## predict genesets associated with survival
# identify set of genes which affects survival ratio
pregeneset <- survGroup(maf = luad, top = 5, geneSetSize = 2, 
                        time = "days_to_last_followup",plot = T,
                        Status = "vital_status",verbose = T)

mafSurvGroup(maf = luad, geneSet = c("KRAS","BRAF"), 
             time = "days_to_last_followup", Status = "vital_status")

