library(TCGAbiolinks)
library(dplyr)
library(DT)
#☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻
#Masked Somatic Mutation = MAF (file)
#****** Changes in DNA => not inherited 
#* Cancer research => Tumor development and progression 
#* Masking => VCF file
#****** 
#only one tumor-normal pair => *Aliquot Selection* => TCGA aliquot barcode
#1♣ Removing: 
Mutation_Status != "Somatic" -> GDC_FILTER
!duplicated()
BadSeq
#2♣ Remaining:
GDC_Valid_Somatic = TRUE 
FILTER != "panel_normals" 
#****** pipeline for *Raw mutation data = MuTect2
###################

luad_maf <- GDCquery(project = "TCGA-LUAD", 
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     access = "open")
GDCdownload(luad_maf)
luadmaf <- GDCprepare(luad_maf)
tp <- luadmaf[luadmaf$Hugo_Symbol=="TP53",]

