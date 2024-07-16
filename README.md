---
editor_options: 
  markdown: 
    wrap: 72
output: html_document
---

# **maftools package in R**

![](https://img.shields.io/badge/version-2.20.0-blue)
![](https://img.shields.io/badge/open%20access-100%25-green)
![](https://img.shields.io/badge/installation-Rstudio-orange)
![](https://img.shields.io/badge/source-bioconductor-9cf)
![](https://img.shields.io/badge/dependenties-16-yellow)

## required **Mutation Annotation Format (MAF)** file

### Masked Somatic Mutation analysis

`maftools` provides various functions to perform most commonly used
analyses in cancer genomics. This package attempts to summarize,
analyze, annotate and visualize MAF files in an efficient manner from
either **The Cancer Genome Atlas (TCGA)** sources or any in-house
studies as long as the data is in **MAF** format.

[DOI:10.18129/B9.bioc.maftools](https://www.bioconductor.org/packages/release/bioc/html/maftools.html)

Such cohort-based large-scale characterizations often produce large
amounts of data in the form of *somatic variants* containing
**single-nucleotide variants (SNV)** and small insertion/deletions
(indels). Somatic variants provide baseline data for many analyses, such
as **Single Nucleotide Polimorphism (SNP)** types, driver gene
detection, mutation frequencies, somatic interactions, visualization,
and estimation of tumor heterogeneity, applied by the package
employment.

Thanks to [Genome Resarch. PMID:
30341162](https://doi.org/10.1101/gr.239244.118) and
[PoisonAlien/maftools](https://github.com/PoisonAlien/maftools)

### Installation

Install from `BiocManager` package

```{r}
BiocManager::install("maftools")
```

```{r}
library(maftools)
```

##### **Required files**

-   MAF files - can be obtained from TCGA or gz compressed.

-   here we downloaded the MAF file of LUAD cohort from TCGA (GDC
    portal)

```{r}
library(TCGAbiolinks)
```

```{r}
library(DT)
```

```{r}
library(dplyr)
```

```{r}
luad_maf <- GDCquery(project = "TCGA-LUAD", 
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     access = "open")
```

```{r}
GDCdownload(luad_maf,method = "client")
```

```{r}
luadmaf <- GDCprepare(luad_maf)
```

Now, we have all genes with mutation types in LUAD cohort along with
features data. SNP mutations are detectable in the output table, which
can be separated based on *Variant_Classification* column in *luadmaf*
table resulted from `GDCprepare` function.

```{r}
head (luadmaf)
```

```         
A tibble: 6 × 140
  Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position End_Position Strand Variant_Classification
  <chr>                <int> <chr>  <chr>      <chr>               <int>        <int> <chr>  <chr>                 
1 ELAPOR1              57535 BI     GRCh38     chr1            109197542    109197542 +      Silent                
2 HIPK1               204851 BI     GRCh38     chr1            113962403    113962403 +      Missense_Mutation     
3 FDPS                  2224 BI     GRCh38     chr1            155319833    155319833 +      Missense_Mutation     
4 YY1AP1               55249 BI     GRCh38     chr1            155668660    155668660 +      Silent                
5 FCRL3               115352 BI     GRCh38     chr1            157697402    157697402 +      Silent                
6 OR6F1               343169 BI     GRCh38     chr1            247712205    247712205 +      Missense_Mutation     
```

We also access to altered alleles and substitutions' position in each
gene, indication nucleotides and amino acid exchanges in columns:
*HGVSc*, *HGVSp*, and *HGVSp_Short*.

```         
 HGVSc     HGVSp       HGVSp_Short
  <chr>     <chr>       <chr>      
1 c.2190C>T p.Phe730=   p.F730=    
2 c.2068C>G p.Gln690Glu p.Q690E    
3 c.964G>C  p.Gly322Arg p.G322R    
4 c.984C>T  p.Leu328=   p.L328=    
5 c.582G>T  p.Leu194=   p.L194=    
6 c.551G>T  p.Trp184Leu p.W184L  
```

##### **Reading MAF file**

-   `read.maf` function reads and summarizes MAF files in various ways,
    and stores as a MAF object.

-   It is recommended to provide annotations associated with samples in
    MAF.

-   Here I read MAF file of LUAD cohort based on an object resulted from
    `GDCprepare` function because there are sample barcodes annotation
    in it.

```{r}
LUAD <- read.maf(maf = luadmaf)
```

```         
#Validating
-Silent variants: 47915 
-Summarizing
--Possible FLAGS among top ten genes:
  TTN
  MUC16
  USH2A
  FLG
-Processing clinical data
-Processing clinical data
-Finished in 30.9s elapsed (25.6s cpu)
```

```         
#Typing laml shows basic summary of MAF file.

LUAD
```

```         
An object of class  MAF 
                        ID summary    Mean Median
                    <char>  <char>   <num>  <num>
 1:             NCBI_Build  GRCh38      NA     NA
 2:                 Center      BI      NA     NA
 3:                Samples     616      NA     NA
 4:                 nGenes   16642      NA     NA
 5:        Frame_Shift_Del    4256   6.909    4.0
 6:        Frame_Shift_Ins    1290   2.094    1.0
 7:           In_Frame_Del     406   0.659    0.0
 8:           In_Frame_Ins      49   0.080    0.0
 9:      Missense_Mutation  126254 204.958  137.5
10:      Nonsense_Mutation   10454  16.971   10.0
11:       Nonstop_Mutation     175   0.284    0.0
12:            Splice_Site    3719   6.037    3.0
13: Translation_Start_Site     211   0.343    0.0
14:                  total  146814 238.334  158.5
```

> sample summary

```{r}
sample <- getSampleSummary(LUAD)
```

> Gene summary

```{r}
gene <- getGeneSummary(LUAD) 
```

> Recommended clinical data associated with each
> sample/Tumor_Sample_Barcode in MAF.

-   Here I get the linicalData from MAF file resulted from `read.maf`
    function. Another way is to obtain clinical data frame from GDC
    portal using `GDCquery` function.

```{r}
clin <- getClinicalData(LUAD)
```

### **Analysis**

#### Somatic Interactions

Interaction between pair genes using `somaticInteractions` function
based on:

-   **Co-occurance**- co-occurance frequency - both genes are often
    mutated together.

-   **Mutually exlusive** - One gene is mutated, the other gene is not.

Useful for identification of relationships between genes based on their
mutation patterns across a popularion of samples. Highlights potential
functional connections between the genes.

-   Significance customization -- *P.Value \< 0.05 or 0.01*

```{r}
somaticinter <- somaticInteractions(maf = LUAD, top = 20, pvalue = 0.0
```

![](C:/Users/Lenovo/Desktop/maftools package/Rplot.TIFF)

