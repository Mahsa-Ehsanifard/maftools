
# **maftools package in R**

![](https://img.shields.io/badge/version-2.20.0-blue)
![](https://img.shields.io/badge/open%20access-100%25-green)
![](https://img.shields.io/badge/installation-Rstudio-orange)
![](https://img.shields.io/badge/source-bioconductor-9cf)
![](https://img.shields.io/badge/dependenties-16-yellow)

## required **Mutation Annotation Format (MAF)** file

### Masked Somatic Mutation analysis


```maftools``` provides various functions to perform most commonly used analyses in cancer genomics. This package attempts to summarize, analyze, annotate and visualize MAF files in an efficient manner from either **The Cancer Genome Atlas (TCGA)** sources or any in-house studies as long as the data is in **MAF** format.

[DOI:10.18129/B9.bioc.maftools](https://www.bioconductor.org/packages/release/bioc/html/maftools.html)

Such cohort-based large-scale characterizations often produce large amounts of data in the form of *somatic variants* containing **single-nucleotide variants (SNV)** and small insertion/deletions (indels).  Somatic variants provide baseline data for many analyses, such as **Single Nucleotide Polimorphism (SNP)** types, driver gene detection, mutation frequencies, somatic interactions, visualization, and estimation of tumor heterogeneity, applied by the package employment. 

Thanks to 
[Genome Resarch. PMID: 30341162](https://doi.org/10.1101/gr.239244.118) and [PoisonAlien/maftools](https://github.com/PoisonAlien/maftools)



### Installation

Install from `BiocManager` package

```{r}
BiocManager::install("maftools")
```

``` {r}
library(maftools)
```


##### **Required files**

* MAF files - can be obtained from TCGA or gz compressed.

+ here we downloaded the MAF file of LUAD cohort from TCGA (GDC portal)

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

Now, we have all genes with mutation types in LUAD cohort along with features data. SNP mutations are detectable in the output table, which can be separated based on *Variant_Classification* column in *luadmaf* table resulted from `GDCprepare` function.

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

We also access to altered alleles and substitutions' position in each gene, indication nucleotides  and amino acid exchanges in columns: *HGVSc*, *HGVSp*, and *HGVSp_Short*.

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


