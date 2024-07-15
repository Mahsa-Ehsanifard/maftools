
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


###### **Required files**

* MAF files - can be obtained from TCGA or gz compressed.

+ here we downloaded the MAF file of LUAD cohort from TCGA (GDC portal)

```{r}
luad_maf <- GDCquery(project = "TCGA-LUAD", 
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     access = "open")
```

```{r}
GDCdownload(luad_maf)
```

```{r}
luadmaf <- GDCprepare(luad_maf)
```




