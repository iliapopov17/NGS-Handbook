# Crohn's disease
>This folder contains a training data for 16S amplicon analysis and an example of analysis.
- `Crohns_disease.Rmd` - `R Markdown` journal with this whole pipeline done

## Introduction

In this study, the gut microbiome composition in individuals with Crohn’s Disease (CD) and Healthy Controls (HC) was investigated. The goal was to identify taxonomic shifts, potential functional implications, and associations with disease. High-throughput sequencing data and various bioinformatics tools were employed to analyze microbial diversity and abundance.

Data obtained from the article [«_A microbial signature for Crohn's disease_»](https://pubmed.ncbi.nlm.nih.gov/28179361/) was used. These data include information on microbial samples from healthy individuals and patients with Crohn's disease.

## Instruction
- Language: `R` v.4.2+
- IDE: `RStudio`

You can run commands below in your `RStudio` in `R script`.<br>
Or if you want to write a beautiful & convenient to read laboratory journal you can use `R Markdown`.

📝 Also see pre-prepared [laboratory journal](https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_01_Introduction/16S_analysis_pipeline.Rmd)

### Loading libraries and data

In this 16S amplicon analysis pipeline we will use `MicrobeR`, `balance`, `NearestBalance` & `selbal` packages. Their installation is a bit difficult. Please follow the code below:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("philr")
BiocManager::install("DECIPHER")

devtools::install_github("jbisanz/MicrobeR")
devtools::install_github("tpq/balance")
devtools::install_bitbucket("knomics/nearestbalance")
devtools::install_github(repo = "malucalle/selbal")
```

First, load all the libraries needed for the data analysis

**_Input_**

```r
library(data.table)
library(openxlsx)
library(MicrobeR)
library(ggplot2)
library(zCompositions)
library(NearestBalance)
library(GUniFrac)
library(vegan)
library(ape)
library(selbal)
```
 
Then set the working directory
 
**_Input_**

```r
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)
```

Load metadata and sort it by participant number

**_Input_**

```r
metadata = fread("data_Crohns_disease/metadata.csv")
metadata[, .N, by = diagnosis_full]
```

**_Output_**

```
diagnosis_full  N
<chr> <int>
CD  34
HC  34
2 rows
```

Load the `Counts` table

**_Input_**

```r
counts <- read.csv("data_Crohns_disease/counts.csv",row.names = 1)
counts <- counts[metadata$sample, colSums(counts) >0]
```

### Check the data

How many samples and microbial samples?

**_Input_**

```r
dim(counts)
```

**_Output_**

```
[1]  68 210
```

What's the coverage?

**_Input_**

```r
range(rowSums(counts))
```

**_Output_**

```
[1]  19414 121234
```

Composition of samples

**_Input_**

```r
raw_abund <- counts/rowSums(counts)
```

**_Input_**

```r
metadata$SampleID <- metadata$sample
```

**_Input_**

```r
Microbiome.Barplot(t(raw_abund), metadata, CATEGORY = "diagnosis_full")
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/microbiome_barplot.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/microbiome_barplot.jpg", width = 11, height = 3.5, dpi = 300)
```

Is there enough coverage?

**_Input_**

```r
metadata$coverage <- rowSums(counts)
```

**_Input_**

```r
ggplot(metadata) + 
  geom_histogram(aes(coverage)) + 
  theme_minimal() + 
  xlab("N samples")
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/coverage_quality.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/coverage_quality.jpg", width = 11, height = 3.5, dpi = 300)
```

Sequencing quality is sufficient

### Filtration from rare and under-represented taxa

Keeping the microbes that occur in >30% of samples

**_Input_**

```r
filt_counts <- counts[, colSums(counts>0)>0.3*nrow(counts)]
metadata[, filt_coverage := rowSums(filt_counts)]
metadata[, proportion_of_prevalent_taxa := 100*filt_coverage/coverage]
```

Coverage of samples after filtration

**_Input_**

```r
ggplot(metadata)+
  geom_histogram(aes(filt_coverage)) + 
  theme_minimal()+
  xlab("Post-filtration coverage") + 
  ylab("N samples")
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/post-filtration_coverage.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/post-filtration_coverage.jpg", width = 11, height = 3.5, dpi = 300)
```

What proportion of the microbes remained in the analysis

**_Input_**

```r
ggplot(metadata)+
  geom_histogram(aes(proportion_of_prevalent_taxa)) + 
  theme_minimal()+ 
  xlab("Proportion of microbes remaining in the assay") + 
  ylab("N samples")
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/remaining_proportion.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/remaining_proportion.jpg", width = 11, height = 3.5, dpi = 300)
```

How many samples and microbial samples?

**_Input_**

```r
dim(filt_counts)
```

**_Output_**

```
[1] 68 89
```

What's the coverage?

**_Input_**

```r
range(rowSums(filt_counts))
```

**_Output_**

```
[1]  12981 120882
```

**_Input_**

```r
matrix_data <- matrix(c("Before", 210, 19412,
                         "After", 89, 12981), 
                      nrow = 2, byrow = TRUE)

data <- as.data.frame(matrix_data)

colnames(data) <- c("", "N microbes", "Minimum coverage")

print(data)
```

**_Output_**

```
  N microbes  Minimum coverage
<chr> <chr> <chr>
Before	210	19412
After	89	12981
```

Calculating relative abundances

**_Input_**

```r
abundance <- cmultRepl(filt_counts)
```

**_Output_**

```
No. adjusted imputations:  694 
```

**_Input_**

```r
heatmap_with_split(abundance,
                   metadata, ~ diagnosis_full,
                   show_samp_names = F) + 
  theme(axis.text.y = element_text(size =5))
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/relative_abundance.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/relative_abundance.jpg", width = 8, height = 10, dpi = 300)
```

### Counting alpha diversity

**_Input_**

```r
alpha_div <- rowMeans(sapply(1:5, function(i){
  counts_rar_i = Rarefy(counts, 19000)$otu.tab.rff
  alpha_div_i = vegan::diversity(counts_rar_i)
}))
metadata$Shannon.index <- alpha_div[metadata$sample]
```

**_Input_**

```r
ggplot(metadata) + 
  geom_boxplot(aes(diagnosis_full, Shannon.index, fill = diagnosis_full)) + 
  theme_minimal() +
  theme(legend.position = 'none') + 
  xlab("")
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/alpha_diversity.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/alpha_diversity.jpg", width = 8, height = 8, dpi = 300)
```

Is it different?

**_Input_**

```r
wilcox.test(Shannon.index ~ diagnosis_full, metadata)$p.value
```

**_Output_**

```
[1] 1.314773e-09
```

### Aitchison's beta diversity: is there a difference in proportions?

**_Input_**

```r
clr <- log(abundance) - rowMeans(log(abundance))
beta_div <- dist(clr)
```

**_Input_**

```r
pcoa_res <- pcoa(beta_div)$vectors
var <- apply(pcoa_res, 2, var)
var_rel <- round(var*100/sum(var), 1)
```

**_Input_**

```r
ggplot(cbind(metadata,pcoa_res)) + 
  geom_point(aes(Axis.1, Axis.2, col=diagnosis_full)) +
  coord_fixed() + 
  theme_minimal() + 
  labs(col="") + 
  xlab(paste0("Axis.1 (",var_rel[1], "%)")) + 
  ylab(paste0("Axis.2 (",var_rel[2], "%)"))
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/beta_diversity.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/beta_diversity.jpg", width = 8, height = 8, dpi = 300)
```

Is the difference statistically significant?

**_Input_**

```r
adonis2(beta_div ~ metadata$diagnosis_full)
```

**_Output_**

```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = beta_div ~ metadata$diagnosis_full)
                        Df SumOfSqs      R2      F Pr(>F)    
metadata$diagnosis_full  1     5400 0.16562 13.101  0.001 ***
Residual                66    27205 0.83438                  
Total                   67    32605 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### What exactly is the difference

**_Input_**

```r
nb <- nb_lm(abundance = abundance,
            metadata = metadata,
            pred = "diagnosis_full")
```

**_Input_**

```r
heatmap_with_split(abundance = abundance[, unlist(nb$nb$b1)],
                   metadata = metadata,
                   formula = ~ diagnosis_full,
                   num_name = "health-related",
                   den_name = "disease-related",
                   show_samp_names = F,
                   balance = nb$nb$b1)
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/heatmap_w_split.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/heatmap_w_split.jpg", width = 10, height = 10, dpi = 300)
```

Illustrating the difference between the average microbiota of healthy and sick people

**_Input_**

```r
psi <- make_psi_from_sbp(nb$coord$sbp)
mean_diff_clr <- drop(nb$lm_res$coefficients[2,] %*% psi)
bal_unit <- balance_to_clr(nb$nb$b1, colnames(abundance))
bal_diff_clr <- drop(bal_unit %*% mean_diff_clr) * bal_unit
tab <- data.table(taxon = names(mean_diff_clr),
                  clr_diff = mean_diff_clr,
                  bal = bal_diff_clr)
setorderv(tab, "clr_diff")
tab$taxon <- factor(tab$taxon, levels = tab$taxon)
```

Mean difference between sick and healthy individuals

**_Input_**

```r
ggplot(tab) + 
  geom_col(aes(clr_diff, taxon), fill = "darkblue") + 
  theme_minimal() + xlab("CLR(v)") + ylab("") + 
  theme(axis.text.y = element_text(size =5))
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/mean_difference.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/mean_difference.jpg", width = 10, height = 10, dpi = 300)
```

Approximate, simplified difference

**_Input_**

```r
ggplot(tab) + 
  geom_col(aes(bal, taxon), fill = "darkblue") + 
  theme_minimal() + xlab("CLR(b)") + ylab("") +
  theme(axis.text.y = element_text(size =5))
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/approximate_difference.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/approximate_difference.jpg", width = 10, height = 10, dpi = 300)
```

Balance value in each sample

**_Input_**

```r
metadata$balance <- balance.fromSBP(abundance, nb$nb$sbp)
```

**_Input_**

```r
ggplot(metadata) + 
  geom_boxplot(aes(diagnosis_full, balance, fill = diagnosis_full)) + 
  theme_minimal() +
  theme(legend.position = 'none') + 
  xlab("") 
```

**_Output_**

<img src="https://github.com/iliapopov17/NGS-Handbook/blob/main/05_16S_amplicon_analysis/05_02_Crohns_disease/imgs/balance_value.jpg" width="100%">

**_Input_**

```r
ggsave("imgs/balance_value.jpg", width = 8, height = 8, dpi = 300)
```

Has it changed?

**_Input_**

```r
wilcox.test(balance ~ diagnosis_full, metadata)$p.value
```

**_Output_**

```
[1] 2.621864e-17
```
