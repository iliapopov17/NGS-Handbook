# NGS Data Analysis Manuals

> The materials in this repository are based on educational courses I have taken<br>
> It can be used as a helpful repository with cheat-sheets for NGS studies.

<div style='justify-content: center'>
<img src="https://github.com/iliapopov17/NGS-Data-Analysis-Manual/blob/main/imgs/NGS_workflow.png" align='center', width="100%">
</div>

_Typical workflow of NGS data analysis_

To recreate any of the steps of this manual please install:

```bash
conda env create -f ngs-manual.yml
```

And of cource do not forget to activate the envinronment!

```bash
conda activate ngs-manual
```

## Much more to be disclosured soon:
- Whole Genome and Pangenome Analyses
- Merging Phylogenetics pipeline into this handbook repository
- 16S Amplicon Analysis

## Genomic Variation Analysis

In the [Genomic Variation Analysis folder](02_Genomic_Variation_Analysis) there is a detailed guide how to conduct studies on Variant Calling using `fastqc`, `trimmomatic`, `bwa`, `samtools`, `abra2`, `bcftools`, `snpEff` & `SnpSift`

## Quality Control of raw data

In the [Quality Control folder](01_Quality_Control) there is a detailed guide how to conduct quality control of raw data using `fastqc` and `trimmomatic`.
