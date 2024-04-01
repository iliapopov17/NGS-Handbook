#!/bin/bash

set -e

# Check the reads quality
fastqc SRR17909485_*.fastq.gz

# Trim low-quality bases
trimmomatic PE -threads 2 SRR17909485_1.fastq.gz SRR17909485_2.fastq.gz SRR17909485_R1.trim.paired.fastq.gz SRR17909485_R1.trim.unpaired.fastq.gz SRR17909485_R2.trim.paired.fastq.gz SRR17909485_R2.trim.unpaired.fastq.gz LEADING:22 TRAILING:22 SLIDINGWINDOW:6:22 MINLEN:32

# Index reference
bwa index EcoliK12MG1655.fa

# Align to reference
bwa mem -t 2 -R '@RG\tID:1\tPL:illumina\tPU:barcode1\tSM:SRR17909485\tLB:Library1\tCN:BioschoolVarCall' EcoliK12MG1655.fa SRR17909485_R1.trim.paired.fastq.gz SRR17909485_R2.trim.paired.fastq.gz | samtools view -b > EcoliK12MG1655.SRR17909485.unsorted.bam

# Sort alignment
samtools sort --threads 2 EcoliK12MG1655.SRR17909485.unsorted.bam > EcoliK12MG1655.SRR17909485.sorted.bam

# Make bam index
samtools index EcoliK12MG1655.SRR17909485.sorted.bam

# Realign indels
abra2 --threads 2 --mad 100 --mbq 24 --ref EcoliK12MG1655.fa --in EcoliK12MG1655.SRR17909485.sorted.bam --out EcoliK12MG1655.SRR17909485.final.bam

# Index final bam
samtools index EcoliK12MG1655.SRR17909485.final.bam

# Call variants
bcftools mpileup -Ou --max-depth 5000 -f EcoliK12MG1655.fa EcoliK12MG1655.SRR17909485.final.bam | bcftools call -mv --ploidy 1 -Ov -o EcoliK12MG1655.SRR17909485.called.bcftools.vcf

# Make snpEff database
echo "EcoliK12MG1655.genome : EcoliK12MG1655
EcoliK12MG1655.chromosomes : EcoliK12MG1655.gb
EcoliK12MG1655.codonTable : Standard" > snpEff.config
mkdir -p data/EcoliK12MG1655
cp EcoliK12MG1655.gb data/EcoliK12MG1655/genes.gbk
snpEff build -c snpEff.config -genbank EcoliK12MG1655

# Annotate variants
snpEff ann -v EcoliK12MG1655  EcoliK12MG1655.SRR17909485.called.bcftools.vcf > EcoliK12MG1655.SRR17909485.annotated.vcf
mv snpEff_genes.txt EcoliK12MG1655.SRR17909485.snpEff_genes.txt
mv snpEff_summary.html EcoliK12MG1655.SRR17909485.snpEff_summary.html

# Filter high-quality variants
SnpSift filter "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')" EcoliK12MG1655.SRR17909485.annotated.vcf > EcoliK12MG1655.SRR17909485.higheffect.vcf
