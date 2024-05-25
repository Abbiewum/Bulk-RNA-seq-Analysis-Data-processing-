#!/bin/bash

# Download datasets
fastq-dump --split-files SRR2842079 SRR28420796 SRR28420797 SRR28420798

# Quality control using fastqc
install fastqc
fastqc SRR28420795_1.fastq SRR28420795_2.fastq

# Trimming using sickle
sudo install sickle
sickle pe -f SRR28420795_1.fastq -r SRR28420795_2.fastq -t sanger -l 50 -q 20 -o output_R1_trimmed_SRR28420795_1.fastq -p output_R2_trimmed_SRR28420795_2.fastq -s sample1_singleton.fastq

# Downloading index file
HISAT2 website index for human

# Downloading Annotation file and Unzip
wget ftp://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
gunzip -d genome Homo_sapiens.GRCh38.111.gtf.gz

# Alignment with HISAT2
install hisat2
hisat2 -p 3 -x ~/ref_genome/grch38_genome/grch38/genome -1 output_R1_trimmed_SRR28420795_1.fastq -2 output_R2_trimmed_SRR28420795_2.fastq -S Mapping/SRR28420795.sam

# Convert to a BAM file
samtools view -b -s Mapping/SRR28420795.sam -o Mapping/SRR28420795.bam

# Sorting the BAM file
samtools sort -@ 3 -o SRR28420795.sorted.bam Mapping/SRR28420795.bam

# Quantification using featureCounts
featureCounts -T 3 -p -a ~/Homo_sapiens.GRCh38.111.gtf -o counts_SRR28420795.txt /home/abby/SRR28420795.sorted.bam


