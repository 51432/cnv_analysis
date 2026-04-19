#!/usr/bin/env bash

# CNVkit WES pipeline shared configuration.
REFERENCE="/data/person/wup/public/liusy_files/reference_genomes/hg38/reference/Homo_sapiens_assembly38.fasta"
CNVKIT_TARGET_BED="/data/person/wup/public/liusy_files/reference_genomes/hg38/intervals/UCSC_Agilent_exome_regions.hg38.bed"
CNVKIT_ACCESS_BED="/data/person/wup/public/liusy_files/reference_genomes/hg38/resources/access.hg38.bed"

# Recommended project-level CNVkit pooled-normal reference path.
CNVKIT_REFERENCE_DEFAULT="/data/person/wup/public/liusy_files/reference_genomes/hg38/resources/cnvkit/wes_reference.cnn"

# Default runtime options.
CNVKIT_MODE_DEFAULT="wes"
CNVKIT_STAGE_DEFAULT="build-reference"
CNVKIT_THREADS_DEFAULT=8
CNVKIT_MAX_PARALLEL_DEFAULT=8

# Slurm defaults (can be overridden via environment variables before submission).
# Partition is restricted to cpu1/cpu2 in this project.
SLURM_PARTITION_DEFAULT="cpu1"
SLURM_PARTITION_ALLOWED="cpu1 cpu2"
SLURM_TIME_DEFAULT="24:00:00"
SLURM_MEM_DEFAULT="16G"
