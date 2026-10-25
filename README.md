# Single-Cell RNA-seq Analysis: A Reproducible Workflow

## Overview

This repository documents a reproducible analysis of publicly available
single-cell RNA sequencing (scRNA-seq) datasets using R and Python.

The project aims to develop practical skills in scRNA-seq data processing,
cellular heterogeneity analysis, cell type annotation, and the investigation
of cellular dynamics and cell–cell interactions.

The workflow is primarily based on the
[scRNA-seq Analysis Vignette](https://github.com/quadbio/scRNAseq_analysis_vignette),
with additional implementations and modifications using Seurat, Scanpy,
scVelo, CellRank, CellPhoneDB, and COMUNET.

## Analysis Workflow

The analysis is organized into four main parts:

1. **Basic scRNA-seq Analysis**
   - Quality control and filtering
   - Data normalization and highly variable gene selection
   - Principal component analysis (PCA)
   - Dimensionality reduction using UMAP and t-SNE
   - Cell clustering and marker gene identification
   - Cell type annotation
   - Pseudotime analysis

2. **Data Integration and Batch Effect Correction**
   - Integration of multiple scRNA-seq datasets
   - Seurat anchor-based integration
   - Harmony integration
   - Comparison of integration approaches

3. **Reference-Based Cell Type Annotation**
   - Cluster-level transcriptomic similarity analysis
   - Cell-level transcriptomic similarity analysis
   - Reference-based label transfer using Seurat

4. **Advanced scRNA-seq Analysis**
   - Cluster connectivity analysis using PAGA
   - Diffusion pseudotime analysis
   - RNA velocity analysis using scVelo
   - Cell fate and trajectory inference using CellRank
   - Cell–cell communication analysis using CellPhoneDB
   - Ligand–receptor interaction network analysis using COMUNET

## Data Sources

The analyses primarily use publicly available datasets from the following studies:

1. [Organoid single-cell genomic atlas uncovers human-specific features
   of brain development](https://www.nature.com/articles/s41586-019-1654-9)

2. [Charting human development using a multi-endodermal organ atlas
   and organoid models](https://www.sciencedirect.com/science/article/pii/S0092867421005316)

These datasets are reanalyzed for educational purposes and to explore
computational approaches for investigating cellular heterogeneity,
developmental trajectories, and intercellular communication.

## Scope

This project represents a computational reanalysis of existing public
datasets. It does not involve the generation of new experimental data.

The repository provides analysis scripts, methodological explanations,
and visualization workflows to support reproducibility and further
exploration of single-cell transcriptomic data.
