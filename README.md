# Single-Cell RNA-seq Analysis: A Reproducible Workflow

## Overview

This repository presents a reproducible computational workflow for analyzing
publicly available single-cell RNA sequencing (scRNA-seq) datasets using
R and Python.

The analysis focuses on characterizing cellular heterogeneity, identifying
cell populations and their molecular signatures, investigating developmental
trajectories, and exploring cell–cell communication.

The workflow covers essential scRNA-seq processing and analysis steps,
including quality control, normalization, dimensionality reduction,
clustering, marker gene identification, and cell type annotation.
It also incorporates advanced computational approaches for data integration,
reference-based label transfer, RNA velocity, cell fate inference,
and ligand–receptor interaction analysis.

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

## Full Analysis

The complete documented analysis, including code,
figures, results, and biological interpretations,
is available here:

[View the full scRNA-seq analysis](analysis.md)

