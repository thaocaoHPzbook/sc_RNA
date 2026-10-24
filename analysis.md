# Single-cell RNA-seq Data Analysis in R
<a id="top"></a>
A reproducible workflow for single-cell RNA sequencing (scRNA-seq) data analysis using R, including quality control, clustering, cell annotation, trajectory inference, data integration, and advanced downstream analyses.

## Table of Contents

- [Preparation](#preparation)

### Part 1. Basic scRNA-seq Data Analysis
- [Step 0. Import Seurat Package](#step-0-import-seurat-package)
- [Step 1. Create a Seurat Object](#step-1-create-a-seurat-object)
- [Step 2. Quality Control](#step-2-quality-control)
- [Step 3. Normalization](#step-3-normalization)
- [Step 4. Highly Variable Feature Selection](#step-4-highly-variable-feature-selection)
- [Step 5. Data Scaling](#step-5-data-scaling)
- [Optional: SCTransform](#optional-sctransform)
- [Step 6. Principal Component Analysis (PCA)](#step-6-principal-component-analysis-pca)
- [Step 7. Non-linear Dimensionality Reduction](#step-7-non-linear-dimensionality-reduction)
- [Step 8. Cell Clustering](#step-8-cell-clustering)
- [Step 9. Cell Type Annotation](#step-9-cell-type-annotation)
- [Step 10. Pseudotime Analysis](#step-10-pseudotime-analysis)
- [Step 11. Save Results](#step-11-save-results)

### Part 2. Integration of Multiple scRNA-seq Datasets
- [Step 0. Load Datasets](#step-0-load-datasets)
- [Step 1. Merge Datasets](#step-1-merge-datasets)
- [Step 2. Data Integration Methods](#step-2-data-integration-methods)
  - [2.1 Seurat Integration](#21-seurat-integration)
  - [2.2 Harmony](#22-harmony)
  - [2.3 LIGER](#23-liger)
  - [2.4 Mutual Nearest Neighbors (MNN)](#24-mutual-nearest-neighbors-mnn)
  - [2.5 RSS Integration with BrainSpan](#25-rss-integration-with-brainspan)
  - [2.6 Cluster Similarity Spectrum (CSS)](#26-cluster-similarity-spectrum-css)
- [Step 3. Comparison of Integration Methods](#step-3-comparison-of-integration-methods)

### Part 3. Reference-based Cell Type Annotation
- [Step 0. Load Reference and Query Datasets](#step-0-load-reference-and-query-datasets)
- [Method 1. Transcriptome Similarity](#method-1-transcriptome-similarity)
  - [1.1 Cluster-level Similarity](#11-cluster-level-similarity)
  - [1.2 Cell-level Similarity](#12-cell-level-similarity)
- [Method 2. Seurat Label Transfer](#method-2-seurat-label-transfer)

### Part 4. Advanced scRNA-seq Analysis
- [4.1 Cluster Connectivity Analysis with PAGA](#41-cluster-connectivity-analysis-with-paga)
- [4.2 Branched Trajectory and Pseudotime Analysis](#42-branched-trajectory-and-pseudotime-analysis)
- [4.3 RNA Velocity Analysis](#43-rna-velocity-analysis)
- [4.4 Trajectory Analysis with CellRank](#44-trajectory-analysis-with-cellrank)
- [4.5 Cell–Cell Communication Analysis](#45-cellcell-communication-analysis)

---


## Preparation
<a id="top"></a>
The analyses begin with preprocessed gene expression count matrices generated using Cell Ranger. The original tutorial uses publicly available scRNA-seq datasets from human cerebral organoids, published by [Velasco et al. (2019)](https://www.nature.com/articles/s41586-019-1654-9). These datasets provide a foundation for reproducing established workflows, exploring downstream analyses, and interpreting cellular heterogeneity and developmental trajectories.

## Part 1. Basic scRNA-seq Data Analysis
<a id="top"></a>
The gene expression count matrix from dataset DS1 was imported using `Read10X()`. A Seurat object was then created using `CreateSeuratObject()` to store the expression data and associated cell metadata for downstream analysis.

```r
library(Seurat)

### Step 0. Import Seurat Package
<a id="top"></a>

### Step 1. Create a Seurat Object
<a id="top"></a>

### Step 2. Quality Control
<a id="top"></a>

### Step 3. Normalization
<a id="top"></a>

### Step 4. Highly Variable Feature Selection
<a id="top"></a>

### Step 5. Data Scaling
<a id="top"></a>

### Optional: SCTransform
<a id="top"></a>

### Step 6. Principal Component Analysis (PCA)
<a id="top"></a>

### Step 7. Non-linear Dimensionality Reduction
<a id="top"></a>

### Step 8. Cell Clustering
<a id="top"></a>

### Step 9. Cell Type Annotation
<a id="top"></a>

### Step 10. Pseudotime Analysis
<a id="top"></a>

### Step 11. Save Results
<a id="top"></a>

## Part 2. Integration of Multiple scRNA-seq Datasets
<a id="top"></a>

### Step 0. Load Datasets
<a id="top"></a>

### Step 1. Merge Datasets
<a id="top"></a>

### Step 2. Data Integration Methods
<a id="top"></a>

#### 2.1 Seurat Integration
<a id="top"></a>

#### 2.2 Harmony
<a id="top"></a>

#### 2.3 LIGER
<a id="top"></a>

#### 2.4 Mutual Nearest Neighbors (MNN)
<a id="top"></a>

#### 2.5 RSS Integration with BrainSpan
<a id="top"></a>

#### 2.6 Cluster Similarity Spectrum (CSS)
<a id="top"></a>

### Step 3. Comparison of Integration Methods
<a id="top"></a>

## Part 3. Reference-based Cell Type Annotation
<a id="top"></a>

### Step 0. Load Reference and Query Datasets
<a id="top"></a>

### Method 1. Transcriptome Similarity
<a id="top"></a>

#### 1.1 Cluster-level Similarity
<a id="top"></a>

#### 1.2 Cell-level Similarity
<a id="top"></a>

### Method 2. Seurat Label Transfer
<a id="top"></a>

## Part 4. Advanced scRNA-seq Analysis
<a id="top"></a>

### 4.1 Cluster Connectivity Analysis with PAGA
<a id="top"></a>

### 4.2 Branched Trajectory and Pseudotime Analysis
<a id="top"></a>

### 4.3 RNA Velocity Analysis
<a id="top"></a>

### 4.4 Trajectory Analysis with CellRank
<a id="top"></a>

### 4.5 Cell–Cell Communication Analysis
<a id="top"></a>


