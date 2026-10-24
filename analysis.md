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


### Step 0. Import Seurat Package
<a id="top"></a>
```r
library(Seurat)
```

### Step 1. Create a Seurat Object
<a id="top"></a>

```r
# Load the 10x Genomics count matrix
counts <- Read10X(data.dir = "data/DS1/")

# Create a Seurat object
seurat <- CreateSeuratObject(
  counts = counts,
  project = "DS1"
)

# Inspect the Seurat object
seurat
```

<img width="1094" height="186" alt="image" src="https://github.com/user-attachments/assets/006419a2-1dcc-4d2f-9591-f4d2301b4d51" />

**Results**
The DS1 dataset was successfully imported into Seurat, containing **33,694 genes across 4,672 cells** within a single RNA assay. The object retains the raw gene expression count matrix, which will be used for subsequent quality control, normalization, and downstream analyses.


### Inspect Cell Metadata

The metadata of the first cells can be examined using:

```r
head(seurat)
```

<img width="996" height="666" alt="image" src="https://github.com/user-attachments/assets/271c7465-3d7a-4d76-842d-e36da86d4535" />

**Results**

The output displays cell-level metadata, including:

- **Cell barcode:** A unique identifier for each cell.
- **`orig.ident`:** The dataset or sample identifier (DS1).
- **`nCount_RNA`:** The total number of UMI counts detected in each cell.
- **`nFeature_RNA`:** The number of unique genes detected in each cell (genes with non-zero counts).

For example, the first cell contains **16,752 UMI counts across 3,287 detected genes**.

These metrics provide an initial overview of sequencing depth and gene detection per cell and are used in subsequent quality control to identify potentially low-quality cells.

### Step 2. Quality Control
<a id="top"></a>

Quality control (QC) was performed to identify potentially low-quality cells and doublets before downstream analysis. Three key metrics were examined:

- **`nFeature_RNA`:** Number of genes detected per cell. Very low values may indicate poor-quality cells, while unusually high values may suggest doublets.
- **`nCount_RNA`:** Total UMI counts per cell, reflecting the amount of captured RNA.
- **`percent.mt`:** Percentage of mitochondrial gene counts. Elevated values may indicate cellular stress or damage.

QC thresholds were determined based on the distribution of these metrics rather than applying fixed cutoffs.

#### 2.1. Identify Gene Categories for QC

```r
# Extract gene names
gene_list <- rownames(seurat)

# Identify ERCC spike-in genes
is_spike <- grep("^ERCC", gene_list, value = TRUE)

# Identify mitochondrial genes
is_mito <- grep("^MT-", gene_list, value = TRUE)

# Identify ribosomal protein genes
is_ribo <- grep("^RP[LS]", gene_list, value = TRUE)

# Preview identified genes
head(is_spike)
head(is_mito)
head(is_ribo)
```
<img width="1026" height="218" alt="image" src="https://github.com/user-attachments/assets/260c56e0-2da6-48b0-8c6b-51818dfe1a96" />

#### 2.2. Calculate Mitochondrial, Spike, and Ribosome RNA Percentages

The proportion of mitochondrial gene counts was calculated for each cell using `PercentageFeatureSet()`.

```r
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.ercc"]] <- PercentageFeatureSet(seurat, pattern = "^ERCC")
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RP[LS]")
```

<img width="1418" height="432" alt="image" src="https://github.com/user-attachments/assets/c4b57639-b41d-4259-ba09-ca1440b6c351" />


#### 2.3. Visualize QC Metrics

Violin plots were generated to examine the distributions of gene counts, detected genes, and mitochondrial percentages across cells.

```r
VlnPlot(
  seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)
```
<img width="1440" height="1200" alt="image" src="https://github.com/user-attachments/assets/a222af21-8cd4-4c8d-a31d-fe5cd0005397" />


#### 2.4. Examine Relationships Between QC Metrics

Scatter plots were used to investigate relationships between total UMI counts, detected genes, and mitochondrial percentages.

```r
library(patchwork)

plot1 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

plot2 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

plot1 + plot2
```
<img width="1440" height="600" alt="image" src="https://github.com/user-attachments/assets/fe6ab457-0abc-4373-b211-3ce03479b078" />

**Interpretation**

The relationship between `nCount_RNA` and `nFeature_RNA` reflects the association between sequencing depth and gene detection. Mitochondrial percentages provide an additional indicator of potential cellular stress or damage.

These distributions will guide the selection of appropriate filtering thresholds to remove potentially low-quality cells while retaining biologically relevant populations.


#### 2.5. Filter Low-Quality Cells

Based on the QC metrics, cells were filtered using the following criteria:

- **`nFeature_RNA > 500`:** Remove cells with very few detected genes.
- **`nFeature_RNA < 5000`:** Exclude cells with unusually high gene counts, which may represent doublets or multiplets.
- **`percent.mt < 5`:** Remove cells with elevated mitochondrial gene expression, potentially indicating cellular stress or damage.

```r
# Number of cells before filtering
n_cells_before <- ncol(seurat)

# Filter cells based on QC thresholds
seurat <- subset(
  seurat,
  subset = nFeature_RNA > 500 &
           nFeature_RNA < 5000 &
           percent.mt < 5
)
```

**Interpretation**

QC filtering removes cells that do not meet the selected quality thresholds. The remaining cells are retained for downstream analyses, including normalization, dimensionality reduction, and clustering.

**Additional QC considerations**

Filtering based on gene counts and mitochondrial percentages may not identify all doublets or damaged cells. Additional methods, such as [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder), can be applied to detect potential doublets when necessary.

### Step 3. Normalization
<a id="top"></a>
Gene expression data were normalized to account for differences in sequencing depth and total UMI counts across cells. Seurat's `NormalizeData()` function was applied using the default **LogNormalize** method.

```r
# Normalize gene expression data
seurat <- NormalizeData(
  seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
```

**Normalization procedure:**

1. Divide each gene's UMI count by the total UMI counts in the corresponding cell.
2. Multiply the normalized values by a scaling factor of 10,000.
3. Apply a natural log transformation using `log1p()` to reduce the skewness of expression values.

**Results**

The normalized expression matrix is stored in the `data` layer of the RNA assay, while the original raw counts are retained in the `counts` layer. This allows gene expression levels to be compared more appropriately across cells in subsequent analyses.


### Step 4. Highly Variable Feature Selection
<a id="top"></a>


Highly variable genes (HVGs) were identified to capture biological heterogeneity across cells. Genes with substantial variation in expression are generally more informative for distinguishing cell populations than genes with relatively uniform expression.

The top **3,000 highly variable genes** were selected using Seurat's `FindVariableFeatures()` function with the default `vst` method.

```r
# Identify the top 3,000 highly variable genes
seurat <- FindVariableFeatures(
  seurat,
  selection.method = "vst",
  nfeatures = 3000
)
```

```bash
# Extract the top 20 highly variable genes
top_features <- head(VariableFeatures(seurat), 20)

# Visualize highly variable genes
plot1 <- VariableFeaturePlot(seurat)

# Label the top 20 genes
plot2 <- LabelPoints(
  plot = plot1,
  points = top_features,
  repel = TRUE
)

# Display plots
plot1 + plot2
```
<img width="1440" height="600" alt="image" src="https://github.com/user-attachments/assets/d0555706-426a-4454-9762-062a81bf50c4" />

**Results and Interpretation**

The variable feature plot illustrates the relationship between average gene expression and standardized variance across cells.

- **Highly variable genes (HVGs):** Genes showing greater expression variability than expected for their average expression level.
- **Non-variable genes:** Genes with relatively lower standardized variability.
- **Top 20 HVGs:** The most highly variable genes, highlighted and labeled in the plot.

The selected 3,000 HVGs will be used in subsequent dimensionality reduction and clustering analyses to characterize cellular heterogeneity.


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


