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

Data scaling was performed to standardize gene expression values across cells, preventing highly expressed genes from disproportionately influencing downstream analyses.

Seurat's `ScaleData()` function was applied to the selected highly variable genes (HVGs), centering and scaling each gene to approximately zero mean and unit variance.

```r
# Scale the selected highly variable genes
seurat <- ScaleData(seurat)
```

**Results and Interpretation**

The scaled expression values are stored in the `scale.data` layer of the RNA assay and will be used for principal component analysis (PCA).

Scaling ensures that genes with different expression levels contribute more comparably to dimensionality reduction and subsequent clustering.

**Optional: Regression of Unwanted Variation**

Technical factors, such as the number of detected genes (`nFeature_RNA`) and mitochondrial percentage (`percent.mt`), can be regressed out if they strongly influence the observed cellular variation.

```r
# Optional: Scale data with regression
seurat <- ScaleData(
  seurat,
  vars.to.regress = c("nFeature_RNA", "percent.mt")
)
```

However, regression was not applied in the initial analysis to avoid potentially removing biologically relevant variation.

---

### Optional: SCTransform Normalization
<a id="top"></a>

[SCTransform](https://satijalab.org/seurat/articles/sctransform_vignette) provides an alternative preprocessing approach based on regularized negative binomial regression.

Unlike the standard LogNormalize workflow, SCTransform combines normalization, variance stabilization, and highly variable feature selection into a single procedure.

```r
# Alternative preprocessing using SCTransform
seurat_sct <- SCTransform(
  seurat,
  assay = "RNA",
  variable.features.n = 3000,
  verbose = FALSE
)
```

**Note:** SCTransform is an alternative to Steps 3–5 and is not required when following the standard LogNormalize workflow. It uses raw RNA counts and generates an SCT assay for downstream analysis.



### Step 6. Principal Component Analysis (PCA)
<a id="top"></a>

Principal component analysis (PCA) was performed on the scaled highly variable genes to reduce data dimensionality while preserving major patterns of variation across cells.

PCA summarizes correlated gene expression patterns into principal components (PCs), providing a lower-dimensional representation of cellular heterogeneity for downstream clustering and visualization.

#### 6.1. Run PCA

The first 50 principal components were calculated using Seurat's `RunPCA()` function.

```r
# Perform PCA using highly variable genes
seurat <- RunPCA(
  seurat,
  npcs = 50
)


# Examine genes contributing to the first PCs
print(
  seurat[["pca"]],
  dims = 1:5,
  nfeatures = 5
)
```
<img width="1462" height="978" alt="image" src="https://github.com/user-attachments/assets/d65f4b99-7c9a-4c80-a8fe-fd7679103127" />

#### 6.2. Determine the Number of Principal Components

An elbow plot was generated to examine the standard deviation associated with each principal component and guide the selection of PCs for downstream analyses.

```r
# Visualize the first 50 principal components
ElbowPlot(
  seurat,
  ndims = 50
)
```
<img width="960" height="960" alt="image" src="https://github.com/user-attachments/assets/2c8e7803-9c5a-4b1a-ac73-aa32bd7a0966" />

**Interpretation**

The elbow plot shows how the standard deviation decreases across successive PCs.

- **Early PCs:** Capture the largest sources of variation in gene expression.
- **Later PCs:** Explain progressively smaller amounts of variation.
- **Elbow point:** Indicates where the decrease becomes less pronounced, providing a heuristic for selecting PCs.

The final number of PCs should be determined by considering both the elbow plot and the biological relevance of the variation captured by individual components.

#### 6.3. Visualize PCA Gene Contributions

PCA heatmaps were generated to examine genes contributing to each principal component and explore the underlying sources of cellular heterogeneity.

The first 20 PCs were visualized in groups of four, displaying 20 contributing genes and 500 selected cells per heatmap.

```r
library(ggplot2)
library(patchwork)

options(repr.plot.width = 8, repr.plot.height = 8)

# Visualize PC1–PC20 in groups of four
for (start in seq(1, 17, by = 4)) {

  p <- PCHeatmap(
    seurat,
    dims = start:(start + 3),
    nfeatures = 20,
    cells = 500,
    balanced = TRUE,
    ncol = 2,
    fast = FALSE
  )

  # Customize heatmap colors and appearance
  p <- p &
    scale_fill_gradient2(
      low = "#58a0cc",
      mid = "#FFFFFF",
      high = "#db1738",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      name = "Expression"
    ) &
    theme(
      plot.title = element_text(
        size = 14,
        face = "bold"
      ),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  print(p)
}
```

**Interpretation**

PCA heatmaps display patterns of scaled gene expression among cells ordered by their PC scores.

- **Rows:** Genes contributing strongly to each PC.
- **Columns:** Cells ordered according to their PCA scores.
- **Colors:** Relative scaled gene expression levels.
- **Gene patterns:** Groups of genes with coordinated expression may reflect distinct biological processes or cellular populations.

Examining the gene contributions helps assess whether individual PCs capture biologically meaningful variation and supports the selection of PCs for downstream clustering.


<img width="960" height="960" alt="image" src="https://github.com/user-attachments/assets/d518d6bf-fa8b-4dcc-b3d6-be7b29791d93" />
<img width="960" height="960" alt="image" src="https://github.com/user-attachments/assets/335d69a5-1274-42df-99f7-4eebb2c336ff" />
<img width="960" height="960" alt="image" src="https://github.com/user-attachments/assets/33ce75b2-f1da-4cea-b8ab-844034ddb116" />
<img width="960" height="960" alt="image" src="https://github.com/user-attachments/assets/f97869bd-a289-4877-81d7-7bd92616b638" />
<img width="960" height="960" alt="image" src="https://github.com/user-attachments/assets/3107b7ea-02b8-44ee-aff2-ba47fed2770b" />


### Step 7. Non-linear Dimensionality Reduction
<a id="top"></a>

Non-linear dimensionality reduction was performed using **t-distributed Stochastic Neighbor Embedding (t-SNE)** and **Uniform Manifold Approximation and Projection (UMAP)** to visualize cellular heterogeneity in two-dimensional space.

Both methods use the principal components obtained from PCA to generate low-dimensional representations of cells while aiming to preserve local neighborhood relationships.

#### 7.1. Run t-SNE and UMAP

The first 20 principal components were used to construct t-SNE and UMAP embeddings.

```r
# Run t-SNE using the first 20 PCs
seurat <- RunTSNE(
  seurat,
  dims = 1:20
)

# Run UMAP using the first 20 PCs
seurat <- RunUMAP(
  seurat,
  dims = 1:20
)
```

#### 7.2. Visualize t-SNE and UMAP

The resulting embeddings were visualized side by side to compare the spatial organization of cells.

```r
library(ggplot2)
library(patchwork)

# Set figure dimensions
options(repr.plot.width = 12, repr.plot.height = 6)

# t-SNE visualization
plot1 <- TSNEPlot(seurat) +
  ggtitle("t-SNE") +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  )

# UMAP visualization
plot2 <- UMAPPlot(seurat) +
  ggtitle("UMAP") +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  )

# Display plots
plot1 + plot2
```
<img width="1800" height="2880" alt="image" src="https://github.com/user-attachments/assets/755ba11e-e061-4899-b0ae-5a6f7cf7280c" />



**Results and Interpretation**

Both t-SNE and UMAP provide two-dimensional representations of cellular heterogeneity.

- **Each point** represents an individual cell.
- **Nearby cells** generally exhibit similar gene expression profiles in the reduced PCA space.
- **Separated groups** may indicate transcriptionally distinct cell populations.
- **t-SNE** emphasizes local neighborhood relationships.
- **UMAP** also preserves local structure and may reveal broader patterns of cellular organization.

These visualizations provide an initial overview of cellular heterogeneity before graph-based clustering and cell type annotation.

#### 7.3. Visualization of Cell Type Markers

To explore potential cell identities and developmental states, the expression patterns of selected canonical marker genes were visualized on both t-SNE and UMAP embeddings.

These markers represent cell proliferation, neural progenitors, immature neurons, and different regions of the developing human brain.

```r
library(patchwork)

# Set figure dimensions
options(repr.plot.width = 15, repr.plot.height = 24)

# Selected neural and developmental markers
genes <- c(
  "MKI67", "NES", "DCX",
  "FOXG1", "DLX2", "EMX1",
  "OTX2", "LHX9", "TFAP2A"
)

<img width="1800" height="2880" alt="image" src="https://github.com/user-attachments/assets/6a53e866-fa3c-4e64-9aea-3e71218fa3a8" />
# Visualize marker expression on t-SNE
plot1 <- FeaturePlot(
  seurat,
  features = genes,
  ncol = 3,
  reduction = "tsne"
)

# Visualize marker expression on UMAP
plot2 <- FeaturePlot(
  seurat,
  features = genes,
  ncol = 3,
  reduction = "umap"
)

# Display both embeddings
plot1 / plot2
```

**Selected Marker Genes**

| Gene | Biological Association |
|------|------------------------|
| MKI67 | Cell proliferation |
| NES | Neural progenitor cells |
| DCX | Immature neurons |
| FOXG1 | Telencephalon |
| DLX2 | Ventral telencephalon |
| EMX1 | Dorsal telencephalon (cortex) |
| OTX2 | Midbrain and diencephalon development |
| LHX9 | Diencephalic and midbrain neuronal populations |
| TFAP2A | Midbrain–hindbrain boundary-associated populations |

**Results and Interpretation**

Feature plots reveal the spatial distribution of marker gene expression across cells in the two-dimensional embeddings.

- **MKI67 and NES:** Help identify proliferating cells and neural progenitor populations.
- **DCX:** Highlights populations associated with immature neurons.
- **FOXG1, DLX2, and EMX1:** Provide information about telencephalic regional identities.
- **OTX2, LHX9, and TFAP2A:** Help explore additional regional identities during brain development.

Comparing marker expression patterns across t-SNE and UMAP provides preliminary insights into cellular heterogeneity and potential developmental states.

These observations will support subsequent cell clustering and biological annotation.


### Step 8. Cell Clustering
<a id="top"></a>

Graph-based clustering was performed to identify transcriptionally distinct cell populations based on their gene expression profiles.

Seurat constructs a Shared Nearest Neighbor (SNN) graph using PCA embeddings and applies a community detection algorithm to identify groups of cells with similar transcriptional characteristics.

#### 8.1. Construct the Nearest Neighbor Graph

The first 20 principal components were used to construct a cell–cell nearest neighbor graph.

```r
# Construct the nearest neighbor graph
seurat <- FindNeighbors(
  seurat,
  dims = 1:20
)
```

#### 8.2. Identify Cell Clusters

Graph-based clustering was performed using the Louvain algorithm with a resolution of 1.0.

```r
# Perform graph-based clustering
seurat <- FindClusters(
  seurat,
  resolution = 1.0
)
```

The `resolution` parameter controls clustering granularity:

- **Lower resolution:** Produces fewer, larger clusters.
- **Higher resolution:** Produces more, smaller clusters.

A resolution of 1.0 was used to explore finer cellular heterogeneity.

#### 8.3. Visualize Cell Clusters

The identified clusters were visualized using the previously generated t-SNE and UMAP embeddings.

```r
library(ggplot2)
library(patchwork)

# Set figure dimensions
options(repr.plot.width = 12, repr.plot.height = 6)

# t-SNE clustering visualization
plot1 <- DimPlot(
  seurat,
  reduction = "tsne",
  label = TRUE
) + ggtitle("t-SNE")

# UMAP clustering visualization
plot2 <- DimPlot(
  seurat,
  reduction = "umap",
  label = TRUE
) + ggtitle("UMAP")

# Display plots side by side
plot1 + plot2
```
<img width="1440" height="720" alt="image" src="https://github.com/user-attachments/assets/2ef322ff-fa0f-41ae-9a37-f342f6f6b9fd" />


#### 8.4. Examine Cluster Distribution

The number of cells assigned to each cluster was calculated to examine cluster sizes.

```r
# Count cells in each cluster
table(Idents(seurat))
```
<img width="1194" height="112" alt="image" src="https://github.com/user-attachments/assets/33b823a0-e234-4145-8369-ce23d63cc59d" />

**Results and Interpretation**

Graph-based clustering identifies groups of cells with similar transcriptomic profiles.

- **Each point:** Represents an individual cell.
- **Colors:** Indicate different clusters.
- **Cluster labels:** Represent cluster identities assigned by Seurat.
- **Cluster sizes:** Indicate the number of cells assigned to each cluster.

The resulting clusters provide a basis for identifying distinct cell populations and exploring cellular heterogeneity.

However, computational clusters do not necessarily correspond directly to biological cell types. Further characterization using differentially expressed genes and canonical cell type markers is required for biological annotation.


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


