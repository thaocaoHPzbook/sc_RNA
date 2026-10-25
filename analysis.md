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

Cell type annotation was performed by examining the expression of canonical marker genes across the identified clusters.

Marker genes associated with neural progenitors, neuronal differentiation, and regional identities of the developing brain were selected to characterize cellular populations in the human cerebral organoid dataset.

#### 9.1. Select Canonical Marker Genes

A panel of established neural and developmental markers was selected to examine potential cell identities across clusters.

```r
# Define canonical cell type markers
ct_markers <- c(
  # Proliferation, neural progenitors, and neurons
  "MKI67", "NES", "DCX", "FOXG1",

  # Ventral telencephalon
  "DLX2", "DLX5", "ISL1", "SIX3",
  "NKX2-1", "SOX6", "NR2F2",

  # Dorsal telencephalon
  "EMX1", "PAX6", "GLI3",
  "EOMES", "NEUROD6",

  # Other brain regions and neuronal populations
  "RSPO3", "OTX2", "LHX9",
  "TFAP2A", "RELN", "HOXB2", "HOXB5"
)

# Check marker availability
available_markers <- intersect(
  ct_markers,
  rownames(seurat)
)

# Identify missing markers
missing_markers <- setdiff(
  ct_markers,
  available_markers
)

print(missing_markers)
```

#### 9.2. Visualize Marker Expression Across Clusters

A heatmap was generated to examine the expression patterns of canonical markers across the identified cell clusters.

```r
library(ggplot2)

# Ensure selected markers are scaled
seurat <- ScaleData(
  seurat,
  features = unique(c(
    VariableFeatures(seurat),
    available_markers
  )),
  verbose = FALSE
)

# Set figure dimensions
options(repr.plot.width = 16, repr.plot.height = 9)

# Visualize marker expression
DoHeatmap(
  seurat,
  features = available_markers,
  group.by = "seurat_clusters",
  size = 4
) +
  scale_fill_gradient2(
    low = "#0973b5",
    mid = "#f6f9f6",
    high = "#d32a46",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    name = "Expression"
  ) +
  theme(
    axis.text.y = element_text(
      size = 11,
      face = "italic"
    ),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
```
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/89f1426a-d9d2-4ea4-9554-96b71b1f3807" />


#### 9.3. Biological Interpretation

The heatmap displays the relative expression of selected marker genes across cell clusters.

| Marker genes | Potential biological identity |
|---|---|
| MKI67 | Proliferating cells |
| NES, PAX6 | Neural progenitor populations |
| DCX | Immature neurons |
| FOXG1 | Telencephalic populations |
| DLX2, DLX5, NKX2-1 | Ventral telencephalic populations |
| EMX1, EOMES, NEUROD6 | Dorsal telencephalic and cortical populations |
| RSPO3, OTX2, LHX9 | Diencephalic or midbrain-associated populations |
| RELN | Cajal–Retzius-associated populations |
| HOXB2, HOXB5 | Hindbrain-associated populations |

**Results and Interpretation**

The expression patterns of canonical markers provide preliminary evidence for the biological identities of individual cell clusters.

Clusters showing coordinated expression of neural progenitor markers may represent progenitor populations, whereas clusters enriched for neuronal differentiation markers may correspond to more differentiated neuronal states.

Cell type annotation should be based on combinations of markers rather than individual genes, as many markers are shared across multiple cell populations.

These results provide a foundation for further annotation using cluster-specific differentially expressed genes and reference-based approaches.

#### 9.4. Identify Cluster-Specific Marker Genes

Differential expression analysis was performed to identify marker genes associated with each cell cluster.

Seurat's `FindAllMarkers()` function was used to compare gene expression in each cluster against all remaining cells, allowing the identification of genes preferentially expressed in individual clusters.

```r
# Identify positive marker genes for each cluster
cl_markers <- FindAllMarkers(
  seurat,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = log2(1.2)
)

# Preview marker results
head(cl_markers)
```

**Parameters**

- **`only.pos = TRUE`:** Identify genes with higher expression in the target cluster.
- **`min.pct = 0.25`:** Consider genes detected in at least 25% of cells in either comparison group.
- **`logfc.threshold = log2(1.2)`:** Apply a minimum log2 fold-change threshold corresponding to a 1.2-fold difference.

#### 9.5. Filter Significant Markers and Summarize Cluster Distribution

Marker genes were filtered using an adjusted p-value threshold of 0.05. The number of significant marker genes was then calculated for each cluster.

```r
library(dplyr)

# Retain statistically significant markers
markers_filtered <- cl_markers %>%
  filter(p_val_adj < 0.05)

# Count significant markers in each cluster
marker_counts <- tibble(
  cluster = levels(Idents(seurat))
) %>%
  left_join(
    markers_filtered %>%
      mutate(cluster = as.character(cluster)) %>%
      count(cluster, name = "n_markers"),
    by = "cluster"
  ) %>%
  mutate(
    n_markers = coalesce(n_markers, 0L)
  )

# Display marker counts
print(marker_counts, n = Inf)
```
<img width="412" height="730" alt="image" src="https://github.com/user-attachments/assets/cdc00232-ae72-4957-bc28-29b01b282892" />

**Interpretation**

The number of significant markers provides an overview of transcriptional differences between each cluster and the remaining cells.

Clusters with numerous significant markers may exhibit more distinct transcriptional signatures. However, marker counts alone do not determine cell type identity or clustering quality.

#### 9.6. Visualize the Top 10 Marker Genes per Cluster

The top 10 significant marker genes from each cluster were selected based on their average log2 fold change.

A heatmap was generated to visualize their expression patterns across the identified cell populations.

```r
library(dplyr)
library(ggplot2)

# Select the top 10 markers per cluster
top10_markers <- markers_filtered %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

# Extract unique genes for visualization
genes_plot <- unique(top10_markers$gene)

# Scale selected genes in a separate Seurat object
seurat_heatmap <- ScaleData(
  seurat,
  assay = "RNA",
  features = genes_plot,
  verbose = FALSE
)

# Set figure dimensions
options(
  repr.plot.width = 16,
  repr.plot.height = max(
    8,
    length(genes_plot) * 0.18
  )
)

# Visualize cluster-specific marker expression
DoHeatmap(
  seurat_heatmap,
  assay = "RNA",
  features = genes_plot,
  slot = "scale.data",
  size = 4
) +
  scale_fill_gradient2(
    low = "#0785e0",
    mid = "#FFFFFF",
    high = "#d81639",
    midpoint = 0,
    name = "Scaled\nexpression"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
```
<img width="1920" height="3175" alt="image" src="https://github.com/user-attachments/assets/3760146b-ea1d-4902-95f4-62c45bc42bbc" />


**Results and Interpretation**

The heatmap illustrates the expression patterns of the top marker genes across cell clusters.

- **Rows:** Significant marker genes identified by differential expression analysis.
- **Columns:** Individual cells grouped by cluster identity.
- **Red:** Relatively higher scaled gene expression.
- **Blue:** Relatively lower scaled gene expression.
- **Cluster-specific patterns:** Groups of genes preferentially expressed in particular cell populations.

The identified marker genes provide a data-driven basis for characterizing transcriptionally distinct cell populations.

Combining these markers with established neural and developmental markers can support more reliable biological annotation of the identified clusters.


#### 9.7. Alternative Marker Identification Using Presto

As an alternative to Seurat's `FindAllMarkers()`, cluster-specific marker genes were identified using the [Presto](https://github.com/immunogenomics/presto) package.

Presto provides an efficient implementation of the Wilcoxon rank-sum test for differential expression analysis in single-cell datasets. Each cluster was compared against all remaining cells using log-normalized RNA expression data.

**1. Identify Marker Genes**

```r
library(presto)
library(dplyr)

# Identify markers using log-normalized RNA data
cl_markers_presto <- wilcoxauc(
  seurat,
  assay = "data",
  seurat_assay = "RNA"
)

# Preview results
head(cl_markers_presto)
```


**2. Filter Significant Marker Genes**

Marker genes were filtered using expression differences, detection percentages, and adjusted p-values.

```r
# Filter significant markers
markers_filtered_presto <- cl_markers_presto %>%
  filter(
    logFC > log(1.2),
    pct_in > 20,
    padj < 0.05
  )
```

**Filtering criteria**

- **`logFC > log(1.2)`:** Retain genes showing a positive difference in mean log-normalized expression above the selected threshold.
- **`pct_in > 20`:** Retain genes detected in more than 20% of cells within the target cluster.
- **`padj < 0.05`:** Retain statistically significant genes after Benjamini–Hochberg correction.

**3. Count Significant Markers per Cluster**

The number of significant marker genes was calculated for each cluster, including clusters without markers meeting the filtering criteria.

```r
# Count significant markers per cluster
marker_counts_presto <- tibble(
  group = levels(Idents(seurat))
) %>%
  left_join(
    markers_filtered_presto %>%
      mutate(group = as.character(group)) %>%
      count(group, name = "n_markers"),
    by = "group"
  ) %>%
  mutate(
    n_markers = coalesce(n_markers, 0L)
  )

# Display marker counts
print(marker_counts_presto, n = Inf)
```
<img width="390" height="730" alt="image" src="https://github.com/user-attachments/assets/8f8df4f8-5289-49d7-92a1-9f7213d43118" />

**4. Select the Top 10 Markers per Cluster**

The top 10 significant markers were selected from each cluster based on their log fold-change values.

```r
# Select top 10 markers per cluster
top10_markers_presto <- markers_filtered_presto %>%
  group_by(group) %>%
  slice_max(
    order_by = logFC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

# Display selected markers
print(
  top10_markers_presto,
  n = Inf,
  width = Inf
)
```
<img width="1234" height="994" alt="image" src="https://github.com/user-attachments/assets/a7c7cd5f-e0b1-4614-82cb-3444ea3fb918" />

**Results and Interpretation**

Presto identifies genes with differential expression patterns across cell clusters, providing an alternative approach to marker discovery.

The resulting marker lists can be used to:

- Identify transcriptional signatures associated with individual clusters.
- Compare marker genes identified using different analytical implementations.
- Support biological annotation through canonical markers and literature-based interpretation.

The selected markers provide additional evidence for characterizing cellular heterogeneity in the cerebral organoid dataset.


#### 9.8. Visualize Presto-Derived Marker Genes

A heatmap was generated to visualize the expression patterns of the top 10 marker genes identified by Presto for each cell cluster.

Marker genes were selected based on their log fold-change values and scaled separately for visualization while preserving the original Seurat object.

```r
library(ggplot2)

# Extract unique marker genes from Presto results
genes_plot_presto <- unique(
  top10_markers_presto$feature
)

# Scale selected genes in a separate Seurat object
seurat_heatmap_presto <- ScaleData(
  seurat,
  assay = "RNA",
  features = genes_plot_presto,
  verbose = FALSE
)

# Set figure dimensions
options(
  repr.plot.width = 16,
  repr.plot.height = max(
    8,
    length(genes_plot_presto) * 0.18
  )
)

# Generate marker expression heatmap
heatmap_presto <- DoHeatmap(
  seurat_heatmap_presto,
  assay = "RNA",
  features = genes_plot_presto,
  slot = "scale.data",
  size = 4
) +
  scale_fill_gradient2(
    low = "#0887e1",
    mid = "#FFFFFF",
    high = "#e01339",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Scaled\nexpression"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Display heatmap
print(heatmap_presto)
```
<img width="1920" height="1987" alt="image" src="https://github.com/user-attachments/assets/dce7fa3f-7ad1-414b-bd83-a48bafbfb79c" />


**Results and Interpretation**

The heatmap illustrates the expression patterns of the top marker genes identified using Presto across cell clusters.

- **Rows:** Marker genes selected from the Presto differential expression results.
- **Columns:** Individual cells grouped by cluster identity.
- **Red:** Higher relative scaled expression.
- **Blue:** Lower relative scaled expression.
- **Cluster-specific patterns:** Genes preferentially expressed in particular cell populations.

These expression patterns provide additional information for identifying transcriptionally distinct cell populations and supporting subsequent biological annotation.

The Presto-derived markers can also be compared with those identified using Seurat's `FindAllMarkers()` to examine similarities and differences between the two marker discovery approaches.

#### 9.9. Detailed Visualization of Selected Marker Genes

The expression patterns of **NEUROD2** and **NEUROD6**, two neuronal differentiation-associated genes, were further examined using UMAP feature plots and violin plots.

These genes were selected to investigate their potential enrichment in Cluster 2 and assess their contribution to the transcriptional identity of this cell population.

```r
library(ggplot2)
library(patchwork)

# Set figure dimensions
options(repr.plot.width = 16, repr.plot.height = 9)

# Generate UMAP feature plots
fp <- FeaturePlot(
  seurat,
  features = c("NEUROD2", "NEUROD6"),
  reduction = "umap",
  cols = c("#DADADA", "#8773e0"),
  pt.size = 0.5,
  order = TRUE,
  combine = FALSE
)

# Preserve UMAP aspect ratio
fp[[1]] <- fp[[1]] + coord_fixed()
fp[[2]] <- fp[[2]] + coord_fixed()

# Generate violin plots
vp <- VlnPlot(
  seurat,
  features = c("NEUROD2", "NEUROD6"),
  pt.size = 0,
  combine = FALSE
)

# Remove redundant legends
vp[[1]] <- vp[[1]] + NoLegend()
vp[[2]] <- vp[[2]] + NoLegend()

# Arrange UMAP and violin plots
p <- wrap_plots(
  A = fp[[1]],
  B = vp[[1]],
  C = vp[[2]],
  D = fp[[2]],
  design = "ABC\nDBC"
) &
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold"
    ),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.text = element_text(size = 11)
  )

# Display combined visualization
print(p)
```

**Results and Interpretation**

The combined visualization provides complementary information about the expression of NEUROD2 and NEUROD6.

- **UMAP FeaturePlots:** Display the spatial distribution of gene expression across individual cells.
- **Violin Plots:** Show the distribution of gene expression levels across the identified clusters.
- **NEUROD2 and NEUROD6:** Genes associated with neuronal differentiation, particularly in excitatory neuronal populations.

Enrichment of these markers in Cluster 2 would support its potential neuronal identity. However, additional canonical markers should be examined before assigning a definitive cell type annotation.
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/68d43f9d-7654-4182-a9da-a673b324e5a3" />


#### 8.5. Cluster Distribution and Visualization

The clustering results were further examined to determine the number of identified clusters, their cell distributions, and the clustering resolutions stored in the Seurat object.

To facilitate comparison between dimensionality reduction methods, consistent cluster colors were applied to both t-SNE and UMAP visualizations.

**1. Inspect Cluster Assignments**

```r
# Number of cells in each cluster
table(Idents(seurat))

# Total number of clusters containing cells
length(unique(Idents(seurat)))

# List stored clustering resolutions
grep(
  "snn_res",
  colnames(seurat[[]]),
  value = TRUE
)
```

<img width="1138" height="280" alt="image" src="https://github.com/user-attachments/assets/84a06ce5-6de8-47b0-9f33-b0e9c82d3017" />


**2. Visualize Clusters on t-SNE and UMAP**

```r
library(ggplot2)
library(patchwork)

# Set figure dimensions
options(repr.plot.width = 14, repr.plot.height = 6)

# Assign consistent colors to clusters
cluster_cols <- scales::hue_pal()(
  length(levels(Idents(seurat)))
)

names(cluster_cols) <- levels(Idents(seurat))

# t-SNE visualization
plot_tsne <- DimPlot(
  seurat,
  reduction = "tsne",
  cols = cluster_cols,
  pt.size = 0.3,
  label = TRUE,
  label.size = 3,
  repel = TRUE
) +
  ggtitle(NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

# UMAP visualization
plot_umap <- DimPlot(
  seurat,
  reduction = "umap",
  cols = cluster_cols,
  pt.size = 0.3,
  label = TRUE,
  label.size = 3,
  repel = TRUE
) +
  ggtitle(NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

# Display both plots
plot_tsne + plot_umap
```
<img width="1680" height="720" alt="image" src="https://github.com/user-attachments/assets/ba3ce7db-c9e1-49aa-abc2-14917e68c7bc" />


**Results and Interpretation**

The clustering results were visualized using both t-SNE and UMAP embeddings, with consistent colors representing the same clusters across both plots.

- **Cluster size:** Indicates the number of cells assigned to each cluster.
- **Cluster number:** Represents the total number of computationally identified cell populations.
- **Cluster distribution:** Shows the organization of transcriptionally similar cells in the reduced-dimensional space.
- **Consistent colors:** Facilitate direct comparison of cluster assignments between t-SNE and UMAP.

Although t-SNE and UMAP may display different spatial arrangements, the cluster identities remain unchanged because both visualizations use the same clustering results.

These findings provide an overview of cellular heterogeneity and establish a foundation for subsequent marker gene analysis and cell type annotation.

#### 9.10. Reference-Based Annotation Using VoxHunt

To further investigate the regional identities of cerebral organoid cell populations, transcriptomic profiles were compared with spatial gene expression patterns from the developing mouse brain using [VoxHunt](https://github.com/quadbio/VoxHunt).

VoxHunt calculates transcriptomic similarities between single-cell or cluster-level expression profiles and in situ hybridization data from the Allen Brain Atlas.

**1. Load Reference Data**

```r
library(voxhunt)

# Load the downloaded Allen Brain Atlas data
load_aba_data("path/to/ABA_data/voxhunt_data")

# Select 300 highly variable genes from the E13.5 atlas
genes_use <- variable_genes("E13", 300)$gene
```

**2. Prepare Expression Data**

The log-normalized RNA expression matrix and cluster assignments were extracted from the Seurat object.

```r
# Extract normalized RNA expression
expr <- SeuratObject::GetAssayData(
  seurat,
  assay = "RNA",
  layer = "data"
)

# Extract cluster assignments in cell order
groups <- seurat[[]][
  colnames(expr),
  "seurat_clusters"
]
```

**3. Map Cell Clusters to the Reference Atlas**

```r
# Perform reference-based spatial mapping
vox_map <- voxhunt::voxel_map(
  object = Matrix::t(expr),
  stage = "E13",
  groups = groups,
  genes_use = genes_use,
  pseudobulk_groups = TRUE
)

# Visualize spatial similarity maps
options(repr.plot.width = 14, repr.plot.height = 10)

voxhunt::plot_map(vox_map)
```
<img width="1680" height="1200" alt="image" src="https://github.com/user-attachments/assets/8adbd53d-d355-448a-b825-d599c9acc95b" />

**Results and Interpretation**

VoxHunt generates spatial similarity maps between cerebral organoid cell clusters and developing mouse brain regions.

- **Reference:** Allen Brain Atlas, embryonic mouse brain (E13.5).
- **Gene selection:** 300 variable genes from the reference atlas.
- **Mapping approach:** Transcriptomic correlation between cluster-level expression profiles and spatial reference patterns.
- **Biological interpretation:** Regions showing greater transcriptomic similarity provide supporting evidence for potential regional identities of organoid cell populations.

These results complement canonical marker analysis and help refine regional cell identity annotations.

Because the reference is derived from mouse brain tissue, spatial similarities should be interpreted cautiously and considered alongside human neural marker expression.


#### 9.11. Assign Cell Type Annotations

Based on canonical marker expression, cluster-specific genes, and reference-based analysis, biological identities were assigned to the identified cell clusters.

The original cluster identifiers were preserved in `seurat_clusters`, while the annotated cell identities were stored separately for downstream analyses.

```r
library(Seurat)

# Restore original cluster identities
Idents(seurat) <- "seurat_clusters"

# Define cell type annotations
new_ident <- c(
  "0"  = "Dorsal telen. NPC",
  "1"  = "Midbrain-hindbrain boundary neuron",
  "2"  = "Dorsal telen. neuron",
  "3"  = "Dien. and midbrain excitatory neuron",
  "4"  = "MGE-like neuron",
  "5"  = "G2M dorsal telen. NPC",
  "6"  = "Dorsal telen. IP",
  "7"  = "Dien. and midbrain NPC",
  "8"  = "Dien. and midbrain IP and excitatory early neuron",
  "9"  = "G2M Dien. and midbrain NPC",
  "10" = "G2M dorsal telen. NPC",
  "11" = "Dien. and midbrain inhibitory neuron",
  "12" = "Dien. and midbrain IP and early inhibitory neuron",
  "13" = "Ventral telen. neuron",
  "14" = "Unknown 1",
  "15" = "Unknown 2"
)

# Verify that all clusters have annotations
stopifnot(
  all(levels(seurat) %in% names(new_ident))
)

# Rename identities
seurat <- do.call(
  RenameIdents,
  c(list(object = seurat), as.list(new_ident))
)

# Store annotations in metadata
seurat$annotation_tutorial <- Idents(seurat)

# Verify annotations
table(seurat$annotation_tutorial)

# Visualize annotated cell populations
options(repr.plot.width = 16, repr.plot.height = 10)

DimPlot(
  seurat,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.size = 3.5,
  pt.size = 0.3
) +
  NoLegend()
```
<img width="1920" height="1200" alt="image" src="https://github.com/user-attachments/assets/dc616d31-7abb-47e9-8d9a-7d444445b396" />


**Results and Interpretation**

The UMAP visualization displays cell populations labeled according to their proposed biological identities rather than numerical cluster identifiers.

These annotations provide an overview of the cellular composition of the cerebral organoid dataset, including neural progenitors, intermediate progenitors, and differentiated neuronal populations.

The annotated identities will be used for downstream analyses of cellular differentiation and developmental trajectories.

### Step 10. Pseudotemporal Cell Ordering
<a id="top"></a>

Pseudotime analysis was performed to investigate potential developmental trajectories among dorsal telencephalic cell populations.

Based on the previous cell type annotation, clusters 0, 2, 5, 6, and 10 were considered potential dorsal telencephalic populations, including neural progenitors, intermediate progenitors, and differentiated neurons.

These populations may represent different stages of neuronal differentiation. To explore this possibility, diffusion-based trajectory analysis was performed using the [destiny](https://bioconductor.org/packages/destiny/) R package.

#### 10.1. Subset Dorsal Telencephalic Cells

Cells belonging to the selected dorsal telencephalic clusters were extracted for subsequent trajectory analysis.

```r
# Subset dorsal telencephalic cell populations
seurat_dorsal <- subset(
  seurat,
  subset = seurat_clusters %in% c(0, 2, 5, 6, 10)
)

# Inspect the subset
seurat_dorsal

# Count cells in each original cluster
table(seurat_dorsal$seurat_clusters)
```

#### 10.2. Re-identify Highly Variable Genes

Highly variable genes were re-identified within the selected dorsal telencephalic populations.

This step focuses feature selection on transcriptional variation among the selected cells rather than differences between unrelated cell populations.

```r
# Identify highly variable genes in the subset
seurat_dorsal <- FindVariableFeatures(
  seurat_dorsal,
  selection.method = "vst",
  nfeatures = 3000
)

# Examine the selected features
head(VariableFeatures(seurat_dorsal), 20)
```

**Results and Interpretation**

The subset contains dorsal telencephalic populations selected for investigating potential neuronal differentiation trajectories.

Re-identifying highly variable genes allows subsequent analyses to focus on transcriptional differences within these populations, potentially capturing variation associated with progenitor states, neuronal differentiation, and maturation.

These selected cells will be used for diffusion map construction and pseudotime inference.

#### 10.3. Reduce Cell Cycle Effects

Cell cycle-related variation may influence the separation of neural progenitor populations and potentially obscure transcriptional changes associated with neuronal differentiation.

To reduce this effect, known cell cycle genes were excluded from the highly variable gene list before repeating dimensionality reduction.

**1. Exclude Cell Cycle Genes**

```r
# Number of HVGs before filtering
n_hvg_before <- length(
  VariableFeatures(seurat_dorsal)
)

# Remove cell cycle-related genes from HVGs
VariableFeatures(seurat_dorsal) <- setdiff(
  VariableFeatures(seurat_dorsal),
  unlist(cc.genes)
)

# Number of HVGs after filtering
n_hvg_after <- length(
  VariableFeatures(seurat_dorsal)
)

# Check the number of excluded genes
n_hvg_before - n_hvg_after
```

**2. Recompute Scaling and Dimensionality Reduction**

Following the exclusion of cell cycle genes, the selected features were rescaled, and PCA and UMAP were recalculated.

```r
# Scale the remaining highly variable genes
seurat_dorsal <- ScaleData(
  seurat_dorsal,
  features = VariableFeatures(seurat_dorsal)
)

# Recompute PCA
seurat_dorsal <- RunPCA(
  seurat_dorsal,
  npcs = 50
)

# Recompute UMAP using the first 20 PCs
seurat_dorsal <- RunUMAP(
  seurat_dorsal,
  dims = 1:20
)
```
<img width="1474" height="876" alt="image" src="https://github.com/user-attachments/assets/f64b4b51-cb59-4e30-9ae7-69c46e0fd1d7" />

**Results and Interpretation**

Cell cycle-associated genes were excluded from the selected highly variable features to reduce their contribution to subsequent dimensionality reduction.

The updated PCA and UMAP embeddings provide an alternative representation of transcriptional variation among dorsal telencephalic populations.

This approach may help emphasize expression differences associated with neuronal differentiation and maturation rather than cell cycle-related variation.

The resulting embeddings will be used to explore potential developmental trajectories through diffusion-based pseudotime analysis.


#### 10.4. Visualization of Developmental Markers

The expression patterns of selected developmental markers were visualized on the updated UMAP embedding to explore potential transcriptional transitions among dorsal telencephalic cell populations.

Four genes representing different cellular and developmental states were examined.

```r
# Set figure dimensions
options(repr.plot.width = 18, repr.plot.height = 4)

# Visualize developmental markers on the updated UMAP
FeaturePlot(
  seurat_dorsal,
  features = c("MKI67", "GLI3", "EOMES", "NEUROD6"),
  reduction = "umap",
  ncol = 4
)
```
<img width="1500" height="300" alt="image" src="https://github.com/user-attachments/assets/52e6c162-4e96-4ade-bc32-537e159ac33e" />

**Selected Marker Genes**

| Gene | Biological Association |
|---|---|
| MKI67 | Cell proliferation |
| GLI3 | Dorsal forebrain development and patterning |
| EOMES | Cortical intermediate progenitors |
| NEUROD6 | Neuronal differentiation and maturation |

**Results and Interpretation**

The feature plots illustrate the distribution of developmental marker expression across dorsal telencephalic populations.

- **MKI67:** Helps identify proliferating cells.
- **GLI3:** Provides information about dorsal forebrain developmental identity.
- **EOMES:** Highlights intermediate progenitor populations.
- **NEUROD6:** Indicates populations associated with neuronal differentiation.

Distinct or gradually changing expression patterns across the UMAP may suggest transcriptional transitions between progenitor and neuronal states.

These observations provide a preliminary basis for subsequent diffusion map and pseudotime analyses.


#### 10.5. Cell Cycle Regression

Although cell cycle-associated genes were excluded from the highly variable features, residual cell cycle-related variation may still influence the organization of dorsal telencephalic populations.

To further reduce this effect, cell cycle scores were calculated for individual cells and regressed out during data scaling.

**1. Calculate Cell Cycle Scores**

Seurat's `CellCycleScoring()` function was used to estimate S-phase and G2/M-phase scores based on established cell cycle marker genes.

```r
# Calculate cell cycle scores
seurat_dorsal <- CellCycleScoring(
  seurat_dorsal,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

# Examine cell cycle scores
head(
  seurat_dorsal[[]][
    , c("S.Score", "G2M.Score", "Phase")
  ]
)
```

**2. Regress Out Cell Cycle Effects**

The calculated cell cycle scores were included as covariates during scaling to reduce their linear contribution to gene expression variation.

```r
# Scale data while regressing out cell cycle scores
seurat_dorsal <- ScaleData(
  seurat_dorsal,
  vars.to.regress = c("S.Score", "G2M.Score")
)
```

**3. Recompute PCA and UMAP**

PCA and UMAP were recalculated using the cell cycle-adjusted expression data.

```r
# Recompute PCA
seurat_dorsal <- RunPCA(
  seurat_dorsal,
  npcs = 50
)

# Recompute UMAP
seurat_dorsal <- RunUMAP(
  seurat_dorsal,
  dims = 1:20
)
```

**4. Examine Developmental Markers**

```r
# Set figure dimensions
options(repr.plot.width = 18, repr.plot.height = 4)

# Visualize developmental marker expression
FeaturePlot(
  seurat_dorsal,
  features = c("MKI67", "GLI3", "EOMES", "NEUROD6"),
  reduction = "umap",
  ncol = 4
)
```
<img width="1500" height="300" alt="image" src="https://github.com/user-attachments/assets/26ea2c02-9d65-4027-b4e0-373d91f0e2a9" />


**Results and Interpretation**

Cell cycle regression was applied to reduce variation associated with S-phase and G2/M-phase activity.

The updated UMAP provides a representation of dorsal telencephalic populations after accounting for these cell cycle-related effects.

Developmental marker expression was examined to assess whether the resulting organization is consistent with potential transitions between neural progenitors, intermediate progenitors, and differentiated neurons.

---

#### 10.6. Diffusion Map and Pseudotime Analysis

Diffusion map analysis was performed using the [destiny](https://bioconductor.org/packages/destiny/) R package to investigate potential developmental relationships among dorsal telencephalic cell populations.

Diffusion pseudotime (DPT) was then used to estimate the relative ordering of cells along the inferred transcriptional trajectory.

**1. Construct the Diffusion Map**

The first 20 principal components from the cell cycle-adjusted analysis were used as input.

```r
library(destiny)

# Extract PCA embeddings
pca_embeddings <- Embeddings(
  seurat_dorsal,
  "pca"
)[, 1:20, drop = FALSE]

# Construct diffusion map
dm <- DiffusionMap(pca_embeddings)
```

**2. Calculate Diffusion Pseudotime**

```r
# Estimate diffusion pseudotime
dpt <- DPT(dm)

# Convert pseudotime to ranks
dpt_rank <- rank(dpt$dpt)

# Match pseudotime values to cell identifiers
names(dpt_rank) <- rownames(pca_embeddings)

# Store pseudotime in Seurat metadata
seurat_dorsal$dpt <- dpt_rank[
  colnames(seurat_dorsal)
]
```

**3. Visualize Pseudotime and Developmental Markers**

```r
# Set figure dimensions
options(repr.plot.width = 18, repr.plot.height = 4)

# Visualize pseudotime and marker expression
FeaturePlot(
  seurat_dorsal,
  features = c("dpt", "GLI3", "EOMES", "NEUROD6"),
  reduction = "umap",
  ncol = 4
)
```
<img width="1500" height="300" alt="image" src="https://github.com/user-attachments/assets/7c109271-5408-4684-a941-3441a771309f" />

**Results and Interpretation**

Diffusion map analysis provides a representation of transcriptional relationships among the selected dorsal telencephalic populations.

Diffusion pseudotime assigns a relative position to each cell along the inferred trajectory.

- **Early pseudotime:** Cells closer to the selected computational starting point.
- **Intermediate pseudotime:** Cells occupying intermediate positions along the trajectory.
- **Late pseudotime:** Cells farther from the starting point.

Expression patterns of GLI3, EOMES, and NEUROD6 were examined to evaluate whether the inferred ordering is consistent with potential transitions from neural progenitors to intermediate progenitors and differentiated neurons.

The resulting pseudotime represents an inferred transcriptional ordering rather than actual developmental time.
#### 10.7. Gene Expression Dynamics Along Pseudotime

To investigate transcriptional changes along the inferred developmental trajectory, the expression patterns of three representative genes (**GLI3, EOMES, and NEUROD6**) were examined across diffusion pseudotime.

Log-normalized gene expression values were extracted from the Seurat object and plotted against pseudotime ranks. LOESS smoothing was applied to visualize overall expression trends.

**Gene Expression Visualization**

```r
# Extract log-normalized RNA expression
expr <- LayerData(
  seurat_dorsal,
  assay = "RNA",
  layer = "data"
)

# Combine gene expression and pseudotime
df <- data.frame(
  dpt = seurat_dorsal[[]][colnames(expr), "dpt"],
  GLI3 = as.numeric(expr["GLI3", ]),
  EOMES = as.numeric(expr["EOMES", ]),
  NEUROD6 = as.numeric(expr["NEUROD6", ])
)

# Convert to long format
df_long <- tidyr::pivot_longer(
  df,
  cols = c(GLI3, EOMES, NEUROD6),
  names_to = "gene",
  values_to = "expression"
)

df_long$gene <- factor(
  df_long$gene,
  levels = c("GLI3", "EOMES", "NEUROD6")
)

# Plot gene expression dynamics
options(repr.plot.width = 15, repr.plot.height = 5)

ggplot(df_long, aes(x = dpt, y = expression)) +
  geom_point(color = "#0a77ca", size = 0.8, alpha = 0.4) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    color = "#d118ea",
    linewidth = 1.2
  ) +
  facet_wrap(~ gene, ncol = 3) +
  labs(
    x = "Pseudotime (DPT rank)",
    y = "Normalized expression"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#EEE7F2", color = NA),
    strip.text = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )
```
<img width="1800" height="600" alt="image" src="https://github.com/user-attachments/assets/569967fa-20bd-4476-b023-f7a097c80aaf" />

**Results and Interpretation**

Gene expression dynamics were examined to assess potential transcriptional changes associated with neuronal differentiation.

- **GLI3:** Associated with neural progenitor identity and developmental patterning.
- **EOMES:** Marks intermediate progenitor populations during cortical development.
- **NEUROD6:** Associated with neuronal differentiation and maturation.

The smoothed expression profiles can help evaluate whether the inferred pseudotime is consistent with the expected transition from neural progenitors through intermediate progenitors to differentiated neurons.

**Note:** Pseudotime direction depends on the selected root cell. Therefore, the inferred ordering should be validated against established developmental markers before assigning early and late differentiation states.


### Step 11. Save the Results
<a id="top"></a>

The final Seurat objects were saved in RDS format to preserve the processed expression data, metadata, dimensionality reductions, clustering results, and cell annotations for future analysis.

Two objects were saved:

- **seurat_obj_all.rds:** Complete scRNA-seq dataset after preprocessing, clustering, and cell type annotation.
- **seurat_obj_dorsal.rds:** Subset of dorsal telencephalic cells used for cell cycle regression and pseudotime analysis.

**Save Seurat Objects**

```r
# Create output directory
dir.create("DS1", showWarnings = FALSE)

# Save processed Seurat objects
saveRDS(seurat, "DS1/seurat_obj_all.rds")
saveRDS(seurat_dorsal, "DS1/seurat_obj_dorsal.rds")
```

**Load Saved Objects**

The saved objects can be reloaded without repeating the entire analysis.

```r
# Load complete dataset
seurat <- readRDS("DS1/seurat_obj_all.rds")

# Load dorsal telencephalic subset
seurat_dorsal <- readRDS("DS1/seurat_obj_dorsal.rds")
```

The saved objects allow downstream analyses and visualizations to be reproduced using the previously processed data.

## Part 2. Integration of Multiple scRNA-seq Datasets
<a id="top"></a>
When analyzing multiple scRNA-seq datasets, technical differences between experiments may introduce batch effects that obscure biological similarities between cell populations.

In this section, two human cerebral organoid datasets (DS1 and DS2) were combined to examine their initial distribution and evaluate the need for data integration.

### Step 1. Load and Merge Datasets

**1. Load Processed Seurat Objects**

The processed Seurat objects from DS1 and DS2 were loaded for joint analysis.

```r
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load processed datasets
seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")
seurat_DS2 <- readRDS("DS2/seurat_obj_all.rds")

# Examine datasets
seurat_DS1
seurat_DS2
```

**2. Merge the Datasets**

The two datasets were merged into a single Seurat object, while preserving their original sample identities.

```r
# Merge DS1 and DS2
seurat_merged <- merge(
  x = seurat_DS1,
  y = seurat_DS2,
  add.cell.ids = c("DS1", "DS2"),
  project = "DS1_DS2"
)

# Examine the merged object
seurat_merged

# Number of cells from each dataset
table(seurat_merged$orig.ident)
```

**3. Preprocessing of the Merged Dataset**

The merged dataset was normalized, followed by highly variable gene selection, scaling, PCA, and UMAP.

```r
# Join layers for the initial unintegrated analysis
seurat <- JoinLayers(seurat_merged)

# Preprocessing and dimensionality reduction
seurat <- seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)
```
<img width="1416" height="756" alt="image" src="https://github.com/user-attachments/assets/8b2904ff-711d-453a-b7bd-bf90b687ab94" />


**4. Visualize Dataset Distribution and Marker Expression**

UMAP visualization was used to examine the distribution of cells from the two datasets. Expression patterns of representative neural developmental markers were also investigated.

```r
options(repr.plot.width = 16, repr.plot.height = 9)

# Visualize cells by dataset
plot1 <- DimPlot(
  seurat,
  group.by = "orig.ident"
)

# Visualize neural developmental markers
plot2 <- FeaturePlot(
  seurat,
  features = c("FOXG1", "EMX1", "DLX2", "LHX9"),
  ncol = 2,
  pt.size = 0.1
)

# Combine plots
plot1 + plot2 +
  plot_layout(widths = c(1.5, 2))
```
<img width="1680" height="720" alt="image" src="https://github.com/user-attachments/assets/7b183d5a-aadb-4666-989f-8037ccf46d08" />

**Results and Interpretation**

The merged UMAP provides an initial overview of the transcriptional relationships between DS1 and DS2 before batch correction.

- **Dataset distribution:** Examining the overlap between DS1 and DS2 helps identify potential dataset-specific separation.
- **FOXG1:** Associated with telencephalic identity.
- **EMX1:** Associated with dorsal telencephalic development.
- **DLX2:** Associated with ventral telencephalic development.
- **LHX9:** Associated with diencephalic and midbrain development.

Differences in the distribution of cells from DS1 and DS2 may reflect technical batch effects, biological differences, or a combination of both.

This initial analysis serves as a reference for evaluating the effects of subsequent data integration methods.

**Note:** Merging datasets does not automatically correct batch effects. The integrated results should be evaluated to determine whether technical differences have been reduced while preserving meaningful biological variation. 
### Step 2-1. Data Integration Using Seurat (CCA)

Seurat provides an anchor-based integration framework to identify corresponding cell populations across multiple datasets and reduce technical batch effects.

The method uses **Canonical Correlation Analysis (CCA)** to identify shared transcriptional patterns between datasets. Integration anchors are then identified between corresponding cells and used to construct a batch-corrected expression representation.

The workflow consists of four main steps:

1. Normalize each dataset and identify highly variable genes.
2. Identify integration anchors using CCA.
3. Integrate the datasets using the identified anchors.
4. Perform dimensionality reduction and clustering on the integrated data.

**1. Prepare Individual Datasets**

Each dataset was independently normalized, and the top 3,000 highly variable genes were identified.

```r
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set RNA as the active assay
DefaultAssay(seurat_DS1) <- "RNA"
DefaultAssay(seurat_DS2) <- "RNA"

# Normalize and identify highly variable genes
seurat_DS1 <- seurat_DS1 %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000)

seurat_DS2 <- seurat_DS2 %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000)
```

**2. Identify Integration Anchors**

The datasets were organized into a list and processed using `FindIntegrationAnchors()`.

```r
# Create a list of Seurat objects
seurat_objs <- list(
  DS1 = seurat_DS1,
  DS2 = seurat_DS2
)

# Identify integration anchors using CCA
anchors <- FindIntegrationAnchors(
  object.list = seurat_objs,
  reduction = "cca",
  dims = 1:30
)
```

The first 30 canonical correlation dimensions were used to identify corresponding cell populations across DS1 and DS2.

**3. Integrate Datasets**

Integration was performed using the previously identified anchors.

```r
# Integrate datasets
seurat_integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

# Examine the integrated object
seurat_integrated
```
<img width="1174" height="312" alt="image" src="https://github.com/user-attachments/assets/cc8ae0c0-4431-4a78-ac24-1f90e0a2a5d4" />


The resulting Seurat object contains two important assays:

- **RNA:** Original, uncorrected expression data.
- **integrated:** Batch-corrected expression values generated during integration.

**4. Dimensionality Reduction and Clustering**

The integrated assay was used for PCA, UMAP, and graph-based clustering.

```r
# Use the integrated assay
DefaultAssay(seurat_integrated) <- "integrated"

# Scale integrated expression data
seurat_integrated <- ScaleData(seurat_integrated)

# PCA and UMAP
seurat_integrated <- RunPCA(
  seurat_integrated,
  npcs = 50
)

seurat_integrated <- RunUMAP(
  seurat_integrated,
  dims = 1:20
)

# Graph-based clustering
seurat_integrated <- seurat_integrated %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.6)
```
<img width="1494" height="912" alt="image" src="https://github.com/user-attachments/assets/3fcb8cc7-feac-4cd5-84cb-340a96ffcba4" />
**5. Visualize Integration Results**

UMAP plots were generated to examine the distribution of cells across datasets and the resulting clusters.

Expression patterns of representative neural developmental markers were also examined using the original RNA assay.

```r
options(repr.plot.width = 14, repr.plot.height = 9)

# Dataset distribution
plot1 <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  group.by = "orig.ident"
)

# Integrated cell clusters
plot2 <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
)

# Switch to RNA assay for marker expression
DefaultAssay(seurat_integrated) <- "RNA"

# Developmental marker expression
plot3 <- FeaturePlot(
  seurat_integrated,
  features = c("FOXG1", "EMX1", "DLX2", "LHX9"),
  reduction = "umap",
  ncol = 2,
  pt.size = 0.1
)

# Combine visualizations
((plot1 / plot2) | plot3) +
  plot_layout(widths = c(1, 2))
```
<img width="1680" height="1080" alt="image" src="https://github.com/user-attachments/assets/3657cd4f-43fe-4687-9a97-78749e3005e1" />

**6. Save Integrated Results**

```r
# Save the integrated Seurat object
saveRDS(
  seurat_integrated,
  file = "integrated_seurat_cca.rds"
)
```

**Results and Interpretation**

CCA-based integration was performed to reduce technical differences between DS1 and DS2 while preserving shared biological variation.

The resulting UMAP was examined to evaluate:

- **Dataset mixing:** Whether cells from DS1 and DS2 with similar transcriptional profiles occupy overlapping regions.
- **Cluster organization:** Whether distinct transcriptional populations remain identifiable after integration.
- **Marker expression:** Whether the integrated cell populations retain biologically meaningful neural developmental signatures.

The integrated UMAP can be compared with the unintegrated UMAP to assess changes in dataset separation and cellular organization.

**Note:** Improved dataset mixing alone does not demonstrate successful integration. Biological differences between cell populations should also be preserved, and potential overcorrection should be considered. 

### Step 2-2. Data Integration Using Harmony

[Harmony](https://github.com/immunogenomics/harmony) is a dimensionality reduction-based integration method designed to reduce batch effects while preserving biological variation across scRNA-seq datasets.

Unlike Seurat's anchor-based integration, Harmony operates directly on PCA embeddings. It iteratively applies soft clustering and batch-specific corrections to align corresponding cell populations across datasets.

The corrected embeddings are subsequently used for dimensionality reduction, visualization, and clustering.

**1. Prepare the Merged Dataset**

DS1 and DS2 were merged into a single Seurat object, followed by normalization, highly variable gene selection, scaling, and PCA.

```r
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

# Merge datasets
seurat_harmony <- merge(
  x = seurat_DS1,
  y = seurat_DS2,
  add.cell.ids = c("DS1", "DS2")
)

# Set RNA as the active assay
DefaultAssay(seurat_harmony) <- "RNA"

# Preprocessing
seurat_harmony <- NormalizeData(seurat_harmony)

seurat_harmony <- FindVariableFeatures(
  seurat_harmony,
  nfeatures = 3000
)

seurat_harmony <- ScaleData(seurat_harmony)

# PCA
seurat_harmony <- RunPCA(
  seurat_harmony,
  npcs = 50
)
```


**2. Perform Harmony Integration**

Harmony was applied to the PCA embeddings using `orig.ident` as the batch variable.

```r
# Perform Harmony integration
seurat_harmony <- RunHarmony(
  seurat_harmony,
  group.by.vars = "orig.ident",
  dims.use = 1:20,
  max_iter = 50
)
```

The resulting batch-corrected embeddings are stored in the `harmony` dimensional reduction.

```r
# Examine Harmony embeddings
head(Embeddings(seurat_harmony, "harmony"))
```

**3. UMAP and Clustering**

The Harmony-corrected embeddings were used for UMAP visualization and graph-based clustering.

```r
# UMAP using Harmony embeddings
seurat_harmony <- RunUMAP(
  seurat_harmony,
  reduction = "harmony",
  dims = 1:20
)

# Construct nearest-neighbor graph
seurat_harmony <- FindNeighbors(
  seurat_harmony,
  reduction = "harmony",
  dims = 1:20
)

# Identify cell clusters
seurat_harmony <- FindClusters(
  seurat_harmony,
  resolution = 0.6
)
```
<img width="1386" height="510" alt="image" src="https://github.com/user-attachments/assets/1b1505b6-34de-44c4-afd8-bcf46dc6f76a" />
<img width="1120" height="284" alt="image" src="https://github.com/user-attachments/assets/c0e43eac-be5a-4c5e-857d-941171a4ceaa" />

**4. Visualize Integration Results**

The integrated UMAP was examined to assess dataset mixing, cluster organization, and the expression of representative neural developmental markers.

```r
options(repr.plot.width = 14, repr.plot.height = 9)

# Dataset distribution
plot1 <- DimPlot(
  seurat_harmony,
  reduction = "umap",
  group.by = "orig.ident"
)

# Cell clusters
plot2 <- DimPlot(
  seurat_harmony,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
)

# Developmental marker expression
plot3 <- FeaturePlot(
  seurat_harmony,
  features = c("FOXG1", "EMX1", "DLX2", "LHX9"),
  reduction = "umap",
  ncol = 2,
  pt.size = 0.1
)

# Combine visualizations
((plot1 / plot2) | plot3) +
  plot_layout(widths = c(1, 2))
```
<img width="1680" height="1080" alt="image" src="https://github.com/user-attachments/assets/93448936-f57e-4a8a-b548-9dd2f8d24e60" />

**5. Save Integrated Results**

```r
# Save Harmony-integrated Seurat object
saveRDS(
  seurat_harmony,
  file = "integrated_harmony.rds"
)
```

**Results and Interpretation**

Harmony integration was performed to reduce dataset-associated variation in the PCA embeddings of DS1 and DS2.

The resulting UMAP was examined to evaluate:

- **Dataset mixing:** The extent to which transcriptionally similar cells from DS1 and DS2 occupy overlapping regions.
- **Cluster organization:** Whether distinct cell populations remain identifiable after integration.
- **Biological preservation:** Whether developmental marker expression remains consistent with expected neural cell identities.

Harmony-corrected embeddings provide an alternative representation of the merged datasets, enabling comparison with the unintegrated data and Seurat CCA integration results.

**Note:** Harmony corrects the low-dimensional embeddings rather than directly modifying the RNA expression matrix. Therefore, the original normalized RNA expression values remain available for downstream gene expression analysis. 

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

## Part 3. Reference-Based Cell Type Annotation

Cell type annotation is an essential step in scRNA-seq analysis. Although marker-based annotation can identify major cell populations, assigning accurate cell identities may be challenging when closely related cell types share similar transcriptional signatures.

Reference-based annotation provides an alternative approach by comparing the gene expression profiles of a query dataset with an existing, well-annotated reference.

In this section, the **Seurat CCA-integrated dataset (DS1 + DS2)** was used as the query, while a previously annotated scRNA-seq dataset was used as the reference.

The objective is to investigate how existing reference annotations can facilitate the identification of cell types in the query dataset.

### Step 0. Load and Explore Reference Data

**1. Load Query and Reference Datasets**

The Seurat-integrated dataset was loaded as the query, together with the annotated reference dataset provided in the original tutorial.

```r
library(Seurat)
library(ggplot2)
library(patchwork)

# Load query dataset (Seurat CCA integration)
seurat_query <- readRDS("integrated_seurat.rds")

# Load annotated reference dataset
seurat_ref <- readRDS("ref_seurat_obj.rds")

# Examine Seurat objects
seurat_query
seurat_ref
```

**2. Examine Reference Annotations**

The metadata of the reference dataset was examined to identify existing cell type annotations and developmental classifications.

```r
# Available metadata columns
colnames(seurat_ref[[]])

# Cell counts by annotated cell type
table(seurat_ref$celltype)

# Cell counts by developmental branch
table(seurat_ref$branch)
```
<img width="1618" height="334" alt="image" src="https://github.com/user-attachments/assets/c3c093a7-fcba-4cf6-817b-47853c07285c" />
<img width="752" height="862" alt="image" src="https://github.com/user-attachments/assets/3845edc3-15b1-4771-b4ec-b6ba71ae8b7c" />

The `celltype` and `branch` metadata provide existing biological annotations that can be used to characterize the reference populations.

**3. Visualize Reference Cell Populations**

UMAP visualizations were generated to examine the distribution of cells according to developmental branches and annotated cell types.

```r
# Reference UMAP colored by developmental branch
plot1 <- DimPlot(
  seurat_ref,
  reduction = "umap",
  group.by = "branch"
) +
  ggtitle("Reference: branch")

# Reference UMAP colored by cell type
plot2 <- DimPlot(
  seurat_ref,
  reduction = "umap",
  group.by = "celltype"
) +
  ggtitle("Reference: cell type")

# Display annotation plots
options(repr.plot.width = 14, repr.plot.height = 5)

print(plot1 | plot2)
```
<img width="1680" height="600" alt="image" src="https://github.com/user-attachments/assets/be0311c6-52b6-42d8-9515-d990a0697974" />

**4. Examine Developmental Marker Expression**

The expression of representative neural developmental markers was examined to assess the biological characteristics of the annotated reference populations.

```r
# Use original RNA expression
DefaultAssay(seurat_ref) <- "RNA"

# Visualize developmental markers
plot3 <- FeaturePlot(
  seurat_ref,
  features = c(
    "SOX2", "DCX", "FOXG1",
    "EMX1", "DLX2", "LHX9"
  ),
  reduction = "umap",
  ncol = 3,
  pt.size = 0.1,
  order = TRUE,
  cols = c("#E5E5E5", "#8B75C9")
)

# Display marker expression
options(repr.plot.width = 14, repr.plot.height = 8)

print(plot3)
```
<img width="1680" height="960" alt="image" src="https://github.com/user-attachments/assets/c3d92ffe-b7ce-4d27-a02e-b5d9d7757321" />

**Results and Interpretation**

The reference dataset was examined to characterize its annotated cellular populations and developmental organization.

UMAP visualizations were used to investigate:

- **Developmental branches:** The distribution of cells across annotated developmental lineages.
- **Cell type annotations:** The organization of previously identified cell populations.
- **Marker expression:** The correspondence between annotated populations and established neural developmental markers.

The expression of SOX2, DCX, FOXG1, EMX1, DLX2, and LHX9 provides additional biological context for evaluating reference cell identities.

These analyses establish the reference dataset for subsequent transcriptomic similarity analysis and reference-based cell type annotation. 

### Method 1-1. Transcriptome Similarity at the Cell Cluster Level

Reference-based cell type annotation can be performed by comparing the transcriptomic profiles of query cell clusters with those of annotated cell populations in a reference dataset.

In this analysis, **DS1 was used as the query dataset**, while the annotated reference dataset provided the known cell type identities.

Average gene expression profiles were calculated for each group, and **Spearman correlation** was used to evaluate transcriptomic similarity.

The analysis consisted of four steps:

1. Calculate average gene expression profiles for DS1 clusters and reference cell types.
2. Select informative genes shared between DS1 and the reference.
3. Calculate Spearman correlations between transcriptional profiles.
4. Visualize transcriptomic similarities using a hierarchically clustered heatmap.

#### 1. Calculate Average Gene Expression Profiles

Log-normalized RNA expression values were extracted from DS1 and the reference dataset. Average expression was calculated for each DS1 cluster and annotated reference cell type.

```r
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

# Load datasets
seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")
seurat_ref <- readRDS("ref_seurat_obj.rds")

# Use RNA assays
DefaultAssay(seurat_DS1) <- "RNA"
DefaultAssay(seurat_ref) <- "RNA"

# Join layers if using Seurat v5
if (inherits(seurat_DS1[["RNA"]], "Assay5")) {
  seurat_DS1 <- JoinLayers(seurat_DS1, assay = "RNA")
}

if (inherits(seurat_ref[["RNA"]], "Assay5")) {
  seurat_ref <- JoinLayers(seurat_ref, assay = "RNA")
}

# Extract normalized RNA expression
expr_ds1 <- GetAssayData(
  seurat_DS1,
  assay = "RNA",
  layer = "data"
)

expr_ref <- GetAssayData(
  seurat_ref,
  assay = "RNA",
  layer = "data"
)

# Extract DS1 cluster identities
ds1_labels <- Idents(seurat_DS1)[colnames(expr_ds1)]

# Extract reference cell type annotations
ref_labels <- seurat_ref[[]][
  colnames(expr_ref), "celltype"
]

# Identify groups
ds1_groups <- levels(droplevels(ds1_labels))
ref_groups <- sort(unique(as.character(ref_labels)))

# Average expression for DS1 clusters
avg_expr_ds1 <- sapply(ds1_groups, function(ct) {
  Matrix::rowMeans(
    expr_ds1[, ds1_labels == ct, drop = FALSE]
  )
})

# Average expression for reference cell types
avg_expr_ref <- sapply(ref_groups, function(ct) {
  Matrix::rowMeans(
    expr_ref[, ref_labels == ct, drop = FALSE]
  )
})

# Assign column names
colnames(avg_expr_ds1) <- ds1_groups
colnames(avg_expr_ref) <- ref_groups
```

#### 2. Select Informative Genes

Highly variable genes from the reference dataset that were also present in DS1 were selected for transcriptomic similarity analysis.

```r
# Select reference HVGs shared with DS1
genes2cor <- Reduce(
  intersect,
  list(
    VariableFeatures(seurat_ref),
    rownames(expr_ds1),
    rownames(expr_ref)
  )
)

# Number of selected genes
length(genes2cor)

# Check that sufficient genes are available
stopifnot(length(genes2cor) >= 50)
```

This approach focuses the comparison on genes that capture transcriptional heterogeneity within the reference dataset.

#### 3. Calculate Transcriptomic Similarity

Spearman correlation was used to measure the similarity between the average gene expression profiles of DS1 clusters and annotated reference cell types.

```r
# Calculate Spearman correlation
corr2ref_cl <- cor(
  avg_expr_ds1[genes2cor, , drop = FALSE],
  avg_expr_ref[genes2cor, , drop = FALSE],
  method = "spearman"
)

# Examine correlation matrix
round(corr2ref_cl, 3)
```



The resulting correlation matrix contains:

- **Rows:** DS1 clusters.
- **Columns:** Annotated reference cell types.
- **Values:** Spearman correlation coefficients.

Higher positive correlations indicate greater similarity between the corresponding transcriptional profiles.

#### 4. Visualize Transcriptomic Similarity

Hierarchical clustering was applied to the correlation matrix to organize DS1 clusters and reference cell types according to their similarity patterns.

```r
library(ggplot2)
library(patchwork)

# Hierarchical clustering
hc_row <- hclust(dist(corr2ref_cl))
hc_col <- hclust(dist(t(corr2ref_cl)))

# Convert dendrogram branches into ggplot segments
tree_segments <- function(hc) {

  n <- length(hc$order)

  positions <- numeric(n)
  positions[hc$order] <- seq_len(n)

  centers <- heights <- numeric(n - 1)
  segments <- list()

  for (i in seq_len(n - 1)) {

    children <- hc$merge[i, ]

    x <- sapply(children, function(j) {
      if (j < 0) positions[-j] else centers[j]
    })

    y <- sapply(children, function(j) {
      if (j < 0) 0 else heights[j]
    })

    centers[i] <- mean(x)
    heights[i] <- hc$height[i]

    segments[[i]] <- data.frame(
      x = c(x[1], x[2], x[1]),
      y = c(y[1], y[2], heights[i]),
      xend = c(x[1], x[2], x[2]),
      yend = rep(heights[i], 3)
    )
  }

  do.call(rbind, segments)
}

# Reorder correlation matrix
mat <- corr2ref_cl[
  hc_row$order,
  hc_col$order,
  drop = FALSE
]

nr <- nrow(mat)
nc <- ncol(mat)

# Prepare heatmap data
heat_data <- expand.grid(
  row = seq_len(nr),
  col = seq_len(nc)
)

heat_data$value <- as.vector(mat)

# Construct heatmap
p_heat <- ggplot(
  heat_data,
  aes(col, nr + 1 - row, fill = value)
) +
  geom_tile() +
  scale_x_continuous(
    breaks = seq_len(nc),
    labels = colnames(mat),
    limits = c(0.5, nc + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(nr),
    labels = rev(rownames(mat)),
    limits = c(0.5, nr + 0.5),
    expand = c(0, 0),
    position = "right"
  ) +
  scale_fill_gradientn(
    colours = c("#0390d7", "#FAF7F2", "#db1d36"),
    name = "Spearman"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(hjust = 0),
    legend.position = "bottom"
  )

# Column dendrogram
p_top <- ggplot(tree_segments(hc_col)) +
  geom_segment(
    aes(x, y, xend = xend, yend = yend),
    linewidth = 0.4
  ) +
  scale_x_continuous(
    limits = c(0.5, nc + 0.5),
    expand = c(0, 0)
  ) +
  theme_void()

# Row dendrogram
p_left <- ggplot(tree_segments(hc_row)) +
  geom_segment(
    aes(
      -y,
      nr + 1 - x,
      xend = -yend,
      yend = nr + 1 - xend
    ),
    linewidth = 0.4
  ) +
  scale_y_continuous(
    limits = c(0.5, nr + 0.5),
    expand = c(0, 0)
  ) +
  theme_void()

# Combine heatmap and dendrograms
options(
  repr.plot.width = 12,
  repr.plot.height = 12
)

wrap_plots(
  A = p_top,
  B = p_left,
  C = p_heat,
  design = "#A
            BC",
  widths = c(1, 9),
  heights = c(1, 8)
)
```
<img width="1440" height="1440" alt="image" src="https://github.com/user-attachments/assets/738f8da9-8ce1-4165-8340-7dab649b2020" />

#### 5. Identify Candidate Cell Type Annotations

For each DS1 cluster, the reference cell type with the highest Spearman correlation was identified as a potential annotation.

```r
# Identify closest reference cell types
best_matches <- data.frame(
  DS1_cluster = rownames(corr2ref_cl),

  closest_celltype = colnames(corr2ref_cl)[
    max.col(corr2ref_cl, ties.method = "first")
  ],

  spearman_rho = apply(
    corr2ref_cl,
    1,
    max,
    na.rm = TRUE
  )
)

# Display candidate annotations
print(best_matches)
```

**Results and Interpretation**

Transcriptomic similarity analysis was performed to investigate potential correspondences between DS1 cell populations and annotated reference cell types.

The correlation heatmap provides an overview of the transcriptional relationships between the two datasets.

- **Higher positive correlations:** Indicate greater transcriptomic similarity between DS1 clusters and reference cell types.
- **Lower correlations:** Suggest more distinct transcriptional profiles.
- **Hierarchical clustering:** Groups populations with similar correlation patterns.

The reference cell type showing the highest correlation with each DS1 cluster was considered a candidate annotation.

These results can support the interpretation of cellular identities previously inferred from canonical marker expression.

**Note:** Transcriptomic similarity does not establish definitive cell type identity. Candidate annotations should be validated using known marker genes and additional biological evidence.

### Method 1-2. Transcriptome Similarity at the Cell Level

While cluster-level transcriptomic similarity analysis compares average gene expression profiles, cell-level analysis allows reference-based annotation of individual cells.

In this approach, the expression profile of each cell in DS1 was compared with the average expression profiles of annotated reference cell types.

**Spearman correlation** was used to measure transcriptomic similarity. To improve computational efficiency, gene expression values were first converted to ranks, followed by sparse matrix correlation.

#### 1. Calculate Cell-Level Transcriptomic Similarity

The same informative genes selected in Method 1-1 were used for cell-level analysis.

```r
library(presto)
library(qlcMatrix)
library(Seurat)

# Rank average expression profiles of reference cell types
ranked_expr_ref <- apply(
  avg_expr_ref[genes2cor, , drop = FALSE],
  2,
  rank
)

# Rank gene expression within individual DS1 cells
ranked_expr_ds1 <- presto::rank_matrix(
  expr_ds1[genes2cor, , drop = FALSE]
)$X_ranked

# Calculate correlations
corr2ref_cell <- as.matrix(
  qlcMatrix::corSparse(
    ranked_expr_ds1,
    ranked_expr_ref
  )
)

# Assign cell and reference identifiers
rownames(corr2ref_cell) <- colnames(expr_ds1)
colnames(corr2ref_cell) <- colnames(avg_expr_ref)
```

The resulting correlation matrix contains:

- **Rows:** Individual DS1 cells.
- **Columns:** Annotated reference cell types.
- **Values:** Spearman correlation coefficients.

#### 2. Predict Cell Types Based on Maximum Correlation

Each DS1 cell was assigned the reference cell type with the highest correlation coefficient.

```r
# Identify the most similar reference cell type
ct_maxcor <- colnames(corr2ref_cell)[
  apply(corr2ref_cell, 1, which.max)
]

# Match predictions to cell barcodes
names(ct_maxcor) <- rownames(corr2ref_cell)

# Store predictions in Seurat metadata
seurat_DS1$celltype_maxcor <- ct_maxcor[
  colnames(seurat_DS1)
]

# Examine predicted cell types
table(seurat_DS1$celltype_maxcor)
```

These predictions provide candidate cell type annotations based on transcriptomic similarity to the reference dataset.

#### 3. Visualize Predicted Cell Types

The reference-based predictions were compared with the original cell annotations using UMAP.

```r
library(ggplot2)
library(patchwork)

options(repr.plot.width = 18, repr.plot.height = 10)

# Original cell annotations
plot1 <- DimPlot(
  seurat_DS1,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  pt.size = 0.2
) +
  ggtitle("Original annotation") +
  NoLegend()

# Reference-based predictions
plot2 <- DimPlot(
  seurat_DS1,
  reduction = "umap",
  group.by = "celltype_maxcor",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  pt.size = 0.2
) +
  ggtitle("Prediction from reference") +
  NoLegend()

# Compare annotations
(plot1 | plot2) &
  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )
```
<img width="2160" height="1200" alt="image" src="https://github.com/user-attachments/assets/4c3cc8ed-4feb-452f-ae29-b80253f964b9" />

**Results and Interpretation**

Cell-level transcriptomic similarity analysis was performed to predict the identities of individual cells in DS1 using an annotated reference dataset.

Each cell was assigned a candidate identity based on the highest Spearman correlation with the reference cell type expression profiles.

The predicted annotations were visualized alongside the original DS1 annotations to assess their consistency.

This approach allows cell-level comparisons without relying exclusively on cluster-level average expression profiles.

**Note:** Maximum-correlation assignment does not account for uncertainty in cell type predictions. Closely related cell types may exhibit similar transcriptional signatures, potentially leading to ambiguous assignments. Therefore, the resulting annotations should be interpreted alongside canonical marker expression and other supporting evidence.

#### 4. Summarize Cell-Level Similarity Across DS1 Cell Populations

To investigate transcriptomic similarities at the population level, cell-level correlation coefficients were standardized and summarized according to the existing DS1 cell annotations.

For each DS1 cell, correlation coefficients across reference cell types were standardized using Z-score transformation.

The standardized scores were then averaged within each DS1 cell population.

**1. Standardize Cell-Level Correlations**

```r
# Standardize correlations across reference cell types
# separately for each DS1 cell
corr2ref_scaled <- scale(t(corr2ref_cell))

# Extract original DS1 annotations
cell_groups <- Idents(seurat_DS1)[
  colnames(corr2ref_scaled)
]

# Identify cell populations
group_names <- levels(droplevels(cell_groups))
```

The standardized scores describe the relative similarity of each DS1 cell to different reference cell types.

**2. Calculate Average Similarity by Cell Population**

```r
# Average standardized similarity within each DS1 group
corr2ref_sum2cl <- t(
  sapply(group_names, function(cl) {
    rowMeans(
      corr2ref_scaled[
        , cell_groups == cl,
        drop = FALSE
      ]
    )
  })
)

# Assign row and column names
rownames(corr2ref_sum2cl) <- group_names
colnames(corr2ref_sum2cl) <- colnames(corr2ref_cell)

# Examine the resulting matrix
round(corr2ref_sum2cl, 3)
```

The resulting matrix contains:

- **Rows:** Original DS1 cell populations.
- **Columns:** Annotated reference cell types.
- **Values:** Average standardized transcriptomic similarity scores.

**3. Hierarchical Clustering**

Hierarchical clustering was applied to organize DS1 populations and reference cell types according to their similarity patterns.

```r
# Hierarchical clustering
hc_row <- hclust(dist(corr2ref_sum2cl))
hc_col <- hclust(dist(t(corr2ref_sum2cl)))

# Reorder similarity matrix
mat <- corr2ref_sum2cl[
  hc_row$order,
  hc_col$order,
  drop = FALSE
]

nr <- nrow(mat)
nc <- ncol(mat)

# Prepare heatmap data
heat_data <- expand.grid(
  row = seq_len(nr),
  col = seq_len(nc)
)

heat_data$value <- as.vector(mat)
```

**4. Visualize Standardized Transcriptomic Similarity**

A hierarchically clustered heatmap was generated to visualize the relative similarities between DS1 cell populations and reference cell types.

```r
library(ggplot2)
library(patchwork)

# Heatmap
p_heat <- ggplot(
  heat_data,
  aes(col, nr + 1 - row, fill = value)
) +
  geom_tile() +
  scale_x_continuous(
    breaks = seq_len(nc),
    labels = colnames(mat),
    limits = c(0.5, nc + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(nr),
    labels = rev(rownames(mat)),
    limits = c(0.5, nr + 0.5),
    expand = c(0, 0),
    position = "right"
  ) +
  scale_fill_gradient2(
    low = "#098aca",
    mid = "#FAF7F2",
    high = "#c91830",
    midpoint = 0,
    name = "Scaled\nsimilarity"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(hjust = 0),
    legend.position = "bottom"
  )

# Column dendrogram
p_top <- ggplot(tree_segments(hc_col)) +
  geom_segment(
    aes(x, y, xend = xend, yend = yend),
    linewidth = 0.4
  ) +
  scale_x_continuous(
    limits = c(0.5, nc + 0.5),
    expand = c(0, 0)
  ) +
  theme_void()

# Row dendrogram
p_left <- ggplot(tree_segments(hc_row)) +
  geom_segment(
    aes(
      -y,
      nr + 1 - x,
      xend = -yend,
      yend = nr + 1 - xend
    ),
    linewidth = 0.4
  ) +
  scale_y_continuous(
    limits = c(0.5, nr + 0.5),
    expand = c(0, 0)
  ) +
  theme_void()

# Combine heatmap and dendrograms
options(
  repr.plot.width = 12,
  repr.plot.height = 12
)

patchwork::wrap_plots(
  A = p_top,
  B = p_left,
  C = p_heat,
  design = "#A\nBC",
  widths = c(1, 9),
  heights = c(1, 8)
)
```
<img width="1440" height="1440" alt="image" src="https://github.com/user-attachments/assets/0931c2d6-44d7-41f8-a3e0-5a072568341d" />

**Results and Interpretation**

Standardized cell-level transcriptomic similarity scores were summarized across the original DS1 cell populations.

The resulting heatmap provides an overview of the relative transcriptional relationships between DS1 populations and annotated reference cell types.

- **Positive scores:** Indicate higher-than-average similarity to a reference cell type relative to other reference types.
- **Negative scores:** Indicate lower-than-average relative similarity.
- **Hierarchical clustering:** Groups populations with similar reference similarity patterns.

This analysis provides complementary evidence for evaluating the consistency between the original DS1 annotations and reference-based predictions.

**Note:** Standardized similarity scores represent relative associations rather than raw Spearman correlation coefficients. Therefore, they should not be interpreted as direct correlation values.

### Method 2. Seurat Label Transfer
<a id="top"></a>

## Part 4. Advanced scRNA-seq Analysis

Beyond standard preprocessing, dimensionality reduction, clustering, and cell type annotation, additional computational approaches can be applied to investigate biological relationships between cell populations.

Advanced scRNA-seq analyses include:

- Trajectory inference and pseudotime analysis.
- Cluster connectivity and branching analysis.
- RNA velocity and cell fate prediction.
- Cell-cell communication analysis.
- Gene regulatory network inference.

In this section, selected advanced methods were applied to the human cerebral organoid dataset (DS1) to explore potential relationships among neural cell populations.

---

### Part 4-1. Cluster Connectivity Analysis with PAGA

**Partition-based Graph Abstraction (PAGA)** is a graph-based approach for investigating connectivity between cellular populations.

While clustering identifies distinct cell populations, it does not directly describe their relationships. PAGA addresses this limitation by evaluating the connectivity between groups based on the single-cell neighborhood graph.

The resulting graph provides a simplified representation of the relationships among cell populations, which may help identify potential developmental connections.

PAGA was implemented using [Scanpy](https://scanpy.readthedocs.io/) in Python through the `reticulate` interface in R.

#### 1. Configure the Python Environment

The Python environment containing Scanpy and AnnData was configured using `reticulate`.

```r
library(reticulate)

# Select Python environment
reticulate::use_python(
  "/home/admin1/miniconda3/envs/anndata-py311/bin/python",
  required = TRUE
)

# Verify Python configuration
reticulate::py_config()

# Check required Python modules
reticulate::py_module_available("scanpy")
reticulate::py_module_available("anndata")
```

The Python environment must contain the required packages before running the analysis.

#### 2. Load the Processed DS1 Dataset

The previously processed DS1 Seurat object was loaded for PAGA analysis.

```r
library(Seurat)
library(Matrix)

# Load DS1
seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")

# Set RNA assay
DefaultAssay(seurat_DS1) <- "RNA"

# Join layers if using Seurat v5
if (inherits(seurat_DS1[["RNA"]], "Assay5")) {
  seurat_DS1 <- JoinLayers(
    seurat_DS1,
    assay = "RNA"
  )
}

# Extract normalized expression
expr <- GetAssayData(
  seurat_DS1,
  assay = "RNA",
  layer = "data"
)

# Extract raw counts
counts <- GetAssayData(
  seurat_DS1,
  assay = "RNA",
  layer = "counts"
)

# Extract cell and gene identifiers
cells <- colnames(expr)
genes <- rownames(expr)
```

#### 3. Convert Seurat to AnnData

The Seurat object was converted into an AnnData object to enable analysis with Scanpy.

The AnnData object contains:

- **X:** Log-normalized gene expression matrix.
- **obs:** Cell-level metadata and annotations.
- **var:** Gene identifiers.
- **layers:** Raw count matrix.
- **obsm:** PCA and UMAP embeddings.

```r
# Import Python AnnData module
ad <- reticulate::import(
  "anndata",
  convert = FALSE
)

# Create AnnData object
adata <- ad$AnnData(
  X = t(expr),

  obs = data.frame(
    celltype = as.character(
      Idents(seurat_DS1)[cells]
    ),
    row.names = cells
  ),

  var = data.frame(
    row.names = genes
  ),

  layers = list(
    counts = t(
      counts[genes, cells, drop = FALSE]
    )
  ),

  obsm = list(
    X_pca = Embeddings(
      seurat_DS1,
      "pca"
    )[cells, 1:20, drop = FALSE],

    X_umap = Embeddings(
      seurat_DS1,
      "umap"
    )[cells, , drop = FALSE]
  )
)

# Convert annotations to categorical data
adata$obs[["celltype"]] <- adata$obs[["celltype"]]$astype(
  "category"
)

# Save AnnData object
adata$write_h5ad("DS1/anndata_obj.h5ad")
```

#### 4. Construct the Neighborhood Graph

A k-nearest-neighbor graph was constructed using the first 20 principal components.

```r
# Import Scanpy
sc <- reticulate::import(
  "scanpy",
  convert = FALSE
)

# Load AnnData
adata_DS1 <- sc$read_h5ad(
  "DS1/anndata_obj.h5ad"
)

# Construct neighborhood graph
sc$pp$neighbors(
  adata_DS1,
  n_neighbors = 20L,
  use_rep = "X_pca"
)
```

The neighborhood graph represents local transcriptional relationships between individual cells in PCA space.

#### 5. Perform PAGA Analysis

PAGA was applied to evaluate connectivity between the annotated DS1 cell populations.

```r
# Perform PAGA
sc$tl$paga(
  adata_DS1,
  groups = "celltype"
)

# Save PAGA results
adata_DS1$write_h5ad(
  "DS1/anndata_paga.h5ad"
)
```

PAGA summarizes the cell-level neighborhood graph into a population-level connectivity graph.

#### 6. Visualize the PAGA Graph

The PAGA graph was visualized to examine connectivity between neural cell populations.

```r
# Configure non-interactive plotting backend
reticulate::py_run_string("
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
")

# Generate PAGA graph
sc$pl$paga(
  adata_DS1,
  color = "celltype",
  fontsize = 10L,
  frameon = FALSE,
  show = FALSE
)

# Save figure
reticulate::py_run_string("
fig = plt.gcf()
fig.set_size_inches(12, 9)
fig.savefig(
    'DS1/paga_connectivity.png',
    dpi=180,
    bbox_inches='tight'
)
")

# Display figure in Jupyter
IRdisplay::display_png(
  file = "DS1/paga_connectivity.png"
)
```
<img width="1708" height="1110" alt="image" src="https://github.com/user-attachments/assets/316fb7e0-022a-4322-a55d-c0f709e9bdc5" />

**Results and Interpretation**

PAGA was used to investigate transcriptional connectivity among annotated cell populations in DS1.

The resulting graph represents:

- **Nodes:** Annotated cell populations.
- **Node sizes:** Relative numbers of cells within each population.
- **Edges:** Connectivity between populations based on the single-cell neighborhood graph.
- **Edge weights:** Strength of the inferred connectivity between populations.

Stronger connections suggest that two populations share a greater degree of neighborhood connectivity in the analyzed transcriptional space.

These relationships may provide insights into potential developmental associations between neural progenitors, intermediate progenitors, and differentiated neurons.

**Note:** PAGA connectivity does not independently establish developmental direction or lineage ancestry. Additional evidence, such as pseudotime analysis or RNA velocity, is required to support directional interpretations.


### Part 4-2. Pseudotime Reconstruction in Complex Cell Populations

In Part 1, diffusion pseudotime (DPT) analysis was performed on a subset of dorsal telencephalic cells to investigate neuronal differentiation.

However, developmental systems often contain multiple cell populations, branching trajectories, and distinct terminal states. Restricting the analysis to a single unbranched trajectory may exclude potentially informative relationships between populations.

In this section, diffusion map and DPT analyses were applied to the **complete DS1 dataset** to explore transcriptional progression across multiple neural cell populations.

Different root cell selection strategies were also examined to evaluate their influence on pseudotime reconstruction.

#### 1. Construct a Diffusion Map Using All DS1 Cells

The first 20 principal components were used to construct a diffusion map representing transcriptional relationships across the complete dataset.

```r
library(Seurat)
library(destiny)
library(ggplot2)

# Load processed DS1
seurat_DS1 <- readRDS("DS1/seurat_obj_all.rds")

# Extract PCA embeddings
pca_embeddings <- Embeddings(
  seurat_DS1,
  "pca"
)[, 1:20, drop = FALSE]

# Cell identifiers in PCA order
pca_cells <- rownames(pca_embeddings)

# Construct diffusion map
set.seed(42)

dm <- DiffusionMap(pca_embeddings)
```

#### 2. Calculate Diffusion Pseudotime

An initial DPT analysis was performed without explicitly defining a biologically informed root cell.

```r
# Calculate DPT
dpt <- DPT(dm)

# Extract pseudotime
dpt_values <- dpt$dpt

# Convert to pseudotime ranks
dpt_rank <- setNames(
  rank(dpt_values, na.last = "keep"),
  pca_cells
)

# Store pseudotime in Seurat metadata
seurat_DS1$dpt <- dpt_rank[
  colnames(seurat_DS1)
]
```

**3. Visualize Pseudotime and Developmental Markers**

The reconstructed pseudotime was visualized together with representative markers of neural progenitors and differentiating neurons.

```r
# Use RNA assay
DefaultAssay(seurat_DS1) <- "RNA"

# Set figure dimensions
options(
  repr.plot.width = 16,
  repr.plot.height = 5
)

# Visualize pseudotime and marker expression
FeaturePlot(
  seurat_DS1,
  features = c(
    "dpt", "SOX2", "NHLH1", "DCX"
  ),
  reduction = "umap",
  ncol = 4,
  pt.size = 0.2,
  order = TRUE,
  cols = c("#E5E5E5", "#8872BD")
) &
  theme(
    aspect.ratio = 1,
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(size = 10)
  )
```
<img width="1920" height="600" alt="image" src="https://github.com/user-attachments/assets/2dec3f4f-0048-4b2f-9b3c-431eedd9fdcc" />


**Interpretation**

The expression of developmental markers provides biological context for interpreting the inferred pseudotime.

- **SOX2:** Neural progenitor marker.
- **NHLH1:** Associated with early neuronal differentiation.
- **DCX:** Marker of immature neurons.

Their expression patterns were examined to assess whether the inferred pseudotime was consistent with potential transitions from neural progenitors toward differentiating neurons.

#### 4. Reconstruct Pseudotime Using Neural Progenitor Root Cells

Because DPT depends on the selected starting point, an alternative analysis was performed using cells annotated as neural progenitors.

Four neural progenitor populations were considered:

```r
# Neural progenitor populations
npc_types <- c(
  "Dorsal telen. NPC",
  "G2M dorsal telen. NPC",
  "Dien. and midbrain NPC",
  "G2M Dien. and midbrain NPC"
)

# Match cell annotations to PCA order
npc_labels <- Idents(seurat_DS1)[pca_cells]

# Identify neural progenitor indices
idx_NPC <- which(npc_labels %in% npc_types)

# Check candidate cells
length(idx_NPC)
```

Three neural progenitor cells were randomly selected as candidate trajectory tips.

```r
# Ensure sufficient NPC cells are available
stopifnot(length(idx_NPC) >= 3)

# Select candidate tips
set.seed(12345)

idx <- sample(idx_NPC, 3)

# Reconstruct pseudotime
dpt2 <- destiny::DPT(
  dm,
  tips = idx
)

# Extract pseudotime
dpt2_values <- dpt2$dpt

# Convert to ranks and preserve cell identifiers
dpt2_rank <- setNames(
  rank(dpt2_values, na.last = "keep"),
  pca_cells
)

# Store in Seurat metadata
seurat_DS1$dpt2 <- dpt2_rank[
  colnames(seurat_DS1)
]
```

#### 5. Refine Root Cell Selection

An additional DPT reconstruction was performed using a neural progenitor cell selected from the beginning of the previously inferred pseudotime.

```r
# Identify NPC cells with valid pseudotime values
idx_valid <- idx_NPC[
  is.finite(
    seurat_DS1$dpt2[pca_cells[idx_NPC]]
  )
]

# Check candidate cells
if (length(idx_valid) == 0L) {
  stop("No neural progenitor cells with valid DPT values.")
}

# Select an early NPC according to dpt2
root_idx <- idx_valid[
  which.min(
    seurat_DS1$dpt2[pca_cells[idx_valid]]
  )
]

# Recalculate DPT
dpt3 <- destiny::DPT(
  dm,
  tips = root_idx
)

# Extract pseudotime
dpt3_values <- dpt3$dpt

# Convert to ranks
dpt3_rank <- setNames(
  rank(dpt3_values, na.last = "keep"),
  pca_cells
)

# Store in Seurat metadata
seurat_DS1$dpt3 <- dpt3_rank[
  colnames(seurat_DS1)
]
```

#### 6. Compare Pseudotime Reconstructions

The three pseudotime estimates were visualized on the original DS1 UMAP embedding.

```r
# Set figure dimensions
options(
  repr.plot.width = 18,
  repr.plot.height = 5
)

# Compare pseudotime estimates
FeaturePlot(
  seurat_DS1,
  features = c("dpt", "dpt2", "dpt3"),
  reduction = "umap",
  ncol = 3,
  cols = c("#E6E6E6", "#8B75C9"),
  order = TRUE
)
```
<img width="2160" height="600" alt="image" src="https://github.com/user-attachments/assets/3d39e1e2-3d6d-4413-96a2-dab95db7e380" />

**Results and Interpretation**

Diffusion pseudotime analysis was applied to the complete DS1 dataset to investigate potential transcriptional progression across multiple neural cell populations.

Three pseudotime reconstructions were examined:

- **DPT:** Initial reconstruction using the default computational root selection.
- **DPT2:** Reconstruction using selected neural progenitor cells as candidate tips.
- **DPT3:** Reconstruction using a neural progenitor root selected according to DPT2.

The reconstructed pseudotime distributions were compared with established developmental marker expression to evaluate their biological consistency.

**Note:** Pseudotime represents inferred transcriptional progression rather than chronological developmental time. In complex datasets containing multiple developmental lineages, a single pseudotime ordering may not adequately describe all differentiation trajectories.

The direction and biological interpretation of each inferred trajectory depend on root cell selection and should be evaluated using independent biological evidence.


### 4.3 RNA Velocity Analysis

### 4.4 Trajectory Analysis with CellRank


### 4.5 Cell–Cell Communication Analysis



