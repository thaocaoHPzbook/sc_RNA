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


