# Genomic Sequence Clustering

Unsupervised identification of cancer subtypes from high-dimensional gene expression data using Singular Value Decomposition (SVD) and Agglomerative Hierarchical Clustering.

## Objective

Cluster patients based on their gene expression profiles to discover distinct cancer subtypes without using any labeled data. The pipeline addresses the Curse of Dimensionality (5,000 genes, 300 patients) through SVD-based dimensionality reduction before applying hierarchical clustering.

## Dataset

**Gene Expression Cancer RNA-Seq** from the UCI Machine Learning Repository (ID: 401). This is a random extraction from the TCGA Pan-Cancer (HiSeq) project.

- **801 patients** (samples) across 5 cancer types
- **20,531 genes** (features) per patient, measured by Illumina HiSeq RNA-Seq
- No missing values
- Cancer types: BRCA (Breast), KIRC (Kidney), COAD (Colon), LUAD (Lung), PRAD (Prostate)

Source: https://archive.ics.uci.edu/dataset/401/gene+expression+cancer+rna+seq

## Pipeline

### 1. Data Preprocessing

- Load the gene expression matrix (patients x genes)
- Compute descriptive statistics and visualize expression distributions
- Standardize each gene to zero mean and unit variance:
  ```
  z = (x - mean) / std
  ```
  This prevents highly-expressed genes from dominating the decomposition.

### 2. Dimensionality Reduction with TruncatedSVD

The gene expression matrix X (n x p) is decomposed as:

```
X = U * S * V^T
```

Where:

- U: left singular vectors (patient loadings)
- S: diagonal matrix of singular values (component importance)
- V^T: right singular vectors (gene loadings)

TruncatedSVD retains only the top k=50 components, compressing 20,531 dimensions to 50 while preserving the majority of cancer-type-discriminating variance. A 410x compression ratio.

### 3. Distance Computation

Three pairwise distance metrics computed on the reduced space:

- **Euclidean**: straight-line distance, sensitive to magnitude
- **Cosine**: angular separation, invariant to expression magnitude
- **Correlation**: Pearson correlation-based, captures co-expression patterns

### 4. Agglomerative Hierarchical Clustering

Bottom-up algorithm with three linkage methods:

- **Ward**: minimizes within-cluster variance at each merge
  ```
  delta(A, B) = (|A| * |B|) / (|A| + |B|) * ||mu_A - mu_B||^2
  ```
- **Complete**: maximum pairwise distance between clusters
- **Average**: mean of all pairwise distances

### 5. Dendrogram Visualization

- Full dendrograms for all three linkage methods
- Truncated dendrogram highlighting major cluster structure
- Automatic cut-height suggestion based on largest gap in merge distances

### 6. Model Evaluation

**Internal Metrics** (no ground truth needed):

- Silhouette Score: cohesion vs separation, range [-1, 1]
- Calinski-Harabasz Index: between/within cluster variance ratio
- Davies-Bouldin Index: inter-cluster similarity, lower is better

**External Metrics** (vs ground truth):

- Adjusted Rand Index (ARI): cluster-to-label agreement, range [-1, 1]
- Normalized Mutual Information (NMI): information-theoretic agreement
- Confusion matrix with Hungarian alignment for per-cancer-type accuracy

## Results

| Metric           | Value                          |
| ---------------- | ------------------------------ |
| Dataset          | TCGA Pan-Cancer RNA-Seq (UCI)  |
| Patients         | 801                            |
| Genes            | 20,531                         |
| SVD Components   | 50                             |
| Compression      | 410x (20,531 to 50 dimensions) |
| Clusters         | 5                              |
| Linkage          | Ward                           |
| Silhouette Score | Computed at runtime            |
| ARI              | Computed at runtime            |
| NMI              | Computed at runtime            |

## Tech Stack

- Python, NumPy, pandas
- scikit-learn (TruncatedSVD, AgglomerativeClustering, evaluation metrics)
- SciPy (linkage, dendrogram, pdist, Hungarian algorithm)
- ucimlrepo (dataset fetching)
- Matplotlib, Seaborn

## How to Run

1. Clone the repository
2. Install dependencies:
   ```
   pip install numpy pandas scikit-learn scipy matplotlib seaborn ucimlrepo
   ```
3. Open `genomic_clustering.ipynb` in Jupyter Notebook or VS Code
4. Run all cells sequentially (the dataset is fetched automatically from UCI)

## Project Structure

```
genomic-sequence-clustering/
    genomic_clustering.ipynb    # Complete analysis pipeline
    README.md                   # Project documentation
```
