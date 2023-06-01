# API

Import PROST as:

```python
import PROST
```

## PROST Index (PI)

PI is a quantifiable index that allows spatial feature gene selection, and the results can be combined with other methods. This is done by importing the data into AnnData and running the function below.

|  Functions                                     | Descriptions                                |
| ---------------------------------------------  | ------------------------------------------- |
|  `PROST.prepare_for_PI(adata[, platform...])`  | Process the raw AnnData for PI calculation. |
|  `PROST.cal_prost_index(adata[, platform...])` | Compute PI.                                 |



### `PROST.prepare_for_PI(adata, percentage = 0.1, platform="visium")`

Process raw AnnData, including gene filtering and data refinement.

**Parameters** **:**
- **adata（AnnData）**: Annotated data object.
- **percentage（float, [0-1]）**: A threshold value for gene expression positivity, and genes smaller than this threshold are removed during processing.
- **platform（str）**: Data platforms for spatial transcriptome, data will be processed accordingly according to the platform. Currently supported platforms are: `Visium`, `Stereo-seq`, `SeqFISH`.

**Returns** **:**
- **adata（AnnData）**: Annotated data object with result.

---

### `PROST.cal_prost_index(adata, connect_kernel_size=5, neighbors=8, del_rate=0.01, platform = "visium")`

The process of PI calculation, including pre-processing, subregion segmentation and calculating PROST Index.

**Parameters** **:**
- **adata（AnnData）**: Pre-processed Annotated data object.
- **connect_kernel_size（int）**: The size of the morphological operator, which controls the size of the rejected orphan points(cells).
- **neighbors（int, [4 or 8]）**: Mode of identify connectivity domain neighbor, usually `4 connected` and `8 connected`.
- **del_rate（float, [0-1]）**: Threshold for the connected domain proportion (`del_rate = connected domain size / total number of cells`), with `default 0.01`. This value removes small plaques of featureless genes.
- **platform（str）**: Data platforms for spatial transcriptome, data will be processed accordingly according to the platform. Currently supported platforms are: `Visium`, `Stereo-seq`, `SeqFISH`.

**Returns** **:**
- **adata（AnnData）**: Annotated data object with result(`adata.var["PI"]`).




## PROST Neural Network (PNN)

PNN is a neural network model that can perform unsupervised clustering of organizational space domains, and it can receive any gene expression information, but we suggest that the genes with the highest PI scores be entered.

| Functions                                     | Descriptions                                 |
| --------------------------------------------- | -------------------------------------------- |
| `PROST.feature_selection(adata[,by...])`      | Selection of genes of interest.              |
| `PROST.run_prost_clust(adata[, platform...])` | Run PNN for clustering.                      |
|


### `PROST.feature_selection(adata, by = 'manual', selected_gene_name=None, n_top_genes = 3000)`

The function to select genes, you can choose by top PI genes or your gene list, or choose randomly.

**Parameters** **:**
- **adata（AnnData）**: Annotated data object.
- **by（str）**: The way to select genes, and there are three types of received parameters: `"prost"`, `"scanpy"`, `"maunal" or any words`. 
- **selected_gene_name（list([str])）**: When `by = "manual" or any words`, this parameter is available. You can input the gene list which you interest.
- **n_top_genes（int）**: When `by = "prost" or "scanpy"`, this parameter is used to set the number you want use. The `"prost"` will get the top SVGs gene by PI, and the `"scanpy"` will get the top HVGs by scanpy.

**Returns** **:**
- **adata（AnnData）**: Annotated data object with result(`adata.var_names` have the selected gene name, and `adata` only have the selected gene data.).

---

### `PROST.run_prost_clust(adata, SEED, n_clusters=None, platform=None, k_neighbors = 7,  min_distance = 50, init="leiden", res = None, key_added = "PROST", tol=5e-3, gnnlayers = 2, num_pcs = 50, lr = 0.1, dropout=0.05,  leaky_alpha = 0.15, max_epochs = 500, laplacin_filter=True, post_processing = False, pp_run_times = 3, cuda=False)`

The process of PNN, including adjacency matrix, PCA embedding, laplacin filter, neural network.

**Parameters** **:**
- **adata（AnnData）**: Annotated data object.
- **SEED（int）**: Random seed values to ensure reproducible results.
- **n_clusters（int）**: Number of targets for clustering.
- **platform（str）**: Data platforms for spatial transcriptome, data will be processed accordingly according to the platform. Currently supported platforms are: `Visium`, `Stereo-seq`, `SeqFISH`.
- **k_neighbors（int）**: If `platform = "Visium"`, you need to set the number of neighbors to calculate the adjacency matrix.
- **min_distance（int）**: If `platform != "Visium"`, such as `"stereo-seq"`, you need to set the minimum distance of your neighbors to build the adjacency matrix.
- **init（str）**: Cluster initialization method, you have four choose that `"kmeans"`, `"mclust"`, `"louvain"`, and `"leiden"`.
- **res（float）**: If `init = "louvain" or "leiden"`, you must set this parameter as the clustering resolution.
- **key_added（str）**: The key will store the embeddig of input genes in input AnnData.
- **tol（float）**: Error tolerance for neural network optimization, `default is 1e-4`.
- **gnnlayers（int）**: Number of layers of graph neural network, `default is 2`.
- **num_pcs（int）**: Number of PCs used for PCA, `default is 50`.
- **lr（float, [0-1]）**: Learning rate of neural network, `default is 0.1`.
- **dropout（float, [0-1]）**: Dropout of neural network, `default is 0.05`.
- **leaky_alpha（float）**: Parameters of the LeakyReLU activation function for the attention layer, `default is 0.15`.
- **max_epochs（int）**: The maximum number of iterations of PNN, `default is 500`.
- **laplacin_filter（bool）**: Whether to use Laplacian filter, `default is True`.
- **post_processing（bool）**: Whether to perform PNN clustering post-processing, `default is True`.
- **pp_run_times（int）**: The maximum number of iterations of post_processing, `default is 3`.
- **cuda（bool）**: If you have gpu(N), you can input `True`, and `False` is default.

**Returns** **:**
- **adata（AnnData）**: Annotated data object with result(`adata.obs["clustering"]`).



## Other

Here are some other functions, including random seed setting, plotting, and calculation of indicators.

| Functions                                     | Descriptions                                 |
| --------------------------------------------- | -------------------------------------------- |
| `PROST.setup_seed(SEED)`                      | Set global random seed.                       |
| `PROST.plot_gene(adata[, platform...])`       | The function of plotting genes.               |
| `PROST.cal_moran_I_and_geary_C_for_PI_SVGs(adata[, save_path...])`| Calculation of moran_I and geary_C of genes.|
| `PROST.cal_metrics_for_DLPFC(labels_pred, labels_true_path[,...])` | Calculating DLPFC dataset metrics.|



### `PROST.setup_seed(SEED)`

Set random seeds in pytorch to make results repeatable.

**Parameters** **:**
- **SEED（int）**: Random seed value.

---

### `PROST.plot_gene(adata, platform, save_path, input_data = None, size=2 ,sorted_by = "PI", top_n = 50, ncols_each_sheet = 5, nrows_each_sheet = 5)`

Batch mapping of gene expression data.

**Parameters** **:**

- **adata（AnnData）**: Annotated data object.
- **platform（str）**: Data platforms for spatial transcriptome, data will be processed accordingly according to the platform. Currently supported platforms are: `Visium`, `Stereo-seq`, `SeqFISH`.
- **save_path（str）**: The path of results.
- **input_data（data）**: If you need to plot other data, you need to enter the data object to `input_data` (the data needs to have the input_data["genename"] field, which stores the name of the gene to be plotted)
- **size（int）**: Plot size, `default is 2`.
- **sorted_by（str）**: Field names for gene sorting, `defalut is PI`.
- **top_n（int）**: Number of genes plotted, default is top 50 of PI or HVGs.
- **ncols_each_sheet（int）**: Number of cols of the picture frame.
- **nrows_each_sheet（int）**: Number of rows of the picture frame.

---

### `PROST.cal_moran_I_and_geary_C_for_PI_SVGs(adata, PI_top_n, save_path=None)`

Calculate moran_I and geary_C, and save PI, moran_I, and geary_C as a table for output.

**Parameters** **:**

- **adata（AnnData）**: Annotated data object.
- **PI_top_n（int）**: Number of genes to be counted, `default is 50`.
- **save_path（str）**: The path of results.

---

### `PROST.cal_metrics_for_DLPFC(labels_pred, labels_true_path, print_result = True)`

Evaluating the clustering results of DLPFC.

**Parameters** **:**
- **labels_pred（list[int]）**: Annotated data object.
- **labels_true_path（str）**: Reference file path
- **print_result（bool）**: Whether to print the results or not, enter `True` if you need to print the results to the command line.

**Returns** **:**
- **ARI（int）**: Adjusted_rand_score.
- **NMI（int）**: Normalized_mutual_info_score.
- **silhouette_score（int）**: Silhouette_score.
