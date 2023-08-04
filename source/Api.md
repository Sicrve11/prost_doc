# API

Import PROST as:

```python
import PROST
```

## PROST Index (PI)

PI is a quantifiable index that allows spatial feature gene selection, and the results can be combined with other methods. This is done by importing the data into AnnData and running the function below.

|  Functions                          | Descriptions                                  |
| ----------------------------------- | -------------------------------------------   |
|  `PROST.pre_process`                | Pre-process gene count. |
|  `PROST.make_image`                 | Convert one-dimensional gene count into two-dimensional interpolated gene image. |
|  `PROST.gene_img_flatten`           | Convert two-dimensional interpolated gene image into one-dimensional gene count. |
|  `PROST.gau_filter_for_single_gene` | Gaussian filter for two-dimensional gene spatial expression images displayed in a regular pixels.  |
|  `PROST.cal_prost_index`            | Use PI to identify spatially variable genes for ST data. |
|


### `pre_process(adata, percentage = 0.1, var_stabilization = True)`

Pre-process gene count.

**Parameters** **:**
- **adata : Anndata**:   
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
- **percentage : float (default: 0.1)**:  
    For each gene, count the number of spots (cells) with expression greater than 0, if number is less than a `percentage` of the total spots (cells) number, remove this gene.
- **var_stabilization : bool (default: True)**:   
    Var-stabilize transformation.

**Returns** **:**
- **gene_use**:  
    Index of genes that `percentage` greater than threshold.
- **rawcount**:  
    Expression matrix of genes that `percentage` greater than threshold.

---

### `make_image(genecount, locates, platform = "visium", get_image_idx = False, grid_size = 20, interpolation_method='linear')`

Convert one-dimensional gene count into two-dimensional interpolated gene image.

**Parameters** **:**
- **genecount : pandas.DataFrame**:    
    The matrix of gene count expression. Rows correspond to genes and columns to cells. 
- **locates : matrix of shape (n_samples, 2)**:  
    The matrix of gene expression locates. Rows correspond to cells and columns to X-coordinate  and Y-coordinateof the position.
- **platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')**:   
    Sequencing platforms for generating ST data.
- **get_image_idx : bool (default: False)**:   
    If `get_image_idx=True`, calculate `image_idx_1d`. 
- **grid_size : int (default: 20)**:   
    The size of grid for interpolating irregular spatial gene expression to regular grids.
- **interpolation_method : str ['nearest','linear',cubic'] (default: linear)**:   
    The method for interpolating irregular spatial gene expression to regular grids. Same as `scipy.interpolate.griddata`

**Returns** **:**
- **image : ndarray**:  
    2-D gene spatial expression images displayed in a regular pixels.
- **image_idx_1d**:  
    If `get_image_idx=True`, which could be input to function `PROST.gene_img_flatten()`.

---

### `gene_img_flatten(I, image_idx_1d)`

Convert two-dimensional interpolated gene image into one-dimensional gene count.

**Parameters** **:**
- **I**:    
    The 2-D gene interpolated image. 
- **image_idx_1d**:  
    The 2-D index for 1-D gene count. Calculated by function `PROST. make_image()` with setting `get_image_idx = True`

**Returns** **:**
- **output**:  
    One-dimensional gene count.

---

### `gau_filter_for_single_gene(gene_data, locates, platform = "visium", image_idx_1d = None)`

Gaussian filter for two-dimensional gene spatial expression images displayed in a regular pixels.

**Parameters** **:**
- **gene_data : pandas.DataFrame**:    
    The matrix of gene count expression. Rows correspond to genes and columns to cells. 
- **locates : matrix of shape (n_samples, 2)**:  
    The matrix of gene expression locates. Rows correspond to cells and columns to X-coordinate  and Y-coordinateof the position. 
- **platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')**:     
    Sequencing platforms for generating ST data. 
- **image_idx_1d**:  
    The 2-D index for 1-D gene count. Calculated by function `PROST.make_image()` with setting `get_image_idx = True` 

**Returns** **:**
- **output**:  
    One-dimensional gene count.

---

### `cal_prost_index(adata, connect_kernel_size=5, del_rate=0.01, platform = "visium")`

Use PI to identify spatially variable genes for ST data.

**Parameters** **:**
- **adata : AnnData**:  
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
- **connect_kernel_size : int (default: 5)**:   
    Define the size of kernel for implementing morphological closing operation on the binary image. The larger size of kernel could eliminate a larger minutiae section inside the foreground.
- **del_rate : float (default: 0.01)**:     
    For a gene, if the largest foreground is less than `del_rate` of the entire spatial expression, it will be recognized as not having a significant spatial pattern.
- **del_rate（float, [0-1]）**:   
    Threshold for the connected domain proportion (`del_rate = connected domain size / total number of cells`), with `default 0.01`. This value removes small plaques of featureless genes.
- **platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')**:   
    Sequencing platforms for generating ST data.

**Returns** **:**
- **adata : Anndata**:  
    adata.obs.PI : PI score for each genes.  
    adata.obs.SEP : Seperability score for each genes.  
    adata.obs.SIG : Significance score for each genes. 

---



## PROST Neural Network (PNN)

PNN is a neural network model that can perform unsupervised clustering of organizational space domains, and it can receive any gene expression information, but we suggest that the genes with the highest PI scores be entered.

| Functions                    | Descriptions                                 |
| -----------------------------| -------------------------------------------- |
| `PROST.feature_selection`    | A feature selection tool for ST data. |
| `PROST.get_adj`              | Calculate adjacency matrix for ST data. |
| `PROST.preprocess_graph`     | Preprocess adj matrix. |
| `PROST.cluster_post_process` | Post_processing tool for cluster label that integrates neighborhood information. |
| `PROST.refine_clusters`      | Reassigning Cluster Labels Using Spatial Domain Information. |
| `PROST.run_prost_clust`      | Use PNN to identify spatial domains for ST data. |
|


### `feature_selection(adata, selected_gene_name = None, by = 'prost', n_top_genes = 3000)`

A feature selection tool for ST data.

**Parameters** **:**
- **adata : Anndata**:   
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
- **selected_gene_name : pd.Series (default: None)**:   
    Manually set the genes' name to seclect. Input as type [pandas.Series].
- **by : str ["prost", "scanpy"] (default: None)**:  
    Method for feature selection. 
    If `by=="prost"`, feature will be selected by PI;
    If `by=="scanpy"`, feature will be selected by Seurat.
- **n_top_genes : int (default: 3000)**:   
    Number of features (spatially variable genes) to select.

**Returns** **:**
- **genes**:   
    adata that include only selected genes.

---

### `get_adj(adata, mode = 'neighbour', k_neighbors = 7, min_distance = 150, self_loop = True)`

Calculate adjacency matrix for ST data.

**Parameters** **:**
- **mode : str ['neighbour','distance'] (default: 'neighbour')**:  
    The way to define neighbourhood. 
    If `mode='neighbour'`: Calculate adjacency matrix with specified number of nearest neighbors;
    If `mode='distance'`: Calculate adjacency matrix with neighbors within the specified distance.
- **k_neighbors : int (default: 7)**:  
    For `mode = 'neighbour'`, set the number of nearest neighbors if `mode='neighbour'`.
- **min_distance : int (default: 150)**:   
    For `mode = 'distance'`, set the distance of nearest neighbors if `mode='distance'`.
- **self_loop : bool (default: True)**:   
    Whether to add selfloop to the adjacency matrix.

**Returns** **:**
- **adj : matrix of shape (n_samples, n_samples)**:   
    Adjacency matrix where adj[i, j] is assigned the weight of edge that connects i to j.

---

### `preprocess_graph(adj, layer = 2, norm = 'sym', renorm = True, k = 2/3)`

Preprocess adj matrix.

**Parameters** **:**
- **adata : Anndata**:  
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
- **selected_gene_name : pd.Series (default: None)**:  
    Manually set the genes' name to seclect. Input as type [pandas.Series]. 
- **by : str ["prost", "scanpy"] (default: None)**:  
    Method for feature selection. 
    If `by=="prost"`, feature will be selected by PI;
    If `by=="scanpy"`, feature will be selected by Seurat.
- **n_top_genes : int (default: 3000)**:  
    Number of features (spatially variable genes) to select.

**Returns** **:**
- **genes**:  
    adata that include only selected genes.

---

### `cluster_post_process(adata, platform, k_neighbors = None, min_distance = None, key_added = "pp_clustering", p = 0.5, run_times = 3)`

Post_processing tool for cluster label that integrates neighborhood information.

**Parameters** **:**
- **adata : Anndata**:  
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
- **platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')**:  
    Sequencing platforms for generating ST data.
- **k_neighbors : int (default: None)**:  
    Same as `PROST.get_adj()`.
- **min_distance : int (default: None)**:   
    Same as `PROST.get_adj()`.
- **key_added : str (default: 'pp_clustering')**:  
    `adata.obs` key under which to add the cluster labels.
- **p : float (default: 0.5)**:  
    Rate of label changes in terms of neighbors.
- **run_times : int (default: 3)**:   
    Number of post-process runs. If the label does not change in two consecutive processes, the run is also terminated.

**Returns** **:**
- **adata.obs[key_added]**:  
    Array of dim (number of samples) that stores the post-processed cluster label for each cell.

---

### `refine_clusters(result, adj, p=0.5)`

Reassigning Cluster Labels Using Spatial Domain Information.

**Parameters** **:**
- **result**:   
    Clustering result to refine.
- **adj**:  
    Adjcency matrix.
- **k_neighbors or min_distance**:  
    Different way to calculate adj.
- **p : float (default: 0.5)**:  
    Rate of label changes in terms of neighbors.
- **run_times**:  
    Number of post-process runs. If the label does not change in two consecutive processes, program will also be terminated.

**Returns** **:**
- **pred_after**:   
    Check post_processed cluster label.

---

### `run_prost_clust(adata, SEED, platform=None, init="leiden", n_clusters=5, res=0.5, k_neighbors=7, min_distance=50, key_added="PROST", gnnlayers=2, lr=0.1, tol=5e-3, max_epochs=500, post_processing=False, pp_run_times=3, cuda=False)`

Use PNN to identify spatial domains for ST data.

**Parameters** **:**
- **adata : Anndata**:   
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
- **SEED : int**:  
    Random seed.
- **platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')**:  
    Sequencing platforms for generating ST data.
- **init : str ["kmeans","mclust","louvain","leiden"] (default: leiden)**:  
    Methods for initializing cluster centroids.
- **n_clusters : int (default: 5)**:  
    If the number of spatial domains is know, set cluster numbers for `init='kmeans'` or `init='mclust'`.
- **res : float (default: 0.5)**:  
    If the number of spatial domains is unknown, set resolutions parameter for `init='kmeans'` or `init='mclust'`.
- **k_neighbors : int (default: 7)**:  
    For `mode = 'neighbour'`, set the number of nearest neighbors if `mode='neighbour'`.
- **min_distance : int (default: 50)**:  
    For `mode = 'distance'`, set the distance of nearest neighbors if `mode='distance'`.  
- **key_added : str (default: 'PROST')**:  
    `adata.obsm` key under which to add the embedding representation generated by PROST.
- **gnnlayers : int (default: 2)**:  
    Number of stacked laplacian filter.
- **lr : float (default: 0.1)**:   
    Learning rate.
- **tol : float (default: 5e-3)**:   
    Stop criterion. The procedure stops when the change of clustering assignment between two consecutive iterations less than `tol`.
- **max_epochs : int (default: 500)**:   
    Number of epoch to train model.
- **post_processing : bool (default: False)**:   
    Whether to post-process oringal cluster result.
- **pp_run_times : int (default: 3)**:   
    For `post_processing=True`, set the number of post-processing run times.
- **cuda : bool (default: False)**:  
    Whether to use cuda acceleration.

**Returns** **:**
- **adata : Anndata**:  
    adata.obs['clustering'] : Original cluster label for each spot(cell).    
    adata.obs['pp_clustering'] : Post-processed cluster label for each spot (cell).    
    adata.obsm['PROST'] : Embedding representation generated by PROST.     

---

## Other

Here are some other functions, including random seed setting, plotting, and calculation of indicators.

| Functions                                     | Descriptions                                 |
| --------------------------------------------- | -------------------------------------------- |
| `PROST.setup_seed`                      | Set global random seed.   |
| `PROST.plot_gene`       | The function of plotting genes.               |
| `PROST.spatial_autocorrelation`| Statistical test of spatial autocorrelation for each gene. |
| `PROST.simulateH5Data` | Get simulated data. |
|


### `PROST.setup_seed(SEED)`

Set random seeds in pytorch to make results repeatable.

**Parameters** **:**
- **SEED（int）**: Random seed value.

---

### `PROST.plot_gene(adata, platform, save_path, input_data = None, size=2 ,sorted_by = "PI", top_n = 50, ncols_each_sheet = 5, nrows_each_sheet = 5)`

Batch mapping of gene expression data.

**Parameters** **:**

- **adata : AnnData**:  
    Annotated data object.
- **platform : str ['visium','Slide-seq','Stereo-seq','osmFISH','SeqFISH'] (default: 'visium')**:   
    Sequencing platforms for generating ST data.
- **save_path : str**:   
    The path of results.
- **input_data : data**:  
    If you need to plot other data, you need to enter the data object to `input_data` (the data needs to have the input_data["genename"] field, which stores the name of the gene to be plotted)
- **size : int**:    
    Plot size, `default is 2`.
- **sorted_by : str**:   
    Field names for gene sorting, `defalut is PI`.
- **top_n : int**:   
    Number of genes plotted, default is top 50 of PI or HVGs.
- **ncols_each_sheet : int**:   
    Number of cols of the picture frame.
- **nrows_each_sheet : int**:       
    Number of rows of the picture frame.

---

### `spatial_autocorrelation(adata, k = 10, permutations = None)`

Statistical test of spatial autocorrelation for each gene.

**Parameters** **:**
- **adata : Anndata**:    
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
- **k : int (default: 10)**:      
    Number of neighbors to define neighborhood.
- **permutations : int (default: None)**:    
    Number of random permutations for calculating pseudo p-values. 
    Default is 'none' to skip this step.

**Returns** **:**
- **adata : Anndata**:    
    adata.var["Moran_I"] : Moran's I   
    adata.var["Geary_C"] : Moran's C   
    adata.var["p_norm"] : p-value under normality assumption   
    adata.var["p_rand"] : p-value under randomization assumption   
    adata.var["fdr_norm"] : FDR under normality assumption   
    adata.var["fdr_rand"] : FDR under randomization assumption   
    
    if set `permutations`:   
    adata.var["p_sim"] : p-value based on permutation test   
    adata.var["fdr_sim"] : FDR based on permutation test   

---

### `simulateH5Data(adata, rr = 0.05)`

Get simulated data by droping part of the gene randomly.

**Parameters** **:**
- **adata : Anndata**:   
    H5 object. The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
- **rr : float (default: 0,05)**:   
    Dropout rate.

**Returns** **:**
- **adata : Anndata**:   
    Adata with specified dropout rates.

---