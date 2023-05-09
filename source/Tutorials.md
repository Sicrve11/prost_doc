# Tutorial 1: DLPFC Analysis
In this vignette, we analyzed tissue section from the human dorsolateral prefrontal cortex (DLPFC) 10x Visium ST dataset, which was manually annotated as the cortical layers and white matter (WM). The manual annotations were used as the ground truth to evaluate the accuracy of spatial domain segmentation.

---
## Identify SVGs
### 1.Load PROST and its dependent packages

    import pandas as pd 
    import numpy as np 
    import scanpy as sc 
    import os 
    import warnings 
    warnings.filterwarnings("ignore") 
    import matplotlib as mpl 
    import matplotlib.pyplot as plt 
    import PROST 
    PROST.__version__ 


    >>> ' 1.1.2 '

### 2.Set up the working environment and import data 

    # the location of R (used for the mclust clustering)
    ENVpath = "your path of PROST_ENV"            # refer to 'How to use    PROST' section
    os.environ['R_HOME'] = f'{ENVpath}/lib/R'
    os.environ['R_USER'] = f'{ENVpath}/lib/python3.7/site-packages/rpy2'
    
    # Set seed
    SEED = 818
    PROST.setup_seed(SEED)
    
    #%% Read in data
    section_num = 151672
    
    # Set directory (If you want to use additional data, please     change the file path)
    rootdir = 'datasets/DLPFC'
    
    input_dir = os.path.join(f'{rootdir}', str(section_num))
    spatial_dir = os.path.join(f'{rootdir}', str(section_num),  'spatial')
    output_dir = os.path.join(f'{rootdir}', str(section_num),   'results')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Read data from input_dir
    adata = sc.read_visium(path=input_dir, count_file='{}_filtered_feature_bc_matrix.h5'.format(section_num))
    adata.var_names_make_unique()

### 3.Calculate and save PI

    adata = PROST.prepare_for_PI(adata, platform="visium")
    adata = PROST.cal_prost_index(adata, platform="visium")
    adata.write_h5ad(output_dir+"/PI_result.h5")

    >>> Variable names are not unique. To make them unique, call `.var_names_make_unique`.
    >>> Variable names are not unique. To make them unique, call `.var_names_make_unique`.
    >>> Filtering genes ...
    >>> Calculating image index 1D:
    >>> 100%|██████████| 4015/4015 [00:00<00:00, 70423.08it/s]
    >>> Trying to set attribute `.var` of view, copying.
    >>> Normalization to each gene:
    >>> 100%|██████████| 5083/5083 [00:00<00:00, 13624.30it/s]
    >>> Gaussian filtering for each gene:
    >>> 100%|██████████| 5083/5083 [01:07<00:00, 74.99it/s]
    >>> Binary segmentation for each gene:
    >>> 100%|██████████| 5083/5083 [03:44<00:00, 22.60it/s]
    >>> Spliting subregions for each gene:
    >>> 100%|██████████| 5083/5083 [01:14<00:00, 68.52it/s]
    >>> Computing PROST Index for each gene:
    >>> 100%|██████████| 5083/5083 [00:03<00:00, 1478.57it/s]
    >>> PROST Index calculation completed !!
    

### 4.Draw SVGs detected by PI
    PROST.plot_gene(adata, platform="visium",size = 2, sorted_by = "PI", top_n = 25,save_path = output_dir)


    >>> ... storing 'feature_types' as categorical
    >>> ... storing 'genome' as categorical
    >>> Drawing pictures:
    >>> 100%|██████████| 1/1 [00:09<00:00,  9.55s/it]
    >>> Drawing completed !!
![DLPFC_pi_output](./_images/DLPFC/DLPFC_pi_output.png "Draw SVGs detected by PI")

### 5.Calculate Moran'I and Geary'C for SVGs dected by PI  
To assess the credibility of SVGs detected by these methods, we respectively used the spatial information of SVGs to calculate Moran’s I and Geary’s C statistics. 

    PROST.cal_moran_I_and_geary_C_for_PI_SVGs(adata, PI_top_n=50, save_path = output_dir)


    >>> 100%|██████████| 50/50 [00:28<00:00,  1.73it/s]
    >>> Average Moran'I of SVGs detected by PI = 0.4530851882366425 
    >>> Median Moran'I of SVGs detected by PI = 0.3949756315601886 
    >>> Average Geary'C of SVGs detected by PI = 0.5463279557091569 
    >>> Median Geary'C of SVGs detected by PI = 0.603879980983322

|    | geneID    | PI       | Moran_I  | Geary_C  |
|----|-----------|----------|----------|----------|
| 0  | MT-ND2    | 0.876836 | 0.710244 | 0.286052 |
| 1  | MT-ATP6   | 0.834577 | 0.718307 | 0.279474 |
| 2  | MT-CO2    | 0.816857 | 0.735137 | 0.265031 |
| 3  | MT-ND1    | 0.811519 | 0.719740 | 0.276575 |
| 4  | MT-ND4    | 0.778537 | 0.709382 | 0.286116 |
| 5  | MT-CO1    | 0.765785 | 0.718532 | 0.278538 |
| 6  | MT-CYB    | 0.742643 | 0.685484 | 0.314671 |
| 7  | MT-ND3    | 0.729128 | 0.694250 | 0.305809 |
| 8  | MT-CO3    | 0.720694 | 0.732460 | 0.268270 |
| 9  | MT-ND5    | 0.521871 | 0.504742 | 0.494602 |
| 10 | MT3       | 0.493878 | 0.503795 | 0.495233 |
| 11 | SCGB2A2   | 0.478912 | 0.680801 | 0.326584 |
| 12 | MTRNR2L8  | 0.473016 | 0.594114 | 0.403804 |
| 13 | FTH1      | 0.371137 | 0.464734 | 0.536722 |
| 14 | CST3      | 0.364869 | 0.484264 | 0.514394 |
| 15 | RPL41     | 0.357543 | 0.430153 | 0.565999 |
| 16 | NDUFA4    | 0.330432 | 0.423861 | 0.579269 |
| 17 | MBP       | 0.326337 | 0.494592 | 0.505158 |
| 18 | MGP       | 0.296306 | 0.468846 | 0.522271 |
| 19 | SCGB1D2   | 0.282661 | 0.495585 | 0.504803 |
| 20 | TUBA1B    | 0.274488 | 0.455331 | 0.543853 |
| 21 | MTRNR2L12 | 0.259394 | 0.409218 | 0.592105 |
| 22 | TMSB10    | 0.252886 | 0.492274 | 0.509731 |
| 23 | COX6A1    | 0.251861 | 0.379592 | 0.615655 |
| 24 | RPS21     | 0.243671 | 0.345257 | 0.655808 |
| 25 | CLU       | 0.243390 | 0.360654 | 0.640380 |
| 26 | CALM1     | 0.241974 | 0.341966 | 0.657869 |
| 27 | RPL34     | 0.241642 | 0.377486 | 0.626673 |
| 28 | RPL37A    | 0.235148 | 0.365629 | 0.630671 |
| 29 | GAPDH     | 0.234806 | 0.362507 | 0.637499 |
| 30 | SELENOW   | 0.232361 | 0.338193 | 0.663821 |
| 31 | COX6C     | 0.232221 | 0.457555 | 0.541740 |
| 32 | ATP5F1E   | 0.232024 | 0.341070 | 0.657950 |
| 33 | GNAS      | 0.223289 | 0.380733 | 0.617523 |
| 34 | COX4I1    | 0.222081 | 0.336472 | 0.662381 |
| 35 | SNAP25    | 0.214184 | 0.410525 | 0.590272 |
| 36 | MT-ND4L   | 0.214032 | 0.322941 | 0.677345 |
| 37 | CKB       | 0.213443 | 0.331925 | 0.666499 |
| 38 | FTL       | 0.212927 | 0.357595 | 0.643404 |
| 39 | NRGN      | 0.210670 | 0.344608 | 0.652806 |
| 40 | RPS28     | 0.210552 | 0.339370 | 0.658813 |
| 41 | PPIA      | 0.210144 | 0.306437 | 0.696177 |
| 42 | SAA1      | 0.208795 | 0.357856 | 0.638305 |
| 43 | CALM2     | 0.206207 | 0.322927 | 0.677944 |
| 44 | OLFM1     | 0.200591 | 0.290383 | 0.708072 |
| 45 | RPL32     | 0.192678 | 0.302734 | 0.696434 |
| 46 | TMSB4X    | 0.191757 | 0.317970 | 0.681999 |
| 47 | PCSK1N    | 0.191165 | 0.297974 | 0.704676 |
| 48 | RPS27A    | 0.189860 | 0.313465 | 0.684624 |
| 49 | COX7C     | 0.188659 | 0.324590 | 0.675992 |

--- 
## Clustering 
    # Set the number of clusters
    n_clusters = 5
    
### 1.Read PI result and Expression data preprocessing
    adata = sc.read(output_dir+"/PI_result.h5")

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata = PROST.feature_selection(adata, save_path = output_dir, by = "prost", n_top_genes = 3000)

### 2.Run PROST clustering
    PROST.run_prost_clust(adata,
                        platform="visium",
                        key_added = "PROST",
                        init="mclust",                         
                        n_clusters = n_clusters,                        
                        gnnlayers = 2,              
                        laplacin_filter = True,                        
                        lr = 0.1,                         
                        SEED=SEED,                          
                        max_epochs = 500,                        
                        tol = 5e-3,                        
                        post_processing = True,                        
                        pp_run_times = 3)

### 3.Save result
    adata.write_h5ad(output_dir+"/PNN_result.h5")   
    clustering = adata.obs["clustering"]
    clustering.to_csv(output_dir+"/clusters.csv",header = False)
    pp_clustering = adata.obs["pp_clustering"] 
    pp_clustering.to_csv(output_dir+"/pp_clusters.csv",header = False)
    embedding = adata.obsm["PROST"]
    np.savetxt(output_dir+"/embedding.txt",embedding)

    
    >>> Calculating adjacency matrix ...
    >>> Running PCA ...
    >>> Laplacian Smoothing ...
    >>> Initializing cluster centers with mclust, n_clusters known
    >>> Epoch: : 501it [09:07,  1.09s/it, loss=0.093866244]                       
    >>> Clustering completed !!
    >>> Post-processing for clustering result ...
    >>> Refining clusters, run times: 1/3
    >>> Refining clusters, run times: 2/3
    >>> Refining clusters, run times: 3/3

### 4.Plot clustering results 

    plt.rcParams["figure.figsize"] = (4,4)
    sc.pl.spatial(adata, 
                    img_key = "hires", 
                    color = ["annotation","clustering","pp_clustering"],
                    title = ["Manual annotation",'clustering','post-processed clustering'],                
                    na_in_legend = False,
                    ncols = 3,
                    size = 1)

    
    >>> storing 'annotation' as categorical
![clustering results](./_images/DLPFC/DLPFC_clusterresult_output.png "clustering results")

### 5.Calculate ARI and NMI 
To compare the domain segmentation performance quantitatively, we used the Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI) to measure the similarity between the predicted domains and the manual annotations across all twelve sections of the DLPFC dataset.

    ARI, NMI, silhouette_score = PROST.cal_metrics_for_DLPFC(adata.obs["pp_clustering"], labels_true_path = input_dir+'/cluster_labels.csv')

    
    >>> ARI = 0.5910397042708356 
    >>> AMI = 0.6813238415316797 
    >>> NMI = 0.6818348825641031 
    >>> v_measure_score = 0.6818348825641031 
    >>> silhouette_score = 0.3681630775671734 
    ==================================================

### 6.Plot UMAP and PAGA graph 
Next, the embeddings generated by PROST was applied to UMAP for visualization and PAGA for inferring trajectory.

    adata = sc.read_visium(path=input_dir, count_file='{}   _filtered_feature_bc_matrix.h5'.format(section_num))
    adata.var_names_make_unique()
    # Read annotation
    labels_true = pd.read_csv(input_dir+'/cluster_labels.csv')
    labels_true.index = labels_true["key"].str[7:]
    adata.obs["annotation"] = labels_true["ground_truth"]
    adata.obs["annotation"] = adata.obs["annotation"].astype('category').   astype('str')
    used_adata = adata[adata.obs["annotation"]!='nan']
    prost_embed = pd.read_csv(output_dir+"/embedding.txt",header = None,    delim_whitespace=True)
    prost_embed.index = labels_true.index
    adata.obsm["PROST"] = prost_embed
    # Plot
    plt.rcParams["figure.figsize"] = (6,5)
    sc.pp.neighbors(used_adata, use_rep="PROST")
    sc.tl.umap(used_adata)
    sc.tl.paga(used_adata,groups='annotation')
    sc.pl.paga_compare(used_adata, color="annotation",random_state=1,
                                 size = 35,legend_fontsize=25,  node_size_scale=4,
                                 frameon=False, show = False,fontoutline = 2)


    >>> Variable names are not unique. To make them unique, call `.var_names_make_unique`.
    >>> Variable names are not unique. To make them unique, call `.var_names_make_unique`.
    >>> ... storing 'annotation' as categorical
    >>> ... storing 'feature_types' as categorical
    >>> ... storing 'genome' as categorical
    >>> [<Axes:xlabel='UMAP1', ylabel='UMAP2'>, <Axes:>]

![DLPFC_annotation_output](./_images/DLPFC/DLPFC_annotation_output.png "UMAP and PAGA graph")
===

# Tutorial 2: Stereo-seq Analysis
In this vignette, we analysis an ST dataset with cellular resolution (~14 μm in diameter per spot) generated by the Stereo-seq platform from mouse olfactory bulb tissue (add citation) to evaluate the performance of PROST on ST datasets with single-cell resolution.

---
## Identify SVGs
### 1.Load PROST and its dependent packages

    import pandas as pd
    import numpy as np
    import scanpy as sc
    import os
    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import sys
    import PROST
    PROST.__version__


    >>> ' 1.1.2 '

### 2.Set up the working environment and import data 

    # the location of R (used for the mclust clustering)
    ENVpath = "your path of PROST_ENV"            # refer to 'How to use PROST' section 
    os.environ['R_HOME'] = f'{ENVpath}/lib/R'
    os.environ['R_USER'] = f'{ENVpath}/lib/python3.7/site-packages/rpy2'
    
    # init
    SEED = 818
    PROST.setup_seed(SEED)
    
    # Set directory (If you want to use additional data, please change the file path)
    rootdir = 'datasets/Stereo-seq/'
    
    input_dir = os.path.join(rootdir)
    output_dir = os.path.join(rootdir, 'results/')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Read data from input_dir
    adata=sc.read(input_dir+"/used_data.h5")


### 3.Calculate and save PI

    adata = PROST.prepare_for_PI(adata, percentage = 0.01, platform="stereo-seq")
    adata = PROST.cal_prost_index(adata, connect_kernel_size=6, neighbors=8,    platform="stereo-seq")
    adata.write_h5ad(output_dir+"/PI_result.h5")


    >>> Filtering genes ...
    >>> Trying to set attribute `.var` of view, copying.
    >>> Normalization to each gene:
    >>> 100%|██████████| 8520/8520 [00:02<00:00, 2871.91it/s]
    >>> Gaussian filtering for each gene:
    >>> 100%|██████████| 8520/8520 [18:07<00:00,  7.83it/s]
    >>> Binary segmentation for each gene:
    >>> 100%|██████████| 8520/8520 [00:29<00:00, 285.26it/s]
    >>> Spliting subregions for each gene:
    >>> 100%|██████████| 8520/8520 [01:39<00:00, 85.66it/s]
    >>> Computing PROST Index for each gene:
    >>> 100%|██████████| 8520/8520 [18:21<00:00,  7.74it/s]
    >>> PROST Index calculation completed !!

    

### 4.Draw SVGs detected by PI
    PROST.plot_gene(adata, platform="stereo-seq", size = 0.3, top_n = 25, ncols_each_sheet = 5, nrows_each_sheet = 5,save_path = output_dir)    


    >>> Drawing pictures:
    >>> 100%|██████████| 1/1 [00:15<00:00, 15.58s/it]
    >>> Drawing completed !!
![Stereo_seq_SVGs](./_images/Stereo-seq/Stereo_seq_SVGs.png "Stereo_seq_SVGs")

--- 
## Clustering 
    # Set the number of clusters
    n_clusters = 11


### 1.Read PI result and Expression data preprocessing
    PROST.setup_seed(SEED)
    adata = sc.read(output_dir+"/PI_result.h5")

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata = PROST.feature_selection(adata, save_path = output_dir, by = "prost",    n_top_genes = 3000)

    adata

    >>> View of AnnData object with n_obs × n_vars = 19109 × 2527
         obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes',    'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes'
         var: 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts',    'n_cells', 'SEP', 'SIG', 'PI', 'selected'
         uns: 'binary_image', 'del_index', 'gau_fea', 'locates', 'nor_counts', 'shape', 'subregions', 'log1p'
         obsm: 'spatial'

### 2.Run PROST clustering
    PROST.run_prost_clust(adata, 
                        platform="stereo-seq", 
                        min_distance = 50,
                        init="mclust",
                        n_clusters = n_clusters,                     
                        tol = 5e-3,
                        laplacin_filter = True,
                        SEED=SEED,
                        max_epochs = 500,
                        post_processing = False)
    

    >>> Calculating adjacency matrix ...
    >>> Running PCA ...
    >>> Laplacian Smoothing ...
    >>> Initializing cluster centers with mclust, n_clusters known
    >>> Epoch: : 501it [3:07:58, 22.51s/it, loss=0. 24661717]                         
    >>> Clustering completed !!

### 3.Save result
    adata.write_h5ad(output_dir + "/PNN_result.h5")
    clustering = adata.obs["clustering"]
    clustering.to_csv(output_dir + "/clusters.csv",header = False)

### 4.Plot clustering results 
    color_list = ['#1f77b4', '#ff7f0e', '#aec7e8', '#d62728', '#aa40fc', '#8c564b',
                  '#e377c2', '#b5bd61', '#17becf', '#279e68', '#ffbb78']
    adata.uns["clustering"+"_colors"] = color_list
    plt.rcParams["figure.figsize"] = (5,5)
    sc.pl.embedding(adata, basis="spatial", color="clustering",size = 7,s=6, show=False, title='clustering')
    plt.axis('off')
    plt.savefig(output_dir+"/clustering.png", dpi=600, bbox_inches='tight')


![Stereo_seq_Clustering](./_images/Stereo-seq/Stereo_seq_Clustering.png "Stereo_seq_Clustering")

### 5.Plot UMAP
Next, the embeddings generated by PROST was applied to UMAP for visualization.

    plt.rcParams["figure.figsize"] = (4,4)
    sc.pp.neighbors(adata, use_rep="PROST")
    sc.tl.umap(adata)
    ax = sc.pl.umap(adata, color="clustering", frameon=False, size=8,
                            show = False,legend_loc='on data',legend_fontoutline=2,legend_fontsize=11,
                            )
    plt.axis('off')
    plt.subplots_adjust()
    plt.savefig(output_dir+"/umap.png", dpi=600,bbox_inches='tight')

![Stereo_seq_umap](./_images/Stereo-seq/Stereo_seq_umap.png "Stereo_seq_umap")
===

# Tutorial 3: SeqFISH mouse embryo Analysis 
In this vignette, We applied PROST onto a SeqFISH-profiled dataset to evaluate its general applicability. 

---
## Identify SVGs
### 1.Load PROST and its dependent packages

    import pandas as pd 
    import numpy as np 
    import scanpy as sc 
    import os 
    import warnings 
    warnings.filterwarnings("ignore") 
    import matplotlib as mpl 
    import matplotlib.pyplot as plt 
    import PROST 
    PROST.__version__ 


    >>> ' 1.1.2 '

### 2.Set up the working environment and import data 

    # the location of R (used for the mclust clustering)
    ENVpath = "your path of PROST_ENV"            # refer to 'How to use PROST' section  
    os.environ['R_HOME'] = f'{ENVpath}/lib/R'
    os.environ['R_USER'] = f'{ENVpath}/lib/python3.7/site-packages/rpy2'
    
    # init
    SEED = 818
    PROST.setup_seed(SEED)
    
    # Set directory (If you want to use additional data, please change the file path)
    rootdir = 'datasets/SeqFISH/'
    
    input_dir = os.path.join(rootdir)
    output_dir = os.path.join(rootdir,'results/')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Read counts and metadata
    counts = pd.read_csv(input_dir + "counts.txt", sep = "\t")
    metadata = pd.read_csv(input_dir + "metadata.txt", sep = "\t")
    gene_name = counts.index

    # Create anndata for embryo1 (embryo2 or embryo3)
    '''Embryo1'''
    metadata_embryo1 = metadata[metadata["embryo"]=="embryo1"]
    counts_embryo1 = counts.loc[:,metadata_embryo1["uniqueID"]]
    spatial_embryo1 = metadata_embryo1[["x_global","y_global"]]
    spatial_embryo1.index = metadata_embryo1["uniqueID"]

    # Create anndata
    adata = sc.AnnData(counts_embryo1.T)
    adata.var_names_make_unique()
    # read spatial
    adata.obsm["spatial"] = spatial_embryo1.to_numpy()

    # read annotation
    annotation = metadata_embryo1["celltype_mapped_refined"]
    annotation.index = metadata_embryo1["uniqueID"]
    adata.obs["annotation"] = annotation
    adata.write_h5ad(output_dir+"/used_data1.h5")


    >>> ... storing 'annotation' as categorical
    >>> 'Embryo3'

### 3.Calculate and save PI

    adata=sc.read(output_dir+"/used_data1.h5")
    adata = PROST.prepare_for_PI(adata, percentage = 0.01, platform="SeqFISH")
    adata = PROST.cal_prost_index(adata, connect_kernel_size=8, neighbors=8,    platform="SeqFISH",del_rate=0.05)
    adata.write_h5ad(output_dir+"/PI_result.h5")


    >>> Filtering genes ...
    >>> Trying to set attribute `.var` of view, copying.
    >>> Normalization to each gene:
    >>> 100%|██████████| 351/351 [00:00<00:00, 5237.45it/s]
    >>> Gaussian filtering for each gene:
    >>> 100%|██████████| 351/351 [00:40<00:00,  8.67it/s]
    >>> Binary segmentation for each gene:
    >>> 100%|██████████| 351/351 [00:00<00:00, 18470.16it/s]
    >>> Spliting subregions for each gene:
    >>> 100%|██████████| 351/351 [00:00<00:00, 8355.38it/s]
    >>> Computing PROST Index for each gene:
    >>> 100%|██████████| 351/351 [00:39<00:00,  8.99it/s]
    >>> PROST Index calculation completed !!



### 4.Draw SVGs detected by PI
    PROST.plot_gene(adata, platform="SeqFISH", size = 0.3, top_n = 25, ncols_each_sheet = 5, nrows_each_sheet = 5,save_path = output_dir)


    >>> Drawing pictures:
    >>> 100%|██████████| 1/1 [00:15<00:00, 15.74s/it]
    >>> Drawing completed !!
![SeqFish_mouse_embryo_svgs](./_images/SeqFish/SeqFish_mouse_embryo_svgs.png "Draw SVGs detected by PI")

### 5.Calculate Moran'I and Geary'C for SVGs dected by PI  
To assess the credibility of SVGs detected by these methods, we respectively used the spatial information of SVGs to calculate Moran’s I and Geary’s C statistics. 

    PROST.cal_moran_I_and_geary_C_for_PI_SVGs(adata, PI_top_n=50, save_path = output_dir)


    >>> 100%|██████████| 50/50 [20:36<00:00, 24.73s/it]
    >>> Average Moran'I of SVGs detected by PI = 0.34560184671132244 
    >>> Median Moran'I of SVGs detected by PI = 0.36483518066319803 
    >>> Average Geary'C of SVGs detected by PI = 0.6171409870242501 
    >>> Median Geary'C of SVGs detected by PI = 0.5998499527724659

|    |  geneID |       PI |  Moran_I |  Geary_C |
|---:|--------:|---------:|---------:|---------:|
|  0 |   Hoxb9 | 1.000000 | 0.488931 | 0.465319 |
|  1 |    Cdx2 | 0.566876 | 0.411493 | 0.548865 |
|  2 |   Hoxc8 | 0.558326 | 0.302270 | 0.661814 |
|  3 |   Wnt5a | 0.554483 | 0.429061 | 0.535140 |
|  4 |   Bambi | 0.483233 | 0.424679 | 0.537970 |
|  5 |   Hoxa9 | 0.479378 | 0.348462 | 0.628304 |
|  6 |   Hoxb4 | 0.465979 | 0.291391 | 0.673090 |
|  7 |   Hoxb3 | 0.453696 | 0.275454 | 0.681236 |
|  8 | Tmem119 | 0.449840 | 0.420086 | 0.532817 |
|  9 |   Fgfr2 | 0.440052 | 0.296640 | 0.670719 |
| 10 |   Dusp6 | 0.436945 | 0.409998 | 0.546796 |
| 11 |  Tfap2b | 0.393255 | 0.370571 | 0.585511 |
| 12 |  Tfap2a | 0.388913 | 0.325753 | 0.623145 |
| 13 |   Smim1 | 0.384021 | 0.259597 | 0.701088 |
| 14 | Aldh1a2 | 0.379967 | 0.441316 | 0.526007 |
| 15 |  Hoxa11 | 0.350988 | 0.294576 | 0.655135 |
| 16 |    Tgm1 | 0.346816 | 0.055714 | 0.915164 |
| 17 |   Snai1 | 0.341068 | 0.361320 | 0.599254 |
| 18 |    Apln | 0.340328 | 0.126965 | 0.830317 |
| 19 |     Ttn | 0.338863 | 0.628871 | 0.343572 |
| 20 |   Cldn4 | 0.330571 | 0.427771 | 0.530761 |
| 21 |    Tbx4 | 0.321771 | 0.442577 | 0.498330 |
| 22 |   Sox10 | 0.318814 | 0.335579 | 0.630862 |
| 23 |   Hoxd9 | 0.318049 | 0.225928 | 0.741588 |
| 24 |   Hemgn | 0.309839 | 0.128971 | 0.838806 |
| 25 |    Tbx5 | 0.308528 | 0.441978 | 0.515955 |
| 26 |    Msx1 | 0.304036 | 0.459791 | 0.498687 |
| 27 |    Hcn4 | 0.297837 | 0.530591 | 0.472810 |
| 28 |   Suz12 | 0.296183 | 0.143624 | 0.823324 |
| 29 |    Evx1 | 0.290166 | 0.313486 | 0.646443 |
| 30 |    Wnt2 | 0.289947 | 0.494649 | 0.471269 |
| 31 |   Cntfr | 0.288735 | 0.378909 | 0.587733 |
| 32 |   Kmt2d | 0.283215 | 0.199363 | 0.765696 |
| 33 |     Afp | 0.273933 | 0.133562 | 0.840113 |
| 34 |   Hand1 | 0.271916 | 0.457340 | 0.501484 |
| 35 |   Hoxb1 | 0.271296 | 0.439619 | 0.527100 |
| 36 |   Sfrp5 | 0.271183 | 0.368350 | 0.600446 |
| 37 |   Podxl | 0.267092 | 0.420093 | 0.544644 |
| 38 |    Rgl1 | 0.265776 | 0.168780 | 0.792695 |
| 39 |  Popdc2 | 0.263662 | 0.552754 | 0.417445 |
| 40 |   Nanog | 0.260469 | 0.080936 | 0.879910 |
| 41 |  Dnmt3a | 0.259922 | 0.285143 | 0.681970 |
| 42 |  Pdgfra | 0.257661 | 0.416672 | 0.540798 |
| 43 |   Gata4 | 0.256379 | 0.544601 | 0.425685 |
| 44 | Cbfa2t3 | 0.255569 | 0.174157 | 0.785340 |
| 45 |  Tcf7l1 | 0.254504 | 0.186455 | 0.777835 |
| 46 |  Slc4a1 | 0.254302 | 0.321962 | 0.624245 |
| 47 |   Pdgfa | 0.252928 | 0.308816 | 0.644064 |
| 48 |   Gata5 | 0.252056 | 0.499571 | 0.458818 |
| 49 |   Gata6 | 0.251892 | 0.434917 | 0.530933 |

--- 
## Clustering 
    # Set the number of clusters
    n_clusters = 24
    

### 1.Read PI result and Expression data preprocessing
    PROST.setup_seed(SEED)
    # Read PI result
    adata = sc.read(output_dir+"/PI_result.h5")

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)


### 2.Run PROST clustering
    PROST.run_prost_clust(adata, 
                        platform="SeqFISH", 
                        min_distance = 3,
                        init="mclust",
                        n_clusters = n_clusters,                      
                        tol = 5e-3,
                        laplacin_filter = True,
                        lr = 0.1, 
                        SEED=SEED,
                        max_epochs = 500,
                        post_processing = False)

### 3.Save result
    adata.write_h5ad(output_dir + "/PNN_result.h5")
    clustering = adata.obs["clustering"]
    clustering.to_csv(output_dir + "/clusters.csv",header = False)
    embedding = adata.obsm["PROST"]
    np.savetxt(output_dir + "/embedding.txt",embedding)

    
    >>> Calculating adjacency matrix ...
    >>> Running PCA ...
    >>> Laplacian Smoothing ...
    >>> Initializing cluster centers with mclust, n_clusters known
    >>> Epoch: : 501it [3:17:42, 23.68s/it, loss=0.28359604]                         
    >>> Clustering completed !!

### 4.Plot annotation

    plt.rcParams["figure.figsize"] = (5,5)
    ax = sc.pl.embedding(adata, basis="spatial", color="annotation",size = 7,s=6, show=False, title='annotation')
    ax.invert_yaxis()
    plt.axis('off')
    plt.savefig(output_dir+"/annotation.png", dpi=600, bbox_inches='tight')

![annotation results](./_images/SeqFish/SeqFish_mouse_embryo_annotation.png "annotation results")

### 5.Plot clustering result
    plt.rcParams["figure.figsize"] = (5,5)
    ax = sc.pl.embedding(adata, basis="spatial", color="clustering",size = 7,s=6, show=False, title='clustering')
    ax.invert_yaxis()
    plt.axis('off')
    plt.savefig(output_dir+"/clustering.png", dpi=600, bbox_inches='tight')


### 6.Plot UMAP
    plt.rcParams["figure.figsize"] = (4,4)
    sc.pp.neighbors(adata, use_rep="PROST")
    sc.tl.umap(adata)
    ax = sc.pl.umap(adata, color="clustering", frameon=False, size=8,show = False)
    plt.axis('off')
    plt.subplots_adjust()
    plt.savefig(output_dir+"/umap.png", dpi=600,bbox_inches='tight')

![SeqFish_mouse_embryo_umapresult](./_images/SeqFish/SeqFish_mouse_embryo_umapresult.png "SeqFish_mouse_embryo_umap_result")


