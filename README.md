Name: Meiheng Liang
Programming Language: [Python]
Date: [Dec 2023]
Description:
This script is a demonstration the application of scRNA analysis using seurat pipeline referencing Seurat's [guided clustering tutorial]: (http://satijalab.org/seurat/pbmc3k_tutorial.html) ([Satija et al., 2015](https://doi.org/10.1038/nbt.3192)). 
The data consist of '3k PBMCs from a Healthy Donor' which available from https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz

Required files:

pbmc4k_raw_gene_bc_matrices.tar.gz

############################################################################
Required packages:
numpy
pandas
scanpy

############################################################################
#unpack data via linux-based system: uncomment and run the following to download and unpack the data. The last line creates a directory for writing processed data.

#Step wise Following script for preprocessing, normalization, PCA analysis, clustering, gene ranking, and annotation

import numpy as np
import pandas as pd
import scanpy as sc

import numpy as np
import pandas as pd
import scanpy as sc


adata = sc.read_10x_mtx(
    'path/to/your/file',  # read .mtx files from the directory
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
adata #check data structure
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white') 
#setting resolution for the graphs

Preprocessing
#Show those genes that yield the highest fraction of counts in each single cell, across all cells.

sc.pl.highest_expr_genes(adata, n_top=20, ) 
#checking top 20 genes that are highly expressed in all cells counts
#Basic filtering:

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata
#filter out mitochondria gene
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) # pp.calculate_qc_metrics, is a function for qc_metrics calculation
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True) 
#ncount=number of expressed genes, total counts= toal expression of each cell(based on UMI), and pct_count_mt reflect the percentage of expression from mitochondria (setting for mt filtration)
#Remove cells that have too many mitochondrial genes expressed or too many total counts:

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt') # checking the percent of expression from mitocondria versus total  (settig for filtration)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')# checking the number of genes versus total genes 

adata = adata[adata.obs.n_genes_by_counts < 3500, :] 
# setting range from 2500-4500 depending on the distribution of scatters 
adata = adata[adata.obs.pct_counts_mt < 5, :] #common setting for mt filtration
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


## normalization expressions levels for comparison across different cells 
sc.pp.normalize_total(adata, target_sum=1e4) 
sc.pp.log1p(adata)
##Logarithmize the data:

adata.raw = adata 
#Identify highly-variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5) 
sc.pl.highly_variable_genes(adata) #ploting highly variable genes
#listing top expressed genes
highly_variable_genes = adata.var[adata.var['highly_variable']].index
print(highly_variable_genes[:10]) 
#retains highly variable genes for further analysis
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10) #sclae gene to unit variance
adata # checking filteration result

# Principal component analysis
Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
sc.tl.pca(adata, svd_solver='arpack')
#sc.pl.pca(adata, color='CST3')

#inspect the contribution of single PCs to the total variance in the data. e.g. used in the clustering function sc.tl.louvain() or tSNE sc.tl.tsne(). 
sc.pl.pca_variance_ratio(adata, log=True) # checking the covariance of each PC for total 
Save the result.

adata.write(results_file) #store file 

#Computing the neighborhood graph
#neighborhood graph calculation and exmaine via biomarkers
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40) # parameters varied by datasets, but here is taking the default values.
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])


#Embedding the neighborhood graph
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')
#"raw" (normalized, logarithmized, but uncorrected) gene expression

sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False) #using normalize, logarithmized, and corrected datasets by setting raw as False

#Clustering the neighborhood graph
#Leiden graph-clustering method (community detection based on optimizing modularity) by Traag et al. (2018): directly clusters the neighborhood graph of cells, which was computed in the previous section.

#using Leiden graph-based clustering method
sc.tl.leiden(adata) #calculate 
sc.pl.umap(adata, color=['leiden']) #plot
Plot the clusters, which agree quite well with the result of Seurat.

#Save the result.
adata.write(results_file)

#Finding marker genes
#ranking highly differential genes in each cluster via t-test.
#ranking top_25 genes via t-test for efficiency
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
sc.settings.verbosity = 3  # reduce the verbosity

##Finding Marker genes
#checking for clusters and p-values
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(10)


Save the result.
adata.write(results_file)
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(10)

marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
adata = sc.read(results_file)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
#examine subclusters by leiden
for i in adata.obs['leiden'].cat.categories:
  number = len(adata.obs[adata.obs['leiden']==i])
  print('the number of category {} is {}'.format(i,number))

#remove genes with minimal cell counts
adata = adata[adata.obs[adata.obs['leiden'].astype(int)<11].index]
adata
# examine clustering result after the removal
for i in adata.obs['leiden'].cat.categories:
  number = len(adata.obs[adata.obs['leiden']==i])
  print('the number of category {} is {}'.format(i,number))
sc.pl.umap(adata, color=['leiden'], wspace=0.4, show=False)
adata #checking status
#Annotation:
marker_genes = ['CD3D', 'CD3E', 'CD3G','NKG7', 'KLRB1', 'MS4A1', 'CD79A', 'CD79B', 'CD68', 'CD1C']
marker_genes1 = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
With the exceptions of IL7R, which is only found by the t-test and FCER1A, which is only found by the other two appraoches, all marker genes are recovered in all approaches.

Louvain Group	Markers	Cell Type
0	IL7R	CD4 T cells
1	CD14, LYZ	CD14+ Monocytes
2	MS4A1	B cells
3	CD8A	CD8 T cells
4	GNLY, NKG7	NK cells
5	FCGR3A, MS4A7	FCGR3A+ Monocytes
6	FCER1A, CST3	Dendritic Cells
7	PPBP	Megakaryocytes
Let us also define a list of marker genes for later reference.

sc.pl.dotplot(adata, marker_genes, groupby='leiden');
sc.pl.dotplot(adata, marker_genes1, groupby='leiden');
cluster2annotation = {
    '0': 'CD4 T',
    '1': 'CD14 Monocytes',
    '2': 'CD8T',
    '3': 'CD4T',
    '4': 'B cells',
    '5': 'NK ',
    '6': 'NK cells',
    '7': 'B cells',
    '8': 'NK cells',
    '9': 'Dendritic cells',
    '10': 'FCGR3A+ Monocytes'
}
adata.obs['major_celltype'] = adata.obs['leiden'].map(cluster2annotation).astype('category')
#new_cluster_names = ['CD4 T', 'CD 14 Monocytes', 'CD8T', 'CD4T','B cell','CD8T', 'NK cell', 'B cell', 'NK cell', 'Dendritic', 'FCGR3A Monocytes']

#adata = sc.read(results_file)
sc.tl.dendrogram(adata,groupby='major_celltype')
sc.pl.dotplot(
    adata,
    marker_genes_dict,
    groupby='major_celltype',
    dendrogram=True,
    color_map="Blues",
    swap_axes=False,
    use_raw=True,
    standard_scale="var",
)
# Get a table with the scores and groups.

ax=sc.pl.embedding(
    adata,
    basis="X_umap",
    color='major_celltype',
    title='RNA-seq',
    frameon=False,
    ncols=3,
    #save='_figure1_celltype.png',
    return_fig=True,
    show=False,
)

sc.pl.umap(adata, color='major_celltype', legend_loc='on data', title='PMC4K_RNA-seq', frameon=False, save='.pdf')
If we want a more detailed view for a certain group, use sc.pl.rank_genes_groups_violin.

sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90);
Reload the object with the computed differential expression (i.e. DE via a comparison with the rest of the groups):

results_file = 'write/pbmc4k_Raw_2nd.h5ad'  # naming file for analysis results storage
adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading
adata.raw.to_adata().write('./write/pbmc4k_2nd_run.h5ad')
If you want to compare a certain gene across groups, use the following.




############################################################################
# Output files: 
resultfile.gz
