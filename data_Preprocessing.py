import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from full_spot import spatial_coordinates_knn
from full_spot import spot_neighbor
from full_spot import Modifying_gene_low_expression
from full_spot import Modifying_gene_high_expression



def read_data(file_way):
    """
    Here, we assume that the data is a .csv format file, the first two or three columns of the 
    data are spatial information, and the rest are gene expression information. The third 
    column name should be z for three-dimensional spatial transcriptome data.

    Returns data in anndata format created from the input data.
    """
    raw_data = pd.read_csv(file_way)
    if raw_data.columns[2] == 'z' or raw_data.columns[2] == 'Z':
        coordinates = raw_data.iloc[:,:3].values.T
        gene_expr = raw_data.iloc[:, 3:].values
        gene_name = raw_data.columns[3:]
        obs = pd.DataFrame(
        {
        'spatial_coordinates_x' : coordinates[0], 
        'spatial_coordinates_y' : coordinates[1],
        'spatial_coordinates_z' : coordinates[2]
        },
        index = ['sport' + str(i + 1) for i in range(len(coordinates[0]))]
                    )
    else: 
        coordinates = raw_data.iloc[:, :2].values.T
        gene_expr = raw_data.iloc[:, 2:].values
        gene_name = raw_data.columns[2:]
        obs = pd.DataFrame(
        {
        'spatial_coordinates_x' : coordinates[0], 
        'spatial_coordinates_y' : coordinates[1], 
        },
        index = ['sport' + str(i + 1) for i in range(len(coordinates[0]))]
                    )
    var = pd.DataFrame(
        {
        },
        index = gene_name
                    )
    adata = ad.AnnData(gene_expr, obs=obs, var=var)
    return adata

def data_pre(adata):
    '''
    # 'total_counts' = 'the total counts per cell'
    # 'pct_counts_mt' = 'the percentage of counts in mitochondrial genes'
    # 'n_genes_by_counts' = 'the number of genes expressed in the count matrix'
    '''
    if len(adata.obs.columns) == 3:
        coordinates = pd.concat([adata.obs['spatial_coordinates_x'], 
                                adata.obs['spatial_coordinates_y'], 
                                adata.obs['spatial_coordinates_z']], axis = 1)
        D = 3
    elif len(adata.obs.columns) == 2:
        coordinates = pd.concat([adata.obs['spatial_coordinates_x'], 
                                adata.obs['spatial_coordinates_y']], axis = 1)
        D = 2
    sc.pp.filter_genes(adata, min_cells = 5) # Remove low-quality genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True) # qc Quality control
    # sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    # sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    # neighbor_low_spot_expression = spot_neighbor(coordinates.T, 
    #                                             adata[adata.obs.total_counts <= 10000, :].obs['spot_order'].values)
    # Modifying_gene_low_expression(neighbor_low_spot_expression, adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
    # sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save = 'scatter_plot2.svg')
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    gene_expr = adata.X
    if D == 3:
        coordinates = pd.concat([adata.obs['spatial_coordinates_x'], 
                                adata.obs['spatial_coordinates_y'], 
                                adata.obs['spatial_coordinates_z']], axis = 1)
        coordinates.columns = ['x', 'y', 'z']
    elif D == 2:
        coordinates = pd.concat([adata.obs['spatial_coordinates_x'], 
                                adata.obs['spatial_coordinates_y']], axis = 1)
        coordinates.columns = ['x', 'y']
    gene_name = np.array([i for i in adata.var.index])
    new_data = pd.DataFrame(gene_expr, columns = gene_name, index = coordinates.index)
    new_data = pd.concat([coordinates, new_data], axis= 1)
    
    return new_data

# example
'''
file_way = 'dataset1.csv'
adata = read_data(file_way)
new_data = data_pre(adata)
'''