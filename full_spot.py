import numpy as np
import scipy.spatial as spt

def spatial_coordinates_knn(spatial_coordinates, target_spot):
    '''
    Returns the space neighbor of a spot
    spatial_coordinates: the coordinates of all spot
    target_spot: the coordinates of target spot
    '''
    kt = spt.KDTree(data = spatial_coordinates, leafsize=10)
    d, x = kt.query(target_spot, 4)  # Returns the distance d of the nearest neighbor points and the order x in the array
    return d, x
def spot_neighbor(spatial_coordinates, target_spots):
    '''
    Calculate the nearest neighbor cells of all low-expressed cells and return the order of the nearest neighbor cells
    spatial_coordinates: the coordinates of all spot
    target_spots: the coordinates of target spot
    '''
    dist = []
    spot_order = np.empty([len(target_spots), 4], dtype='int')
    i = 0
    for item in target_spots:
        d, s_n = spatial_coordinates_knn(spatial_coordinates, spatial_coordinates[item])
        dist.append(d)
        spot_order[i] = s_n
        i += 1
    return spot_order
def Modifying_gene_low_expression(neighbor_low_spot_expression, adata):
    '''To modify the genes of the cells with low expression, the expression value of 
    the genes shared by the three nearest cells and the average value of their own gene 
    expression were selected to modify the gene expression level
    neighbor_low_spot_expression: The order of the nearest 3 neighbor cells
    adata: type is anndata
    '''
    for item in neighbor_low_spot_expression:
        transition = np.empty([len(item), len(adata.X[0])])
        transition = np.copy(adata.X[item, :])
        zero = np.where(transition[1:, :] == 0)[1] # Pick out the genes that don't share
        transition[1:, zero] = 0 # sets the expression value of genes that are not shared to 0
        adata.X[item[0], :] = np.mean(transition, axis= 0)
    return None
def Modifying_gene_high_expression(neighbor_hig_spot_expression, adata):
    '''To modify the genes of the cells with low expression, the expression value of the 
    genes shared by the three nearest cells and the average value of their own gene expression 
    were selected to modify the gene expression level
    neighbor_low_spot_expression: The order of the nearest 3 neighbor cells
    adata: type is anndata
    '''
    for item in neighbor_hig_spot_expression:
        transition = np.empty([len(item), len(adata.X[0])])
        transition = np.copy(adata.X[item, :])
        zero = np.where(transition == 0)[1]
        transition[:, zero] = 0
        adata.X[item[0], :] = np.mean(transition, axis= 0)
    return None