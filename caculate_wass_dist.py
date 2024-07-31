import numpy as np
from scipy import  interpolate
from scipy.stats import wasserstein_distance
import pandas as pd
import os


def sc_fit(x, expr):
    '''
    x: the coordinates of spot(2D or 3D)
    expr: expression of a gene
    This function fits the expression function of a certain gene in the spatial domain by linear interpolation.
    '''
    if len(x) == 2:
        func = interpolate.Rbf(x[0], x[1], expr, soomth = 0.01, function='linear')
    elif len(x) == 3:
        func = interpolate.Rbf(x[0], x[1], x[2], expr, soomth = 0.01, function='linear')
    else:
        func = 0
    return func

def sc_Wasserstein_Dist(x, expr, func):
    '''
    x: the coordinates of spot(2D or 3D)
    expr: expression of a gene
    func: expression function of a certain gene
    '''
    s = 50 # scale parameter
    mean_expr = np.mean(expr)
    if len(x) == 2:
        x1, y1 = np.mgrid[min(x[0]):max(x[0]):(max(x[0]) - min(x[0])) / s,
                min(x[1]):max(x[1]):(max(x[1]) - min(x[1])) / s]
        mean_func = np.ones(len(x1) * 2) * mean_expr
        random_exper1 = np.random.normal(loc=0, scale = mean_expr / 4.674506, size=len(x1) * 2) # Generating random noise
        random_exper2 = np.random.normal(loc=0, scale = mean_expr / 4.674506, size=len(x1) * 2)
        random_exper3 = np.random.normal(loc=0, scale = mean_expr / 4.674506, size=len(x1) * 2)
        mean_func1 = mean_func + random_exper1
        mean_func2 = mean_func + random_exper2
        mean_func3 = mean_func + random_exper3
        mean_func = (mean_func1+mean_func2+mean_func3)
        curve_value = func(x1, y1).flatten()
    elif len(x) == 3:
        x1, y1, z1 = np.mgrid[
            min(x[0]):max(x[0]):(max(x[0]) - min(x[0])) / s,
            min(x[1]):max(x[1]):(max(x[1]) - min(x[1])) / s, 
            min(x[2]):max(x[2]):(max(x[2]) - min(x[2])) / s] # The definition domain is generated based on the spatial position information in the data
        mean_func = np.ones(len(x1) * 3) * mean_expr
        random_exper1 = np.random.normal(loc=0, scale = mean_expr / 4.674506, size=len(x1) * 3)
        random_exper2 = np.random.normal(loc=0, scale = mean_expr / 4.674506, size=len(x1) * 3)
        random_exper3 = np.random.normal(loc=0, scale = mean_expr / 4.674506, size=len(x1) * 3)
        mean_func1 = mean_func + random_exper1
        mean_func2 = mean_func + random_exper2
        mean_func3 = mean_func + random_exper3
        mean_func = (mean_func1+mean_func2+mean_func3)
        curve_value = func(x1, y1, z1).flatten()
    else:
        print('spatial data error')
    curve_value[np.where(curve_value < 0 )[0]] = 0 # de-negative
    # wass_dist1 = wasserstein_distance(curve_value, mean_func1)
    # wass_dist2 = wasserstein_distance(curve_value, mean_func2)
    # wass_dist3 = wasserstein_distance(curve_value, mean_func3)
    return wasserstein_distance(curve_value, mean_func) # Take the mean of the distances as the final distance
def run(data, dataset_name):
    """
    data: Input data in dataframe format, with the first two or three columns of 
    spatial coordinate data and the following columns of gene expression data. 
    The first two columns of the column name are required to be xy, and the third 
    column name should be z for three-dimensional spatial transcriptome data.
    dataset_name: The data type is a str
    """
    if data.columns[2] == 'z' or data.columns[2] == 'Z':
        expr = data.iloc[:, 3:].values.T
        coordinates = data.iloc[:, :3].values.T
        gene_name = np.array([i for i in data.columns[3:]])
    else:
        expr = data.iloc[:, 2:].values.T
        coordinates = data.iloc[:, :2].values.T
        gene_name = np.array([i for i in data.columns[2:]])
    sorted = list()
    i = 1
    for item in expr:
        item_Wasserstein_Dist = sc_Wasserstein_Dist(coordinates, item, sc_fit(coordinates, item))
        sorted.append(item_Wasserstein_Dist)
        print('The calculation of {} genes was completed'.format(i))
        i += 1
    sorted_arr = np.array(sorted)##Wasserstein_Dist的array
    # sorted_arr = sorted_arr / max(sorted_arr)
    gene_wd = pd.DataFrame({'gene name': gene_name, 'wass dist': sorted_arr})
    gene_wd = gene_wd.sort_values(by='wass dist')
    if not os.path.exists('result'):
        os.makedirs('result')
    gene_wd.to_csv('result/'+dataset_name + '_wd.csv', index = False)
    return gene_wd
