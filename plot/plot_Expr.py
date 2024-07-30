from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
def plt_3Dsvg(gene_name, coordinates, data_gene):
    fig = plt.figure()
    ax1 = plt.axes(projection='3d')
    #ax = fig.add_subplot(111,projection='3d')
    test_gene = data_gene.values[np.where(data_gene['gene'].values == gene_name)[0][0], 1:]
    mean_test_gene = np.mean(test_gene)
    std_test_gene = np.std(test_gene)
    boolean = []
    for item in test_gene:
        if item > mean_test_gene:
            boolean.append(item)
        else:
            boolean.append(0)
    ax1.scatter3D(coordinates[0], coordinates[1], coordinates[2], c = boolean, cmap = 'Greens')
    plt.show()
###############plot high wass dist gene#########
def selected_gene_fig(pro_svg_array, coordinates, data_gene1):
    for item in pro_svg_array[-10:, 0]:
        # plt_svg(item, coordinates, data_gene1).savefig('result/fig/experssion ' + str(item) + '.svg', 
        #                                                    format = 'svg')
        plt_3Dsvg(item, coordinates, data_gene1)

def plt_2Dsvg(gene_name, spatial_coordinates, data_gene): 
    #set_figsize(figsize=(3.5, 2.5))
    fig = plt.figure(figsize=(3.5, 2.5))
    test_gene = data_gene.values[np.where(data_gene['gene'].values == gene_name)[0][0], 2:]
    plt.scatter(spatial_coordinates[0], spatial_coordinates[1], c = test_gene, s = 0.1, cmap = 'PiYG')
    plt.title('The expression of '+ gene_name)
    plt.xticks(())
    plt.yticks(())
    plt.colorbar()
    plt.show()
    return fig

def selected_gene_fig(pro_svg_array, coordinates, data_gene1):
    for item in pro_svg_array[-20:, 0]:
        plt_2Dsvg(item, coordinates, data_gene1).savefig('result/fig/experssion ' + str(item) + '.svg', 
                                                        format = 'svg')