import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad


# --------------------------clustering-----------------------
def sc_cluster(adata, save_way1, save_way2):
    color_list = [(53, 183, 119), (145, 211, 192), (254, 168, 9), (188, 62, 3), (114, 170, 207), 
                  (255, 178, 208), (9, 147, 150), (235, 215, 165), (238, 155, 0), (204, 102, 2), 
                  (174, 32, 18), (155, 34, 39)]
    color_list = [(r / 255, g / 255, b / 255) for r, g, b in color_list]
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    # sc.pl.umap(adata)
    sc.tl.leiden(adata)
    sc.tl.louvain(adata = adata)
    adata.uns['leiden_colors'] = color_list[:len(adata.uns['leiden_colors'])]
    adata.uns['louvain_colors'] = color_list[:len(adata.uns['louvain_colors'])]
    sc.pl.umap(adata, color=["leiden"], save=save_way1)
    sc.pl.umap(adata, color=["louvain"], save=save_way2)

    return adata
from matplotlib.lines import Line2D
def cluster_fig(adata, coordinates, title_txt, save_path):
    colors = {
    '0': '#35b777',
    '1': '#91d3c0',
    '2': '#fea809',
    '3': '#bc3e03',
    '4': '#72aacf',
    '5': '#ffb2d0',
    '6': '#099396',
    '7': '#ebd7a5',
    '8': '#ee9b00',
    '9': '#cc6602',
    '10': '#ae2012',
    '11': '#9b2227'}
    fig = plt.figure()
    if len(coordinates) == 2:
        # plt.scatter(coordinates[0], coordinates[1], c = adata.obs['leiden'].map(colors), s = 10)
        plt.scatter(coordinates[0], coordinates[1], c = adata.obs['louvain'].map(colors), s = 80)
        plt.xticks([])
        plt.yticks([])
    if len(coordinates) == 3:
        ax = plt.axes(projection='3d')
        # ax.scatter3D(coordinates[0], coordinates[1], coordinates[2], coordinates[2], c = adata.obs['leiden'].map(colors))
        ax.scatter3D(coordinates[0], coordinates[1], coordinates[2], coordinates[2], c = adata.obs['louvain'].map(colors))
    plt.title(title_txt)
    # 隐藏 x 和 y 坐标轴的刻度
    # leiden_labels = [str(i) for i in range(len(adata.uns['leiden_colors']))]
    leiden_labels = [str(i) for i in range(len(adata.uns['louvain_colors']))]
    legend_labels = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[label], markersize=10, label=f'Cluster {label}') for label in leiden_labels]
    plt.legend(handles=legend_labels)
    plt.legend(handles=legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(save_path, format='svg', bbox_inches='tight')
    # 显示图形
    plt.show()
# ----------------------sownstream analysis-----------------------
from sklearn.metrics import adjusted_rand_score
def Rand_Index(y_true, y_pred):
    # caculate Rand_Index
    rand_index = adjusted_rand_score(y_true, y_pred)
    print("Adjusted Rand Index:", rand_index)
def to_adata(nnSVG_svg_Rep11_MOB, data_gene1, coordinates, adata, title_txt, save_way):
    """
    nnSVG_svg_Rep11_MOB: The P-value dataframe of the resulting gene
    data_gene1: raw data
    coordinates: 
    adata: nndata type of raw data
    """
    nnSVG_selected = data_gene1.T.iloc[1:, :]
    nnSVG_selected.columns = data_gene1['gene'].values
    # nnSVG_svg_Rep11_MOB = nnSVG_svg_Rep11_MOB.sort_values(by='p_values')
    # nnSVG_svg_name = nnSVG_svg_Rep11_MOB.iloc[:50, 0].values # for seqfish_3D dataset
    nnSVG_svg_name = [item for item in nnSVG_svg_Rep11_MOB.values[np.where(nnSVG_svg_Rep11_MOB.values[:, 1] < 0.05)[0], 0]]
    nnSVG_svg_name_1 = nnSVG_svg_name.copy()
    for item in nnSVG_svg_name_1:
        if item not in nnSVG_selected.columns:
            nnSVG_svg_name.remove(item)
    nnSVG_selected = nnSVG_selected[nnSVG_svg_name]
    if len(coordinates) == 2:
        obs = pd.DataFrame(
            {
            'spot_order' : np.array([i for i in range(len(coordinates.T))]), 
            'spatial_coordinates_x' : coordinates[0], 
            'spatial_coordinates_y' : coordinates[1]
            },
            index = np.array(nnSVG_selected.index)
                        )
    else:
        obs = pd.DataFrame(
            {
            'spot_order' : np.array([i for i in range(len(coordinates.T))]), 
            'spatial_coordinates_x' : coordinates[0], 
            'spatial_coordinates_y' : coordinates[1],
            'spatial_coordinates_z' : coordinates[2]
            },
            index = np.array(nnSVG_selected.index)
                        )
    var = pd.DataFrame(
        {
        },
        index = nnSVG_selected.columns.values
                    )
    nnSVG_adata = ad.AnnData(nnSVG_selected.values,obs=obs,var=var)
    sc.pp.neighbors(nnSVG_adata)
    sc.tl.umap(nnSVG_adata)
    sc.pl.umap(nnSVG_adata)
    sc.tl.leiden(nnSVG_adata)
    sc.pl.umap(nnSVG_adata, color=["leiden"])
    sc.tl.louvain(adata = nnSVG_adata)
    sc.pl.umap(nnSVG_adata, color=["louvain"])
    Rand_Index(adata.obs['leiden'].values, nnSVG_adata.obs['leiden'].values)
    Rand_Index(adata.obs['louvain'].values, nnSVG_adata.obs['louvain'].values)
    cluster_fig(nnSVG_adata, coordinates, title_txt, save_way)