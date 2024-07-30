import data_Preprocessing
import caculate_wass_dist
import fitter_dist
import pandas as pd



# -----------------------------real datasets-----------------------------
file_way = 'data/2D_Rep11_MOB.csv'
adata = data_Preprocessing.read_data(file_way)
new_data = data_Preprocessing.data_pre(adata)
wass_dist = caculate_wass_dist.run(new_data, '2D_Rep11_MOB')
threshold_gene_site, SVG_numbers, SVGs = fitter_dist.calculate_threshold_gene_site(wass_dist)

file_way = 'data/2D_DLPFC.csv'
adata = data_Preprocessing.read_data(file_way)
new_data = data_Preprocessing.data_pre(adata)
wass_dist = caculate_wass_dist.run(new_data, '2D_DLPFC')
threshold_gene_site, SVG_numbers, SVGs = fitter_dist.calculate_threshold_gene_site(wass_dist)

file_way = 'data/3D_RA2.csv'
adata = data_Preprocessing.read_data(file_way)
new_data = data_Preprocessing.data_pre(adata)
wass_dist = caculate_wass_dist.run(new_data, '3D_RA2')
threshold_gene_site, SVG_numbers, SVGs = fitter_dist.calculate_threshold_gene_site(wass_dist)

file_way = 'data/3D_seqfish_27764670.csv'
adata = data_Preprocessing.read_data(file_way)
new_data = data_Preprocessing.data_pre(adata)
wass_dist = caculate_wass_dist.run(new_data, '3D_seqfish_27764670')
threshold_gene_site, SVG_numbers, SVGs = fitter_dist.calculate_threshold_gene_site(wass_dist)

# -----------------------------simulation datasets-----------------------------
file_way = 'data/3DhighExpr.csv'
new_data = pd.read_csv(file_way)
wass_dist = caculate_wass_dist.run(new_data.iloc, '3DhighExpr')
threshold_gene_site, SVG_numbers, SVGs = fitter_dist.calculate_threshold_gene_site(wass_dist)

file_way = 'data/3DlowExpr.csv'
new_data = pd.read_csv(file_way)
wass_dist = caculate_wass_dist.run(new_data, '3DlowExpr')
threshold_gene_site, SVG_numbers, SVGs = fitter_dist.calculate_threshold_gene_site(wass_dist)

file_way = 'data/2DhighExpr.csv'
new_data = pd.read_csv(file_way)
wass_dist = caculate_wass_dist.run(new_data, '2DhighExpr')
threshold_gene_site, SVG_numbers, SVGs = fitter_dist.calculate_threshold_gene_site(wass_dist)

file_way = 'data/2DlowExpr.csv'
new_data = pd.read_csv(file_way)
wass_dist = caculate_wass_dist.run(new_data, '2DlowExpr')
threshold_gene_site, SVG_numbers, SVGs = fitter_dist.calculate_threshold_gene_site(wass_dist)