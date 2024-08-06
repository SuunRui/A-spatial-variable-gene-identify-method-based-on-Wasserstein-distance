# 拟合分布
from fitter import Fitter
from scipy import stats
from statsmodels.stats.multitest import multipletests
import numpy as np

def calculate_threshold_gene_site(store_sorted):
    '''
    计算阈值基因位点的函数，通过对wass_dist进行拟合，得到拟合出来的函数的95%分位点，将分位点右侧的基因都看成可变基因
    dist_names：备选分布
    fitted_param：备选分布的名字和参数
    best_distri_name：最优分布的名字
    option_distri：求分位点的命令
    norm_a_right：右侧分位点
    threshold_gene_site：阈值基因所在的顺序
    '''
    dist_names = ['beta',
                'expon',
                'gamma',
                'lognorm',
                'norm',
                'geninvgauss',
                'erlang']
    wass_dist = store_sorted.iloc[:, 1]
    f = Fitter(wass_dist, distributions= dist_names, timeout =100)  # 创建Fitter类
    f.fit()  # 调用fit函数拟合分布
    f.summary()   #返回排序好的分布拟合质量（拟合效果从好到坏）,并绘制数据分布
    print(f.df_errors) #返回这些分布的拟合质量（均方根误差的和）
    fitted_param = f.fitted_param #返回拟合分布的参数
    f.fitted_pdf #使用最适合数据分布的分布参数生成的概率密度
    best_distri = f.get_best(method='sumsquare_error') #返回最佳拟合分布及其参数
    for item in best_distri.keys():
        best_distri_name = item
    fitted_param = f.fitted_param[best_distri_name]
    print('best_distri: {}'.format(best_distri_name))
    cdf_value = stats.__getattribute__(best_distri_name).cdf(wass_dist, *fitted_param)
    p_value = 1 - cdf_value
    SVG_numbers = len(np.where(p_value<0.05)[0])
    threshold_gene_site = len(wass_dist) - SVG_numbers
    SVGs = store_sorted.iloc[threshold_gene_site:, :]
    print(SVG_numbers, SVGs)
    return threshold_gene_site, SVG_numbers, SVGs
# threshold_gene_site, SVG_numbers = calculate_threshold_gene_site(svg_dist)
# calculate_threshold_gene_site(svg_dist)