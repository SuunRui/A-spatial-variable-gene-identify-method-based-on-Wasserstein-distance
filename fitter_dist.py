
from fitter import Fitter
from scipy import stats
import numpy as np

def calculate_threshold_gene_site(store_sorted):
    '''
    The function of threshold gene loci was calculated, and 95% of the 
    fitting function was obtained by fitting wass_dist, and the genes on 
    the right side of the loci were regarded as variable genes.
    dist_names: Name of the alternative distribution
    fitted_param: Name and parameters of the alternative distribution
    best_distri_name: The name of the optimal distribution
    option_distri: The command to find the segmentation point
    norm_a_right: Right side loci
    threshold_gene_site: The order in which the threshold genes are located
    '''
    svg_dist = store_sorted.iloc[:, 1]
    dist_names = ['beta',
            'expon',
            'gamma',
            'lognorm',
            'norm',
            'erlang']
    f = Fitter(svg_dist, distributions= dist_names, timeout =100)  # 创建Fitter类
    f.fit()  # Call the fit function to fit the distribution
    f.summary()  # Returns the ordered distribution fit quality (from good to bad fit) and plots the data distribution
    print(f.df_errors) # Returns the fitting mass of these distributions (sum of root-mean-square errors)
    fitted_param = f.fitted_param # Returns the parameters that fit the distribution
    f.fitted_pdf # The probability density generated using the distribution parameter that best fits the data distribution
    best_distri = f.get_best(method='sumsquare_error') # Returns the best-fit distribution and its parameters
    for item in best_distri.keys():
        best_distri_name = item
    print('best_distri: {}'.format(best_distri_name))
    # x = np.linspace(0, 1, 1000)
    if len(fitted_param[best_distri_name]) == 2:
        option_distri = 'stats.' + best_distri_name + '.ppf(0.95, loc=fitted_param[best_distri_name][0], scale=fitted_param[best_distri_name][1])'
    elif len(fitted_param[best_distri_name]) == 3:
        option_distri = 'stats.' + best_distri_name + '.ppf(0.95, fitted_param[best_distri_name][0], loc=fitted_param[best_distri_name][1], scale=fitted_param[best_distri_name][2])'
    else:
        option_distri = 'stats.' + best_distri_name + '.ppf(0.95, fitted_param[best_distri_name][0], fitted_param[best_distri_name][1], loc=fitted_param[best_distri_name][2], scale=fitted_param[best_distri_name][3])'
    norm_a_right = eval(option_distri)
    if np.argwhere(svg_dist >= norm_a_right).size != 0:
        threshold_gene_site = np.argwhere(svg_dist >= norm_a_right)[0][0]
        SVG_numbers = len(svg_dist) - threshold_gene_site
    else:
        threshold_gene_site = np.argwhere(svg_dist >= 0.6)[0][0]
        SVG_numbers = len(svg_dist) - threshold_gene_site
    print('gene_number: {}\nthreshold_gene_site: {}\nSVG_numbers: {}'.format(len(svg_dist), threshold_gene_site, SVG_numbers))
    SVGs = store_sorted.iloc[threshold_gene_site:, :]
    return threshold_gene_site, SVG_numbers, SVGs
# threshold_gene_site, SVG_numbers = calculate_threshold_gene_site(svg_dist)