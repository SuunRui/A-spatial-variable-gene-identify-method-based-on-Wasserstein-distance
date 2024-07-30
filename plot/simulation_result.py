import matplotlib.pyplot as plt
import numpy as np
import matplotlib
def simulation_result_plt(data, title):
    data = np.array(data)
    fig, ax = plt.subplots()
    if len(data) == 3:
        methods = ['SWDG', 'BSP', 'Spark-X']
        custom_colors = [(130, 178, 154),(244, 241, 222),(21, 151, 165)]
        custom_colors = [(r / 255, g / 255, b / 255) for r, g, b in custom_colors]
    else:
        methods = [ 'SWDG','BSP', 'nnSVG', 'spatialDE']
        custom_colors = [(130, 178, 154),(244, 241, 222),(223, 122, 94),(60, 64, 91)]
        custom_colors = [(r / 255, g / 255, b / 255) for r, g, b in custom_colors]

    metric1_scores = data[:, 0]  # TP
    metric2_scores = data[:, 1]  # FP
    
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.serif'] = 'Times New Roman'
    matplotlib.rcParams['font.size'] = 14 
    
    bar_width = 0.35
    index = np.arange(len(methods))
    plt.axhline(y=100, color='red', linestyle='--', label='Quantity is 100')
    
    plt.bar(index, metric1_scores, width=bar_width, label='True Positive Numbers')
    plt.bar(index + bar_width, metric2_scores, width=bar_width, label='False Positive Numbers')

    
    plt.xlabel('Methods')
    plt.ylabel('Numbers')
    plt.title('The resulting graph of the '+title+' dataset')

    
    plt.xticks(index + bar_width / 2, methods)
    # plt.show()
    # plt.legend()
    # legend = ax.legend(loc=(0.3, -0.25), frameon=False, ncol=3)
    plt.savefig('legend.svg', format='svg', bbox_inches='tight')
    plt.show()
    return fig