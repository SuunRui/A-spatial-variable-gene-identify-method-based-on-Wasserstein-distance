import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 读取数据
df = pd.read_csv("D:/A_study/A_study/ST_cluster/code_and_data/data/simulation/our_simulation_data/our_2DfullExpr.csv")

# 创建一个包含9个子图的大图
fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True, sharey=True)

# 循环遍历每个子图
for i, ax in enumerate(axes.flatten()):
    # 绘制散点图
    scatter = ax.scatter(df['x'].values, df['y'].values, c=df[f'svgs00{i + 1}'].values, cmap='PiYG', s=10)
    # 隐藏刻度标签
    ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
    # 添加子图的边框
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.5)
    # 添加子图标题
    if i+1 in [1, 2, 3, 8, 9]:
        ax.set_title('disconnect')
    else:
        ax.set_title('connect')
# # 在第一行上方添加 "Disconnect Set" 标签
# fig.text(0.46, 0.9, 'Disconnect Set', ha='center', va='center', fontsize=12, color='black')

# # 在第二行上方添加 "Connect Set" 标签
# fig.text(0.46, 0.63, 'Connect Set', ha='center', va='center', fontsize=12, color='black')

# 显示颜色条
fig.colorbar(scatter, ax=axes, orientation='vertical', fraction=0.03, pad=0.1)
fig.text(0.95, 0.5, 'logcounts', ha='center', va='bottom', rotation='vertical', fontsize=12, color='black')
# 保存图形
fig.savefig('D:/A_study/A_study/ST_cluster/code_and_data/result/our/fig/simulation.svg', format='svg')
