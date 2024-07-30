import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np

# 读取数据
df = pd.read_csv("D:/A_study/A_study/ST_cluster/code_and_data/data/simulation/our_simulation_data/our_3DfullExpr.csv")

# 创建图形
fig = plt.figure(figsize=(15, 12))

# 列名
# svgs_columns = df.columns[[1005, 1006, 1007, 1008,]]
svgs_columns = np.array(['svgs001', 'svgs002', 'svgs003', 'svgs004', 'svgs011', 'svgs012', 'svgs013', 'svgs014', 'svgs015'])

# 创建颜色条
cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.8])
# 在颜色条的侧边位置添加 'logcounts' 文字
fig.text(0.92, 0.5, 'logcounts', ha='center', va='bottom', rotation='vertical', fontsize=12, color='black')

# 循环绘制3D子图
for i, column in enumerate(svgs_columns, start=1):
    ax = fig.add_subplot(3, 3, i, projection='3d')
    scatter = ax.scatter(df['x'].values, df['y'].values, zs=df['z'].values, c=df[column].values, cmap='Greens', s=3)
    # 隐藏坐标轴刻度
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    if i in [1, 5, 6, 7]:
        ax.set_title('connect')
    else:
        ax.set_title('disconnect')

# 调整布局
plt.tight_layout()

# 显示颜色条
fig.colorbar(scatter, cax=cbar_ax)

# 保存图形为SVG格式
fig.savefig('D:/A_study/A_study/ST_cluster/code_and_data/result/our/fig/simulation_3D.svg', format='svg')

# 显示图形
plt.show()
