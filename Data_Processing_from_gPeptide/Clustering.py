
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Union
from matplotlib.axes import Axes


def plot_expression_matrix(data: pd.DataFrame, fig_name: str, value_name: str, format: str = 'png', dpi: int = 500,
                            method: str = 'average', metric: str = 'euclidean',
                            row_cluster: bool = True, col_cluster: bool = True,
                            row_colors: Union[None, pd.DataFrame, pd.Series, list] = None,
                            col_colors: Union[None, pd.DataFrame, pd.Series, list] = None,
                            yticklabels = 1, yticklabel_size: float = 7.5, xticklabel_size=10.0 , xrotated=False, linewidth: float = 1.62, linewidths: float = 0.0,
                            z_score: int = None, standard_scale: int = None,
                            mask: Union[pd.DataFrame, np.ndarray, None] = None,
                            cmap: str = 'coolwarm', annot: bool = False, width_scale: float = 1.0, height_scale: float = 1.0, colors_ratio=(0.024, 0.0114)) -> None:
    plt.rcParams['savefig.dpi'] = dpi
    plt.rcParams['figure.dpi'] = dpi
    cg = sns.clustermap(data, method=method, metric=metric, figsize=(16 * width_scale, 20 * height_scale),
                        xticklabels =1, yticklabels=yticklabels, linewidths=linewidths,
                        z_score=z_score, standard_scale=standard_scale,
                        row_cluster=row_cluster, col_cluster=col_cluster,
                        row_colors=row_colors, col_colors=col_colors,
                        colors_ratio=colors_ratio, cbar_pos=(.02, .8, .05, .18),
                        cmap=plt.cm.get_cmap(cmap), tree_kws={'linewidths': linewidth, 'colors': 'black'},
                        annot=annot, mask=mask,
                        cbar_kws={'orientation':'vertical', 'extend':'both'})
    fig = cg._figure

    ax: Axes = fig.axes[-2]
    xlabels = ax.get_xticklabels()
    ax.set_xticklabels(xlabels, fontdict={'font': 'Times New Roman', 'size': 18})
    if type(yticklabels) is int and yticklabels < 3:
        ylabels = ax.get_yticklabels()
        ax.set_yticklabels(ylabels, fontdict={'font': 'Times New Roman', 'size': yticklabel_size},  rotation='horizontal')
    if xrotated:
        xlabels = ax.get_xticklabels()
        ax.set_xticklabels(xlabels, fontdict={'font': 'Times New Roman', 'size': xticklabel_size},
                           rotation='vertical')
    else:
        xlabels = ax.get_xticklabels()
        ax.set_xticklabels(xlabels, fontdict={'font': 'Times New Roman', 'size': xticklabel_size})

    fig.axes[-1].set_title(value_name, fontdict={'font': 'Times New Roman', 'size': 24})
    plt.yticks(family='Times New Roman', size=20)
    row_order = cg.dendrogram_row.reordered_ind
    if col_cluster:
        col_order = cg.dendrogram_col.reordered_ind
        data = data.iloc[row_order].T.iloc[col_order].T
    else:
        col_order = data.columns.tolist()
        data = data.iloc[row_order].T.loc[col_order].T


    if not col_colors is None:
        ax = cg.ax_col_colors
        ax.tick_params(axis='y', which='both', labelsize=20)

    if not row_colors is None:
        ax = cg.ax_row_colors
        ax.tick_params(axis='x', which='both', labelsize=20)

    if row_colors is not None and type(row_colors) == pd.DataFrame:
        row_colors = row_colors.iloc[row_order]
        for category in reversed(row_colors.columns): data.insert(0, category, row_colors[category], )
    data.to_csv('Expression_Matrix_Clustered.csv', index=True)
    plt.savefig(f'{fig_name}.{format}', bbox_inches='tight')
    os.startfile(f'{fig_name}.{format}')
    
    
if __name__ == '__main__':
    
    data = pd.read_csv('Expression_Matrix_For_Cluster.csv', index_col=0)
    plot_expression_matrix(data, 'cluster_map', 'Scaled intensity', method='ward',
                           row_cluster=True, col_cluster=True, row_colors=None, col_colors=None,
                           linewidth=1.14, linewidths=0.0,
                           yticklabels=1, yticklabel_size=14.0, xticklabel_size=30.0, xrotated=True, cmap='bwr',
                           width_scale=1.14, height_scale=0.9,
                           colors_ratio=(0.03121, 0.0162))