# -*- coding: utf-8 -*-

from pkg_resources import resource_filename
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import q2templates
import qiime2

TEMPLATES = resource_filename("q2_metnet._clustermap", "assets")

clustermap_choices = {
    'metric': {'braycurtis', 'canberra', 'chebyshev', 'cityblock',
               'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
               'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski',
               'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener',
               'sokalsneath', 'sqeuclidean', 'yule'},
    'method': {'single', 'complete', 'average', 'weighted', 'centroid',
               'median', 'ward'},
    'cluster': {'samples', 'features', 'both', 'none'},
    'color_scheme': {'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG',
                     'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap',
                     'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r',
                     'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd',
                     'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r',
                     'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2',
                     'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn',
                     'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r',
                     'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy',
                     'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r',
                     'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r',
                     'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral',
                     'Spectral_r', 'Vega10', 'Vega10_r', 'Vega20', 'Vega20_r',
                     'Vega20b', 'Vega20b_r', 'Vega20c', 'Vega20c_r', 'Wistia',
                     'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r',
                     'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot',
                     'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r',
                     'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cool',
                     'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r',
                     'cubehelix', 'cubehelix_r', 'flag', 'flag_r',
                     'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r',
                     'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r',
                     'gist_rainbow', 'gist_rainbow_r', 'gist_stern',
                     'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot',
                     'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r',
                     'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r',
                     'inferno', 'inferno_r', 'jet', 'jet_r', 'magma',
                     'magma_r', 'mako', 'mako_r', 'nipy_spectral',
                     'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r',
                     'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow',
                     'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r',
                     'spectral', 'spectral_r', 'spring', 'spring_r', 'summer',
                     'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r',
                     'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain',
                     'terrain_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r',
                     'winter', 'winter_r', 'cividis', 'cividis_r'}
}


_clustering_map = {'both': {'col_cluster': True, 'row_cluster': True},
                   'samples': {'col_cluster': False, 'row_cluster': True},
                   'features': {'col_cluster': True, 'row_cluster': False},
                   'none': {'col_cluster': False, 'row_cluster': False}}

def plotClusteMap(output_dir: str, table: pd.DataFrame,
                sample_metadata: qiime2.CategoricalMetadataColumn = None,
                feature_metadata: qiime2.CategoricalMetadataColumn = None, 
                title: str = None, metric: str = 'euclidean', 
                method: str = 'average', cluster: str = 'both', 
                color_scheme: str = 'rocket', xlabels: bool = True, ylabels: bool = True) -> None:
    
    if table.empty:
        raise ValueError('Empty table.')

    # filter out depending on selected samples
    if sample_metadata is not None:
        sample_metadata = sample_metadata.to_dataframe()
        table = table.loc[:,sample_metadata.index.values]

    # filter out depending on selected reactions, exchanges or subsystems
    if feature_metadata is not None:
        feature_metadata = feature_metadata.to_dataframe()
        table = table.loc[feature_metadata.index.values,:]
    
    clustermap = sns.clustermap(table, xticklabels = table.columns,
                                method=method, metric = metric,
                                **_clustering_map[cluster],
                                cmap=color_scheme)
    
    if title is not None:
        clustermap.fig.suptitle(title)
   
    if xlabels:
        plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    if ylabels:
        plt.setp(clustermap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    img_fp = os.path.join(output_dir, 'clustermap.png')
    clustermap.savefig(img_fp, dpi=100, bbox_inches='tight')

    index_fp = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index_fp, output_dir, context={'metric': metric, 'method': method})
    