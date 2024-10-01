# -*- coding: utf-8 -*-

from pkg_resources import resource_filename
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sns
import os
import q2templates
import qiime2

TEMPLATES = resource_filename("q2_metnet._pca", "assets")

pca_choices = {
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

        
def plotPCA(output_dir: str, table: pd.DataFrame, 
            sample_metadata: qiime2.CategoricalMetadataColumn = None,
            feature_metadata: qiime2.CategoricalMetadataColumn = None, 
            point_label: bool = False, title : str = None,
            color_scheme: str = 'rocket') -> None:
    
    def label_point(x, y, val, ax):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            ax.text(point['x']+.05, point['y']+.05, str(point['val']))
    
    if table.empty:
        raise ValueError('Empty table.')

    # filter out depending on selected samples
    if sample_metadata is not None:
        # Transform metadata to dataframe
        sample_metadata = sample_metadata.to_dataframe()
        # Filter metadata
        sample_metadata = sample_metadata.loc[sample_metadata.index.isin(table.columns)]
        # Filter scores table
        table = table.loc[:,sample_metadata.index.values]

    # filter out depending on selected reactions, exchanges or subsystems
    if feature_metadata is not None:
        feature_metadata = feature_metadata.to_dataframe()
        table = table.loc[feature_metadata.index.values,:]

    pca = PCA(n_components = 2)
    pca_fit = pca.fit_transform(table.transpose())
    
    points = pd.DataFrame(data = {'label': table.columns.values})
    
    reduced_data = pd.DataFrame(data = {'0':pca_fit[:,0],
                                        '1':pca_fit[:,1],
                                        'label':[x[0] for x in sample_metadata.values]})
    
    ax = sns.lmplot(x='0',
                y='1',
                data=reduced_data,
                fit_reg=False,
                legend=False,
                height=9,
                hue='label',
                palette=color_scheme)
    
    plt.xlabel("PCA 1 ("+str(round(pca.explained_variance_ratio_[0]*100,2))+'%)',
                weight = "bold", fontsize = 14, labelpad=15)
    plt.ylabel("PCA 2 ("+str(round(pca.explained_variance_ratio_[1]*100,2))+'%)',
                weight = "bold", fontsize = 14, labelpad=15)
    
    plt.legend(fontsize = 14, loc = "best", frameon=False)
    
    if title is not None:
        ax.set(title = title)

    sns.set_context("notebook", font_scale=1.2)
    sns.set_style("ticks")
    
    if point_label:
        label_point(reduced_data['0'], reduced_data['1'], points['label'], plt.gca())

    img_fp = os.path.join(output_dir, 'pca.png')
    ax.savefig(img_fp, dpi=100, bbox_inches='tight')

    index_fp = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index_fp, output_dir)
