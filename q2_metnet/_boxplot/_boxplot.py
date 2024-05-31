# -*- coding: utf-8 -*-

from pkg_resources import resource_filename
import matplotlib.pyplot as plt
import pandas as pd
import os
import q2templates
import qiime2
import re

TEMPLATES = resource_filename("q2_metnet._boxplot", "assets")

def plotBoxplot(output_dir: str,  table: pd.DataFrame, differentialresults: pd.DataFrame, 
                sample_metadata: qiime2.CategoricalMetadataColumn, namefeature: str,
                condition_name: str, control_name: str, title: str = None) -> None:
    
    sample_metadata = sample_metadata.to_dataframe()
    condition = [sample_metadata.index[x] for x in range(len(sample_metadata)) if sample_metadata.values[x] == condition_name]
    control = [sample_metadata.index[x] for x in range(len(sample_metadata)) if sample_metadata.values[x] == control_name]

    if re.findall("S\d+ \| ", namefeature):
        tmp_name = namefeature.split(" | ")[1]
    else:
        tmp_name = namefeature.split(" | ")[0]

    try:
        control_group = table.loc[tmp_name,control]
        condition_group = table.loc[tmp_name,condition]
    except:
        raise AttributeError("The ID for reaction/subsystem or for the samples are not present in the data")
    
    fig, ax = plt.subplots()
    plt.boxplot([control_group,condition_group], positions=[0.5,1])
    plt.xticks([0.5,1], [control_name, condition_name])
    ax.set_xlim(left = 0.2, right=1.3)
    
    temp = differentialresults.index[differentialresults.index == namefeature].values[0]
    name = temp.split(" | ")[1]
    adj_pval = differentialresults.Adjusted_p_Value[differentialresults.index == namefeature].values[0]
    fold_change = differentialresults.FC[differentialresults.index == namefeature].values[0]
    if title is None:
        plt.title(name.capitalize() + ' (adj.pval = %f)' % (adj_pval))
    else:
        plt.title(title)
    plt.ylabel(r"Activity Score", weight = "bold", fontsize = 12, labelpad=15)
    
    plt.xlabel(r"Condition", weight = "bold", fontsize = 12, labelpad=15)
    
    img_fp = os.path.join(output_dir, 'boxplot.png')
    plt.savefig(img_fp, dpi=100, bbox_inches='tight')

    index_fp = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index_fp, output_dir, context={'p_value': adj_pval, 'fold_change': fold_change})