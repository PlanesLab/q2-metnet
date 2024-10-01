# -*- coding: utf-8 -*-

import pkg_resources
from scipy import stats
from statsmodels.stats import multitest
import pandas as pd
import qiime2
import biom
from q2_metnet._generateNetwork import _generateModel

def _extractExchanges(reactions, Model, SampleID, input_interest, stream):
    if input_interest:
        inputs = pd.read_csv(stream, sep = "\t", encoding = "ISO-8859-1")
        exchanges = reactions.loc[inputs.inputs,SampleID]
    else:
        exchanges = reactions.loc[Model.Exchanges(),SampleID]
        
    return exchanges

def differentialExchanges(reactions: biom.Table, metadata: qiime2.MetadataColumn, condition_name: str, control_name: str,
                           selection_model: str = "AGREDA", input_interest: str = True) -> pd.DataFrame:
    
    df_reactions = reactions.to_dataframe().transpose()
    df_metadata = metadata.to_dataframe()
    df_condition_list = metadata.to_series()
    
    df_reactions = df_reactions.loc[:,df_metadata.index.values]
    condition_sample = df_metadata.index[[x == condition_name for x in df_condition_list]]
    control_sample = df_metadata.index[[x == control_name for x in df_condition_list]]
    
    if selection_model == "AGREDA":
        stream_reactions = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)

        stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_Exchange_metabolites.tsv')
        ex_mets = pd.read_csv(stream, sep = "\t", encoding = "ISO-8859-1")

        stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_Input_reactions.tsv')
        exchanges = _extractExchanges(df_reactions, Model, df_metadata.index.values, input_interest, stream)       
    elif selection_model == "AGORAv103":
        stream_reactions = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_Exchange_metabolites.tsv')
        ex_mets = pd.read_csv(stream, sep = "\t", encoding = "ISO-8859-1")
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_Input_reactions.tsv')
        exchanges = _extractExchanges(df_reactions, Model, df_metadata.index.values, input_interest, stream)
    elif selection_model == "AGORAv201":
        stream_reactions = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv201/AGORA_v2.0.1_Exchange_metabolites.tsv')
        ex_mets = pd.read_csv(stream, sep = "\t", encoding = "ISO-8859-1")
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv201/AGORA_v2.0.1_Input_reactions.tsv')
        exchanges = _extractExchanges(df_reactions, Model, df_metadata.index.values, input_interest, stream)
    else:
        raise ValueError("Select a valid metabolic reconstruction among: AGREDA, AGORAv103, AGORAv201")
    
    results = pd.DataFrame(index = exchanges.index, columns = ['p-value', 'FC'])
    
    for idx in exchanges.index:
        
        control_group = exchanges.loc[idx,control_sample]
        condition_group = exchanges.loc[idx,condition_sample]
        try:
            p_value = stats.mannwhitneyu(list(control_group.values),list(condition_group.values)).pvalue
        except ValueError:
            p_value = 1
        results.loc[idx,'p-value'] = p_value
        results.loc[idx,'FC'] = condition_group.mean() - control_group.mean()
    
    p_adj = multitest.fdrcorrection(results['p-value'].values)
    temp = list(Model.rxnID)
    adjusted_results = pd.DataFrame(data = {'FC': results['FC'].values,
                                            "p_Value": results['p-value'].values,
                                            "Adjusted_p_Value": p_adj[1]})

    metnames = [ex_mets.loc[ex_mets.rxnID.isin([x]),'metNames'].values[0] for x in results.index]
    rxnID = results.index.values
    temp = [' | '.join([rxnID[x],metnames[x]]) for x in range(len(rxnID))]
    adjusted_results.index = temp

    # Sort by adjusted p-values and the absolute value of FC
    adjusted_results['Absolute_FC'] = adjusted_results['FC'].abs()
    sorted_adjusted_results = adjusted_results.sort_values(by=['Adjusted_p_Value', 'Absolute_FC'], ascending=[True, False])
    sorted_adjusted_results = sorted_adjusted_results.drop(columns=['Absolute_FC'])
    
    return sorted_adjusted_results

def differentialSubSystems(subsystems: biom.Table, metadata: qiime2.MetadataColumn, condition_name: str, control_name: str) -> pd.DataFrame:

    df_subsystem = subsystems.to_dataframe().transpose()
    df_metadata = metadata.to_dataframe()
    df_condition_list = metadata.to_series()
    
    df_subsystem = df_subsystem.loc[:,df_metadata.index.values]
    condition_sample = df_metadata.index[[x == condition_name for x in df_condition_list]]
    control_sample = df_metadata.index[[x == control_name for x in df_condition_list]]

    results = pd.DataFrame(index = df_subsystem.index, columns = ['p-value', 'FC'])

    for idx in results.index:
        control_group = df_subsystem.loc[idx,control_sample]
        condition_group = df_subsystem.loc[idx,condition_sample]
        try:
            p_value = stats.mannwhitneyu(list(control_group.values),list(condition_group.values)).pvalue
        except ValueError:
            p_value = 1
        results.loc[idx,'p-value'] = p_value
        results.loc[idx,'FC'] = condition_group.mean() - control_group.mean()
    
    p_adj = multitest.fdrcorrection(results['p-value'].values)
    adjusted_results = pd.DataFrame(data = {'FC': results['FC'].values,
                                            "p_Value": results['p-value'].values,
                                            "Adjusted_p_Value":p_adj[1]})
    
    # temp = [' | '.join(['S%d' % x,results.index[x]]) for x in range(len(results.index))]
    # adjusted_results.index = temp
    adjusted_results.index = results.index

    # Sort by adjusted p-values and the absolute value of FC
    adjusted_results['Absolute_FC'] = adjusted_results['FC'].abs()
    sorted_adjusted_results = adjusted_results.sort_values(by=['Adjusted_p_Value', 'Absolute_FC'], ascending=[True, False])
    sorted_adjusted_results = sorted_adjusted_results.drop(columns=['Absolute_FC'])
    
    return sorted_adjusted_results

def differentialReactions(reactions: biom.Table, metadata: qiime2.MetadataColumn, condition_name: str, control_name: str,
                           selection_model: str = "AGREDA") -> pd.DataFrame:

    if selection_model == "AGREDA":
        stream_reactions = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)        
    elif selection_model == "AGORAv103":
        stream_reactions = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)
    elif selection_model == "AGORAv201":
        stream_reactions = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)
    else:
        raise ValueError("Select a valid metabolic reconstruction among: AGREDA, AGORAv103, AGORAv201")

    df_reactions = reactions.to_dataframe().transpose()
    df_metadata = metadata.to_dataframe()
    df_condition_list = metadata.to_series()
    
    df_reactions = df_reactions.loc[:,df_metadata.index.values]
    condition_sample = df_metadata.index[[x == condition_name for x in df_condition_list]]
    control_sample = df_metadata.index[[x == control_name for x in df_condition_list]]

    results = pd.DataFrame(index = df_reactions.index, columns = ['p-value', 'FC'])
    
    for idx in results.index:
        
        control_group = df_reactions.loc[idx,control_sample]
        condition_group = df_reactions.loc[idx,condition_sample]
        try:
            p_value = stats.mannwhitneyu(list(control_group.values),list(condition_group.values)).pvalue
        except ValueError:
            p_value = 1
        results.loc[idx,'p-value'] = p_value
        results.loc[idx,'FC'] = condition_group.mean() - control_group.mean()

    p_adj = multitest.fdrcorrection(results['p-value'].values)
    temp = list(Model.rxnID)
    adjusted_results = pd.DataFrame(data = {'FC': results['FC'].values,
                                            "p_Value": results['p-value'].values,
                                            "Adjusted_p_Value": p_adj[1]})
    
    # rxnnames = [Model.rxns[temp.index(x)] for x in results.index]
    rxnID = results.index.values
    # temp = [' | '.join([rxnID[x],rxnnames[x]]) for x in range(len(rxnID))]
    # adjusted_results.index = temp
    adjusted_results.index = rxnID

    # Sort by adjusted p-values and the absolute value of FC
    adjusted_results['Absolute_FC'] = adjusted_results['FC'].abs()
    sorted_adjusted_results = adjusted_results.sort_values(by=['Adjusted_p_Value', 'Absolute_FC'], ascending=[True, False])
    sorted_adjusted_results = sorted_adjusted_results.drop(columns=['Absolute_FC'])
    
    return sorted_adjusted_results
