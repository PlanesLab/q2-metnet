# Data types
import biom

# Functions
import pandas as pd
import pkg_resources
from q2_metnet._generateNetwork import _generateModel
from q2_metnet._inputFiles import _extractTaxaPresentAGREDA

def _reactionsBetweenSamples(Samples, PresentTaxa, Frequency, Model, rxnTax):
    
    count_rxns = pd.DataFrame(index = Model.rxnID, columns = Frequency.ID.values)
    n_asv = 0
    for keys, values in PresentTaxa.items():
        print("Calculating Reaction scores, ASV: %d/%d" % (n_asv,len(PresentTaxa.keys())-1))
        n_asv += 1
        tmp =rxnTax.iloc[:,[x-1 for x in values["TAXA"].index.values]]
        count_rxns.loc[:,keys] = tmp.T.mean().values
        
    count_rxns.fillna(0, inplace = True)
    
    tmp_freq = Frequency.loc[:,Samples]
    tmp_freq.index = count_rxns.columns
    Reactions = count_rxns.dot(tmp_freq)
        
    return Reactions

def _subsystemsBetweenSamples(Reactions, Model, class_exchange):
    
    tmp_sub = pd.DataFrame(data = {'Sub': Model.subSystems})
    tmp_sub.fillna('', inplace = True)
    tmp_sub = tmp_sub.Sub.drop_duplicates()
    
    all_sub = []
    
    for x in tmp_sub:
        if x == "":
            continue
        y = x.split(";")
        all_sub += y
    
    all_sub = list(set(all_sub))
    all_sub += list(class_exchange.Class.values)
    all_sub = list(set(all_sub))
    
    sub_rxns = pd.DataFrame(index = all_sub, columns = Reactions.index)
    
    for x in range(len(Model.rxnID)):
        print("Reaction: %d/%d" % (x,len(Model.rxnID)-1))
        try:
            tmp = class_exchange.loc[class_exchange.rxnID.isin([Model.rxnID[x]]),'Class'].values[0]
        except IndexError:
            tmp = Model.subSystems[x]
        try:
            tmp = tmp.split(";")
        except AttributeError:
            continue
        for each in tmp:
            sub_rxns.loc[each,Model.rxnID[x]] = 1
    
    sub_rxns.fillna(0,inplace = True)
    
    sub_rxns = sub_rxns.div(sub_rxns.sum(axis=1), axis=0)
    sub_rxns.dropna(inplace = True)    
    SubSystems_Sample = sub_rxns.dot(Reactions)
    
    temp = [' | '.join(['S%d' % x,SubSystems_Sample.index[x]]) for x in range(len(SubSystems_Sample.index))]
    SubSystems_Sample.index = temp
    
    return SubSystems_Sample

def generateFeatures(frequency: biom.Table, taxa: pd.DataFrame, 
                     selection: str = 'AGREDA', level: str = "s", input_interest: str = True) -> (pd.DataFrame,pd.DataFrame,pd.DataFrame):
    if selection == "AGREDA":
        stream_reactions = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_spInfo.tsv')
        
        PresentTaxa, newFrequency, Samples = _extractTaxaPresentAGREDA(frequency, taxa, stream, level)

        stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_rxnTaxMat.csv')
        rxnTax = pd.read_csv(stream, sep = ",")
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_Exchange_metabolites.tsv')
        class_exchange = pd.read_csv(stream, sep = "\t")
    elif selection == "AGORAv103":
        stream_reactions = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv103/AGORA_v1.0.3_spInfo.tsv')
        
        PresentTaxa, newFrequency, Samples = _extractTaxaPresentAGREDA(frequency, taxa, stream, level)

        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_rxnTaxMat.csv')
        rxnTax = pd.read_csv(stream, sep = ",")

        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_Exchange_metabolites.tsv')
        class_exchange = pd.read_csv(stream, sep = "\t")
    elif selection == "AGORAv201":
        stream_reactions = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__, 'data/AGORAv201/AGORA_v2.0.1_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv201/AGORA_v2.0.1_spInfo.tsv')
        
        PresentTaxa, newFrequency, Samples = _extractTaxaPresentAGREDA(frequency, taxa, stream, level)

        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv201/AGORA_v2.0.1_rxnTaxMat.csv')
        rxnTax = pd.read_csv(stream, sep = ",")

        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv201/AGORA_v2.0.1_Exchange_metabolites.tsv')
        class_exchange = pd.read_csv(stream, sep = "\t")
    else:
        raise ValueError("Select a valid metabolic reconstruction among: AGREDA, AGORAv103, AGORAv201")

    Reactions = _reactionsBetweenSamples(Samples, PresentTaxa, newFrequency, Model, rxnTax)
    Subsystems = _subsystemsBetweenSamples(Reactions, Model, class_exchange)
    
    Xmatrix = newFrequency.copy()
    new_index = []
    for idx in Xmatrix.index:
        tmp_asv = Xmatrix.ID[idx]
        tmp_taxa = list(PresentTaxa[tmp_asv]['TAXA'].loc[:,'AGORA.NAMES'].values)
        new_index.append(';'.join([Xmatrix.ID[idx]] + tmp_taxa))
    Xmatrix.drop(columns = ['LINEAGE','ID'], inplace = True)
    Xmatrix.index = new_index
    
    return (Reactions, Subsystems, Xmatrix)
