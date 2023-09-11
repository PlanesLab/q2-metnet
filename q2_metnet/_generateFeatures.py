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
        print("Calculate Reaction scores, ASV: %d/%d" % (n_asv,len(PresentTaxa.keys())-1))
        n_asv += 1
        tmp =rxnTax.iloc[:,[x-1 for x in values["TAXA"].index.values]]
        count_rxns.loc[:,keys] = tmp.T.mean().values
        
    count_rxns.fillna(0, inplace = True)
    
    tmp_freq = Frequency.loc[:,Samples]
    tmp_freq.index = count_rxns.columns
    Reactions = count_rxns.dot(tmp_freq)
        
    return Reactions

def _subsystemsBetweenSamples(Reactions, Model):
    
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
    
    sub_rxns = pd.DataFrame(index = all_sub, columns = Reactions.index)

    for x in range(len(Model.rxnID)):
        print("Calculate Subsystem scores, Reaction: %d/%d" % (x,len(Model.rxnID)-1))
        tmp = Model.subSystems[x]
        try:
            tmp = tmp.split(";")
        except AttributeError:
            continue
        for each in tmp:
            sub_rxns.loc[each,Model.rxnID[x]] = 1
        
    sub_rxns.fillna(0,inplace = True)
    
    sub_rxns = sub_rxns.div(sub_rxns.sum(axis=1), axis=0)
    
    SubSystems_Sample = sub_rxns.dot(Reactions)
    
    return SubSystems_Sample

def _subclassesBetweenSamples(Reactions, class_exchange, selection, input_interest):
    
    if input_interest:
        if selection == "AGREDA":
            stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_Input_reactions.tsv')
        elif selection_model == "AGORAv103":
            stream = pkg_resources.resource_stream(__name__, 'data/AGORAv103/AGORA_v1.0.3-M_Input_reactions.tsv')
        else:
            raise ValueError("Select a valid metabolic reconstruction among: AGREDA, AGORAv103")
        
        inputs = pd.read_csv(stream, sep = "\t", encoding = "ISO-8859-1")
        class_exchange = class_exchange.loc[inputs.inputs,:]
            
    sub_ex = pd.DataFrame(columns = class_exchange.rxnID, index = class_exchange.Class.drop_duplicates())

    for each_ex in sub_ex.columns:
        print("Calculate Classes scores, Reaction: %d/%d" % (x,len(Model.rxnID)-1))
        sub_class = class_exchange.loc[class_exchange.rxnID.isin([each_ex]),'Class'].values[0]
        sub_ex.loc[sub_class,each_ex] = 1

    sub_ex.fillna(0, inplace = True)
    ex_samples = Reactions.loc[sub_ex.columns,]

    sub_ex = sub_ex.div(sub_ex.sum(axis=1), axis=0)
    Classes_Exchange_Sample = sub_ex.dot(ex_samples)

    return Classes_Exchange_Sample

def generateFeatures(frequency: biom.Table, taxa: pd.DataFrame, 
                     selection: str = 'AGREDA', level: str = "s", input_interest: str = True) -> (pd.DataFrame,pd.DataFrame, pd.DataFrame):
    if selection == "AGREDA":
        stream_reactions = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_rxnInfo.csv')
        stream_metabolites = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_metInfo.csv')
        stream_taxonomy = pkg_resources.resource_filename(__name__,'data/AGREDA/AGREDA_taxonomy.csv')
        Model = _generateModel(stream_reactions, stream_metabolites, stream_taxonomy)
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_spInfo.tsv')
        
        PresentTaxa, newFrequency, Samples = _extractTaxaPresentAGREDA(frequency, taxa, stream, level)

        stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_rxnTaxMat.csv')
        rxnTax = pd.read_csv(stream, sep = ",")
        
        stream = pkg_resources.resource_stream(__name__, 'data/AGREDA/AGREDA_Exchange_metabolites_class.tsv')
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

        stream = pkg_resources.resource_stream(__name__, 'data/AGORAv103/AGORA103_Exchange_metabolites_class.tsv')
        class_exchange = pd.read_csv(stream, sep = "\t")
    else:
        raise ValueError("Select a valid metabolic reconstruction among: AGREDA, AGORAv103")

    Reactions = _reactionsBetweenSamples(Samples, PresentTaxa, newFrequency, Model, rxnTax)
    Subsystems = _subsystemsBetweenSamples(Reactions, Model)
    Classes = _subclassesBetweenSamples(Reactions, class_exchange, selection, input_interest)
    return (Reactions, Subsystems, Classes)
