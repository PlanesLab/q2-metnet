# -*- coding: utf-8 -*-

import pickle
import pandas as pd
import numpy as np
import re

# Define the class for the supra-organism model
class SupraModel:
    
    def __init__(self, S, reactions, metabolites, taxonomy):
        # self.S = S
        self.rxnID = reactions.rxnID.values
        self.rxns = reactions.rxns.values
        self.rxnNames = reactions.rxnNames.values
        self.subSystems = reactions.subSystems.values
        self.lb = reactions.lb.values
        self.ub = reactions.ub.values
        self.c = reactions.c.values
        #self.rxnConfidenceScores = reactions.rxnConfidenceScores.values
        #self.rxnECNumbers = reactions.rxnECNumbers.values
        #self.rxnKEGGID = reactions.rxnKEGGID.values
        #self.rxnAliases = reactions.rxnAliases.values
        #self.rxnSource = reactions.rxnSource.values
        #self.comments = reactions.comments.values
        #self.citations = reactions.citations.values
        self.eqMet = reactions.eqMet.values
        self.eqS = reactions.eqS.values
        self.rxnTax = reactions.taxonomy.values.astype(str)
        self.metID = metabolites.metID.values
        self.mets = metabolites.mets.values
        self.metNames = metabolites.metNames.values
        self.metFormulas = metabolites.metFormulas.values
        self.metCharges = metabolites.metCharges.values
        self.metKEGGID = metabolites.metKEGGID.values
        self.metHMDBID = metabolites.metHMDBID.values
        #self.metChEBIID = metabolites.metChEBIID.values
        #self.metPubChemID = metabolites.metPubChemID.values
        #self.metInChIString = metabolites.metInChIString.values
        #self.metSmiles = metabolites.metSmiles.values
        #self.metAliases = metabolites.metAliases.values
        #self.metSource = metabolites.metSource.values
        self.b = metabolites.b.values
        #self.csense = metabolites.csense.values
        self.Taxonomy = taxonomy.taxonomy.values
    
    # Define a module to extract the reactions' index present in a list of specific taxa
    # Extract the whole list of indexes and a dictionary with index of reacions as key and the count
    # of appearances across taxa as value
    def indexPresentReactonsAndCounts(self, taxa: list):
        idx = []
        counts = {}
        for idx_rxn in range(len(self.rxnID)):
            present_tax = [self.Taxonomy[int(x)-1]  for x in self.rxnTax[idx_rxn].split(";")]
            if any([x in present_tax for x in taxa]):
                idx.append(idx_rxn)
                counts[idx_rxn] = sum([x in present_tax for x in taxa])
        
        return idx, counts
    
    # Define a module to extract the subsystems related to the reactions present in a list of specific taxa
    # Extract the subsystems and the related counts
    def relatedSubSystems(self, taxa: list):
        idx_rxn, counts = self.indexPresentReactonsAndCounts(taxa)
        whole_list_subsystem = self.subSystems[idx_rxn]
        unique_counts = np.unique(whole_list_subsystem[~pd.isna(whole_list_subsystem)], 
                                  return_counts = True)
        return unique_counts
    
    # Define a module to extract the exchanges can be generated from a list of specific taxa
    # Extract the whole list of output and a dictionary with name of metabolite as key and the count
    # of appearances across taxa as value
    def Exchanges(self):#, taxa:list):
        # idx_rxn, counts = self.indexPresentReactonsAndCounts(taxa)
        exchanges = set()
        # ex_counts = {}
        # for each_idx in idx_rxn:
        for each_idx in range(len(self.rxnID)):
            if pd.isna(self.subSystems[each_idx]):
                continue
            if re.findall("^Exchange/demand reaction",self.subSystems[each_idx]):
                tmp_ex = self.metNames[np.where(self.mets == self.eqMet[each_idx])][0]
                # ex = self.metID[np.where(self.mets == self.eqMet[each_idx])][0]
                if re.findall("mucin", tmp_ex):
                    continue
                exchanges.add(self.rxnID[each_idx])
                # if not ex in exchanges:
                #     exchanges.append(ex)
                    # ex_counts[out] = counts[each_idx]
        return exchanges#, ex_counts
    

# Load the files to crete the model    
def _generateModel(filename_reactions, filename_metabolites, filename_taxonomy):
    # Load the files
    reactions = pd.read_csv(filename_reactions, sep = ",", encoding = "ISO-8859-1")
    metabolites = pd.read_csv(filename_metabolites, sep = ",", encoding = "ISO-8859-1")
    taxonomy = pd.read_csv(filename_taxonomy, sep = ",", encoding = "ISO-8859-1")
    # Create the model class
    return SupraModel([], reactions, metabolites, taxonomy)
    
