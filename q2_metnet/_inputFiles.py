# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

# Extract the corresponding models in the database to the taxonomic level of the samples
# Define the level of interest in the phylogenetic tree

def _fillTaxa(taxon):

    ##########################################################################
    # Fill missing taxonomic levels of the FeatureData[Taxonomy]
    ##########################################################################

    # Split string
    levels = taxon.split(';')
    
    # Number of levels
    n = len(levels)
    
    # Taxonomic levels
    taxLevels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    
    # Append missing levels with placeholders
    if n < len(taxLevels):
        missLevels = taxLevels[n:]
        levels.extend(missLevels)
        
    # Join back the taxonomic levels
    return ';'.join(levels)

def _contextTaxa(frequency, taxa, Reference, level):

    ##########################################################################
    # Filter frequency table, removing duplicates and summing their counts
    ##########################################################################

    correspondent_level = {"k":"KINGDOM",
                           "p":"PHYLUM",
                           "c":"CLASS",
                           "o":"ORDER",
                           "f":"FAMILY",
                           "g":"GENUS"}
    
    depth_level = {"k":0,
                   "p":1,
                   "c":2,
                   "o":3,
                   "f":4,
                   "g":5,
                   "s":6}
    
    if level not in depth_level.keys():
        raise ValueError("Select a valid lineage level among: k, p, c, o, f, g, s (kingdom, phylum, class, order, family, genus, species)")
    
    # Transform feature table to dataframe
    df_frequency = frequency.to_dataframe()
    
    # Extract size of feature table and feature dataframe
    nR1, nC1 = df_frequency.shape
    nR2, nC2 = taxa.shape
    
    # Transpose feature table if necessary
    if nR1 != nR2:
        df_frequency = df_frequency.transpose()
    
    # Extract sample IDs
    samples = df_frequency.columns.values
    
    # Fill all taxonomic levels of FeatureData[Taxonomy]
    taxa['Taxon'] = taxa['Taxon'].apply(_fillTaxa)
    
    # Transfrom FeatureData[Taxonomy] from dataframe to dictionary
    dictionary_otu_taxa = taxa['Taxon'].to_dict()

    df_frequency['LINEAGE'] = [dictionary_otu_taxa[x] for x in df_frequency.index]
    df_frequency['LINEAGE'] = df_frequency['LINEAGE'].str.replace("NA", ""). str.replace("\w__","", regex =True). str.replace("[\[\]]","", regex =True)
    
    to_rem = []
    for idx in df_frequency.index:
        tmp = df_frequency['LINEAGE'][idx].split(";")
        tmp = [x.strip() for x in tmp]
        if len(tmp) < depth_level[level]:
            to_rem.append(idx)
            continue
        if not "" == tmp[depth_level[level]]:
            df_frequency.loc[idx, 'LINEAGE'] = ";".join(tmp[:depth_level[level]+1])
        else:
            to_rem.append(idx)
        
    df_frequency.drop(to_rem, inplace = True)
    
    unique_taxa = df_frequency.LINEAGE.drop_duplicates()
    
    new_Frequency = pd.DataFrame(index = range(len(unique_taxa)), columns = ["LINEAGE", "ID"] + list(samples))
    
    count_asv = 0
    for each_taxa in unique_taxa:
        tmp = df_frequency.loc[df_frequency.LINEAGE.isin([each_taxa]),samples]
        
        tmp_freq = tmp.sum(axis = 0)
        tmp_freq['ID'] = "M_ASV%d" % count_asv
        tmp_freq['LINEAGE'] = each_taxa
        
        new_Frequency.loc[count_asv,:] = tmp_freq
        
        count_asv += 1
    
    for each_sample in samples:
        new_Frequency.loc[:,each_sample] = new_Frequency.loc[:,each_sample].div(new_Frequency.loc[:,each_sample].sum())
    
    # ##########################################################################
    # ##########################################################################
    
    results = {}
    to_rem = []
    for idx in new_Frequency.index:
        
        # Extract the phyla of interest, continuing if empty or unknown
        tmp = new_Frequency.LINEAGE[idx].split(";")
        name_phyla = tmp[depth_level[level]]
            
        # if interested on species level, extract the specific species, 
        # otherwise extract all the species related to the level of interest
        if level.lower() == "s":
            try:
                genus = tmp[depth_level[level]-1].split("__")[1]
            except:
                genus = tmp[depth_level[level]-1]
                
            name_phyla = " ".join([genus,name_phyla])
            
            present = []
            for index in Reference.index:
                if name_phyla in Reference['AGORA.NAMES'][index] or \
                                    name_phyla in Reference['NCBI.NAMES'][index]:
                    present.append(index)
            if len(present) == 0:
                to_rem.append(idx)
                continue
            results[new_Frequency['ID'][idx]] = {"NAME_LEVEL": level.lower()+"__"+name_phyla,
                                                    "TAXA": Reference.loc[present,]}
            
        else:
            present = []
            for index in Reference.index:
                if Reference[correspondent_level[level]][index] is np.nan:
                    continue
                if name_phyla in Reference[correspondent_level[level]][index]:
                    present.append(index)
            if len(present) == 0:
                to_rem.append(idx)
                continue
            results[new_Frequency['ID'][idx]] = {"NAME_LEVEL": level.lower()+"__"+name_phyla,
                                                    "TAXA": Reference.loc[present,]}
            
    new_Frequency.drop(index = to_rem, inplace = True)
    new_Frequency.index = range(len(new_Frequency))
    if len(new_Frequency.index) == 0:
        raise Warning("No samples have a correspondent species in AGREDA")

    return results, new_Frequency, samples
    

def _extractTaxaPresentAGREDA(frequency, taxa, stream, level):
    
    ##########################################################################
    # Load database species and extract species model information
    ##########################################################################
    
    # Load the reference for each database
    reference = pd.read_csv(stream, sep = "\t")
    
    # Extract the samples taxa present in the database
    return _contextTaxa(frequency, taxa, reference, level)
