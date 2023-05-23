import qiime2.plugin
from qiime2.plugin import  MetadataColumn, Categorical, Str, Bool, Choices
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy

import q2_metnet
from q2_metnet._generateFeatures import generateFeatures
from q2_metnet._functional_analysis import differentialSubSystems, differentialReactions, differentialExchanges
from q2_metnet._clustermap import plotClusteMap, clustermap_choices
from q2_metnet._pca import plotPCA, pca_choices
from q2_metnet._boxplot import plotBoxplot

cites = qiime2.plugin.Citations.load('citations.bib',
    package='q2_metnet')

plugin = qiime2.plugin.Plugin(
    name='metnet',
    version=q2_metnet.__version__,
    website='',
    package='q2_metnet',
    citations=[cites['blasco2021extended']],
    description=('Package to contextualize table of taxonomical frequency of samples to metabolic reconstruction (AGORA or AGREDA) and extract features based on score of active reactions or subsystems.'),
    short_description='Plugin for taxonomy contextualization to metabolic reconstruction and feature extraction.',
    user_support_text=('Raise an issue to the corresponding author: fplanes@tecnun.es')
)

# Register GenerateFeatures function
plugin.methods.register_function(
    function=generateFeatures,
    inputs={'frequency': FeatureTable[Frequency],
            'taxa': FeatureData[Taxonomy]
    },
    outputs=[('reactions', FeatureTable[Frequency]),
             ('subsystems', FeatureTable[Frequency])
             ],
    input_descriptions={'frequency': 'table of frequency',
        'taxa': 'table of assigned taxonomy'
    },
    parameters={'selection': Str,
                'level': Str},
    output_descriptions={'reactions': 'Reaction scores based on the samples and the taxonomy present in the selected reconstruction',
                         'subsystems': 'Subsystem scores based on the samples and the taxonomy present in the selected reconstruction'
                         },
    parameter_descriptions={'selection': 'selection metabolic network among AGREDA, AGORAv103, AGORAv201',
                            'level': 'taxonomical level of interest: k (kingdom), p (phylum), c (class), o (order), f (family), g (genus), s (species, default)'},
    name='Reactions and subsystems features extraction',
    description='Extraction of the score related to each reaction and subsystem present in the metabolic reconstruction considering the taxonomy included in the samples'
)

# Register DifferentialSubsystems function
plugin.methods.register_function(
    function=differentialSubSystems,
    inputs={'subsystems': FeatureTable[Frequency]},
    outputs=[('differential_analysis', FeatureTable[Frequency])],
    input_descriptions={'subsystems': 'table of frequency'},
    parameters={'metadata': MetadataColumn[Categorical],
                'condition_name': Str,
                'control_name': Str},
    output_descriptions={'differential_analysis': 'Differential analysis of the subsystems scores'},
    parameter_descriptions={'metadata': 'list of the condition states',
                            'condition_name': 'name of the condition category under analysis, taken from the metadata file',
                            'control_name': 'name of the control category under analysis, taken from the metadata file'},
    name='Differential score analysis of the subsystems',
    description='Differential score analysis of the subsystems'
)

# Register DifferentialReactions function
plugin.methods.register_function(
    function=differentialReactions,
    inputs={'reactions': FeatureTable[Frequency]},
    outputs=[('differential_analysis', FeatureTable[Frequency])],
    input_descriptions={'reactions': 'table of frequency'},
    parameters={'metadata': MetadataColumn[Categorical],
                'condition_name': Str,
                'control_name': Str,
                'selection_model': Str},
    output_descriptions={'differential_analysis': 'Differential analysis of the reactions scores'},
    parameter_descriptions={'metadata': 'list of the condition states',
                            'condition_name': 'name of the condition category under analysis, taken from the metadata file',
                            'control_name': 'name of the control category under analysis, taken from the metadata file',
                            'selection_model': 'selection of the metabolic reconstruction among AGREDA, AGORAv103, AGORAv201'},
    name='Differential score analysis of the subsystems',
    description='Differential score analysis of the subsystems'
)

# Register DifferentialExchanges function
plugin.methods.register_function(
    function=differentialExchanges,
    inputs={'reactions': FeatureTable[Frequency]},
    outputs=[('differential_analysis', FeatureTable[Frequency])],
    input_descriptions={'reactions': 'table of frequency'},
    parameters={'metadata': MetadataColumn[Categorical],
                'condition_name': Str,
                'control_name': Str,
                'selection_model': Str, 
                'input_interest': Bool},
    output_descriptions={'differential_analysis': 'Differential analysis of the exchanges scores'},
    parameter_descriptions={'metadata': 'list of the condition states',
                            'condition_name': 'name of the condition category under analysis, taken from the metadata file',
                            'control_name': 'name of the control category under analysis, taken from the metadata file',
                            'selection_model': 'selection of the metabolic reconstruction among AGREDA, AGORAv103, AGORAv201',
                            'input_interest': 'Boolean to define if focus on the exchanges that can be input (True, default) or all of them (False)'},
    name='Differential score analysis of the subsystems',
    description='Differential score analysis of the subsystems'
)

# Register visualizer for the clustermap
plugin.visualizers.register_function(
    function=plotClusteMap,
    inputs={
        'table': FeatureTable[Frequency]
    },
    parameters={
        'sample_metadata': MetadataColumn[Categorical],
        'feature_metadata': MetadataColumn[Categorical],
        'title': Str,
        'metric': Str % Choices(clustermap_choices['metric']),
        'method': Str % Choices(clustermap_choices['method']),
        'cluster': Str % Choices(clustermap_choices['cluster']),
        'color_scheme': Str % Choices(clustermap_choices['color_scheme']),
        'xlabels': Bool, 
        'ylabels': Bool,
    },
    name='Generate a clustermap representation of reactions or subsystems feature tables',
    description='Generate a clustermap representation of reactions or subsystems feature tables',
    input_descriptions={
        'table': 'The feature table to visualize.'
    },
    parameter_descriptions={
        'sample_metadata': 'list of the condition states of interest. Optional',
        'feature_metadata': 'list of the features (reactions or subsystems) of interest. Optional',
        'title': 'Optional title for the plot.',
        'metric': 'Metrics exposed by seaborn (see http://seaborn.pydata.org/generated/seaborn.clustermap.html#seaborn.clustermap for more detail). Euclidean is the default.',
        'method': 'Clustering methods exposed by seaborn (see http://seaborn.pydata.org/generated/seaborn.clustermap.html#seaborn.clustermap for more detail). Average is the default.',
        'cluster': 'Specify which axes to cluster. Both is the default.',
        'color_scheme': 'The matplotlib colorscheme to generate the clustermap.',
        'xlabels': 'boolean to choose if display xlabels', 
        'ylabels': 'boolean to choose if display ylabels',
    }
)

# Register visualizer for the clustermap
plugin.visualizers.register_function(
    function=plotPCA,
    inputs={
        'table': FeatureTable[Frequency]
    },
    parameters={
        'sample_metadata': MetadataColumn[Categorical],
        'feature_metadata': MetadataColumn[Categorical],
        'title': Str,
        'color_scheme': Str % Choices(pca_choices['color_scheme']),
        'point_label': Bool, 
    },
    name='Generate a PCA representation of reactions or subsystems feature table',
    description='Generate a PCA representation of reactions or subsystems feature table.',
    input_descriptions={
        'table': 'The feature table to visualize.'
    },
    parameter_descriptions={
        'sample_metadata': 'list of the condition states of interest. Optional.',
        'feature_metadata': 'list of the features (reactions or subsystems) of interest. Optional.',
        'title': 'Optional title for the plot.',
        'color_scheme': 'The matplotlib colorscheme to generate the PCA',
        'point_label': 'boolean to choose if display labels for the PCA points', 
    }
)

# Register visualizer for the clustermap
plugin.visualizers.register_function(
    function=plotBoxplot,
    inputs={
        'table': FeatureTable[Frequency],
        'differentialresults': FeatureTable[Frequency]
    },
    parameters={
        'sample_metadata': MetadataColumn[Categorical],
        'namefeature': Str,
        'title': Str,
        'condition_name': Str,
        'control_name': Str, 
    },
    name='Generate a boxplot representation of a reaction, exchange or subsystem',
    description='Generate a boxplot representation of a reaction, exchange or subsystem to compare the scores in the condition to the ones in the control.',
    input_descriptions={
        'table': 'The whole feature table.',
        'differentialresults': 'Results of diffential expression analysis'
    },
    parameter_descriptions={
        'sample_metadata': 'list of the condition states of interest.',
        'namefeature': 'name of the reaction, exchanges or subsystem to plot. The name must be the same as the one present in the differential analysis table.',
        'title': 'Optional title for the plot.',
        'condition_name': 'Name of the condition category, extracted from the metadata.',
        'control_name': 'Name of the control category, extracted from the metadata.', 
    }
)
