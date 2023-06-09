# q2-metnet

QIIME 2 plugin to contextualize taxonomic abundance data into high-quality metabolic reconstruction of the human gut microbiota.
Read more about the method in our [paper] .

# Installing

You can install this plugin by cloning this repo and installing manually.

To clone:

```
git clone https://github.com/francesco-balzerani/q2-metnet.git
```

Before installing, please unzip the data folder within the q2_metnet folder:

```
cd q2_metnet
unzip data.zip
rm data.zip
cd ..
```

To install from this repo, check to be into the main directory (running the `ls` command the setup.py file must be displayed), and run:

```
python setup.py install
```

Then, update the plugin cache by typing:

```
qiime dev refresh-cache
```

You can check that the installation worked by typing `qiime` on the command line.
The `metnet` plugin should show up in the list of available plugins.

# Using the plugin

There are available four methods in this plugin: 
- `generateFeatures`, which creates the tables of normalized scores for any reaction and subsystem present in the selected metabolic reconstruction across the different samples under study;
- `differentialReactions`, which computes a differential activity analysis about any kind of reaction for the different conditions under analysis;
- `differentialExchanges`, which computes a differential activity analysis about exchange reactions for the different conditions under analysis;
- `differentialSubSystems`, which computes a differential activity analysis about subsystems for the different conditions under analysis;
- `plotClusteMap`, which permits to visualize different samples by a hierarchically-clustered heatmap;
- `plotPCA`, which allows to visualize different samples by conducting a Principal Component Analysis;
- `plotBoxplot`, which generates boxplot of both reaction or subsystem scores to visualize differences between conditions.

## Preparing your data

You'll need to prepare your ASV (OTU) table, respective taxonomical assignation and metadata file for use with this plugin.
Your ASV(OTU) table should be imported as a [QIIME 2 artifact](https://docs.qiime2.org/2019.1/concepts/#data-files-qiime-2-artifacts), with **ASVs (OTUs) in rows** and **samples in columns**.
Similarly, the taxonomical assignation should be imported as [QIIME 2 artifact](https://docs.qiime2.org/2019.1/concepts/#data-files-qiime-2-artifacts).

Metadata should be a tab-delimited file with a column that contains samples labeled depending on the type of samples they are (for example: `case` and `control`).
This column need to be defined to permit the plugin to compute the differential activity analysis.

If your OTU table is already a QIIME 2 artifact, you can skip directly to running the code.
Otherwise, follow the instructions on the QIIME 2 webpage to use your own tab-delimited OTU table.

## Generate the normalized scores tables

You then run the `generateFeatures` script from the `metnet` qiime plugin. The `AGREDA` and `s` parameter of `--p-selection` and `--p-level` are the default parameters.
The parameter `--p-selection` can be chosen among `AGREDA` (default), `AGORAv103` and `AGORAv201`.
The parameter `--p-level` correspond to the taxonomic depth, ranging from kingdom (`k`) to species (`s`, default).

```
qiime metnet generateFeatures \
	--i-frequency ../asv_table.asv.qza \
	--i-taxa ../assigned_taxonomy.qza \
	--p-selection AGREDA \
	--p-level s \
	--o-reactions ../output_reactions.qza \
	--o-subsystems ../output_subsystems.qza
```

## Calculate differential analysis scores for exchange reactions

You then run the `differentialExchanges` script from the `metnet` qiime plugin.
The parameter `--p-condition-name` represents the category of the samples related to the condition state in the metadata column.
The parameter `--p-control-name` represents the category of the samples related to the control state in the metadata column.
The parameter `--p-selection-model` corresponds to the same metabolic reconstruction used in the table generation (`AGREDA` (default), `AGORAv103`, `AGORAv201`).
The parameter `--p-input-interest` defines if the focus is just on those exchange that can be input (and may be output as well) or just output (default = True).

```
qiime metnet differentialExchanges \
	--i-reactions ../output_reactions.qza \
	--m-metadata-file ../metadata.tsv \
	--m-metadata-column columnLabel \
	--p-condition-name conditionLabel \
	--p-control-name controlLabel \
	--p-selection-model AGREDA \
	--p-input-interest True \
	--o-differential-analysis ../exchanges_differential.qza
```

## Calculate differential analysis scores for any kind of reactions

You then run the `differentialReactions` script from the `metnet` qiime plugin.
The parameter `--p-condition-name` represents the category of the samples related to the condition state in the metadata column.
The parameter `--p-control-name` represents the category of the samples related to the control state in the metadata column.
The parameter `--p-selection-model` corresponde to the same metabolic reconstruction used in the table generation (`AGREDA` (default), `AGORAv103`, `AGORAv201`).

```
qiime metnet differentialReactions \
	--i-reactions ../output_reactions.qza \
	--m-metadata-file ../metadata.tsv \
	--m-metadata-column columnLabel \
	--p-condition-name conditionLabel \
	--p-control-name controlLabel \
	--p-selection-model AGREDA \
	--o-differential-analysis ../reactions_differential.qza
```

## Calculate differential analysis scores for subsystems

You then run the `differentialSubSystems` script from the `metnet` qiime plugin.
The parameter `--p-condition-name` represents the category of the samples related to the condition state in the metadata column.
The parameter `--p-control-name` represents the category of the samples related to the control state in the metadata column.

```
qiime metnet differentialSubSystems \
	--i-subsystems ../output_subsystems.qza 	
	--m-metadata-file ../metadata.tsv 	
	--m-metadata-column columnLabel\
	--p-condition-name conditionLabel \
	--p-control-name controlLabel \
	--o-differential-analysis ../subsystems_differential.qza
```

## Generate the hierarchically-clustered heatmap

You then run the `plotClusteMap` script from the `metnet` qiime plugin.
The `--i-table` can be the normalized scores table for reactions or subsystems.
The input `--m-sample-metadata-file` and the specification of the `--m-sample-metadata-column` define the label groups.
It is possible to introduce parameters such as metrics, methods, color_scheme, title etc, not present below.

```
qiime metnet plotClusteMap 
	--i-table ../output_reactions.qza \
	--m-sample-metadata-file ../metadata.tsv \
	--m-sample-metadata-column columnLabel \
	--o-visualization ../clustermap.qzv
```

## Compute and visualize the Principal Component Analysis

You then run the `plotPCA` script from the `metnet` qiime plugin.
The `--i-table` can be the normalized scores table for reactions or subsystems.
The input `--m-sample-metadata-file` and the specification of the `--m-sample-metadata-column` define the label groups.
It is possible to introduce parameters such as color_scheme, title and individual point label (boolean, default = False) not present below.
```
qiime metnet plotPCA \
	--i-table ../output_reactions.qza \
	--m-sample-metadata-file ../metadata.tsv \
	--m-sample-metadata-column columnLabel \
	--o-visualization ../pca.qzv
```

## Generate the boxplots

You then run the `plotBoxplot` script from the `metnet` qiime plugin.
The `--i-table` can be the normalized scores table for reactions or subsystems.
The `--i-differentialresults` represents the differential activity analysis corresponding to exchanges, reactions or subsystems.
The input `--m-sample-metadata-file` and the specification of the `--m-sample-metadata-column` define the label groups.
The parameter `--p-condition-name` represents the category of the samples related to the condition state in the metadata column.
The parameter `--p-control-name` represents the category of the samples related to the control state in the metadata column.
The parameter `--p-namefeature` represents the name of the exchange, reaction or subsystem for which we generate the boxplot. The name must be extracted exactly from the differential activity table.
It is possible to define a title.

```
qiime metnet plotBoxplot \
	--i-table ../output_reactions.qza \
	--i-differentialresults ../exchanges_differential.qza \
	--m-sample-metadata-file ../metadata.tsv \
	--m-sample-metadata-column columnLabel \
	--p-condition-name conditionLabel \
	--p-control-name controlLabel \
	--p-namefeature 'nameExchange' \
	--o-visualization ../boxplot.qzv
```
