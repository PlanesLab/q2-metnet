# q2-metnet

QIIME 2 plugin to contextualize taxonomic abundance data into high-quality metabolic reconstruction of the human gut microbiota.
Read more about the method in our [paper](https://academic.oup.com/bioinformatics/article/40/11/btae455/7715876).

# Installing

## Creating a conda environment with QIIME2 and q2-metnet

```
conda env create \
 -n q2-metnet \
 -f https://raw.githubusercontent.com/planesLab/q2-metnet.git/main/environment-files/q2-metnet-qiime2-amplicon-2024.10.yml
```

## Using an existing QIIME2 environment

Once you create your QIIME2 environment, you can install this plugin by cloning this repo and installing manually. To clone the repo:

```
git clone https://github.com/PlanesLab/q2-metnet.git
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
The parameter `--p-selection` can be chosen between `AGREDA` (default) and `AGORAv103`.
The parameter `--p-level` correspond to the taxonomic depth, ranging from kingdom (`k`) to species (`s`, default).

```
qiime metnet generateFeatures \
	--i-frequency ../asv_table.asv.qza \
	--i-taxa ../assigned_taxonomy.qza \
	--p-selection AGREDA \
	--p-level s \
	--o-reactions ../output_reactions.qza \
	--o-subsystems ../output_subsystems.qza \
	--o-xmatrix ./output_X_matrix.qza
```

## Calculate differential analysis scores for exchange reactions

You then run the `differentialExchanges` script from the `metnet` qiime plugin.
The parameter `--p-condition-name` represents the category of the samples related to the condition state in the metadata column.
The parameter `--p-control-name` represents the category of the samples related to the control state in the metadata column.
The parameter `--p-selection-model` corresponds to the same metabolic reconstruction used in the table generation (`AGREDA` (default), `AGORAv103`).
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
The parameter `--p-selection-model` corresponde to the same metabolic reconstruction used in the table generation (`AGREDA` (default), `AGORAv103`).

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

## Test files and application

The **test** folder contains all the input data sources to replicate the manuscript results. To run the following tutorial, we will need the following files:
+ **asv_table.qza:** feature table with microbial abundances, where rows correspond to each sequence and columns to the different samples. The column **ID** stores the identificators of each sequence (*e.g.* ASV_1).
+ **taxa.qza:** contains the taxonomic lineage of each sequence. The column **ID** stores the identificators of each sequence, which are written in the same order as in the feature table rows.
+ **meta.tsv:** file including sample phenotype data. The rows correspond to each sample and are written in the same order as the feature table columns.

In order to use your own data, please refer to [QIIME2 documentation](https://docs.qiime2.org/2023.7/semantic-types/#common-semantic-types) to check allowed column headers. In any case, you can use the same columnn headers as in this tutorial files. Once QIIME2 and q2-metnet have been successfully installed, the first step is to generate the reactions and subsystem normalized activity scores of our data. This step may take a couple of minutes to be performed.

```
qiime metnet generateFeatures \
	--i-frequency ./test/asv_table.qza \
	--i-taxa ./test/taxa.qza \
	--p-selection AGREDA \
	--p-level s \
	--o-reactions ./rxns_scores.qza \
	--o-subsystems ./subs_scores.qza \
	--o-xmatrix ./Xmatrix.qza
```

By the selection of parameters `--p-selection` and `--p-level` we decided to map our different taxonomies to the AGREDA database, until the species level. At this point, the normalized activity scores can be applied to perform a PCA and classify samples. As shown in the **meta.tsv** file, the different clinical conditions are stored in the **Condition** column. Then, we can generate a PCA of our samples and all the reactions scores by the following command. Below is displayed the correspondent figure. Remember that any **qzv** file can be displayed at [QIIME2 View](https://view.qiime2.org/).

```
qiime metnet plotPCA \
	--i-table ./rxns_scores.qza \
	--m-sample-metadata-file ./test/meta.tsv \
	--m-sample-metadata-column Condition \
	--o-visualization ./rxns_pca.qzv
```

<img src="https://raw.githubusercontent.com/PlanesLab/q2-metnet/refs/heads/main/Figure/pca.png">

Additionally, a hierarchically-clustered heatmap can be generated by running the `plotClusteMap`function. In the following example, we are going to generate the heatmap of the reaction scores across samples. Following, the generated figure is shown.

```
qiime metnet plotClusteMap \
	--i-table ./rxns_scores.qza \
	--m-sample-metadata-file ./test/meta.tsv \
	--m-sample-metadata-column Condition \
	--o-visualization ./rxns_hclust.qzv
```

<img src="https://raw.githubusercontent.com/PlanesLab/q2-metnet/refs/heads/main/Figure/clustermap.png">

We are able to show that the normalized activity scores for some of these reactions are altered through the different clinical conditions of our dataset. We are now going to investigate deeper these differences. The different clinical conditions of our data are **Allergic**, **Celiac**, **Lean** and **Obese**, where **Lean** can be considered as the control state. In this example, we are now going to check the differentially expressed reactions and subsystems between **Celiac** and **Lean** samples.

```
qiime metnet differentialReactions \
	--i-reactions ./rxns_scores.qza \
	--m-metadata-file ./test/meta.tsv \
	--m-metadata-column Condition \
	--p-condition-name Celiac \
	--p-control-name Lean \
	--p-selection-model AGREDA \
	--o-differential-analysis ./diff_reactions.qza
```
```
qiime metnet differentialSubSystems \
	--i-subsystems ./subs_scores.qza \
	--m-metadata-file ./test/meta.tsv \
	--m-metadata-column Condition \
	--p-condition-name Celiac \
	--p-control-name Lean \
	--o-differential-analysis ./diff_subs.qza
```
```
qiime metnet differentialExchanges \
	--i-reactions ./rxns_scores.qza \
 	--m-metadata-file ./test/meta.tsv \
  	--m-metadata-column Condition \
   	--p-condition-name Celiac \
    	--p-control-name Lean \
       	--o-differential-analysis ./diff_ex.qza
```

Finally, the differential activity score of a given feature (reaction or subsystem) can be displayed by means of the `plotBoxplot`function. In the following example, we are going to generate the boxplot of the **Bile acid metabolism** subsytem across samples of the **Celiac** and **Lean** condition.
```
qiime metnet plotBoxplot \
	--i-table ./subs_scores.qza \
	--i-differentialresults ./diff_subs.qza \
	--m-sample-metadata-file ./test/meta.tsv \
	--m-sample-metadata-column Condition \
	--p-condition-name Celiac \
	--p-control-name Lean \
	--p-namefeature "S144 | Bile acid metabolism" \
	--o-visualization ./sub_boxplot.qzv
```
<img src="https://raw.githubusercontent.com/PlanesLab/q2-metnet/refs/heads/main/Figure/boxplot.PNG">

Remember that you can run the following command and visualize the output file [here](https://view.qiime2.org/) to select the different values for '--p-name-feature'.
```
qiime metadata tabulate \
    --m-input-file ./diff_subs.qza \
    --o-visualization ./diff_subs.qzv
```