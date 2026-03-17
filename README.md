# RNAseq-Analysis-and-Visualization-Pipeline
This repository and the included tutorial were created to speed up initial bulk RNAseq analysis for the Jaspers lab at UNC Chapel Hill, as no unified pipeline had been established despite several researchers independently conducting RNA sequencing projects. As each researcher tended to generate their own analysis pipeline, different techniques were often employed making direct comparisons between datasets complicated. Additionally, each researcher having to create a pipeline from scratch after their first RNAseq project was slow and represented a wasted effort that could have been spent on more in depth analysis of the dataset. This pipeline seeks to accelerate that initial analysis to quickly move on to the more interesting and novel work.

## 1.	Raw Data Isolation and Cleaning

The Jaspers lab has largely used Genewiz services to conduct sequencing. This initial folder contains code for the ingestion of aligned count data as provided by Genewiz (in multiple .csv files) followed by initial QC checks including MDS and PCA plots as well as correlation dendrograms to identify potential outlier samples. A synthetic dataset was generated within the script to allow for a follow along tutorial here and within the following sections.

<p align="center">
  <img width="700" height="600" alt="image" src="https://github.com/user-attachments/assets/fb9d5f1e-1892-4cb0-a517-f085b8d7d0df" />
</p>

## 2.	DEG and IPA

This folder provides scripts for determining differentially expressed genes (DEGs) using a limma-voom focused workflow followed by data extraction and visualization after Ingenuity Pathway Analysis (IPA). 

i.	As the lab largely utilizes primary samples, our experimental designs frequently demonstrate a repeated measures design where cells from each individual donor are split between toxicologic exposures and appropriate controls. As such, accounting for sample donor within the differential expression analysis is often desired. 

Below is a variancePartition plot which demonstrates the variance explainable by each designated factor. Users are directed to specify their own desired factors to determine which factors should be included in the model during differential expression.

<p align="center">
  <img width="600" height="350" alt="image" src="https://github.com/user-attachments/assets/a55eb97d-afa4-42a8-a394-8d4d6fde5e18" />
</p>

The tutorial compares results generated using duplicateCorrelation() from limma and dream() (doi: 10.1093/bioinformatics/btaa687)  from variancePartition. Dream() is certainly a more appropriate approach as it was specifically designed for this scenario while duplicateCorrelation() is an off-label application originally designed to correct for spot duplication and technical replicates in microarrays. 

The below plot demonstrates the differences in calculated p-values between duplicateCorrelation() and dream() for the synthetic dataset. While this dataset does not demonstrate very large donor-specific effects (see above), differences are clearly visible between the two methods. The important note, though, is that the dream() workflow takes substantially longer to run than duplicateCorrelation(); however, parellelization in dream() is possible and the overall run time is not too bad on a modern computer unless sample sizes are quite large. 

<p align="center">
  <img width="800" height="250" alt="image" src="https://github.com/user-attachments/assets/1146213e-9cf0-43c5-91b3-39ff6b0763fa" />
</p>

Following differential expression, results can be visualized using volcano plots and compared between groups of interest using euler plots.

<p align="center">
  <img width="750" height="230" alt="image" src="https://github.com/user-attachments/assets/aba33c8b-1e3b-4782-96e8-ecbdd7edebb0" />
  <img width="850" height="285" alt="image" src="https://github.com/user-attachments/assets/95db0726-91c8-4a77-aa56-7f607e9cde29" />
</p>

ii.	IPA is a powerful tool which allows for the translation of differential expression results to the prediction of pathway-level activity. However, a massive amount of information can be generated from toxicology studies, requiring significant manual effort to extract, analyze, and visualize these results. This script ingests all data which is exported from an IPA core analysis and sorts data into their specific readouts (canonical pathways, master regulators, etc) and allows for preliminary visualization using a variety of techniques including bar plots and heatmaps.

<p align="center">
  <img width="1000" height="500" alt="image" src="https://github.com/user-attachments/assets/fe4a749a-596a-402f-9544-7a636a2a63d6" />
  <img width="500" height="600" alt="image" src="https://github.com/user-attachments/assets/36b5ee57-5059-44f2-885f-e3775c22273c" />
</p>

## 3.	Heatmaps

Following differential expression and IPA, it is often useful to visualize gene expression patterns and perform dimensionality reduction to look for less organized but more universal patterns. This script identifies genes which demonstrate significantly different expression levels between any treatment groups (ANOVA with BH correction) followed by visualization using heatmaps and principal component analysis (PCA). These approaches often provide a wider range of details regarding overall trends in each treatment group and have proven valuable for hypothesis generation and as an excellent starting point for deeper analysis.

<p align="center">
  <img width="900" height="500" alt="image" src="https://github.com/user-attachments/assets/0f51ab03-fb19-47d7-ad55-9e66c9fc1186" />
</p>

## 4.	Gene Ontology and Gene Set Variation Analysis

In conjunction with IPA, a more general approach to investigating broader patterns in RNAseq data are gene set enrichment approaches. Here, gene ontology (GO) enrichment followed by gene set variation analysis is conducted. GO terms are curated lists of genes linked to specific cellular and biological processes, allowing for the assessment of sets of genes which are known to be related. Differentially expressed genes are used to determine which GO gene sets are enriched before gene set variation analysis is used to test those gene sets for particular enrichment or suppression compared to other groups. Resulting GSVA values can easily be visualized using heatmaps to quickly compare pathway activity levels.
<p align="center">
  <img width="644" height="600" alt="image" src="https://github.com/user-attachments/assets/1a553f38-28bf-4dd5-9ac8-6fc181b08460" />
</p>
