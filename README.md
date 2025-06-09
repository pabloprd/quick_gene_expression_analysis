# R35 Moore Nanostring Analysis Functions

This is a simple and specific R package that is created for ease of use for pairwise and 
multiple differential pairwise gene expression. It streamlines a lot of code into functions
so your R code is neater and easier to read. This package specifically only uses limma for pairwise
and multiple differential pairwise gene expression. There are also some plotting functions
that include a boxplot, barplots, and a good volcano plot. 

## Workflow
- Comparing columns with only 2 variables (Pairwise comparison):
  sumstats_and_dge(normalized_counts, metadata, "Feature column", "group1", "group2", "Character the sample starts with")
  ex: sumstats_and_dge(normalized_counts, mmetadata, 'Sex', 'Male', 'Female', 'NA')

- Comparing columns with >2 variables (Multivariate comparison or ANOVA)
  Limma_mdge_results(normalized_counts, metadata, column of interest)
  ex: Limma_mdge_results(normalized_counts, mmetadata, 'Total_Discrimination')

Both of these output a table with the fold change and the adjusted p value which should give you an idea on any significant differences.
To investigate further with a multivariate comparison do a pairwise comparison between the features with the highest log fold change.



## Input Data

I will describe the input data here. You will only be working with RCC files and an associated metadata.csv file  
- You will need a metadata file with the sample names along with the associated features
- You will need to read in the counts from the RCC file so they look like this:
 
<br> 

![Screenshot of the application](./Raw_counts_table.png)


## Results

Example output figures:

Boxplot gene expression

<br>

![Screenshot of the application](./CCL4_expression.png)

Volcano plot gene expression

<br>

![Screenshot of the application](./Ancestry_gene_expression.png)

## Acknowledgments
This package was created for the Moore Lab at University of Maryland for quicker R analysis and visualization.
<br>

Created by Paul Parodi

<br>

Last updated: 6/9/2025

