% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gene_ontology_GOdatabase_functions.R,
%   R/Summary_stats_and_pairwise_dge.R
\name{Limma_dge_prep}
\alias{Limma_dge_prep}
\title{Pairwise differential gene expression with Limma}
\usage{
Limma_dge_prep(
  normalized_counts,
  metadata_table,
  testing_column,
  column1,
  column2
)

Limma_dge_prep(
  normalized_counts,
  metadata_table,
  testing_column,
  column1,
  column2
)
}
\arguments{
\item{normalized_counts}{Dataframe of your normalized counts with sample names at the top and genes as the index.}

\item{metadata_table}{Metadata table associated with your dataframe}

\item{testing_column}{Column you are conducting the pairwise comparison in}

\item{column1}{One of the categorical group names you are comparing 'Male'}

\item{column2}{Another of the categorical group names you are comparing 'Female'}
}
\value{
A data matrix containing the differential gene expression between two groups. Returns measurements logFC, Average expression, test statistic, P. Value, adjusted P. Value, and B statistic

A data matrix containing the differential gene expression between two groups. Returns measurements logFC, Average expression, test statistic, P. Value, adjusted P. Value, and B statistic
}
\description{
This function runs pairwise differential expression using the Limma package between the two groups you are comparing in a column

This function runs pairwise differential expression using the Limma package between the two groups you are comparing in a column
}
\examples{
# Limma_dge_prep(normalized_counts_dataframe, metadata_table, 'Sex", 'M', 'F')

# Limma_dge_prep(normalized_counts_dataframe, metadata_table, 'Sex", 'M', 'F')

}
