% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visual_plotting_functions.R
\name{Single_gene_boxplot}
\alias{Single_gene_boxplot}
\title{Single gene boxplot}
\usage{
Single_gene_boxplot(
  counts,
  metadata_table,
  column,
  gene,
  sample_names = "Sample Name"
)
}
\arguments{
\item{counts}{Dataframe of normalized counts}

\item{metadata_table}{Metadata table associated with your dataframe}

\item{column}{Metadata table column you are interested in age, sex, etc.}

\item{gene}{Gene you want to visualize the expression of ex: 'CD4'}

\item{sample_names}{column that the sample names are under in the metadata.}
}
\value{
a nice boxplot comparing groups
}
\description{
This function boxplots a column of your choice and a gene of your choice and allows you to compare them.
}
\examples{
# Single_gene_boxplot(normalized_counts, metadata_tablee , 'Ancestry', 'CCL4', "Sample Name")
}
