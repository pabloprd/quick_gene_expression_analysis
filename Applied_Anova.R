#ANOVA FUNCTION for differential gene expression.
#Takes in the counts, the metadata, table, and the column of interest.
#' @export
applied_anova <- function(counts, table, column){
  #ALL THIS IS FOR POTENTIAL LIMMA STUFF DOWN THE LINE
  #################################################################
  row.names(table) <- table$Sample_names #Making the index of the metadata equal the sample names instead of the genes
  print('Is there a gene index and a gene column??')
  table$Sample_names <- NULL #making sure there is not a gene index AND a gene column
  truth <- all(colnames(counts) %in% rownames(metadata)) #
  table <- table[colnames(counts),, drop = FALSE] #the column names of the counts
  group <- factor(table[[paste0(column)]]) #What group each sample belongs too
  design <- model.matrix(~0 + group) #Each column corresponds to
  colnames(design) <- levels(group) #turning groupHD to HD etc.
  truth <- all(colnames(counts) %in% rownames(metadata))
  print(truth) #If it prints false it means the tables are not aligned.
  ##################################################################

  #GETTING ANOVA P
  get_anova_p <- function(x, group) {
    fit <- aov(x ~ group)
    summary(fit)[[1]][["Pr(>F)"]][1]
  }

  #GETTING ANOVA F
  get_anova_F <- function(x, group) {
    fit <- aov(x ~ group)
    summary(fit)[[1]][["F value"]][1]
  }


  anova_pvals <- apply(counts, 1, get_anova_p, group=group) #Getting the p vals of each row
  anova_Fstats <- apply(counts, 1, get_anova_F, group=group) #Getting the F statistics of each row.

  #Total means. Divide by group sizes.
  #COME BACK AND ANNOTATE CODE!
  group_levels <- levels(group) #get the unique group levels
  gene_means <- sapply(group_levels, function(g){
    rowMeans(counts[, group == g, drop = FALSE])
  })

  #For each group g in group levels the function selects the columns of counts where group == g
  #Then it calcualtes the mean expression for each gene across those samples using rowMeans
  #It applies this function to every group which returns a matrix where rows are the genes
  #The columns are groups, and the mean expression of each gene in each group


  rownames(gene_means) <- rownames(counts)

  anova_results <- data.frame(
    means = gene_means,
    p_value = anova_pvals,
    F_statistic = anova_Fstats,
    adj_p_value = p.adjust(anova_pvals, method="BH")
  )
  return(anova_results)
}
