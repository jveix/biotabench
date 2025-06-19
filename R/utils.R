#' Wrapper function for \link[phyloseq]{otu_table}.
#'
#' @param object a feature table containing features at columns and samples at rows
#'
#' @returns a \link[phyloseq]{phyloseq-class} feature table
#' @export
#'
#' @examples
#'data(GlobalPatterns)
#'feature_table(test_phylo_object)
feature_table <- function(object) {
  return(phyloseq::otu_table(object, taxa_are_rows = F))
}
