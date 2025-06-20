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

parse_taxonomy_gtdb <- function (char.vec)
{
  char.vec = phyloseq::parse_taxonomy_default(char.vec)
  Tranks = c(d = "Domain", p = "Phylum", c = "Class", o = "Order",
             f = "Family", g = "Genus", s = "Species")
  ti = grep("[[:alpha:]]{1}\\_\\_", char.vec)
  if (length(ti) == 0L) {
    warning("No GTDB prefixes were found. \n", "Consider using parse_taxonomy_default() instead if true for all OTUs. \n",
            "Dummy ranks may be included among taxonomic ranks now.")
    taxvec = char.vec
  }
  else {
    taxvec = gsub("[[:alpha:]]{1}\\_\\_", "", char.vec)
    repranks = Tranks[substr(char.vec[ti], 1, 1)]
    names(taxvec)[ti[!is.na(repranks)]] = repranks[!is.na(repranks)]
  }
  return(taxvec)
}


parse_taxonomy_helper <- function(char.vec) {
  tax_object <- parse_taxonomy_gtdb(strsplit(char.vec, ";", TRUE)[[1]])
  last_rank <- NULL

  for(rank in tax_object){
    if (rank == "") {
      tax_object[length(tax_object)] <- paste(last_rank, "sp.", sep = "_")
      break
    }
    last_rank <- rank
  }
  return(tax_object)
}
