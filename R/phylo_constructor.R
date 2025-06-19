#' Phyloseq constructor
#'
#' This function creates a \link[phyloseq]{phyloseq-class} object.
#'
#' @param microbiota_data (Required). A data frame containing samples in rows and features with full taxonomy in cols.
#' @param id_col (Required). String for a column in the data frame to be used as sample id.
#' @param metadata (Optional). A data frame containing metadata for the experiment. Should contain the id_col.
#'
#' @returns A \link[phyloseq]{phyloseq-class} object.
#' @export
#'
#' @examples
#' phylo_constructor(test_data, id_col = "Infants_pdp_ID")
phylo_constructor <- function(microbiota_data, id_col, metadata) {

  if (missing(microbiota_data)) {
    stop("No microbiota_data was provided")
  }
  if (missing(id_col)) {
    stop("You need to provide a id_col identifier")
  }

  if (inherits(microbiota_data, "tbl_df")) {
    microbiota_data <- as.data.frame(microbiota_data)
  }

  ft_table <- microbiota_data[, grepl("d__", names(microbiota_data))]
  id_names <- microbiota_data[[id_col]]
  rownames(ft_table) <- id_names

  tax_table <- do.call(rbind, lapply(colnames(ft_table), phyloseq::parse_taxonomy_qiime))
  colnames(tax_table)[1] <- "Domain"

  tax_names <- paste("ASV", 1:nrow(tax_table), sep = "_")
  rownames(tax_table) <- tax_names
  colnames(ft_table) <- tax_names

  if (!missing(metadata)) {
    return(phyloseq::phyloseq(feature_table(ft_table), phyloseq::tax_table(tax_table), phyloseq::sample_data(metadata)))
  }
  else {
    message("metadata not provided, creating a phyloseq object without metadata.")
    return(phyloseq::phyloseq(feature_table(ft_table), phyloseq::tax_table(tax_table)))
  }
}
