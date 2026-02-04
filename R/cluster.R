##############
# Clustering #
##############

#' Read VSEARCH OTU-tab output into dataframe
#'
#' Read OTU table, cleanup colnames and return
#'
#' @param otu_table_path Path to OTU table file (can also be link to Googe Drive)
#'
#' @returns Dataframe with read counts per (z)OTU (rows) per sample (columns)
#' @export
#'
#' @examples
#' #To add
read.otu_table_vsearch = function(otu_table_path){
  otu_df <- read.table_gdrive(otu_table_path, sep='\t', header = T)
  colnames(otu_df) = gsub(pattern = "^X\\.*", replacement = "", colnames(otu_df))
  rownames(otu_df) = otu_df[,1]
  return(otu_df[-1])
}

#' Read OBITools3 MERGED SAMPLE OTU-tab output into dataframe
#'
#' Read OTU table, cleanup colnames and return
#'
#' @param otu_table_path Path to OTU table file (can also be link to Googe Drive)
#'
#' @returns Dataframe with read counts per ASV (rows) per sample (columns)
#' @export
#'
#' @examples
#' #To add
read.otu_table_obitools3 = function(otu_table_path){
  # Read table
  otu_df <- read.table_gdrive(otu_table_path, sep='\t', header = T)
  # Make OTU ID the rownames
  rownames(otu_df) = otu_df$ID
  # Extract only the ASV count data
  sample_slice = grepl("^MERGED_sample\\.*", colnames(otu_df))
  otu_df = otu_df[, sample_slice]
  # Clean up sample names
  colnames(otu_df) = gsub(pattern = "^MERGED_sample\\.*", replacement = "", colnames(otu_df))
  # return in phyloseq::otu_table() format
  # return(phyloseq::otu_table(otu_df, taxa_are_rows = T))
  return(otu_df)
}

#' Read MUMU log and return in dataframe
#'
#' Read mumu log, add column names based on documentation.
#' https://github.com/frederic-mahe/mumu/blob/09a510bca992562405f1904debfbd10ae8d7bdef/man/mumu.1#L124
#'
#' > Output file for OTU merging statistics (18 columns separated by tabulations, first line is a header line with column names)
#'
#' - This is incorrect, no header is present in the file. -> Hence this function.
#'
#' @param mumu_log_path Path to MUMU logfile (can also be link to Googe Drive)
#'
#' @returns Dataframe with MUMU log
#' @export
#'
#' @examples
#' #To add
read.log_mumu = function(mumu_log_path){
  mumu_log = read.table_gdrive(mumu_log_path, header = F)
  colnames(mumu_log) = c("query_otu", "parent_otu", "percent_similarity", "total_abundance_q", "total_abundance_p",
                         "overlap_abundance_q", "overlap_abundance_p", "incidence_q", "incidence_p", "incidence_both",
                         "min_abundance_ratio", "sum_abundance_ratio", "mean_abundance_ratio",
                         "min_abundance_ratio_nonzero", "mean_abundance_ratio_nonzero", "max_abundance_ratio",
                         "relative_cooccurence", "status")
  return(mumu_log)
}
