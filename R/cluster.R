##############
# Clustering #
##############

# Read VSEARCH OTU-tab output into dataframe
read_otu_table = function(otu_table_path){
  otu_df <- read.csv(otu_table_path, sep='\t', header = T)
  colnames(otu_df) = gsub(pattern = "^X\\.*", replacement = "", colnames(otu_df))
  rownames(otu_df) = otu_df$OTU.ID
  return(otu_df[-1])
}

# Read OBITools3 MERGED SAMPLE OTU-tab output into phyloseq
read_merged_otu_table = function(otu_table_path){
  # Read table
  otu_df <- read.csv(otu_table_path, sep='\t', header = T)
  # Make OTU ID the rownames
  rownames(otu_df) = otu_df$ID
  # Extract only the ASV count data
  sample_slice = grepl("^MERGED_sample\\.*", colnames(otu_df))
  otu_df = otu_df[, sample_slice]
  # Clean up sample names
  colnames(otu_df) = gsub(pattern = "^MERGED_sample\\.*", replacement = "", colnames(otu_df))
  # return in phyloseq::otu_table() format
  return(phyloseq::otu_table(otu_df, taxa_are_rows = T))
}

# read mumu log (mainly give colnames)
read_mumu_log = function(mumu_log_path){
  mumu_log = read.table(mumu_log_path)
  colnames(mumu_log) = c("query_otu", "parent_otu", "percent_similarity", "total_abundance_q", "total_abundance_p",
                         "overlap_abundance_q", "overlap_abundance_p", "incidence_q", "incidence_p", "incidence_both",
                         "min_abundance_ratio", "sum_abundance_ratio", "mean_abundance_ratio",
                         "min_abundance_ratio_nonzero", "mean_abundance_ratio_nonzero", "max_abundance_ratio",
                         "relative_cooccurence", "status")
  return(mumu_log)
}
