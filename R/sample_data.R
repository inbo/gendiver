################
# Sample sheet #
################

# Function to try to remove/standardize column names created by humans
# remove any non-letter/number and replace with underscore
clean_colnames = function(x){
  xout = toupper(trimws(gsub("[^a-zA-Z0-9]", " ", x)))
  xout = gsub("\\s+", "_", xout)
  return(xout)
}

# All purpose starting point for reading the lab sample_sheet
read_lab_sample_sheet_basic = function(sample_metadata_path){
  # 1. read xlsx file and convert to dataframe
  if (file.exists(sample_metadata_path)){
    sample_metadata = as.data.frame(readxl::read_xlsx(sample_metadata_path))
  } else {
    # try to interpret as google sheets ID-key
    # need to have run
    ## library(googlesheets4)
    ## googlesheets4::gs4_auth(USER_NAME)
    sample_metadata = as.data.frame(googlesheets4::read_sheet(sample_metadata_path, sheet=1))
  }

  # 2. clean the colnames, i.e. all caps and no whitespace
  colnames(sample_metadata) = clean_colnames(colnames(sample_metadata))

  # Remove nay row where the USI is NA
  sample_metadata = sample_metadata[!is.na(sample_metadata$UNIQUE_SAMPLE_CODE),]
  # Or "special" NA
  sample_metadata = sample_metadata[sample_metadata$UNIQUE_SAMPLE_CODE != "#N/A",]

  # Change rownames to USI for phyloseq
  rownames(sample_metadata) = sample_metadata$UNIQUE_SAMPLE_CODE

  return(sample_metadata)
}

# Function to loop over all multiplexed sample names and write TSV with USI - barcodes
write_barcode_files = function(lab_sample_sheet, out_dir, LIB_COL=22){

  target_col = lab_sample_sheet[,c(LIB_COL)]

  for (lib_code_i in unique(target_col)){
    print(lib_code_i)

    lib_data = lab_sample_sheet[target_col == lib_code_i , ]
    outtable = lib_data[, c("UNIQUE_SAMPLE_CODE",
                            "SEQUENCE_BARCODE_F_PRIMER",
                            "SEQUENCE_BARCODE_R_PRIMER")]

    write.table(outtable, file = file.path(out_dir, paste0(lib_code_i, ".tsv")),
                sep="\t", quote = F, col.names = F, row.names = F, )

  }
}
