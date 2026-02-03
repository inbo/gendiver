################
# Sample sheet #
################

# Function to try to remove/standardize column names created by humans
# Remove any non-letter/number and replace with underscore
clean_colnames = function(x){
  xout = toupper(trimws(gsub("[^a-zA-Z0-9]", " ", x)))
  xout = gsub("\\s+", "_", xout)
  return(xout)
}

#' Read lab-sample-sheet
#'
#' All purpose starting point for reading the lab sample_sheet.
#'
#' @param sample_metadata_path Filepath for \link[readxl]{read_xlsx} or google-sheets id \link[googlesheets4]{range_read}
#'
#' @returns Very basic cleaned and formatted dataframe from lab-sample-sheet
#' @export
#'
#' @examples
#' # If access to INBO GDRIVE, run:
#' my_sample_data_df = read.sample_sheet_lab("1cWkVqk3y7668OVRhIlpa0fUmTh4oJTRWUqp60GRJ5J8")
read.sample_sheet_lab = function(sample_metadata_path){
  # 1. read xlsx file or GSHEET and convert to dataframe
  if (file.exists(sample_metadata_path)){
    sample_metadata = as.data.frame(readxl::read_xlsx(sample_metadata_path))
  } else {
    # try to interpret as google sheets ID-key or URL
    sample_metadata = as.data.frame(googlesheets4::read_sheet(sample_metadata_path, sheet=1))
  }

  # 2. clean the colnames, i.e. all caps and no whitespace
  colnames(sample_metadata) = clean_colnames(colnames(sample_metadata))

  # Remove any row where the USI is NA
  sample_metadata = sample_metadata[!is.na(sample_metadata$UNIQUE_SAMPLE_CODE),]
  # Or "special" NA
  sample_metadata = sample_metadata[sample_metadata$UNIQUE_SAMPLE_CODE != "#N/A",]

  # Change rownames to USI for phyloseq
  rownames(sample_metadata) = sample_metadata$UNIQUE_SAMPLE_CODE

  return(sample_metadata)
}

#' Write barcode-file for demultiplexing
#'
#' Loop over all multiplexed sample names (POOLS) and write TSV with USI - barcodes combinations
#' Check the lab-sample-sheet you use that each USI-barcode combination is UNIQUE!
#'
#' @param lab_sample_sheet Dataframe with lab-sample-data (see \link[gendiver]{read.sample_sheet_lab})
#' @param out_dir Directory path to write barcode-files to
#' @param LIB_COL Column name or ID designating the multiplexed POOLs
#'
#' @returns Nothing, prints POOL-ids and writes files to `outdir`
#' @export
#'
#' @examples
#' my_sample_data_df = read.sample_sheet_lab("1cWkVqk3y7668OVRhIlpa0fUmTh4oJTRWUqp60GRJ5J8")
#' write.barcode_files(my_sample_data_df, out_dir=tempdir())
write.barcode_files = function(lab_sample_sheet, out_dir, LIB_COL=1){

  target_col = lab_sample_sheet[,c(LIB_COL)]

  for (lib_code_i in unique(target_col)){
    print(lib_code_i)

    lib_data = lab_sample_sheet[target_col == lib_code_i , ]
    outtable = lib_data[, c("UNIQUE_SAMPLE_CODE",
                            "SEQUENCE_BARCODE_F_PRIMER",
                            "SEQUENCE_BARCODE_R_PRIMER")]

    utils::write.table(outtable, file = file.path(out_dir, paste0(lib_code_i, ".tsv")),
                sep="\t", quote = F, col.names = F, row.names = F, )

  }
}
