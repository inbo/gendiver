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


#' Get sample sheets URLs from LST002
#'
#' Read sample sheet URLs from \href{https://docs.google.com/spreadsheets/d/1h70fKsQCqPQS5tSQYOwws2zeCJJu9BOlH7mI0xUg5rA/edit?gid=1637141996#gid=1637141996}{LST002}
#'
#' @returns List with URLs
#' @export
#'
#' @examples
#' # todo
get.sample_sheet_urls = function(){
  LST002 = googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1h70fKsQCqPQS5tSQYOwws2zeCJJu9BOlH7mI0xUg5rA/edit?gid=1637141996#gid=1637141996")

  out = list()
  out[as.character(LST002$`Run name`)] = LST002$sample_gsheet_url

  # custom additions
  out$`25004` = list(
    "insekp" = "https://docs.google.com/spreadsheets/d/1SxZwDzzA7mcSlpO5YDMXcZKpUCtxO_c56bNyGxktwE0/edit?usp=drive_link",
    "riaz_teleo" = "https://docs.google.com/spreadsheets/d/1Zlge3srOWzk4jnR5XCJjbOkrmIPPrhMA8y2jsq0-Vfc/edit?usp=drive_link")

  # custom additions
  out$`23005` = list(
    "olig_coll" = "https://docs.google.com/spreadsheets/d/12cIhbTq-NDT4WFMYuqcAeTxJr5RgC5A5AX7TxYLysxU/edit?gid=2119359410#gid=2119359410",
    "riaz" = "https://docs.google.com/spreadsheets/d/1hcgX9Iz4IpA2oCQtez6zby03k9CLDclZNjJknr6b-sE/edit?usp=drive_link")


    return(out)
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
#' # example on public available data
#' # googlesheets4::gs4_deauth()
#' # my_sample_data_df = read.sample_sheet_lab("1cWkVqk3y7668OVRhIlpa0fUmTh4oJTRWUqp60GRJ5J8")
read.sample_sheet_lab = function(sample_metadata_path){
  # 1. read xlsx file or GSHEET and convert to dataframe
  if (file.exists(sample_metadata_path)){
    sample_metadata = as.data.frame(readxl::read_xlsx(sample_metadata_path))
  } else {
    # try to interpret as google sheets ID-key or URL
    # Note the columns are all read as character due to inconsistent col types
    sample_metadata = as.data.frame(googlesheets4::read_sheet(sample_metadata_path,
                                                              sheet=1,
                                                              col_types = 'c'))
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


#' Convert lab-sample-sheet from pre 2024 analyses
#'
#' Change column names from "old" sample sheets to "new" format.
#'
#' @param sample_sheet_df Sample sheet data frame as per \link[gendiver]{read.sample_sheet_lab}
#'
#' @returns Sample sheet data frame with converted colnames/data
#' @export
#'
#' @examples
#' # example to be added
convert.sample_sheet_pre2024 = function(sample_sheet_df){
  ## Convert "Annelies" lab-sample-sheet format to new INBO format
  # Assume sample sheets is the same format always.

  colnames(sample_sheet_df) = c(
     "MATRIX", "MATRIX_CODE", "PROJECT", NA,
     "PROJECT_CODE", "EXTRACTION_DATE_YYMMDD", "EXTRACTION_DATE_CODE", "ANALYSE_TYPE", "ANALYSE_TYPE_CODE", "INDIVIDUAL", NA,
     "INDIVIDUAL_CODE", "LIBRARY_TYPE", "LIBRARY_TYPE_REMARKS", NA,
     "PROTOCOL_CODE", "TECHNICAL_REPLICATE" , "TECHNICAL_REPLICATE_CODE", "RUN", "SEQUENCING_PLATFORM", "RUN_CODE", "REMARKS", "UNIQUE_SAMPLE_CODE",
     rep(NA, 8),  "BARCODE_F_PRIMER", "SEQUENCE_BARCODE_F_PRIMER", "BARCODE_R_PRIMER", "SEQUENCE_BARCODE_R_PRIMER",
     rep(NA, 7), "PLATE_NOTES", "WELL",  "POOL_NAME_SEQ"
  )

  sample_sheet_df$ANALYSE_TYPE[sample_sheet_df$ANALYSE_TYPE_CODE == "999"] = "Metabarcoding"
  sample_sheet_df$ANALYSE_TYPE[sample_sheet_df$ANALYSE_TYPE_CODE == "000"] = "NTC"
  sample_sheet_df$ANALYSE_TYPE[sample_sheet_df$ANALYSE_TYPE_CODE == "111"] = "PAC"



  sample_sheet_df = sample_sheet_df[!is.na(colnames(sample_sheet_df))]
  return(sample_sheet_df)


}


#' Write barcode-file for demultiplexing
#'
#' Loop over all multiplexed sample names (POOLS) and write TSV with USI - barcodes combinations
#' Check the lab-sample-sheet you use that each USI-barcode combination is UNIQUE!
#'
#' @param lab_sample_sheet Dataframe with lab-sample-data (see \link[gendiver]{read.sample_sheet_lab})
#' @param out_dir Directory path to write barcode-files to
#' @param LIB_COL Column name or ID designating the multiplexed POOLs
#' @param nextits Write barcode files to fasta format as expected by NextITS (cutadapt)
#'
#' @returns Nothing, prints POOL-ids and writes files to `outdir`
#' @export
#'
#' @examples
#' # example on public available data
#' # googlesheets4::gs4_deauth()
#' # my_sample_data_df = read.sample_sheet_lab("1cWkVqk3y7668OVRhIlpa0fUmTh4oJTRWUqp60GRJ5J8")
#' # write.barcode_files(my_sample_data_df, out_dir=tempdir())
write.barcode_files = function(lab_sample_sheet, out_dir, LIB_COL=1, nextits=FALSE){

  target_col = lab_sample_sheet[,c(LIB_COL)]

  for (lib_code_i in unique(target_col)){

    lib_data = lab_sample_sheet[target_col == lib_code_i , ]

    cat(paste0("\n",lib_code_i, ": ", nrow(lib_data), " samples"))

    if (nextits){
      # In NextITS, we use the SequencingRunID__SampleID naming convention (please note the double underscore separating RunID and SampleID parts)
      # https://next-its.github.io/usage/#barcodes
      out_data = paste0(
        ">Run", lib_data$RUN_CODE,"", lib_data$LIBRARY_CODE, "__", lib_data$UNIQUE_SAMPLE_CODE, '\n',
        lib_data$SEQUENCE_BARCODE_F_PRIMER, "...", lib_data$SEQUENCE_BARCODE_R_PRIMER )
      fileConn<-file(file.path(out_dir, paste0(lib_code_i, ".fasta")))
      writeLines(out_data, fileConn)
      close(fileConn)

    } else {
      outtable = lib_data[, c("UNIQUE_SAMPLE_CODE",
                              "SEQUENCE_BARCODE_F_PRIMER",
                              "SEQUENCE_BARCODE_R_PRIMER")]
      utils::write.table(outtable, file = file.path(out_dir, paste0(lib_code_i, ".tsv")),
                         sep="\t", quote = F, col.names = F, row.names = F, )
    }


  }
}

