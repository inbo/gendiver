#########################
### Utility functions ###
#########################

extract_sample_name = function(x){
  # note the loose regex (*) for third digit -> extraction date can be missing
  # unique_sample_code_regex = "\\d+_\\d+_\\d*_\\d+_\\d+_\\d+_\\d+_\\d+"
  # unique_sample_code_regex = "\\d*_*\\d*_*\\d*_*\\d*_*\\d*_*\\d*_*\\d*"
  # read_count_track.df$sample_id = stringr::str_extract(read_count_track.df$file_name, unique_sample_code_regex)
  y = gsub("\\..*", "",x)

  return(y)

}

# Look at the per-sample files generated during pre-processing and extract process info
extract_filename_info = function(file_name, sample_id_regex){
  info_part = stringr::str_remove_all(file_name, "\\.fast(q|a)(\\.gz)?$|^_|^\\.|^-|merged_|\\.log" )
  # info_part = sub(pattern = sample_id_regex, replacement = "", info_part)
  info_part = paste0(strsplit(info_part, split = "\\.")[[1]][-1],collapse = "." )

  return(info_part)
}

#' Try to read a table from File Path or google Drive
#'
#' At INBO Google Drive is the main platform for RDM. This functions can read CSV tables from Google Drive link
#' (right-click file, "Share", "Copy link")
#'
#' @param file_id Path or URL (Google Drive) to a tabular file that can be read with \link[utils]{read.table}
#' @param ... See \link[utils]{read.table}
#'
#' @returns Dataframe
#' @export
#'
#' @examples
#' #To add
read.table_gdrive =  function(file_id, ...){
  if (file.exists(file_id)){
    message('Reading from File Path')
    out = utils::read.csv(
      file=file_id, ...)
  } else {
    message('Reading from Google Drive')
    out = utils::read.csv(
      text = googledrive::drive_read_string(file_id),
      ...)
  }
  return(out)
}


#' Try to read a FASTX file from File Path or google Drive
#'
#' At INBO Google Drive is the main platform for RDM. This function reads FASTX files into \link[Biostrings]{DNAStringSet} from Google Drive link
#' (right-click file, "Share", "Copy link")
#'
#' @param file_id Path or URL (Google Drive) to a FASTX file that can be read with \link[Biostrings]{readDNAStringSet}
#' @param tempf File path for downloaded data (default=\link[base]{tempfile})
#' @param ... See \link[Biostrings]{readDNAStringSet}
#'
#' @returns \link[Biostrings]{DNAStringSet}
#' @export
#'
#' @examples
#' #To add
read.fastx_gdrive =  function(file_id, tempf=tempfile() , ...){
  if (file.exists(file_id)){
    message('Reading from File Path')
    f_path = file_id
  } else {
    message('Reading from Google Drive')
    f_path= tempf
    googledrive::drive_download(
      file= file_id,
      path=f_path)

  }
  out = Biostrings::readDNAStringSet(f_path, ...)
  return(out)
}


