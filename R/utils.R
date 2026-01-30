#########################
### Utility functions ###
#########################
# not for pkg-users

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
