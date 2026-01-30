##################
# Demultiplexing #
##################

#' Read per-sample readcounts from summary files after demultiplexing job
#'
#' Find and read auto-generated demultiplexing results from demultiplexing folder
#'
#' @param data_dir File path to the `02_demultiplex/` dir
#' @param summary_file_pattern Pattern to recognize summary files (default="sample.txt)
#'
#' @returns Dataframe with read counts per sample
#' @export
#'
#' @examples
#' #To add
read.demultiplex_summary_files = function(data_dir, summary_file_pattern="sample.txt"){
  # Initialize readcount tracking df
  read_count_track.df = data.frame("file_name"=character(),
                                   "read_count"=numeric(),
                                   "multiplexed_source"=character())

  ## Read for each DEMULTIPLEXED lib the sample counts (from count seqs PBS)
  demultiplex_res_lib_paths = list.dirs(path = file.path(data_dir), full.names = TRUE, recursive = F)

  for (mylib in demultiplex_res_lib_paths){
    readcounts_file = list.files(mylib, pattern = summary_file_pattern, full.names = T)
    if (length(readcounts_file) > 0) {
      count_tab = utils::read.table(readcounts_file)
      colnames(count_tab) = c("file_name", 'read_count')
      count_tab$multiplexed_source = basename(mylib)
      read_count_track.df=rbind(read_count_track.df, count_tab)
    }

  }

  # note the loose regex (*) for third digit -> extraction date can be missing
  # unique_sample_code_regex = "\\d+_\\d+_\\d*_\\d+_\\d+_\\d+_\\d+_\\d+"
  # unique_sample_code_regex = "\\d*_*\\d*_*\\d*_*\\d*_*\\d*_*\\d*_*\\d*"
  # read_count_track.df$sample_id = stringr::str_extract(read_count_track.df$file_name, unique_sample_code_regex)
  read_count_track.df$sample_id = gsub("\\..*", "",read_count_track.df$file_name)
  read_count_track.df$sample_id[is.na(read_count_track.df$sample_id)] = "no_sample"

  read_count_track.df$type = sapply(read_count_track.df$file_name, extract_filename_info, "\\..*")

  read_count_track.df$type[grepl("_F$|_R1$", read_count_track.df$sample_id)] = "demultiplexed_R1"
  read_count_track.df$type[grepl("_R$|_R2$", read_count_track.df$sample_id)] = "demultiplexed_R2"
  read_count_track.df$type[grepl("unknown", read_count_track.df$sample_id)] = paste0("unknown-", read_count_track.df$type[grepl("unknown", read_count_track.df$sample_id)])
  read_count_track.df$sample_id = gsub(pattern = "_F$|_R1$|_R$|_R2$", "", read_count_track.df$sample_id)

  return(read_count_track.df)
}
