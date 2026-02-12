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
read.summary_files_demultiplex = function(data_dir, summary_file_pattern="sample.txt"){
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

  # Add sample ID
  read_count_track.df$sample_id = extract_sample_name(read_count_track.df$file_name)
  read_count_track.df$sample_id[is.na(read_count_track.df$sample_id)] = "no_sample"
  read_count_track.df$sample_id = gsub(pattern = "_F$|_R1$|_R$|_R2$", "", read_count_track.df$sample_id)

  # Add type
  read_count_track.df$type = sapply(read_count_track.df$file_name, extract_filename_info, "\\..*")
  read_count_track.df$type[grepl("_F\\.|_R1\\.", read_count_track.df$file_name)] = "demultiplexed_R1"
  read_count_track.df$type[grepl("_R\\.|_R2\\.", read_count_track.df$file_name)] = "demultiplexed_R2"
  read_count_track.df$type[grepl("unknown", read_count_track.df$sample_id)] = paste0("unknown-", read_count_track.df$type[grepl("unknown", read_count_track.df$sample_id)])
  read_count_track.df$type = factor(df$type, levels = rev(unique(read_count_track.df$type)))

  return(read_count_track.df)
}


#' Summarize demultiplexing readcounts for plotting
#'
#'  demultiplexed readcounts per sample
#'
#' @param df Data frame with read counts per sample (output of \link[gendiver]{read.summary_files_demultiplex})
#'
#' @returns Summary data ready for plotting
#' @export
#'
#' @examples
#' #To add
summary.summary_files_demultiplex = function(df){

  # Convert type to factor (for grouping ggplot)
  df$type = factor(df$type, levels = rev(unique(df$type)))

  #1. Overview of all reads
  demultiplex_summary_df = aggregate(read_count ~ type, data = df, sum)

  # only keep R1 info
  demultiplex_summary_df = demultiplex_summary_df[grepl("1$|FWD$", x = demultiplex_summary_df$type),]

  # add total raw reads:
  # because demultiplexing is done in 2 steps, unk+res is 2x raw
  # For cutadapt, unknown-demultiplexed is not just the discarded reads, iIts all discarded *combined* of the 2 steps
  # ==> so it also includes "good" reads, but in reverse complement from other step
  tot_raw_sum = sum(demultiplex_summary_df$read_count) / 2

  demultiplex_summary_df = rbind(demultiplex_summary_df, data.frame(type ="raw_data_R1", read_count=tot_raw_sum))

  # create use info column
  demultiplex_summary_df$type_use = "used_data"
  demultiplex_summary_df$type_use[grepl(pattern = "unknown|discard|unassembled|5prime|Undetermined", x = demultiplex_summary_df$type )] = "discarted_data"

  tot_sample_count = demultiplex_summary_df$read_count[demultiplex_summary_df$type == "demultiplexed_R1"]

  # Create order for plotting
  demultiplex_summary_df$type = factor(demultiplex_summary_df$type,
                                       levels=rev(c('raw_data_R1',
                                                    "unknown_R1",
                                                    "unknown-demultiplexed_R1",
                                                    "unknown-round1-1",
                                                    "unknown-round1bis-1",
                                                    'unknown_round2_1',
                                                    "not-demultiplexed_R1",
                                                    "demultiplexed_R1"
                                       )))
  return(demultiplex_summary_df)
  }
