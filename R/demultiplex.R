##################
# Demultiplexing #
##################

#' Read per-sample readcounts from summary files after demultiplexing job
#'
#' Find and read auto-generated demultiplexing results from demultiplexing folder
#'
#' @param data_dir File path to the `02_demultiplex/` dir
#' @param summary_file_pattern Pattern to recognize summary files (default="_stats.tsv).
#' @param legacy_mode Enable to read legacy "readcounts_per_sample.txt" files (default=FALSE)
#'
#' @returns Dataframe with read counts per sample
#' @export
#'
#' @examples
#' #To add
read.summary_files_demultiplex = function(data_dir,
                                          summary_file_pattern = "stats.tsv",
                                          legacy_mode=F){
  # Initialize readcount tracking df
  read_count_track.df = data.frame("file_name" = character(),
                                   "read_count" = numeric(),
                                   "multiplexed_source" = character())

  ## Read for each DEMULTIPLEXED lib the sample counts (from count seqs PBS)
  demultiplex_res_lib_paths = list.dirs(path = file.path(data_dir), full.names = TRUE, recursive = F)

  for (mylib in demultiplex_res_lib_paths){
    readcounts_file = list.files(mylib, pattern = summary_file_pattern, full.names = T)
    if (length(readcounts_file) > 0) {
      if (legacy_mode){
        # Old style: readcounts_per_sample.txt -> 2 columns USI readcounts
        count_tab = utils::read.table(readcounts_file)
      } else {
        # New: seqkit stats output -T
        count_tab = utils::read.table(readcounts_file, header = T, fill = T)
        count_tab = count_tab[, c("file", "num_seqs")]
      }
      colnames(count_tab) = c("file_name", 'read_count')
      count_tab$multiplexed_source = basename(mylib)
      read_count_track.df = rbind(read_count_track.df, count_tab)

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
  read_count_track.df$type = factor(read_count_track.df$type, levels = rev(unique(read_count_track.df$type)))

  return(read_count_track.df)
}


#' Summarize demultiplexing readcounts for plotting
#'
#'  demultiplexed readcounts per sample
#'
#' @param df Data frame with read counts per sample (output of \link[gendiver]{read.summary_files_demultiplex})
#' @param formula Formula to use for aggregation (default: `read_count ~ type + multiplexed_source`)
#'
#' @returns Summary data ready for plotting
#' @export
#'
#' @examples
#' #To add
summarize.readcounts_demultiplex = function(df, formula=read_count ~ type + multiplexed_source){

  # Convert type to factor (for grouping ggplot)
  df$type = factor(df$type, levels = rev(unique(df$type)))

  #1. Overview of all reads
  demultiplex_summary_df = stats::aggregate(formula, data = df, sum)

  if (ncol(demultiplex_summary_df) < 3){
    tot_raw_sum = sum(demultiplex_summary_df$read_count) / 2

  } else {

  }

  tot_raw_sum_pp = stats::aggregate(formula, data = demultiplex_summary_df, sum)

  tot_raw_sum_pp$read_count = tot_raw_sum_pp$read_count/2
  tot_raw_sum_pp$type = "raw_data_R1"
  demultiplex_summary_df_pp = rbind(demultiplex_summary_df_pp, tot_raw_sum_pp)

  demultiplex_summary_df = rbind(demultiplex_summary_df,
                                 data.frame(type = "raw_data_R1", read_count = tot_raw_sum))

  # create use info column
  demultiplex_summary_df$type_use = "used_data"
  demultiplex_summary_df$type_use[grepl(pattern = "unknown|discard|unassembled|5prime|Undetermined", x = demultiplex_summary_df$type )] = "discarted_data"

  tot_sample_count = demultiplex_summary_df$read_count[demultiplex_summary_df$type == "demultiplexed_R1"]

  # Create order for plotting
  demultiplex_summary_df$type = factor(demultiplex_summary_df$type,
                                       levels = rev(c('raw_data_R1',
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
