##################
# Pre-processing #
##################

#' Read read-counts from logfile - PEAR merging
#'
#'
#' @param logfile Path to a PEAR log file
#'
#' @returns Integer with readcount of sample **after** PEAR merging (i.e. merged)
#' @export
#'
#' @examples
#' #To add
extract.pear_log.readcount = function(logfile){
  readsl = readLines(logfile)[6]
  return(as.integer(trimws(strsplit(readsl, "Merged")[[1]][1])))
}


#' Read read-counts from logfile - VSEARCH --mergepairs
#'
#'(version 2.22.1)
#'
#' @param logfile Path to a VSEARCH log file
#'
#' @returns Integer with readcount of sample **after** VSEARCH --mergepairs (i.e. merged)
#' @export
#'
#' @examples
#' #To add
extract.mergepairs_log.readcount = function(logfile){
  readsl = readLines(logfile)[6]
  return(as.integer(trimws(strsplit(readsl, "Merged")[[1]][1])))
}


#' Read vsearch --fastq_mergepairs log file to get total input reads
#'
#'(version 2.22.1)
#'
#' @param logfile Path to a VSEARCH log file
#'
#' @returns Integer with total readcount of sample **before** VSEARCH --mergepairs (i.e. paired-end input)
#' @export
#'
#' @examples
#' #To add
extract.mergepairs_log.tot_input = function(logfile){
  fcon = file(logfile)
  merge_input = trimws(readLines(fcon)[5])
  close(fcon)
  return(as.numeric(gsub("[^0-9]", "", merge_input)))
}


#' Read read-counts from logfile - CUTADAPT
#'
#'cutadapt --version 5.0
#'
#' @param logfile Path to a CUTADAPT log file
#'
#' @returns Integer with readcount of sample **after** CUTADAPT (i.e. trimmed)
#' @export
#'
#' @examples
#' #To add
extract.cutadapt_log.readcount = function(logfile){
  if (readLines(logfile)[4] == "No reads processed!" ) {
    return(0)
  }
  readsl = readLines(logfile)[12]
  x = strsplit(readsl, "\\s")[[1]]
  y = x[length(x)-1]

  return(as.integer(gsub(",","",y)))
}


#' Read read-counts from logfile - VSEARCH - fastq_filter
#'
#'(version 2.22.1)
#'
#' @param logfile Path to a VSEARCH log file
#'
#' @returns Integer with total readcount of sample **after** VSEARCH - fastq_filter (i.e. qc-filtered)
#' @export
#'
#' @examples
#' #To add
extract.fastqfilter_log.readcount = function(logfile){
  readsl = readLines(logfile)[5]
  x = strsplit(readsl, "\\s")[[1]][1]
  return(as.integer(x))
}


#' Map pre-processing output files to used tool
#'
#'In the INBO eDNA emtabarcoding workflow files have structured filenames. Translate these filenames to the used tool.
#'
#' @param x Filename (string)
#'
#' @returns One of: "saber-demultiplexing", "pear-merging", "vsearch-merging", "cutadapt-trimming" or "vsearch-quality-filter"
#' @export
#'
#' @examples
#' #To add
map.filetype2tool = function(x){
  # Initiate translation-result list
  x_trans = rep(NA, length(x))
  # demultiplexing
  x_trans[grepl('_round2|-round1|demultiplex', x)] = 'saber-demultiplexing'
  # pre-processing
  x_trans[grepl('assembled|discarded', x)] = 'pear-merging'
  x_trans[grepl('merge', x)] = 'vsearch-merging'
  x_trans[grepl('trimm', x)] = 'cutadapt-trimming'
  x_trans[grepl('filtered', x)] = 'vsearch-quality-filter'
  x_trans[grepl('fastqfilter', x)] = 'vsearch-quality-filter'
  # Output
  return(x_trans)
}


#' Read multiple log files for pre-processing workflow
#'
#' Wrapper to extract the remaining readcounts from .log files for each step of the INBO eDNA metabarcoding pre-processing workflow
#' Implemented for:
#'   - merged (VSEARCH --fastq_mergepairs)
#'   - primer-trimmed (cutadapt)
#'   - QC filtered (VSEARCH --fastq_filter)
#'
#' @param mydir Path to a folder containing .log files from pre-processing
#' @param stepname One of: "merged", "trimmed" or "filter"
#' @param logpattern Pattern to find correct log files, see \link[base]{list.files} (default = ".log")
#'
#' @returns Dataframe with readcounts per file
#' @export
#'
#' @examples
#' #To add
read.logs_preprocess = function(mydir, stepname, logpattern = ".log"){

  # Check if logfiles
  lfiles = list.files(mydir, pattern = logpattern, full.names = T)

  read_count_track.df = data.frame("file_name" = basename(lfiles),
                                   "read_count" = NA)

  if (stepname=="merged"){
    read_count_track.df$read_count = sapply(lfiles, extract.mergepairs_log.readcount)
  }

  if (stepname=="trimmed"){
    read_count_track.df$read_count = sapply(lfiles, extract.cutadapt_log.readcount)
  }

  if (stepname=="filter"){
    read_count_track.df$read_count = sapply(lfiles, extract.fastqfilter_log.readcount)
  }

  # Fancy wrapping functionality
  read_count_track.df$sample_id = extract_sample_name(read_count_track.df$file_name)
  read_count_track.df$type = sapply(read_count_track.df$file_name, extract_filename_info, ".*\\.")
  read_count_track.df$toolname = map.filetype2tool(read_count_track.df$type)

  return(read_count_track.df)
}

