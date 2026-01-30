##################
# Pre-processing #
##################

hpc.read.preprocess_logs = function(mydir, stepname, logpattern = ".log", runShell=T){

  # Check if logfiles
  lfiles = list.files(mydir, pattern = logpattern, full.names = T)

  read_count_track.df = data.frame("file_name" = basename(lfiles),
                                   "read_count" = NA)

  if (stepname=="merged"){
    read_count_track.df$read_count = sapply(lfiles, extract.mergepairs_log.readcount)
    new_pattern=".fastq"
  }

  if (stepname=="trimmed"){
    read_count_track.df$read_count = sapply(lfiles, extract.cutadapt_log.readcount)
    new_pattern=".fastq"
  }

  if (stepname=="filter"){
    read_count_track.df$read_count = sapply(lfiles, extract.fastqfilter_log.readcount)
    new_pattern=".fasta"
  }

  return(read_count_track.df)
}

# # read read-counts from logfile - PEAR merging
# extract_readcount_from_pear_log = function(logfile){
#   readsl = readLines(logfile)[6]
#   return(as.integer(trimws(strsplit(readsl, "Merged")[[1]][1])))
# }

# read read-counts from logfile - VSEARCH --mergepairs
extract.mergepairs_log.readcount = function(logfile){
  readsl = readLines(logfile)[6]
  return(as.integer(trimws(strsplit(readsl, "Merged")[[1]][1])))
}

# read vsearch --fastq_mergepairs log file to get total input reads (version 2.22.1)
extract.mergepairs_log.tot_input = function(log.path){
  fcon = file(log.path)
  merge_input = trimws(readLines(fcon)[5])
  close(fcon)
  return(as.numeric(gsub("[^0-9]", "", merge_input)))
}

# read read-counts from logfile - CUTADAPT
extract.cutadapt_log.readcount = function(logfile){
  if (readLines(logfile)[4] == "No reads processed!" ) {
    return(0)
  }
  readsl = readLines(logfile)[12]
  x = strsplit(readsl, "\\s")[[1]]
  y = x[length(x)-1]

  return(as.integer(gsub(",","",y)))
}

# read read-counts from logfile - VSEARCH - fastq_filter
extract.fastqfilter_log.readcount = function(logfile){
  readsl = readLines(logfile)[5]
  x = strsplit(readsl, "\\s")[[1]][1]
  return(as.integer(x))
}

translate.type2tool = function(x){
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
