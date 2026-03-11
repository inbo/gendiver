## Wrappers for rendering
### Should maybe included in gendiver

### SAMPLE SHEET ###

#' Render Rmarkdown Report - Lab Sample Sheet
#'
#' //////////////////
#'
#' @param sample_sheet_url Filepath for \link[readxl]{read_xlsx} or google-sheets id \link[googlesheets4]{range_read}
#' @param out_dir Output folder (Default="sample-sheet-QC-report")
#'
#' @returns HTML report and barcode-files for demultiplexing, in `out_dir`
#' @export
#'
#' @examples
#' #To add
report.samplesheet = function(sample_sheet_url, out_dir="sample-sheet-QC-report") {
  dir.create(out_dir, showWarnings = F)
  out_dir = normalizePath(out_dir)
  tmp=tempdir()
  rmarkdown::render(input = system.file("rmd_report","00_sample_sheet_report.Rmd", package = "gendiver"),
                    output_file = "samplesheet-QC-report.html",
                    output_format="html_document",
                    output_dir = out_dir,
                    knit_root_dir = tmp, intermediates_dir = tmp,
                    params = list("SAMPLE_SHEET_URL"=sample_sheet_url,
                                  "BC_OUT_DIR"=out_dir),
                    quiet = F
  )
}

### DEMULTIPLEX ###

#' Render Rmarkdown Report - Demultiplexing
#'
#' //////// HPC-ONLY //////////
#'
#' @param sample_sheet_url Filepath for \link[readxl]{read_xlsx} or google-sheets id \link[googlesheets4]{range_read}
#' @param demultiplex_result_path HPC-ONLY: file_path to demultiplexing results
#' @param out_dir Output folder (Default="demultiplex-QC-report")
#'
#' @returns HTML report for demultiplexing results, in `out_dir`
#' @export
#'
#' @examples
#' #To add
report.demultiplex = function(sample_sheet_url, demultiplex_result_path, out_dir="demultiplex-QC-report") {
  dir.create(out_dir, showWarnings = F)
  out_dir = normalizePath(out_dir)
  tmp=tempdir()
  rmarkdown::render(input = system.file("rmd_report","01_demultiplexing_report.Rmd", package = "gendiver"),
                    output_file = "demultiplex-QC-report.html",
                    output_format="html_document",
                    output_dir = out_dir,
                    knit_root_dir = tmp, intermediates_dir = tmp,
                    params = list("SAMPLE_SHEET_URL"=sample_sheet_url,
                                  "DEMULTIPLEX_RESULT_PATH"=demultiplex_result_path,
                                  "OUT_DIR"=out_dir),
                    quiet = F
  )

}

### PRE-PROCESS ###

#' Render Rmarkdown Report - Pre-Processing
#'
#' //////// HPC-ONLY //////////
#'
#' @param sample_sheet_url Filepath for \link[readxl]{read_xlsx} or google-sheets id \link[googlesheets4]{range_read}
#' @param preprocess_result_path HPC-ONLY: file_path to pre-processing results
#' @param out_dir Output folder (Default="pre-process-QC-report")
#'
#' @returns HTML report for demultiplexing results, in `out_dir`
#' @export
#'
#' @examples
#' #To add
report.preprocess = function(sample_sheet_url, preprocess_result_path, out_dir="pre-process-QC-report") {
  dir.create(out_dir, showWarnings = F)
  out_dir = normalizePath(out_dir)
  tmp=tempdir()
  rmarkdown::render(input = system.file("rmd_report","02_preprocessing_report.Rmd", package = "gendiver"),
                    output_file = "preprocess-QC-report.html",
                    output_format="html_document",
                    output_dir = out_dir,
                    knit_root_dir = tmp, intermediates_dir = tmp,
                    params = list("SAMPLE_SHEET_URL"=sample_sheet_url,
                                  "PREPROCESS_RESULT_PATH"=preprocess_result_path,
                                  "OUT_DIR"=out_dir),
                    quiet = F
  )

}

### Basic analysis ###

#' Render Rmarkdown Report - Basic-analysis Clustering and Taxonomy
#'
#' //////// //////////
#'
#' @param sample_sheet_url Filepath for \link[readxl]{read_xlsx} or google-sheets id \link[googlesheets4]{range_read}
#' @param tax_file_path Path to taxonomy table file (can also be link to Googe Drive)
#' @param otu_table_path Path to OTU table file (can also be link to Googe Drive)
#' @param seq_fasta_path Path or URL (Google Drive) to a FASTX file that can be read with \link[Biostrings]{readDNAStringSet}
#' @param out_dir Output folder (Default="operational-data-QC-report")
#'
#' @returns HTML report for Clustering and Taxonomy results, in `out_dir`
#' @export
#'
#' @examples
#' #To add
report.operational_data = function(sample_sheet_url,
                                   tax_file_path,
                                   otu_table_path,
                                   seq_fasta_path,
                                   out_dir="operational-data-QC-report") {
  dir.create(out_dir, showWarnings = F)
  out_dir = normalizePath(out_dir)
  tmp=tempdir()

  rmarkdown::render(input = system.file("rmd_report","03_basic_analysis_report.Rmd", package = "gendiver"),
                    output_file = "operational-data-report.html",
                    output_format="html_document",
                    output_dir = out_dir,
                    knit_root_dir = tmp, intermediates_dir = tmp,
                    params = list("SAMPLE_SHEET_URL"=sample_sheet_url,
                                  "TAX_PATH"=tax_file_path,
                                  "OTU_PATH"=otu_table_path,
                                  "FASTA_PATH"=seq_fasta_path,
                                  "OUT_DIR"=out_dir),
                    quiet = F
  )

}
