######################
### Basic Analysis ###
######################

#' General barplot function for taxa inspection
#'
#' Quickly visualize taxa groups in a barplot
#'
#' @param ps_obj Phyloseq object as INPUT
#' @param taxa_rank Taxonomic rank to visualise the data on
#' @param taxa_select String to select the data to plot based on taxonomy. Multiple taxonomies can be passed using "|". Taxonomies do not have to be of the same level, and are case insensitive. E.g. taxa_select='Amphibia|Homo sapiens`
#' @param taxa_excl String to select the taxonomies to exclude from the data before plotting. Same logic applies as for taxa_select. (Default 'Homo')
#' @param cutoff integer, remove taxa that have less than n percent reads in the dataset.
#' @param RRA Convert reads to relative read abundance (RRA) per sample. (Default FALSE)
#'
#' @returns ggplot2 object
#' @export
#'
#' @examples
#' #To add
make.taxa_barplot = function(ps_obj, taxa_rank="species", taxa_select=NA, taxa_excl="Homo", cutoff=1, RRA=FALSE){
  # Collapse taxonomy to 1 string
  tax_tab_obj = apply(data.frame(phyloseq::tax_table(ps_obj)),  1, paste0, collapse="")

  # Search for any match with taxa to select
  if (!is.na(taxa_select)){
    select_mask = tax_tab_obj[grepl(taxa_select, tax_tab_obj, ignore.case = T)]
    ps_obj = phyloseq::prune_taxa(names(select_mask), ps_obj)
  }

  # Search for any match with taxa to exclude
  if (!is.na(taxa_excl)){
    excl_mask = tax_tab_obj[!grepl(taxa_excl, tax_tab_obj, ignore.case = T)]
    ps_obj = phyloseq::prune_taxa(names(excl_mask), ps_obj)
  }

  # Filter low read taxa
  cutoff_mask = rowSums(phyloseq::otu_table(ps_obj)) / sum(rowSums(phyloseq::otu_table(ps_obj))) * 100
  ps_obj = phyloseq::prune_taxa(cutoff_mask > cutoff, ps_obj)

  # Remove any 0 read taxa (if any)
  ps_obj = phyloseq::prune_taxa(phyloseq::taxa_sums(ps_obj) > 0, ps_obj)

  # Merge taxa
  ps_obj = phyloseq::tax_glom(ps_obj, taxrank = taxa_rank)

  # Remove 0 read samples
  # ps_obj = phyloseq::prune_samples(phyloseq::sample_sums(ps_obj) > 0, ps_obj)

  if (RRA){
    otu_df = data.frame(phyloseq::otu_table(ps_obj))
    otu_df[otu_df == 0] = NA
    otu_df =  phyloseq::otu_table(vegan::decostand(otu_df, method = "total", MARGIN = 2, na.rm = T), taxa_are_rows = T)
    phyloseq::sample_names(otu_df) = gsub("^X", "", phyloseq::sample_names(otu_df))
    phyloseq::otu_table(ps_obj) = otu_df
  }

  # plot
  pl1 = phyloseq::plot_bar(ps_obj, fill = taxa_rank)

  # Return
  return(pl1)
}


#' Create default barplots for project level QC
#'
#' Pass run level data to generate project level barplots with some pre-configured parameters to visualize an overview, a fish and amphibian specific plot and a most abundant species plot.
#'
#' @param ps_obj Phyloseq object as INPUT
#' @param RRA Convert reads to relative read abundance (RRA) per sample. (Default FALSE)
#' @param out_path [RECOMMENDED] Directory to store the plots in. For projects with a lot of samples this should always be used. Otherwise your R-session might crash!
#'
#' @returns List of ggplot2 objects if out_path is not defined.
#' @export
#'
#' @examples
#' #To add
make.default_project_barplots = function(ps_obj, RRA=F, out_path=NA){

  PRJ=phyloseq::sample_data(ps_obj)$PROJECT

  print(table(PRJ, useNA="ifany"))

  prj_l = list()

  for (sub_project_i in unique(PRJ)){

    sub_ps = phyloseq::prune_samples(phyloseq::sample_data(ps_obj)$PROJECT==sub_project_i, ps_obj)

    # Make overview barplot
    cutoff_pct = 0
    title_text1 = paste0(sub_project_i, "_overview")
    pl1 = make.taxa_barplot(
      sub_ps, taxa_rank = "custom_taxon", taxa_excl =NA, RRA = RRA, cutoff = cutoff_pct) +
      ggplot2::ggtitle(title_text1)

    # Make Fish/amphibian barplot
    cutoff_pct = 0.1
    title_text2 = paste0(sub_project_i, "_amphibia_fish", "_cutoff_", cutoff_pct, "_percent")
    pl2 = make.taxa_barplot(
      sub_ps, taxa_rank = "species", taxa_select = 'Actinopteri|Amphibia',
      RRA = RRA, cutoff = cutoff_pct) +
      ggplot2::ggtitle(title_text2)

    # Make most abundant barplot
    cutoff_pct = 2
    title_text3 = paste0(sub_project_i, "_species_no_human", "_cutoff_", cutoff_pct, "_percent")
    pl3 = make.taxa_barplot(
      sub_ps, taxa_rank = "species", taxa_select = NA, RRA = RRA,
      cutoff = cutoff_pct) +
      ggplot2::ggtitle(title_text3)

    if (!is.na(out_path)){
      # save
      out_dir_pproj = file.path(out_path, sub_project_i)
      nsamples = phyloseq::nsamples(sub_ps)
      dir.create(out_dir_pproj, showWarnings = FALSE, recursive = T)

      grDevices::pdf(
        file=file.path(out_dir_pproj, paste0(title_text1, "_barplot.pdf")),
          width=5+ 0.2*nsamples(sub_ps),height=8)
      print(pl1)
      grDevices::dev.off()

      grDevices::pdf(
        file=file.path(out_dir_pproj, paste0(title_text2, "_barplot.pdf")), width=5+0.2*nsamples(sub_ps),height=8)
      print(pl2)
      grDevices::dev.off()

      grDevices::pdf(
        file=file.path(out_dir_pproj, paste0(title_text3, "_barplot.pdf")), width=5+0.2*nsamples(sub_ps),height=8)
      print(pl3)
      grDevices::dev.off()

    } else {
      prj_l[[sub_project_i]]$overview = pl1
      prj_l[[sub_project_i]]$amph_fish = pl2
      prj_l[[sub_project_i]]$all = pl3
    }

  }
  if (is.na(out_path)){
    return(prj_l)
  } else {
    return()
  }

}


combine.obi_otu_tax = function(otu_df, tax_df){
  res_table = merge(tax_df, otu_df, by="row.names")
  colnames(res_table)[1] = "ID"
  rownames(res_table) = res_table$ID

  res_table$COUNT = rowSums(otu_df)

  table_sorted <- res_table[order(res_table$COUNT, decreasing = T),] #sort from largest to smallest number of total counts per ASV
  return(table_sorted)
}


combine.ps_otu_tax = function(ps){
  t_out <- merge(
    as.data.frame(phyloseq::tax_table(ps)),
    as.data.frame(phyloseq::otu_table(ps)), by="row.names")
  row.names(t_out) = t_out$Row.names
  t_out = t_out[, -1]
  return(t_out)
}


#' Export eDNA water datasets
#'
#' Convert "raw" operational data to input for data analysis
#'
#' @param sample_sheet_df Data frame of (lab) sample sheet
#' @param otu_df Data frame of OTU/ASV table
#' @param tax_df Data frame with taxonomic annotation data
#' @param ID_cutoff Reference database similarity needed to filter OTU/ASV data.
#' @param out_path Directory to store the plots in. (Default = ".")
#'
#' @returns List of ggplot2 objects if out_path is not defined.
#' @export
#'
#' @examples
#' #To add
export.data_sets = function(
    sample_sheet_df,
    otu_df,
    tax_df,
    ID_cutoff=1, out_path="."){

  ### PREPARE DATA ###
  PRJ_DIR = out_path
  # Merge and sort all otu+tax data
  table_sorted = combine.obi_otu_tax(otu_df, tax_df)

  # remove sequences that have a database match with less than "cutoff" identity with a database sequence
  unknowns = table_sorted[table_sorted$BEST_IDENTITY < ID_cutoff,]

  #finally, we remove the sequences with less than X% identity to their reference sequence
  table.filt<-table_sorted[table_sorted$BEST_IDENTITY >= ID_cutoff,]

  # Make PS
  taxa_select = tax_df |> select.taxonomy_obitools3() |> as.matrix()

  ref_seqs = Biostrings::DNAStringSet(table.filt$NUC_SEQ)
  names(ref_seqs) = rownames(table.filt)

  myps = phyloseq::phyloseq(
    phyloseq::sample_data(sample_sheet_df),
    phyloseq::otu_table(otu_df, T),
    phyloseq::tax_table(taxa_select),
    ref_seqs
  )

  ssdata = phyloseq::sample_data(myps)
  PROJ_NAME = paste0(paste0(unique(ssdata$RUN_CODE), collapse = "_"), "_", paste0(unique(ssdata$PRIMERS), collapse = "_"))

  # Make RRA version
  ps.RRA = phyloseq::transform_sample_counts(myps, function(OTU) OTU/sum(OTU)*100)
  table.filt.RRA = combine.ps_otu_tax(ps.RRA)[row.names(table.filt),]

  # Tax glom
  ps.filt.taxglom <- phyloseq::tax_glom(myps, taxrank="species")
  table.filt.taxglom = combine.ps_otu_tax(ps.filt.taxglom)

  ps.filt.taxglom.frac <- phyloseq::transform_sample_counts(ps.filt.taxglom, function(OTU) OTU/sum(OTU)*100)
  table.filt.taxglom.RRA = combine.ps_otu_tax(ps.filt.taxglom.frac)

  ### WRITE DATA ###
  #save sorted table to output txt file
  message("Writing 'raw_table_sorted.tsv' ...")
  utils::write.table(
    table_sorted,
    file=file.path(PRJ_DIR, paste0(PROJ_NAME, "_all_reads.tsv")),
    sep="\t", row.names=T, col.names=NA, quote=F)

  message("Writing 'unknowns.tsv' ...")
  #write to output
  utils::write.table(
    unknowns,
    file=file.path(PRJ_DIR, paste0(PROJ_NAME, "_ID", ID_cutoff,"_unknown_reads.tsv")),
    sep="\t",row.names=F,col.names=T,quote=F)

  message("Writing 'filtered_table_sorted.txt' ...")
  utils::write.table(
    table.filt,
    file.path(PRJ_DIR, paste0(PROJ_NAME, "_ID",ID_cutoff, "_reads.tsv")),
    sep="\t",row.names=F,col.names=T,quote=F)

  message("Writing 'filtered_table_sorted_RRA.txt' ...")
  utils::write.table(
    table.filt.RRA,
    file.path(PRJ_DIR, paste0(PROJ_NAME, "_ID",ID_cutoff, "_RRA.tsv")),
    sep="\t",row.names=F,col.names=T,quote=F)

  message("Writing 'agglomerated_species.tsv' ...")
  utils::write.table(
    table.filt.taxglom,
    file.path(PRJ_DIR, paste0(PROJ_NAME, "_ID",ID_cutoff, "_species_reads.tsv")),
    sep="\t",row.names=F,col.names=T,quote=F)

  message("Writing 'agglomerated_species_RRA.tsv' ...")
  utils::write.table(
    table.filt.taxglom.RRA,
    file.path(PRJ_DIR, paste0(PROJ_NAME, "_ID",ID_cutoff, "_species_RRA.tsv")),
    sep="\t",row.names=F,col.names=T,quote=F)

  message("Writing 'phyloseq.rds' ...")
  saveRDS(myps, file.path(PRJ_DIR, paste0(PROJ_NAME, "_ID",ID_cutoff , "_phyloseq", ".rds" )))

  #write number of left reads to file
  message("Writing 'number_of_reads_after_processing.tsv' ...")
  utils::write.table(
    as.data.frame(phyloseq::sample_sums(myps)),
    file.path(PRJ_DIR, paste0(PROJ_NAME, "_ID",ID_cutoff,"sample_readcount", ".tsv")),
    sep="\t",row.names=T,col.names=F,quote=F)

  message("Done!")


}


#' Merge technical replicates into samples
#'
#' Merge the data from (generally 3) replicate amplicons into a single sample. Here, this is done by taking the average of each replicate's relative read abundance (RRA).
#'
#' @param otu_df Data frame of OTU/ASV table
#' @param tax_df Data frame with taxonomic annotation data
#' @param sample_sheet_df Data frame of (lab) sample sheet
#'
#' @returns phyloseq object
#' @export
#'
#' @examples
#' #To add
ps.merge_replicates = function(otu_df, tax_df, sample_sheet_df){
  ## First RRA
  otu_rra = vegan::decostand(otu_df, method = "total", MARGIN = 2)
  ## Then, merge replicate per filter - mean()
  x = merge(t(otu_rra), sample_sheet_df[, "filter_code", drop=F], by=0)
  x = stats::aggregate(data=x[-1], .~filter_code, mean)
  otu_rra_merged = data.frame(t(x[-1]))
  colnames(otu_rra_merged) = x[[1]]

  # Update Sample data
  ss_merged = sample_sheet_df[!duplicated(sample_sheet_df$filter_code),]
  rownames(ss_merged) = ss_merged$filter_code

  ## make PS
  ps = phyloseq::phyloseq(
    phyloseq::otu_table(otu_rra_merged, taxa_are_rows = T),
    phyloseq::tax_table(as.matrix(tax_df)),
    phyloseq::sample_data(ss_merged)
  )

  return(ps)
}

