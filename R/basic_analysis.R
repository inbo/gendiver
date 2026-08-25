######################
### Basic Analysis ###
######################

#' General barplot function for taxa inspection
#'
#' Quickly visualize taxa groups in a barplot
#'
#' @param ps_obj Phyloseq object as INPUT
#' @param taxa_rank Taxonomic rank to visualise the data on
#' @param taxa_select String to select the data to plot based on taxonomy. Multiple taxonomies can be passed using "|". Taxonomies do not have to be of the same level, and are case insensitive. E.g. `taxa_select='Amphibia|Homo sapiens'`
#' @param taxa_excl String to select the taxonomies to exclude from the data before plotting. Same logic applies as for `taxa_select`. (Default 'Homo')
#' @param cutoff integer, remove taxa that have less than n % reads in the dataset.
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
  ps_obj = tax_glom(ps_obj, taxrank = taxa_rank)

  # Remove 0 read samples
  # ps_obj = phyloseq::prune_samples(phyloseq::sample_sums(ps_obj) > 0, ps_obj)

  if (RRA){
    otu_df = data.frame(phyloseq::otu_table(ps_obj))
    otu_df[otu_df == 0] = NA
    otu_df =  phyloseq::otu_table(vegan::decostand(otu_df, method = "total", MARGIN = 2, na.rm = T), taxa_are_rows = T)
    sample_names(otu_df) = gsub("^X", "", sample_names(otu_df))
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
show.default_project_barplots = function(ps_obj, RRA=F, out_path=NA){

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
      dir.create(out_dir_pproj, showWarnings = FALSE)

      pdf(file=file.path(out_dir_pproj, paste0(title_text1, "_barplot.pdf")),
          width=5+ 0.2*nsamples(sub_ps),height=8)
      print(pl1)
      dev.off()

      pdf(file=file.path(out_dir_pproj, paste0(title_text2, "_barplot.pdf")), width=5+0.2*nsamples(sub_ps),height=8)
      print(pl2)
      dev.off()

      pdf(file=file.path(out_dir_pproj, paste0(title_text3, "_barplot.pdf")), width=5+0.2*nsamples(sub_ps),height=8)
      print(pl3)
      dev.off()

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


# Legacy function to make a barplot for a given taxon based on a phyloseq object (code from Annelies Haegeman via e-mail)
# ps_barplot_annelies = function(ps, tax_rank="species"){
#   #reorder taxa according to abundance values, for this we need to convert the phyloseq object to a (long) dataframe
#   ps.df <- psmelt(ps)
#   ps.df$tax_rank <- with(ps.df, reorder(get(tax_rank), Abundance))
#   orderedtaxa<-rev(levels(ps.df$tax_rank))
#
#   #define number of colors based on number of full_names
#   getPalette = colorRampPalette(brewer.pal(12, "Set3"))
#   #taxaList = unique(tax_table(ps.filt.taxglom.fish.filt.frac)[,"full_name"])
#   taxaList = levels(ps.df$tax_rank)
#   taxaPalette = getPalette(length(taxaList))
#   #names(taxaPalette) = taxaList
#   names(taxaPalette) = levels(ps.df$tax_rank) #assign names to the colors in the order of abundance as defined above
#   #put color of level "unknown" (=last element of vector) to dark grey
#   #taxaPalette[length(taxaPalette)] <- "#7E7E7E"
#
#   #barplot
#   barplot<-ggplot(ps.df, aes(x=Sample, y=Abundance, fill=tax_rank)) +
#     theme_bw() +
#     theme(axis.text.y=element_text(size=22),legend.text=element_text(size=22),legend.key.size=unit(0.4, "cm"),
#           axis.title=element_text(size=40),legend.title=element_text(size=40),
#           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
#     geom_bar(aes(color=tax_rank, fill=tax_rank), stat="identity", position="stack", color="white", width=0.6, linetype=0) +
#     #facet_grid(Depth~Season, scales="free", space="free_x") +
#     #facet_wrap(~Timepoint, scales = "free_x", ncol=1) +
#     scale_y_continuous(expand = c(0,0)) +
#     guides(fill = guide_legend(ncol = 1, title = paste(tax_rank))) + #add this line to force the legend in 1 column
#     scale_fill_manual(values= taxaPalette, limits=levels(ps.df$tax_rank))
# }
