################
### QC Plots ###
################


#' Barplot and ordination for QC on possible strip-swaps
#'
#' Taxon-barplot and ordination with data grouped per plate and per column to investigate lab-specific errors/trends
#'
#' @param ps_plate Phyloseq object as INPUT  (Ideally of 1 plate)
#' @param topX integer, topX taxa per COLUMN on the plate (default=25)
#'
#' @returns List with 2 plots: barplot and ordination
#' @export
#'
#' @examples
#' #To add
qcplot.plate_column_asv_barplot = function(ps_plate, topX=25){
  out_list = c()
  # Create sample column to merge on
  phyloseq::sample_data(ps_plate)$X = paste(
    phyloseq::sample_data(ps_plate)$WELL_COLUMN,
    phyloseq::sample_data(ps_plate)$TECHNICAL_REPLICATE, sep = "_")

  # And taxglom for easier interpretation
  ps_plate = phyloseq::tax_glom(ps_plate, taxrank = "species")

  # Create relative abundance ps
  ps_plate.frac = phyloseq::transform_sample_counts(ps_plate, function(OTU) OTU/sum(OTU)*100)

  # Plot bars per column for topX taxa
  order.sum = tapply(phyloseq::taxa_sums(ps_plate),
                     phyloseq::tax_table(ps_plate)[, "species"], sum, na.rm=TRUE)


  ps_plate_topX = phyloseq::prune_taxa(
    (phyloseq::tax_table(ps_plate)[, "species"] %in% names(sort(order.sum, TRUE))[1:topX]),
    ps_plate)

  # Plot data
  phyloseq::plot_bar(ps_plate_topX, x="WELL_COLUMN", fill="species", facet_grid=~TECHNICAL_REPLICATE)

  ## MERGE the columns into samples
  ps_plate_merged = phyloseq::merge_samples(x=ps_plate, "X")

  ## Barplot
  phyloseq::plot_bar(ps_plate_merged, x="TECHNICAL_REPLICATE", fill="custom_taxon", facet_grid=~WELL_COLUMN)

  ps_plate_merged_topX = phyloseq::prune_taxa(
    (phyloseq::tax_table(ps_plate_merged)[, "species"] %in% names(sort(order.sum, TRUE))[1:topX]),
    ps_plate_merged)
  ps_plate_merged_topX.frac = phyloseq::transform_sample_counts(ps_plate_merged_topX, function(OTU) OTU/sum(OTU)*100)

  pb = phyloseq::plot_bar(ps_plate_merged_topX.frac, x="TECHNICAL_REPLICATE", fill="species", facet_grid=~WELL_COLUMN) +
    ggplot2::theme(legend.position="none")

  out_list$barplot = pb

  ## Ordination
  ps_plate_merged.frac.ordination = phyloseq::ordinate(ps_plate_merged_topX.frac, "PCoA", "bray")

  phyloseq::sample_data(ps_plate_merged_topX.frac)$WELL_COLUMN = as.factor(phyloseq::sample_data(ps_plate_merged_topX.frac)$WELL_COLUMN)
  phyloseq::sample_data(ps_plate_merged_topX.frac)$TECHNICAL_REPLICATE = as.factor(phyloseq::sample_data(ps_plate_merged_topX.frac)$TECHNICAL_REPLICATE)

  po = phyloseq::plot_ordination(ps_plate_merged_topX.frac, ps_plate_merged.frac.ordination,
                       color="WELL_COLUMN", shape = "TECHNICAL_REPLICATE") +
    ggplot2::geom_text(
      ggplot2::aes(label=.data$WELL_COLUMN), size = 5, nudge_y = -0.015)

  out_list$ordination_plot = po

  return(out_list)

}

clean.ps_sample_sheet = function(ps_obj){
  sample_metadata = as.data.frame(phyloseq::sample_data(ps_obj))
  sample_metadata$PLATE = factor(sample_metadata$PLATE_NOTES)
  sample_metadata$WELL_COLUMN = factor(sample_metadata$WELL_COLUMN)
  sample_metadata$WELL_ROW = factor(sample_metadata$WELL_ROW)
  sample_metadata$TECHNICAL_REPLICATE = factor(sample_metadata$TECHNICAL_REPLICATE)
  return(sample_metadata)
}


#' Top Taxonomy Density Heatmap for QC
#'
#' Map most abundant ASV per sample to plate layout to investigate lab-specific errors/trends
#'
#' @param ps_obj Phyloseq object as INPUT
#' @param omit_cutoff integer, omit samples with less than n reads (default=100). This is important for color gradient to make sense.
#'
#' @returns ggplot2::geom_tile()
#' @export
#'
#' @examples
#' #To add
qcplot.plate_heatmap_toptaxa = function(ps_obj, omit_cutoff = 100){
  # plot top taxa over the plate layout
  sample_metadata = clean.ps_sample_sheet(ps_obj)

  taxa_df = phyloseq::otu_table(ps_obj)
  # Collect data on top ASV of ps_obj
  top_asv_otu_tab = as.data.frame(apply(taxa_df, MARGIN = 2, FUN =  which.max))
  colnames(top_asv_otu_tab) = c("TOP_ASV")
  # top_asv_otu_tab$TOP_ASV_SPECIES = as.data.frame(tax_table(ps_obj))[top_asv_otu_tab$TOP_ASV,][,"species"]
  top_asv_otu_tab$TOP_ASV_COUNT = apply(taxa_df, MARGIN = 2, FUN = max)
  top_asv_otu_tab$TOP_ASV_PROP = top_asv_otu_tab$TOP_ASV_COUNT / colSums(taxa_df)

  xx = merge(sample_metadata, top_asv_otu_tab, by=0,all.x = T)

  xx$TOP_ASV = as.factor(xx$TOP_ASV)
  # xx$TOP_ASV_ALPHA = (car::logit(xx$TOP_ASV_PROP, adjust = 0.025) + 2) / 6

  # plot(xx$TOP_ASV_PROP, xx$TOP_ASV_ALPHA)

  # remove NA samples and remove samples with low abundant top ASV (noise)
  xx_p = xx[!is.na(xx$TOP_ASV) & xx$TOP_ASV_COUNT > omit_cutoff,]

  plate_layout_plot = ggplot2::ggplot( data=xx_p, ggplot2::aes(
      x = .data$WELL_COLUMN, y=.data$WELL_ROW, fill=.data$TOP_ASV, alpha=.data$TOP_ASV_PROP)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~.data$PLATE, ncol = 3, scales = "free") +
    ggplot2::scale_y_discrete(limits=rev) +
    ggplot2::theme_classic() + ggplot2::theme(legend.position="none") +
    ggplot2::ggtitle("Most abundant ASV per well",
            subtitle = paste("Colors correspond to top ASV, opacity to top ASV proportion in the sample. Samples with top ASV <", omit_cutoff, "reads are omitted")
    )

  return(plate_layout_plot)

}

#' Total Readcount Density Heatmap for QC
#'
#' Map readcounts per sample to plate layout to investigate lab-specific errors/trends
#' Recommended to first make a appropriate phyloseq object of the data you want to investigate
#' Example: PAC bleed-over analysis -> input phyloseq is of all exotic ASVs over all samples.
#'
#' @param ps_obj Phyloseq object as INPUT
#' @param omit_cutoff integer, omit samples with less than n reads (default=100). This is important for color gradient to make sense.
#'
#' @returns ggplot2::geom_tile()
#' @export
#'
#' @examples
#' #To add
qcplot.plate_heatmap_readcount = function(ps_obj, omit_cutoff = 100){
  # plot top taxa over the plate layout
  sample_metadata = clean.ps_sample_sheet(ps_obj)

  tot_count = data.frame("total_reads" = phyloseq::sample_sums(ps_obj))

  xx_p = merge(sample_metadata, tot_count, by=0,all.x = T)

  xx_p[xx_p$total_reads < omit_cutoff] = 0

  plate_layout_plot = ggplot2::ggplot(
    data=xx_p,
    ggplot2::aes(x = .data$WELL_COLUMN, y=.data$WELL_ROW, fill=.data$total_reads)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~.data$PLATE, ncol = 3, scales = "free") +
    ggplot2::scale_fill_gradient(high="red", low="green") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("Number of reads per well",
            subtitle = paste("Colors correspond total read number. Samples with <", omit_cutoff, "reads are omitted")
    )

  return(plate_layout_plot)

}
