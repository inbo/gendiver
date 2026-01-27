################
### QC Plots ###
################

# Taxon-barplot with data grouped per plate and per column to investigate lab-specific errors/trends
# Custom function for run 24009
plate_column_analysis_barplot = function(ps_plate, topX=25){
  # phyloseq object as INPUT to make barplot of topX taxa per COLUMN on the plate
  # Ideally of 1 plate

  out_list = c()
  # Create sample column to merge on
  sample_data(ps_plate)$X = paste(sample_data(ps_plate)$WELL_COLUMN, sample_data(ps_plate)$TECHNICAL_REPLICATE, sep = "_")

  # And taxglom for easier interpretation
  ps_plate = tax_glom(ps_plate, taxrank = "species")

  # Create relative abundance ps
  ps_plate.frac = transform_sample_counts(ps_plate, function(OTU) OTU/sum(OTU)*100)

  # Plot bars per column for topX taxa
  order.sum = tapply(taxa_sums(ps_plate), tax_table(ps_plate)[, "species"], sum, na.rm=TRUE)


  ps_plate_topX = prune_taxa((tax_table(ps_plate)[, "species"] %in% names(sort(order.sum, TRUE))[1:topX]), ps_plate)

  # Plot data
  plot_bar(ps_plate_topX, x="WELL_COLUMN", fill="species", facet_grid=~TECHNICAL_REPLICATE)


  ## MERGE the columns into samples
  ps_plate_merged = phyloseq::merge_samples(x=ps_plate, "X")

  ## Barplot
  plot_bar(ps_plate_merged, x="TECHNICAL_REPLICATE", fill="custom_taxon", facet_grid=~WELL_COLUMN)

  ps_plate_merged_topX = prune_taxa((tax_table(ps_plate_merged)[, "species"] %in% names(sort(order.sum, TRUE))[1:topX]), ps_plate_merged)
  ps_plate_merged_topX.frac = transform_sample_counts(ps_plate_merged_topX, function(OTU) OTU/sum(OTU)*100)

  pb = plot_bar(ps_plate_merged_topX.frac, x="TECHNICAL_REPLICATE", fill="species", facet_grid=~WELL_COLUMN) +
    theme(legend.position="none")

  out_list$barplot = pb

  ## Ordination
  ps_plate_merged.frac.ordination = ordinate(ps_plate_merged_topX.frac, "PCoA", "bray")

  sample_data(ps_plate_merged_topX.frac)$WELL_COLUMN = as.factor(sample_data(ps_plate_merged_topX.frac)$WELL_COLUMN)
  sample_data(ps_plate_merged_topX.frac)$TECHNICAL_REPLICATE = as.factor(sample_data(ps_plate_merged_topX.frac)$TECHNICAL_REPLICATE)

  po = plot_ordination(ps_plate_merged_topX.frac, ps_plate_merged.frac.ordination,
                       color="WELL_COLUMN", shape = "TECHNICAL_REPLICATE") + geom_text( aes(label=WELL_COLUMN), size = 5, nudge_y = -0.015)
  out_list$ordination_plot = po

  return(out_list)

}

# Map most abundant ASV to plate layout to investigate lab-specific errors/trends
# Custom function for run 24009
plot_plate_layout_top_ASV = function(ps_obj, sample_metadata_24009, omit_cutoff = 100){
  # plot top taxa over the plate layout

  taxa_df = otu_table(ps_obj)
  # Collect data on top ASV of ps_obj
  top_asv_otu_tab = as.data.frame(apply(taxa_df, MARGIN = 2, FUN =  which.max))
  colnames(top_asv_otu_tab) = c("TOP_ASV")
  # top_asv_otu_tab$TOP_ASV_SPECIES = as.data.frame(tax_table(ps_obj))[top_asv_otu_tab$TOP_ASV,][,"species"]
  top_asv_otu_tab$TOP_ASV_COUNT = apply(taxa_df, MARGIN = 2, FUN = max)
  top_asv_otu_tab$TOP_ASV_PROP = top_asv_otu_tab$TOP_ASV_COUNT / colSums(taxa_df)

  xx = merge(sample_metadata_24009, top_asv_otu_tab, by=0,all.x = T)

  xx$TOP_ASV = as.factor(xx$TOP_ASV)
  # xx$TOP_ASV_ALPHA = (car::logit(xx$TOP_ASV_PROP, adjust = 0.025) + 2) / 6

  # plot(xx$TOP_ASV_PROP, xx$TOP_ASV_ALPHA)

  # remove NA samples and remove samples with low abundant top ASV (noise)
  xx_p = xx[!is.na(xx$TOP_ASV) & xx$TOP_ASV_COUNT > omit_cutoff,]

  plate_layout_plot = ggplot(xx_p, aes(x = WELL_COLUMN, y=WELL_ROW, fill=TOP_ASV, alpha=TOP_ASV_PROP)) + geom_tile() +
    facet_wrap(~PLATE, ncol = 3, scales = "free") + scale_y_discrete(limits=rev) +
    theme_classic() + theme(legend.position="none") +
    ggtitle("Most abundant ASV per well",
            subtitle = paste("Colors correspond to top ASV, opacity to top ASV proportion in the sample. Samples with top ASV <", omit_cutoff, "reads are omitted")
    )

  return(plate_layout_plot)

}
