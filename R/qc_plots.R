################
### QC Plots ###
################

clean.ps_sample_sheet = function(ps_obj){
  sample_metadata = as.data.frame(phyloseq::sample_data(ps_obj))
  sample_metadata$PLATE = factor(gsub("R.*$", "", sample_metadata$PLATE_NOTES))
  sample_metadata$PLATE_REPLICATE = factor(paste0("R",gsub("^.*R", "", sample_metadata$PLATE_NOTES)))
  sample_metadata$PLATE_NOTES = factor(
    sample_metadata$PLATE_NOTES,
    levels = paste(sort((unique(sample_metadata$PLATE_NOTES)))))
  sample_metadata$WELL_COLUMN = factor(
    sample_metadata$WELL_COLUMN,
    levels = paste(sort(as.integer(unique(sample_metadata$WELL_COLUMN)))))
  sample_metadata$WELL_ROW = factor(
    sample_metadata$WELL_ROW,
    levels = paste(sort((unique(sample_metadata$WELL_ROW)))))
  sample_metadata$TECHNICAL_REPLICATE = factor(
    sample_metadata$TECHNICAL_REPLICATE,
    levels = paste(sort(as.integer(unique(sample_metadata$TECHNICAL_REPLICATE)))))

  return(sample_metadata)
}


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
  phyloseq::sample_data(ps_plate) = clean.ps_sample_sheet(ps_plate)
  out_list = c()
  # Create sample column to merge on
  phyloseq::sample_data(ps_plate)$X = paste(
    phyloseq::sample_data(ps_plate)$WELL_COLUMN,
    phyloseq::sample_data(ps_plate)$PLATE_REPLICATE, sep = "_")

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
  phyloseq::plot_bar(ps_plate_topX, x="WELL_COLUMN", fill="species", facet_grid=~PLATE_REPLICATE)

  ## MERGE the columns into samples
  ps_plate_merged = phyloseq::merge_samples(x=ps_plate, "X")

  ## Barplot
  phyloseq::plot_bar(ps_plate_merged, x="PLATE_REPLICATE", fill="custom_taxon", facet_grid=~WELL_COLUMN)

  ps_plate_merged_topX = phyloseq::prune_taxa(
    (phyloseq::tax_table(ps_plate_merged)[, "species"] %in% names(sort(order.sum, TRUE))[1:topX]),
    ps_plate_merged)
  ps_plate_merged_topX.frac = phyloseq::transform_sample_counts(ps_plate_merged_topX, function(OTU) OTU/sum(OTU)*100)

  pb = phyloseq::plot_bar(ps_plate_merged_topX.frac, x="PLATE_REPLICATE", fill="species", facet_grid=~WELL_COLUMN) +
    ggplot2::theme(legend.position="none")

  out_list$barplot = pb

  ## Ordination
  ps_plate_merged.frac.ordination = phyloseq::ordinate(ps_plate_merged_topX.frac, "PCoA", "bray")

  phyloseq::sample_data(ps_plate_merged_topX.frac)$WELL_COLUMN = as.factor(phyloseq::sample_data(ps_plate_merged_topX.frac)$WELL_COLUMN)
  phyloseq::sample_data(ps_plate_merged_topX.frac)$PLATE_REPLICATE = as.factor(phyloseq::sample_data(ps_plate_merged_topX.frac)$PLATE_REPLICATE)

  po = phyloseq::plot_ordination(ps_plate_merged_topX.frac, ps_plate_merged.frac.ordination,
                                 color="WELL_COLUMN", shape = "PLATE_REPLICATE") +
    ggplot2::geom_text(
      ggplot2::aes(label=.data$WELL_COLUMN), size = 5, nudge_y = -0.015)

  out_list$ordination_plot = po

  return(out_list)

}


#' [WRAPPER] Barplot and ordination for QC on possible strip-swaps
#'
#' Taxon-barplot and ordination with data grouped per plate and per column to investigate lab-specific errors/trends
#'
#' @param ps_obj Phyloseq object as INPUT
#' @param topX integer, topX taxa per COLUMN on the plate (default=25)
#' @param out_path Directory to store the plots in. (If this argument is used the default list will not be returned.)
#'
#' @returns List with n (plates) objects, for each a barplot and ordination
#' @export
#'
#' @examples
#' #To add
qcplot.plate_column_asv_analysis = function(ps_obj, topX=25, out_path=NA){

  sample_metadata = clean.ps_sample_sheet(ps_obj)
  platesx = sample_metadata$PLATE
  pl_list = list()

  for (pl_i in levels(platesx)){
    sub_ps = phyloseq::prune_samples(sample_metadata$PLATE == pl_i, ps_obj)
    sub_ps = phyloseq::prune_taxa(phyloseq::taxa_sums(sub_ps) > 0, sub_ps)

    # make plots
    res = gendiver::qcplot.plate_column_asv_barplot(sub_ps, topX)
    res$barplot = res$barplot + ggplot2::ggtitle(paste0("Plate ", pl_i, " top ", topX, " species per column/replicate" ))
    res$ordination_plot = res$ordination_plot + ggplot2::ggtitle(paste0("Plate ", pl_i, " top ", topX, " species per column/replicate" ))

    if (!is.na(out_path)){
      # save
      out_dir_pproj = file.path(out_path, pl_i)
      dir.create(out_dir_pproj, showWarnings = FALSE, recursive = T)

      grDevices::pdf(file=file.path(out_dir_pproj, paste0(pl_i, "_barplot.pdf")), width=12,height=8)
      print(res$barplot)
      grDevices::dev.off()

      grDevices::pdf(file=file.path(out_dir_pproj, paste0(pl_i, "_ordination.pdf")), width=12,height=8)
      print(res$ordination_plot )
      grDevices::dev.off()
    } else {
      pl_list[[pl_i]] = res
    }

  }

  if (is.na(out_path)){
    return(pl_list)
  } else {
    return()
  }


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
  xx_p = xx[!is.na(xx$TOP_ASV) & xx$TOP_ASV_COUNT >= omit_cutoff,]

  # Check if there is data remaining
  if (nrow(xx_p) == 0){
    errorCondition("No data. Try lowering the omit_cutoff, or check your input object.")
    return()
    }

  plate_layout_plot = ggplot2::ggplot( data=xx_p, ggplot2::aes(
      x = .data$WELL_COLUMN, y=.data$WELL_ROW, fill=.data$TOP_ASV, alpha=.data$TOP_ASV_PROP)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~.data$PLATE_NOTES, ncol = 3, scales = "free") +
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
#' @param transform (default="identity") See \link[ggplot2]{scale_continuous}, e.g. "log2", "log10", "sqrt", ...
#'
#' @returns ggplot2::geom_tile()
#' @export
#'
#' @examples
#' #To add
qcplot.plate_heatmap_readcount = function(ps_obj, omit_cutoff = 100, transform="identity"){
  # plot top taxa over the plate layout
  sample_metadata = clean.ps_sample_sheet(ps_obj)

  tot_count = data.frame("total_reads" = phyloseq::sample_sums(ps_obj))

  xx_p = merge(sample_metadata, tot_count, by=0,all.x = T)
  xx_p = xx_p[xx_p$total_reads >= omit_cutoff,]
  # xx_p$total_reads[xx_p$total_reads < omit_cutoff] = NA
  # Check if there is data remaining
  if (nrow(xx_p) == 0){
    errorCondition("No data. Try lowering the omit_cutoff, or check your input object.")
    return()
  }

  plate_layout_plot = ggplot2::ggplot(
    data=xx_p,
    ggplot2::aes(x = .data$WELL_COLUMN, y=.data$WELL_ROW, fill=.data$total_reads)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~.data$PLATE_NOTES, ncol = 3) +
    ggplot2::scale_fill_gradient(high="red", low="lightblue", transform=transform) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_discrete(limits=rev) +
    ggplot2::ggtitle(paste0("Number of reads (",transform,") per well"),
            subtitle = paste0("Colors correspond total read number (",transform,"). Samples with < ", omit_cutoff, " reads are omitted")
    )

  return(plate_layout_plot)

}

#' Basic exploration plots
#'
#' Plots of basic statistics of th operational data
#'
#' @param ps_obj Phyloseq object as INPUT
#'
#' @returns List with ggplot objects
#' @export
#'
#' @examples
#' #To add
qcplot.basic_exploration = function(ps_obj){
  ssdata = phyloseq::sample_data(ps_obj)
  PROJ_NAME = paste0(paste0(unique(ssdata$RUN_CODE), collapse = "_"), "_", paste0(unique(ssdata$PRIMERS), collapse = "_"))

  ## BASIC EXPLORATION PLOTS
  tax_df = data.frame(phyloseq::tax_table(ps_obj))
  otu_df = data.frame(phyloseq::otu_table(ps_obj))
  seqs = phyloseq::refseq(ps_obj)

  # if unclass propagated, revert
  tax_df = data.frame(apply(tax_df, MARGIN = 2, gsub, pattern=".*unclassif.*", replacement=NA))

  annot_overview = colSums(!is.na(tax_df))
  annot_overview = annot_overview[names(annot_overview) != "custom_taxon"]
  annot_overview = c(nrow(tax_df), annot_overview)
  names(annot_overview)[1] = 'TOTAL'

  b1 = ggplot2::ggplot(data=data.frame("tax"=names(annot_overview), "number_otu" = annot_overview)) +
    ggplot2::geom_col(ggplot2::aes(x=factor(.data$tax, levels=names(annot_overview)), y = .data$number_otu)) +
    ggplot2::ggtitle(paste0(PROJ_NAME,": Taxonomy assignments overview - OTUs")) +
    ggplot2::ylab("Number of OTUs") + ggplot2::xlab("")

  tax_df_rc = tax_df
  otu_tot_rc = rowSums(otu_df)
  tax_df_rc$read_count = NA
  tax_df_rc[names(otu_tot_rc), c("read_count")] = otu_tot_rc
  rm(otu_tot_rc)

  # remove taxonomy not in OTU table
  tax_df_rc = tax_df_rc[!is.na(tax_df_rc$read_count),]

  tb_rc = sapply(colnames(tax_df), function(x) sum(tax_df_rc$read_count[!is.na(tax_df_rc[[x]])]))
  annot_overview2 = tb_rc[names(tb_rc) != "custom_taxon"]
  annot_overview2 = c(sum(tax_df_rc$read_count), annot_overview2)
  rm(tax_df_rc)

  names(annot_overview2)[1] = 'TOTAL'

  b2 = ggplot2::ggplot(data=data.frame("tax"=names(annot_overview2), "number_reads" = annot_overview2)) +
        ggplot2::geom_col(ggplot2::aes(x=factor(.data$tax, levels=names(annot_overview2)), y = .data$number_reads)) +
    ggplot2::ggtitle(paste0(PROJ_NAME,": taxonomy assignments overview - READS")) +
    ggplot2::ylab("Number of reads") + ggplot2::xlab("")

  # print(annot_overview2, row.names = F)

  # Sequence length
  # remove sequences not in OTU table

  seqs = seqs[rownames(otu_df)]
  asv_lengths = Biostrings::width(seqs)

  h1 = ggplot2::ggplot(data=data.frame("asv_lengths" = asv_lengths)) +
    ggplot2::geom_histogram(ggplot2::aes(x=.data$asv_lengths), binwidth = 1) +
    ggplot2::ggtitle(paste0(PROJ_NAME,": OTU seq length")) +
    ggplot2::ylab("Number of OTUs") + ggplot2::xlab("Length (bp)")

  reads_length = data.frame("reads"=rowSums(otu_df), "length"=asv_lengths)
  rm(asv_lengths)

  len_reads.df = stats::aggregate(reads~length, reads_length, sum)
  len_reads.df$cumulative = cumsum(len_reads.df$reads)
  rm(reads_length)

  p1 = ggplot2::ggplot(data=len_reads.df) +
    ggplot2::geom_point(ggplot2::aes(x=.data$length, y=.data$reads)) +
    ggplot2::ggtitle(paste0(PROJ_NAME, ": number of reads per OTU seq length")) +
    ggplot2::ylab("Number of reads") + ggplot2::xlab("OTU seq length")


  p2 = ggplot2::ggplot(data=len_reads.df) +
    ggplot2::geom_point(ggplot2::aes(x=.data$length, y=.data$cumulative)) +
    ggplot2::ggtitle(paste0(PROJ_NAME, ": number of reads per OTU seq length")) +
    ggplot2::ylab("Cumulative sum reads") + ggplot2::xlab("OTU seq length")

  h2 = ggplot2::ggplot(data=data.frame("x" = rowSums(otu_df))) +
    ggplot2::geom_histogram(ggplot2::aes(x=.data$x), bins = 100) +
    ggplot2::ggtitle(paste0(PROJ_NAME, ": number of reads per OTU")) +
    ggplot2::ylab("Frequency") + ggplot2::xlab("Total number of reads")

  h3 = ggplot2::ggplot(data=data.frame("x" = colSums(otu_df))) +
    ggplot2::geom_histogram(ggplot2::aes(x=.data$x), bins = 100) +
    ggplot2::ggtitle(paste0(PROJ_NAME, ": number of reads per sample")) +
    ggplot2::ylab("Frequency") + ggplot2::xlab("Total number of reads")

  out_list = list(
    "readcount.otu" = h2,
    "readcount.sample" = h3,
    "taxonomy.otu" = b1,
    "taxonomy.reads" = b2,
    "length.otu" = h1,
    "length.reads" = p1,
    "length.reads.cumsum" = p2
    )

  return(out_list)

}


#' Check taxonomic assignment levels per sample
#'
#' Per sample the read % attributed to each taxonomic level
#'
#' @param ps_obj Phyloseq object as INPUT
#'
#' @returns ggplot2 object
#' @export
#'
#' @examples
#' #To add
qcplot.sample_barplot_taxonomic_coverage = function(ps_obj){
  tax_df = data.frame(phyloseq::tax_table(ps_obj))
  tax_df_no_custom=tax_df[,!grepl("custom",colnames(tax_df))]
  max_rank = function(x) if (sum(!is.na(x)) > 0) colnames(tax_df_no_custom)[max(which(!is.na(x)))] else "_root"
  base_rank = apply(tax_df_no_custom, FUN=max_rank, MARGIN=1)
  table(base_rank)

  otu_df_rank = data.frame(phyloseq::otu_table(ps_obj))
  otu_df_rank$max_tax_rank = base_rank[row.names(otu_df_rank)]
  rm(base_rank)

  otu_df_rank = reshape2::melt(otu_df_rank, id.vars=c("max_tax_rank"), value.name = "read_count", variable.name="USI")
  otu_df_rank = stats::aggregate(data=otu_df_rank, read_count ~ max_tax_rank + USI, sum)

  # rank on non-lowest level (aka how many reads not species level annot)
  llv = utils::tail(colnames(tax_df),1)
  usilvl = otu_df_rank

  if (length(unique(otu_df_rank$max_tax_rank)) != llv) {
    usilvl = stats::aggregate(data=otu_df_rank[otu_df_rank$max_tax_rank != llv,], read_count ~ USI, sum)
    }
  usilvl = usilvl$USI[order(usilvl$read_count, decreasing = F)]

  otu_df_rank$USI = factor(otu_df_rank$USI, levels=usilvl)
  otu_df_rank$max_tax_rank = factor(otu_df_rank$max_tax_rank, levels = rev(c("_root", colnames(tax_df_no_custom))))
  rm(tax_df_no_custom)

  pss_plot = ggplot2::ggplot(
    data=otu_df_rank, ggplot2::aes(x=.data$USI, y=.data$read_count, fill = .data$max_tax_rank)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::coord_flip() + ggplot2::ggtitle("Total sample reads by maximum taxonomic level") +
    ggplot2::theme(legend.text = ggplot2::element_text(size=12),
          legend.key.size=ggplot2::unit(0.4, "cm"),
          legend.position="top") +
    ggplot2::scale_y_continuous(expand = c(0,0))

  return(pss_plot)
}
