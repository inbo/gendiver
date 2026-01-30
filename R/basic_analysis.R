######################
### Basic Analysis ###
######################

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
