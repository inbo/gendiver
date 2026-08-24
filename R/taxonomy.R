############
# Taxonomy #
############


#' Read OBITools3 ECOTAG taxonomy output into dataframe
#'
#' Reading and general cleaning/pruning of OBITOOLS3 ecotag output
#'
#' @param tax_table_path Path to taxonomy table file (can also be link to Googe Drive)
#'
#' @returns Dataframe with taxonomy information (columns) per ASV (rows)
#' @export
#'
#' @examples
#' #To add
read.taxonomy_obitools3 = function(tax_table_path){
  tax_df = read.table_gdrive(tax_table_path, sep='\t', header = T)
  rownames(tax_df) = gsub(pattern = ";size=\\d*", replacement = "", x = tax_df$ID)
  return(tax_df)
}

#' Filter ASV data for eDNA-water (12S primers)
#'
#' Basic pruning of OBITools3 taxonomy dataframe (for RIAZ primer), formatting, adding custom taxon for team gendiv,
#' adding unclassified labels, cleaning and converting to phyloseq tax_table() compatible dataframe
#'
#' @param tax_df Dataframe of (obitools3) taxonomy table (see \link[gendiver]{read.taxonomy_obitools3})
#' @param min_BEST_ID Minimum identity match for ASV to be retained. Default=1, retain only ASVs with perfect match to RefDB
#' @param max_BEST_ID Maximum identity match for ASV to be retained. Default=1, retain only ASVs with perfect match to RefDB
#' @param fillNA See \link[gendiver]{propagate.unclassified_taxonomy} (default=TRUE)
#'
#' @returns Dataframe with ASV selection
#' @export
#'
#' @examples
#' #To add
select.taxonomy_obitools3 = function(tax_df, min_BEST_ID=1, max_BEST_ID=1, fillNA=TRUE){

  # Filter taxa based on BEST_IDENTITY
  filt_nr_asv = nrow(tax_df[! (max_BEST_ID >= tax_df$BEST_IDENTITY & tax_df$BEST_IDENTITY >= min_BEST_ID),])
  tot_nr_asv = nrow(tax_df)

  message(paste0(
    "[INFO] Using min_BEST_ID = ", min_BEST_ID, " and max_BEST_ID = ", max_BEST_ID,
    "\n[INFO] This removes ", filt_nr_asv," (", round(filt_nr_asv/tot_nr_asv*100, 2),"%) OTUs"))

  table.filt=tax_df[max_BEST_ID >= tax_df$BEST_IDENTITY & tax_df$BEST_IDENTITY >= min_BEST_ID,]

  # only keep species names
  taxonomy<-table.filt[,c("phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")]
  colnames(taxonomy) = c("phylum", "class", "order", "family", "genus", "species")

  # Add custom taxonomy and propagate taxonomy to NA
  taxonomy_tab = add_custom_taxonomy(taxonomy, fillNA=fillNA)

  return(taxonomy_tab)
}

#' Read and Filter ASV data for eDNA-water (12S primers)
#'
#' Wrapper for reading and cleaning obitools3 tax table to phyloseq tax_table() compatible dataframe
#'
#' @param tax_table_path Path to taxonomy table file (can also be link to Googe Drive)
#' @param min_BEST_ID Minimum identity match for ASV to be retained. Default=1, retain only ASVs with perfect match to RefDB
#' @param max_BEST_ID Maximum identity match for ASV to be retained. Default=1, retain only ASVs with perfect match to RefDB
#'
#' @returns Dataframe with ASV selection
#' @export
#'
#' @examples
#' #To add
pull.taxonomy_obitools3 = function(tax_table_path, min_BEST_ID=1, max_BEST_ID=1){
  tax_df = gendiver::read.taxonomy_obitools3(tax_table_path)
  tax_out = gendiver::select.taxonomy_obitools3(tax_df, min_BEST_ID, max_BEST_ID)
  return(tax_out)
}

# Add custom taxonomy
add_custom_taxonomy = function(taxonomy_df, fillNA=F){
  exp_cols = c("phylum", "class", "order", "family", "genus", "species")
  if (sum( ! exp_cols %in% base::colnames(taxonomy_df)) != 0) {
    warning(paste("Colnames not expected taxonomy-format:", paste(exp_cols, collapse = ",")))
    return(taxonomy_df)
  }
  taxonomy_df$custom_taxon = taxonomy_df$class
  taxonomy_df$custom_taxon[taxonomy_df$genus == "Homo"] = "Human"
  taxonomy_df$custom_taxon[taxonomy_df$custom_taxon == "Mammalia"] = "Non-human mammals"
  taxonomy_df$custom_taxon[taxonomy_df$custom_taxon == "Aves"] = "Birds"
  taxonomy_df$custom_taxon[taxonomy_df$custom_taxon == "Actinopteri"] = "Fish"
  taxonomy_df$custom_taxon[taxonomy_df$custom_taxon == "Hyperoartia"] = "Fish"
  taxonomy_df$custom_taxon[taxonomy_df$custom_taxon == "Amphibia"] = "Amphibians"
  taxonomy_df$custom_taxon[is.na(taxonomy_df$custom_taxon)]="Other"

  if (fillNA){
    # change NA names to unclassified
    taxonomy_df = propagate.unclassified_taxonomy(taxonomy_df)
    }

  p_i = which(colnames(taxonomy_df) == "phylum")
  taxonomy_df = taxonomy_df[, c(1:p_i, length(colnames(taxonomy_df)), (p_i+1):(length(colnames(taxonomy_df))-1))]

  return(taxonomy_df)
}

# clean unrsolved
clean_NA_taxonomies = function(tax_df, unk_pattern="unresolv|unknown|unclassif|undefin|uncertain|unassign"){
  mask_unk = sapply(tax_df, function(x) grepl(unk_pattern,x))
  tax_df[mask_unk] = NA

  # Extra
  tax_df[tax_df == "NA"] = NA
  tax_df[tax_df == ""] = NA

  return(tax_df)
}

# (Try to) Clean taxonomies
clean_taxonomies = function(tabraw){
  tabraw = lapply(tabraw, gsub, pattern = "_\\d+", replacement = "")
  tabraw = lapply(tabraw, gsub, pattern = "species_", replacement = "")
  tabraw = lapply(tabraw, gsub, pattern = "genus_", replacement = "")
  tabraw = lapply(tabraw, gsub, pattern = "family_", replacement = "")
  tabraw = lapply(tabraw, gsub, pattern = "order_", replacement = "")
  tabraw = lapply(tabraw, gsub, pattern = "class_", replacement = "")
  tabraw = lapply(tabraw, gsub, pattern = "phylum_", replacement = "")
  tabraw = lapply(tabraw, gsub, pattern = ".* sp\\. .*", replacement = "")

  out = as.data.frame(tabraw)
  out = clean_NA_taxonomies(out)
  return(out)
}

#' Parse VSEARCH SINTAX taxonomy string to dataframe
#'
#'
#'@description
#'
#' Each taxonomic identifier must start with an indication of the rank by one of the letters d (for domain) k (kingdom), p (phylum), c (class), o (order), f (family), g (genus), s (species), or t (strain). The letter is followed by a colon (:) and the name of that rank. Commas and semicolons are not allowed in the name of the rank. Non-ascii characters should be avoided in the names.
#'
#' Example:
#'
#' >X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli,t:str._K-12_substr._MG1655
#'
#'
#' @param x VSEARCH SINTAX taxonomy string
#' @param is_pred (SINTAX output specific) If bootstrap values are part of the tax_string. Default = FALSE.
#'
#' @returns Data.frame with (1) row and (n = tax_levels) columns
#' @export
#'
#' @examples
#' tax_string = "d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli,t:str._K-12_substr._MG1655"
#' tax_string |> parse.utax_string()
parse.utax_string = function(x, is_pred=F){
  # Return empty s to fit in when rbind
  if (x==""){
    return(data.frame("species"=NA))
  }
  # split taxon-levels ,
  xspl = unlist(strsplit(x, ",*[dpcofgst]:"))[-1]

  # name to first letter
  names(xspl) = c("domain", "phylum", "class", "order", "family", "genus", "species", "strain")[1:length(xspl)]
  xcol = as.data.frame.list(xspl)

  if (is_pred){
    mycols = names(xspl)
    xcol_bs = stringr::str_extract(xspl, "\\d+\\.\\d+")
    # remove pred
    xcol_tax = vapply(strsplit(as.character(xcol[1,]), "\\("), `[`, 1, FUN.VALUE=character(1))
    # update xcol
    xcol = c(xcol_tax, xcol_bs)
    names(xcol) = c(mycols, paste0(mycols, "_bs"))
    xcol = as.data.frame.list(xcol)
  }

  return(xcol)

}

#' Read VSEARCH SINTAX taxonomy output into dataframe
#'
#' Reading and general cleaning/pruning of VSEARCH SINTAX output
#'
#' @param sintax_path Path to taxonomy table file
#' @param read_supported_only Only return bootstrap supported annotations (default=TRUE)
#' @param add_source Add basename of `sintax_path` to output dataframe (default=FALSE)
#' @param clean_tax Try to remove non-informative text in taxonomic annotations (default=FALSE)
#' @param add_custom  Add a custom taxonomy-level between Phylum and Class with common names for basic quality control (INBO specific) (default=FALSE)
#'
#' @returns Dataframe with taxonomy information (columns) per OTU (rows)
#' @export
#'
#' @examples
#' #To add
read.taxonomy_sintax = function(sintax_path, read_supported_only=TRUE, add_source=FALSE, clean_tax=FALSE, add_custom=F){
  # Read input
  rtab = read.table_gdrive(
    sintax_path, sep = "\t", header = F,
    col.names = c("QUERY_LABEL","TAX_PRED", "STRAND", "TAX_PRED_CUTOFF")
    )

  # order the input on string length to increase chance of full taxonomy template
  rtab = rtab[order(nchar(rtab$TAX_PRED_CUTOFF), decreasing = T),]
  if (read_supported_only) {
    ## Read supported
    tax_supp = data.table::rbindlist(lapply(rtab$TAX_PRED_CUTOFF, parse.utax_string), fill = T)
  } else {
    ## Read complete prediction
    tax_supp = data.table::rbindlist(lapply(rtab$TAX_PRED, parse.utax_string, is_pred = T), fill = T)
  }

  if (clean_tax) {tax_supp = clean_taxonomies(tax_supp)}
  if (add_source) {tax_supp$SOURCE_FILE = basename(sintax_path)}
  if (add_custom){
    # Add full names
    colnames(tax_supp) = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    tax_supp = add_custom_taxonomy(tax_supp)
  }
  tax_supp = as.data.frame(tax_supp)
  # Add OTU-IDs back to data (and remove possible ';size=' info)
  rownames(tax_supp) = gsub(";size=\\d*", "" , rtab$QUERY_LABEL)

  return(tax_supp)
}



paste_string_listOBI = function(x){
  return(paste0("['",paste0(unique(x), collapse="', '"), "']"))
}

paste_string_list = function(x){
  return(paste0("",paste0(unique(x), collapse=" / "), ""))
}

#' Read VSEARCH BLAST6 format to dataframe
#'
#' Reading and basic formatting of vsearch --blast6 output
#'
#'
#' @param ublast6_path Path to ublast6 output file (tsv). Can be file_path or google-drive URL
#'
#' @returns Dataframe with taxonomy information (columns) of target per query hit
#' @export
#'
#' @examples
#' # to add
read.taxonomy_usearch_blast6 = function(ublast6_path){
  usearch_blast6 = gendiver::read.table_gdrive(ublast6_path, sep="\t", header=F)
  colnames(usearch_blast6) = c("QUERY_LABEL", "TARGET_LABEL", "ID", "ALNLEN", "MISM", "OPENS", "QLO", "QHI", "TLO", "THI", "EVALUE", "BITS")
  usearch_blast6$QUERY_LABEL = gsub(";.*", "", usearch_blast6$QUERY_LABEL)
  usearch_blast6$TARGET_TAX = gsub(".*tax=", "", usearch_blast6$TARGET_LABEL)
  usearch_blast6_tax_parsed = data.table::rbindlist(lapply(usearch_blast6$TARGET_TAX, parse.utax_string), fill = T)
  usearch_blast6 = cbind(usearch_blast6,usearch_blast6_tax_parsed)
  usearch_blast6$TARGET_LABEL = gsub(";tax=.*", "", usearch_blast6$TARGET_LABEL)
  usearch_blast6$TARGET_TAX = NULL
  return(usearch_blast6)
}

#' Read VSEARCH LCA format to dataframe
#'
#' Reading and basic formatting of vsearch --lca_out output
#'
#'
#' @param lca_path Path to lca output file (tsv). Can be file_path or google-drive URL
#'
#' @returns Dataframe with taxonomy information (columns) of target per query hit
#' @export
#'
#' @examples
#' # to add
read.taxonomy_usearch_lca = function(lca_path){
  usearch_lca = read.table_gdrive(lca_path, sep="\t", header=F)
  usearch_lca = cbind(usearch_lca, data.table::rbindlist(lapply(usearch_lca[[2]], parse.utax_string), fill = T))
  usearch_lca$V1 = gsub(";.*", "", usearch_lca$V1)
  usearch_lca$V2 = NULL
  colnames(usearch_lca)[1] = "QUERY_LABEL"
  return(usearch_lca)
}

#' Merge VSEARCH BLAST6 and LCA format
#'
#'@description
#' Merging BLAST6 (multiple RefDB hits per query/OTU) with LCA output (1 taxonomy per query/OTU) to generate OBITools like taxonomy table.
#'
#' BLAST6 data of multiple hits are collapsed and represented as a single element per query/OTU row.
#'
#'
#' @param usearch_blast6 Path/URL to, or data.frame with, the BLAST6 output
#' @param usearch_lca Path/URL to, or data.frame with, the LCA output
#'
#' @returns Dataframe with aggregated BLAST6 taxonomy information (columns) of target OTU
#' @export
#'
#' @examples
#' # to add
merger.taxonomy_usearch = function(usearch_blast6, usearch_lca) {
  # allow path/gdrive-id as input, else assume is dataframe
  if (! is.data.frame(usearch_blast6)) {
    usearch_blast6 = read.taxonomy_usearch_blast6(usearch_blast6)
  }
  if (! is.data.frame(usearch_lca)) {
    usearch_lca = read.taxonomy_usearch_lca(usearch_lca)
  }

  usearch_ID_agg = aggregate(data=usearch_blast6, ID ~ QUERY_LABEL, FUN = max )
  usearch_target_IDs = aggregate(data=usearch_blast6, TARGET_LABEL ~ QUERY_LABEL, FUN = paste_string_listOBI )
  usearch_target_species = aggregate(data=usearch_blast6, species ~ QUERY_LABEL, FUN = paste_string_list)
  usearch_agg_info = cbind(usearch_ID_agg, usearch_target_species[-1], usearch_target_IDs[-1] )
  colnames(usearch_agg_info)= c("QUERY_LABEL", "ID", "TARGET_LABEL_species", "TARGET_LABEL")

  usearch_lca$otu_order = row.names(usearch_lca)

  usearch_combi_res = merge(usearch_lca, usearch_agg_info, by="QUERY_LABEL", all=T)
  usearch_combi_res = usearch_combi_res[order(as.numeric(usearch_combi_res$otu_order)),]
  usearch_combi_res$otu_order = NULL
  row.names(usearch_combi_res) = usearch_combi_res$QUERY_LABEL

  return(usearch_combi_res)
}

#' Read DECIPHER IdTaxa() taxonomy output into dataframe
#'
#' Reading and general cleaning/pruning of DECIPHER IdTaxa() output
#'
#' @param tax_tab_path Path to taxonomy CSV table file (can also be link to Googe Drive)
#' @param read_supported_only Only positive annotations (default=TRUE)
#' @param add_source Add basename of `tax_tab_path` to output dataframe (default=FALSE)
#' @param clean_tax Try to remove non-informative text in taxonomic annotations (default=FALSE)
#'
#' @returns Dataframe with taxonomy information (columns) per OTU (rows)
#' @export
#'
#' @examples
#' #To add
read.taxonomy_idtaxa = function(tax_tab_path, read_supported_only=TRUE, add_source=FALSE, clean_tax=FALSE){

  taxranks = c("d", "p", "c", "o", "f", "g", "s")

  tabraw = read.table_gdrive(tax_tab_path, sep = ";", header = F)
  colnames(tabraw) = c("QUERY_LABEL", taxranks)

  # WHY did we invent a new format (mixed \t and ;sv) for the taxonomy table of IDTAXA ????
  # ADjust this in the code!!!
  # for now, remove everythin after SWARM name

  tabraw$QUERY_LABEL = gsub(pattern = "\t.*", "", tabraw$QUERY_LABEL)

  if (read_supported_only) {
    # Remove predictions
    tabraw[taxranks] = lapply(tabraw[taxranks], gsub, pattern = "\\(\\d+\\.?\\d+%*\\)", replacement = "")
    # remove unclassified taxonomies
    tabraw[taxranks] = lapply(tabraw[taxranks], gsub, pattern = ".*unclassified.*", replacement = "")
    tabraw[grepl("unclassified", tabraw)] = NA
  }

  # Remove taxID
  tabraw[taxranks] = lapply(tabraw[taxranks], gsub, pattern = "_\\d+", replacement = "")
  tabraw[tabraw == ""]=NA

  if (clean_tax){tabraw = clean_taxonomies(tabraw)}
  if (add_source){ tabraw$SOURCE_FILE=basename(tax_tab_path)}

  return(tabraw)
}

#' Read BOLDigger3 taxonomy output into dataframe
#'
#' Read BOLDigger3 identification_result.xlsx output
#'
#' @param tax_tab_path Path to taxonomy XLSX table file
#' @param add_source Add basename of `tax_tab_path` to output dataframe (default=FALSE)
#' @param clean_tax Try to remove non-informative text in taxonomic annotations (default=FALSE)
#'
#' @returns Dataframe with taxonomy information (columns) per OTU (rows)
#' @export
#'
#' @examples
#' #To add
read.taxonomy_boldigger3 = function(tax_tab_path, add_source = FALSE, clean_tax=FALSE){
  x_ = readxl::read_xlsx(tax_tab_path)
  x_ = as.data.frame(x_[, c(1:7)])
  colnames(x_) = c("QUERY_LABEL" ,"p"         ,   "c"       ,     "o"       ,     "f"         ,   "g"        ,    "s"  )
  x_[x_ == "no-match"] = NA
  x_$d = NA
  if (clean_tax){x_ = clean_taxonomies(x_)}
  if (add_source){x_$SOURCE_FILE=basename(tax_tab_path)}
  return(x_)
}

#' Complete NA identifications in taxonomic annotation results
#'
#' Complete (OBITools3) taxonomy by replacing NA values with last known taxon + suffix
#'
#' @param tax Taxonomy table
#' @param uncl_suffix Suffix to use in propagation (default="_unclassified")
#'
#' @returns Dataframe with taxonomy information (columns) per OTU (rows)
#' @export
#'
#' @examples
#' #To add
propagate.unclassified_taxonomy = function(tax, uncl_suffix="_unclassified"){
  tax = as.data.frame(tax)
  # take copy
  new_tax = tax
  #loop columns
  unclassified_tax = paste0(tax[,1], uncl_suffix)
  for (taxrank in colnames(new_tax)){
    # find rows that have no taxonomy
    unclassified_mask = is.na(new_tax[,taxrank])
    # propagate taxonomy these rows with no taxonomy
    new_tax[unclassified_mask , taxrank] = unclassified_tax[unclassified_mask]

    # update taxonomy
    unclassified_tax = new_tax[,taxrank]
    # only do + "_unclassified"to not already unclassified
    unclassified_tax[!unclassified_mask] = paste0(unclassified_tax[!unclassified_mask], uncl_suffix)

  }

  # return(tax_table(as.matrix(new_tax)))
  return(new_tax)
}
