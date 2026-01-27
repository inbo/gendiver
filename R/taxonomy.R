############
# Taxonomy #
############

# Reading and general cleaning/pruning of OBITOOLS3 ecotag output
read_obitools_tax_data = function(tax_table_path){
  tax_df = read.csv(tax_table_path, sep='\t', header = T)
  rownames(tax_df) = gsub(pattern = ";size=\\d*", replacement = "", x = tax_df$ID)
  return(tax_df)
}

# Basic pruning of OBITools3 taxonomy dataframe (for RIAZ primer) and converting to phyloseq tax_table()
clean_and_filter_tax_table =function(tax_df, min_BEST_ID=1, max_BEST_ID=1){
  #only keep taxa with BEST_ID==1
  print(paste("using min_BEST_ID =", min_BEST_ID, "and max_BEST_ID =", max_BEST_ID))
  filt_nr_asv = nrow(tax_df[! (max_BEST_ID >= tax_df$BEST_IDENTITY & tax_df$BEST_IDENTITY >= min_BEST_ID),])
  tot_nr_asv = nrow(tax_df)

  print(paste0("This removes ", filt_nr_asv," (", round(filt_nr_asv/tot_nr_asv*100, 2),"%) ASVs"))

  table.filt=tax_df[max_BEST_ID >= tax_df$BEST_IDENTITY & tax_df$BEST_IDENTITY >= min_BEST_ID,]

  # only keep species names
  taxonomy<-table.filt[,c("phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")]
  colnames(taxonomy) = c("phylum", "class", "order", "family", "genus", "species")

  taxonomy$custom_taxon = taxonomy$class
  taxonomy$custom_taxon[taxonomy$genus == "Homo"] = "Human"
  taxonomy$custom_taxon[taxonomy$custom_taxon == "Mammalia"] = "non-human mammal"
  taxonomy$custom_taxon[taxonomy$custom_taxon == "Aves"] = "Bird"
  taxonomy$custom_taxon[taxonomy$custom_taxon == "Actinopteri"] = "Fish"
  taxonomy$custom_taxon[taxonomy$custom_taxon == "Hyperoartia"] = "Fish"
  taxonomy$custom_taxon[taxonomy$custom_taxon == "Amphibia"] = "Amphibians"

  taxonomy = taxonomy[, c("phylum","custom_taxon", "class", "order", "family", "genus", "species")]

  # change NA names to unclassified
  taxonomy_tab = propagate_unclassified_taxonomy(taxonomy)

  return(taxonomy_tab)
}

# Wrapper for reading and cleaning obitools3 tax table to phyloseq tax_table()
read_and_filter_tax_table = function(tax_table_path, min_BEST_ID=1, max_BEST_ID=1){
  tax_df = read_obitools_tax_data(tax_table_path)
  tax_out = clean_and_filter_tax_table(tax_df, min_BEST_ID, max_BEST_ID)
  return(tax_out)
}

# Clean MIDORI2 raw
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
  out[out == "NA"] = NA

  return(out)
}

# Very basic parse of SINTAX like taxonomy output (x:xxxxx,y:yyyyyy,....)
parse_tax_string = function(x, is_pred=F){
  # Return empty s to fit in when rbind
  if (x==""){
    return(data.frame("s"=NA))
  }
  # split taxon-levels ,
  xspl = unlist(strsplit(x, ","))

  # Per taxon level, split key:value
  xcol = data.frame(strsplit(xspl, ":"))
  # Format and return
  colnames(xcol) = xcol[1,]
  xcol = xcol[-1,, drop=FALSE]
  row.names(xcol) = NULL

  if (is_pred){
    mycols = colnames(xcol)
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

read_sintax_tabbedout = function(sintax_path, read_supported_only=TRUE, add_source=FALSE, clean_tax=FALSE){

  rtab=read.csv(sintax_path, sep = "\t", header = F, col.names = c("QUERY_LABEL","TAX_PRED", "STRAND", "TAX_PRED_CUTOFF" ))

  # order the input on string length to increase chance of full taxonomy template
  rtab = rtab[order(nchar(rtab$TAX_PRED_CUTOFF), decreasing = T),]
  if (read_supported_only){
    ## Read supported
    tax_supp = data.table::rbindlist(lapply(rtab$TAX_PRED_CUTOFF, parse_tax_string), fill = T)

  } else {

    ## Read complete prediction
    tax_supp = data.table::rbindlist(lapply(rtab$TAX_PRED, parse_tax_string, is_pred=T), fill = T)

  }

  tax_supp = cbind(rtab[1], tax_supp)

  if (clean_tax){tax_supp = clean_taxonomies(tax_supp)}

  if (add_source){tax_supp$SOURCE_FILE = basename(sintax_path)}

  return(tax_supp)

}

# Read idtaxa output Cecilia code
read_idtaxa_table = function(tax_tab_path, read_supported_only=TRUE, add_source=FALSE, clean_tax=FALSE){

  taxranks = c("d", "p", "c", "o", "f", "g", "s")

  tabraw = read.csv(tax_tab_path, sep = ";", header = F)
  colnames(tabraw) = c("QUERY_LABEL", taxranks)

  # WHY did we invent a new format (mixed \t and ;sv) for the taxonomy table of IDTAXA ????
  # ADjust this in the code!!!
  # for now, remove everythin after SWARM name

  tabraw$QUERY_LABEL = gsub(pattern = "\t.*", "", tabraw$QUERY_LABEL)

  if (read_supported_only){
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

# Read BOLDigger3 identification_result.xlsx output
read_boldigger3_xlsx = function(tax_tab_path, add_source = FALSE, clean_tax=FALSE){
  x_ = readxl::read_xlsx(tax_tab_path)
  x_ = as.data.frame(x_[, c(1:7)])
  colnames(x_) = c("QUERY_LABEL" ,"p"         ,   "c"       ,     "o"       ,     "f"         ,   "g"        ,    "s"  )
  x_[x_ == "no-match"] = NA
  x_$d = NA
  if (clean_tax){x_ = clean_taxonomies(x_)}
  if (add_source){x_$SOURCE_FILE=basename(tax_tab_path)}
  return(x_)
}

# Complete OBITools3 taxonomy by replacing NA values with last known taxon + "_unclassified" suffix
# from https://github.com/slambrechts/INBO_eDNA_metabarcoding_BODEM/blob/main/R_helpers/make_phyloseq_helper.R
propagate_unclassified_taxonomy = function(tax){
  tax = as.data.frame(tax)
  # take copy
  new_tax = tax
  #loop columns
  unclassified_tax = paste0(tax[,1], " unclassified")
  for (taxrank in colnames(new_tax)){
    # find rows that have no taxonomy
    unclassified_mask = is.na(new_tax[,taxrank])
    # propagate taxonomy these rows with no taxonomy
    new_tax[unclassified_mask , taxrank] = unclassified_tax[unclassified_mask]

    # update taxonomy
    unclassified_tax = new_tax[,taxrank]
    # only do + "_unclassified"to not already unclassified
    unclassified_tax[!unclassified_mask] = paste0(unclassified_tax[!unclassified_mask], " unclassified")

  }

  return(tax_table(as.matrix(new_tax)))
}
