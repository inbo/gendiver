
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gendiver

<!-- badges: start -->
<!-- badges: end -->

-   R utilies for I/O (and working with) bioIT output files
-   Facilitate loading the results of different tools and software in R
-   Generic visualization and analysis functions for QC of MB-runs
-   ….. And maybe more to come!

## Installation

You can install the development version of gendiver from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("inbo/gendiver")
```

## Example

### INPUT / OUTPUT

Read Operational Data from INBO Shared Google Drive (institutional
access required)

``` r
library(gendiver)

my_sample_data = read.sample_sheet_lab("<GSHEET_ID/URL>")

my_otu_tab = read.otu_table_vsearch("<path/url>")

my_tax_tab = read.taxonomy_sintax("<path/url>")
```

OBITools3 taxonomy with filter would be like this:

``` r
my_tax_tab = read.taxonomy_obitools3("<path/url>") %>% select.taxonomy_obitools3() 
```

Preparing these dataframes to use with `phyloseq`:

``` r
library(phyloseq)

ps_sample_data = read.sample_sheet_lab("<GSHEET_ID/URL>") %>% sample_data()

ps_otu_tab = read.otu_table_vsearch("<path/url>") %>% otu_table(taxa_are_rows = T)

ps_tax_tab = read.taxonomy_sintax("<path/url>") %>% as.matrix() %>% tax_table()
```

Combine these inputs into a `phyloseq` object:

``` r
my_ps = phyloseq(
  ps_sample_data,
  ps_otu_tab,
  ps_tax_tab
)
```

### QC plots

Minimal examples, assumes data is loaded in phyloseq object

``` r
# heatmaps
gendiver::qcplot.plate_heatmap_readcount(ps_obj)

gendiver::qcplot.plate_heatmap_toptaxa(ps_obj)

gendiver::qcplot.plate_column_asv_barplot(ps_obj_plate)
```
