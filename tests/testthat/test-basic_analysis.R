test_o = read.otu_table_vsearch("https://drive.google.com/file/d/19nQ9bW69dFkVJ0RYRxgDJVvB3Z0clqfx/view?usp=drive_link")
test_t = read.taxonomy_sintax("https://drive.google.com/file/d/1kQOOScnu8UM9QFg6NwKlsMgKeznIlGG5/view?usp=drive_link")
test_s = read.sample_sheet_lab("https://docs.google.com/spreadsheets/d/11WIP5yQtAfthzgwX2S-vcvSM5ybSPAKk6AHCP2IiGq8/edit?usp=drive_link")

test_that("testdata is good", {

  # samplesheet vs otu_tab
  expect_true(sum( ! colnames(test_o) %in% row.names(test_s)) == 0)

  # otu_tab vs tax_tab
  expect_true(nrow(test_o) == nrow(test_t) )
  expect_true(sum( ! row.names(test_o) %in% row.names(test_t)) == 0)

})


test_that("phyloseq pipe works", {

  test_ps = phyloseq::phyloseq(
    phyloseq::sample_data(test_s),
    phyloseq::otu_table(test_o, taxa_are_rows = T),
    phyloseq::tax_table(as.matrix(test_t))
  )

  expect_equal(class(test_ps)[1], c("phyloseq"))

  # See also: https://testthat.r-lib.org/reference/expect_setequal.html
  # Tax table
  expect_equal(phyloseq::ntaxa(test_ps), nrow(test_t))
  expect_setequal(phyloseq::taxa_names(test_ps), rownames(test_t))
  expect_equal(colnames(phyloseq::tax_table(test_ps)), colnames(test_t))
  expect_equal(as.data.frame(phyloseq::tax_table(test_ps))[rownames(test_t),], test_t)

  # OTU table
  expect_equal(phyloseq::ntaxa(test_ps), nrow(test_o))
  expect_setequal(phyloseq::taxa_names(test_ps), rownames(test_o))
  expect_setequal(phyloseq::sample_names(test_ps), colnames(test_o))
  expect_equal(phyloseq::sample_sums(test_ps), colSums(test_o))
  expect_equal(as.data.frame(phyloseq::otu_table(test_ps))[rownames(test_o),], test_o)

  # Sample sheet
  expect_in(phyloseq::sample_names(test_ps), rownames(test_s))
  expect_equal(data.frame(phyloseq::sample_data(test_ps)), test_s[phyloseq::sample_names(test_ps),])

})
