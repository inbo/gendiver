test_that("read.sample_sheet_lab can read the template gsheet", {
  sample_sheet = read.sample_sheet_lab(
      "11WIP5yQtAfthzgwX2S-vcvSM5ybSPAKk6AHCP2IiGq8"
    )

  expect_true(is.data.frame(sample_sheet))
  expect_equal(ncol(sample_sheet), 39)
  expect_equal(nrow(sample_sheet), 54)
})

