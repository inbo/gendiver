test_that("read.otu_table_vsearch works", {
  mytab = read.otu_table_vsearch(
    test_path("testdata","example.zotutab.table")
    )
  expect_equal(ncol(mytab), 12)
  expect_equal(nrow(mytab), 24)
  expect_true(is.data.frame(mytab))

})


test_that("read.otu_table_obitools3 works", {
  mytab = read.otu_table_obitools3(
    test_path("testdata","example.asv_obitools3.table")
  )
  expect_equal(ncol(mytab), 12)
  expect_equal(nrow(mytab), 31)
  expect_true(is.data.frame(mytab))

})

test_that("read.log_mumu works", {
  mytab = read.log_mumu(
    test_path("testdata","example.mumu.log")
  )
  expect_equal(ncol(mytab), 18)
  expect_equal(nrow(mytab), 17)
  expect_true(is.data.frame(mytab))

})

test_that("read.otu_table_vsearch local v.s. remote", {
  x1c = read.otu_table_vsearch(
    "https://drive.google.com/file/d/19nQ9bW69dFkVJ0RYRxgDJVvB3Z0clqfx/view?usp=drive_link")

  x2c = read.otu_table_vsearch(
    test_path("testdata","example.zotutab.table"))

  expect_true(identical(x1c, x2c))

})
