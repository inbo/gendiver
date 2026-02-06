test_that("read table wrapper works LOCAL", {
  mytab_local = read.table_gdrive(
    test_path("testdata","example.mumu.log"),
    sep="\t")
  expect_equal(ncol(mytab_local), 18)
})

test_that("read table wrapper works GDRIVE", {
  mytab_remote = read.table_gdrive(
    "https://drive.google.com/file/d/19nQ9bW69dFkVJ0RYRxgDJVvB3Z0clqfx/view?usp=drive_link",
    sep="\t")
  expect_true(is.data.frame(mytab_remote))
  expect_equal(ncol(mytab_remote), 13)
  expect_equal(nrow(mytab_remote), 24)

})


test_that("read table wrapper gives same result local v.s. remote", {
  x1 = read.table_gdrive(
    "https://drive.google.com/file/d/19nQ9bW69dFkVJ0RYRxgDJVvB3Z0clqfx/view?usp=drive_link",
    sep="\t")

  x2 = read.table_gdrive(
    test_path("testdata","example.zotutab.table"),
    sep="\t")

  expect_true(identical(x1, x2))

})
