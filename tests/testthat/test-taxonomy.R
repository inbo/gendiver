# read otu table
myotu = read.otu_table_vsearch(
  test_path("testdata","example.zotutab.table")
)

test_that("read.taxonomy_obitools3 works", {
  mytax = read.taxonomy_obitools3(
    test_path("testdata","example.taxonomy_obitools3.table")
  )

  expect_true(is.data.frame(mytax))
  expect_equal(nrow(myotu), nrow(mytax))
  expect_true(sum(!row.names(myotu) %in% row.names(mytax)) == 0 )

})

test_that("read.taxonomy_sintax works", {
  mytax = read.taxonomy_sintax(
    test_path("testdata","example.taxonomy_sintax.table")
  )

  expect_true(is.data.frame(mytax))
  expect_equal(nrow(myotu), nrow(mytax))
  expect_true(sum(!row.names(myotu) %in% row.names(mytax)) == 0 ) # should be same to expect_setequal()
  expect_setequal(row.names(myotu), row.names(mytax))
})

test_that("read.taxonomy_sintax local v.s. remote", {
  x1t = read.taxonomy_sintax(
    "https://drive.google.com/file/d/1kQOOScnu8UM9QFg6NwKlsMgKeznIlGG5/view?usp=drive_link")

  x2t = read.taxonomy_sintax(
    test_path("testdata","example.taxonomy_sintax.table"))

  expect_true(identical(x1t, x2t))

})
