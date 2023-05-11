test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("AKAR4 is findable", {
	expect_vector(uqsa_example(),ptype=character())
})

test_that("AKAR4 tsv files are findable (f)", {
	expect_vector(uqsa_example("AKAR4",f="tsv"),ptype=character())
})

test_that("AKAR4 c files are findable ", {
	expect_vector(uqsa_example("AKAR4",f="c"),ptype=character(),size=1)
})

test_that("non-existing example fails",{
	expect_error(uqsa_example("ABCDEF"))
})
