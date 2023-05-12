test_that("multiplication works", {
  expect_equal(2 * 2, 4)
	x <- 9
	expect_equal(3 * 3, x)
})

test_that("Examples are findable", {
	expect_vector(uqsa_example(),ptype=character())
})

test_that("AKAR4 tsv files are findable (f) and valid", {
	expect_vector(AKAR4.tsv <- uqsa_example("AKAR4",f="tsv"),ptype=character())
	expect_vector(tab <- SBtabVFGEN::sbtab_from_tsv(AKAR4.tsv,verbose=FALSE))
	expect_null(suppressMessages(SBtabVFGEN::sbtab.valid(tab)))
})

test_that("AKAP79 tsv files are findable (f) and valid", {
	expect_vector(AKAP79.tsv <- uqsa_example("AKAP79",f="tsv"),ptype=character())
	expect_vector(tab <- SBtabVFGEN::sbtab_from_tsv(AKAP79.tsv,verbose=FALSE))
	expect_warning(suppressMessages(SBtabVFGEN::sbtab.valid(tab)))
})

test_that("AKAR4 c files are findable ", {
	expect_vector(uqsa_example("AKAR4",f="c"),ptype=character(),size=1)
})

test_that("non-existing example fails",{
	expect_error(uqsa_example("ABCDEF"))
})

