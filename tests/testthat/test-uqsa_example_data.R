test_that("Examples are findable", {
	expect_vector(uqsa_example(),ptype=character())
})

test_that("AKAR4 tsv files are findable (f) and valid", {
	expect_vector(AKAR4.tsv <- uqsa_example("AKAR4",f="tsv"),ptype=character())
	expect_vector(tab <- model_from_tsv(AKAR4.tsv))
	expect_true(suppressWarnings(suppressMessages(SBtabVFGEN::sbtab.valid(tab))))
})

test_that("non-existing example fails",{
	expect_error(uqsa_example("ABCDEF"))
})

