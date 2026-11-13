test_that("column works", {
	# select a column by number, names are preserved
	x <- column(mtcars,1)
	y <- mtcars[[1]]
	names(y) <- rownames(mtcars)
	expect_equal(x,y)
	# several option, first match is accepted
	x <- column(mtcars,c("gear","carb"))
	y <- mtcars[["gear"]]
	names(y) <- rownames(mtcars)
	expect_equal(x,y)
	# no match: NULL
	expect_null(column(mtcars,c("GEAR","CARB")))
})
