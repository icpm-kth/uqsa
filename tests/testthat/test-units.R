test_that("%as% works on this system", {
	skip_if(nchar(Sys.which('units'))==0)
	u <- "11 pound"
	y <- u %as% "kg"
	expect_named(y,u)
	expect_lt(y,5)
})

test_that("%as% fails gracefully on Linux", {
	skip_if(nchar(Sys.which('units'))==0)
	platform <- Sys.info()['sysname']
	skip_if_not(platform == 'Linux')
	expect_warning(z <- "pounds" %as% "meters")
	expect_equal(attr(z,"status"),1)
	expect_true(startsWith(attr(z,"stderr"),"conformability error"))
	expect_equal(z,NA,ignore_attr = TRUE)
})

test_that("%as% fails gracefully on MacOS and FreeBSD", {
	skip_if(nchar(Sys.which('units'))==0)
	platform <- Sys.info()['sysname']
	skip_if(platform == 'Linux')
	z <- "pounds" %as% "meters"            # no warning on FreeBSD and MacOS
	## no status is returned on Mac and FreeBSD, which is their fault
	expect_true(startsWith(attr(z,"stderr"),"conformability error"))
	expect_equal(z,NA,ignore_attr = TRUE)
})
