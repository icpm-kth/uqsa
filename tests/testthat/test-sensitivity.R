test_that("sum.of.bin.variance functions work", {
	X <- c(rnorm(100,sd=1),7) # 7>3*sd to have at least a few NA values in binMeansY
	ny <- 2
	Y <- matrix(c(sin(X),cos(X)),ncol=ny)
	hst<-hist(X,plot=FALSE,breaks=30)
	expect_true(any(hst$counts==0))
	interval <- findInterval(X,hst$breaks,all.inside=TRUE)
	expect_lt(max(interval),30)
	expect_gt(min(interval),0)
	binMeansY <- observable.mean.in.bin(interval,Y)
	expect_equal(nrow(binMeansY),length(hst$counts))
	V <- sum.of.bin.variance(hst,binMeansY,colMeans(Y))
	expect_type(V,"double")
	expect_vector(V, ptype = numeric(0), size = ny)
	expect_lt(max(V),1.0)
	expect_gt(min(V),0.0)
})

test_that("sensitivity function works",{
	X <- c(rnorm(100,sd=1),7) # 7>3*sd to have at least a few NA values in binMeansY
	ny <- 2
	Y <- matrix(c(sin(X),cos(X)),ncol=ny)
	S <- globalSensitivity(as.matrix(X),Y,nBins=30)
	expect_type(S,"double")
	expect_true(is.matrix(S))
	expect_equal(dim(S),c(2,1))
	expect_lt(max(S),1.0)
	expect_gt(min(S),0.0)
})

