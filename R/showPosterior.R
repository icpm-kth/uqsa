# upper and lower panel defaults
lp <- function(x,y,num=round(min(NROW(x)/2,NROW(y)/2)),...){
	rx <- range(x)
	ry <- range(y)
	C <- cor(head(x,num),head(y,num)) # in [-1,1]
	i <- 1+round((1+C)*50)            # 1+C in [0,2], (1+C)*50 in [0,100], 1+round((1+C)*50) in [1,101]
	col <- colorspace::diverging_hcl(101,rev=FALSE,palette="Blue Red") # maps to possible correlation values
	polygon(c(rx,rev(rx)),rep(ry,each=2),col=col[i])
	text(mean(rx),mean(ry),sprintf("%.3f",C),cex=2+round(abs(C)*2))
}

up <- function(x,y,num=round(min(NROW(x)/2,NROW(y)/2)),subscripts,...){
	m <- mean(x)
	dx <- diff(range(x))
	dy <- diff(range(y))
	LIM <- c(min(x)-dx,max(x)+dx,min(y)-dy,max(y)+dy)
	if (min(y) < min(x)){
		poly_X <- LIM[c(1,1,2,2)]
		poly_Y <- LIM[c(1,3,3,2)]
	} else {
		poly_X <- LIM[c(1,2,2)]
		poly_Y <- LIM[c(1,1,2)]
	}
	k1 <- MASS::kde2d(head(x,num),head(y,num),n=500,lims=LIM)
	k2 <- MASS::kde2d(tail(x,num),tail(y,num),n=500,lims=LIM)
	C <- contourLines(k2)
	colPatch <- c("#FFFFFF",colorspace::sequential_hcl(length(C),rev=TRUE,palette="Blues 3"))
	colContour <- c("#FFFFFF",colorspace::sequential_hcl(length(C),rev=TRUE,palette="Blues 3"))
	image(k1$x,k1$y,k1$z,add=TRUE,col=colPatch)
	for (i in seq_along(C)){
		lines(C[[i]],col=colContour[i])
	}
	polygon(poly_X,poly_Y,col=rgb(0.2,1,0.2,0.05),border=NA)
}

#' showPosterior makes a pairs plot for a sample
#'
#' This function will display the difference between the posterior and
#' prior by plotting the posterior as shaded density plots and the
#' prior as contour lines of level sets. If the two are identical, the
#' lines will be invisible as they blend into the density
#' plot. Otherwise the contour lines will show up as a distinct
#' feature.
#'
#' @param posterior a matrix, each row is a sample member
#' @param prior a matrix of the same size as the posterior
#' @param lower.panel a lower panel plot function with a default that
#'     shows correlation values (color coded)
#' @param upper.panel an upper panel plot function which defaults to
#'     shaded density plots of the posterior with contour lines for
#'     the prior distribution.
#' @return pairs plot object
#' @export
showPosterior <- function(posterior, prior, lower.panel=lp, upper.panel=up,...){
	return(
		pairs(
			rbind(posterior,prior),
			upper.panel=up,
			lower.panel=lp,
			num=NROW(posterior),
			...
		)
	)
}


#' plots a sample in parallel coordinates
#'
#' This function makes a plot that is quite similar to parallel
#' coordinates. It includes information about the prior as error-bars,
#' centered around th eprior's median.
#'
#' @param x a matrix, with N rows (sample-members), and M columns
#'     (different model parameters). The columns must be named.
#' @param prior a data.frame with at least $median, and $stdv
#'     columns. This data.frame may also include the fields: color,
#'     and colorOutline to change the prior error-bars.
#' @param color the color of the sample lines, should have some
#'     transparency.
#' @param ... parameters are passed to matplot.
#' @export
#' @return produces a plot
pcDist <- function(x,prior,color=rgb(0.5,0.5,0.5,0.05),...){
	pI <- seq(NCOL(x))
	stopifnot(all(c("median","stdv") %in% names(prior)))
	m <- prior$median
	sd <- prior$stdv
	u <- m+sd
	l <- m-sd
	names(pI) <- colnames(x)

	matplot(pI,t(x),type="l",lty=1,col=color,xlab=NA,axes = FALSE,lw=3,...)
	if ("color" %in% names(prior)){
		pColor <- prior$color
	} else {
		pColor <- "green2"
	}
	if ("colorOutline" %in% names(prior)){
		colorOutline <- prior$colorOutline
	} else {
		pColorOutline <- "black"
	}

	arrows(pI,m,pI,u,angle=90,col=pColorOutline,lw=2)
	arrows(pI,m,pI,l,angle=90,col=pColorOutline,lw=2)

	lines(pI,x[1,],col="red3",lw=2)
	arrows(pI,m,pI,u,angle=90,col=pColor)
	arrows(pI,m,pI,l,angle=90,col=pColor)
	axis(2)
	axis(1,at=pI,labels=colnames(x),xpd=TRUE)
}
