#' Plots a Diagram of the Gillespie algorithm
#'
#' This produces an animation, recreating the Gillespie algorithm
#' visually, without trying to be numerically correct.
#'
#' We use only the default values for every model constituent.
#'
#' The C backend is not used here.
#'
#' @param x a stochastic model
#' @param y final time of animation
#' @param arrows use the arrows function to draw reaction arrows
#' @param ... ignored
#' @return x, invisibly
#' @export
#' @examples
#' m <- uqsa_example("AKAR4") |> model_from_tsv()
#' G <- as_cme(m)
#' G$initialCount["C"] <- 100
#' G$initialCount["AKAR4_C"] <- 100
#' if (interactive()) plot(G,100)
plot.cme <- function(x,y=100,arrows=FALSE,...){
	on.exit(
		{
			while (dev.cur() > 1 && dev.hold(0) > 0) dev.flush()
		},
		add = TRUE
	)
	t_final <- y
	n <- length(x$initialCount)
	r <- 2*NROW(x$kinetic.law)
	reactants <- c(x$left,x$right)
	products <- c(x$right,x$left)
	## x-Positioning (y-Positioning depends on number of compounds and reactions)
	xPOS <- c(-1,0,1)
	yPOS <- seq((n-1)/2,by=-1,length.out=n)
	rPOS <- seq((r-1)/2,by=-1,length.out=r)
	names(yPOS) <- names(x$initialCount)
	## thermometer template
	tht <- c(0.2,0.5) # width and height
	A <- str2lang(
		paste(
			c(
				"{",
				if (!is.null(x$expression)) paste0(names(x$expression),"=",x$expression),
				paste0("c(",paste(x$kinetic.law,collapse=", "),")"),
				"}"
			),
			collapse="\n"
		)
	)
	state <- x$initialCount
	ENV <- list2env(as.list(c(x$par,state)))
	a_hat <- max(eval(A, ENV), 1e-12)
	t_ <- 0
	while (t_ < t_final){
		a <- eval(A,ENV)
		a0 <- sum(a)
		if (a0<=0) break
		tau <- rexp(1,rate=a0)
		t_ <- t_ + tau
		j <- sample.int(length(a),1,prob=a/a0)
		L <- reactants[[j]]
		R <- products[[j]]
		rP <- rPOS[j]
		## redraw
		dev.hold()
		frame()
		plot.window(xlim=c(-2,2),ylim=c(-max(n,r)/2,max(n,r)/2))
		title(main=x$name)
		mtext(sprintf("time: %g",t_),side=1)
		## 1. Reactants
		cBars <- cbind(rep(tht[1],n),rep(tht[2],n),state/sum(state,1))
		symbols(rep(xPOS[1],n),yPOS,thermometers=cBars,fg='red4',inches=FALSE,add=TRUE)
		text(
			rep(xPOS[1]-tht[1]/2-0.01,n),
			yPOS,
			labels=names(x$initialCount),
			pos=2
		)
		text(xPOS[1],max(yPOS)+tht[2],labels="REACTANTS",pos=3)
		text(xPOS[3],max(yPOS)+tht[2],labels="PRODUCTS",pos=3)
		## 2. Reactions
		rBars <- cbind(rep(tht[1],r),rep(tht[2],r),a/max(a_hat,a0))
		symbols(rep(xPOS[2],r),rPOS,thermometers=rBars,fg='yellow4',inches=FALSE,add=TRUE)
		text(
			rep(xPOS[2],r),
			rPOS+tht[2]/2,
			labels=paste0("r",seq(r)),
			pos=3
		)
		## 3. Products
		if (length(L)>0) state[names(L)] <- state[names(L)]-L
		if (length(R)>0) state[names(R)] <- state[names(R)]+R
		list2env(as.list(state), envir = ENV)
		cBars <- cbind(rep(tht[1],n),rep(tht[2],n),state/sum(state,1))
		symbols(rep(xPOS[3],n),yPOS,thermometers=cBars,fg='red4',inches=FALSE,add=TRUE)
		text(
			rep(xPOS[3]+tht[1]/2+0.01,n),
			yPOS,
			labels=names(x$initialCount),
			pos=4
		)
		COL <- "blue4"
		if (isTRUE(arrows)){
			if (length(L)>0){
				arrows(rep(xPOS[1],length(L)),yPOS[names(L)],rep(xPOS[2],length(L)),rep(rP,length(L)),col=COL,lwd=2)
			}
			if (length(R)>0){
				arrows(rep(xPOS[2],length(R)),rep(rP,length(R)),rep(xPOS[3],length(R)),yPOS[names(R)],col=COL,lwd=2)
			}
		} else { # bezier curves
			dx <- tht[1]*0.6
			DX <- 4*dx
			for (i in seq_along(L)){
				xspline(
					x=c(xPOS[1]+dx,xPOS[1]+DX,xPOS[2]-DX,xPOS[2]-dx),
					y=c(yPOS[names(L)[i]],yPOS[names(L)[i]],rP,rP),
					shape=1,
					border=COL,
					lwd=2
				)
			}
			for (i in seq_along(R)){
				xspline(
					x=c(xPOS[2]+dx,xPOS[2]+DX,xPOS[3]-DX,xPOS[3]-dx),
					y=c(rP,rP,yPOS[names(R)[i]],yPOS[names(R)[i]]),
					shape=1,
					border=COL,
					lwd=2
				)
				symbols(
					xPOS[3]-dx,
					yPOS[names(R)[i]],
					stars=matrix(c(2,1,2,0.5,2,1)*4e-2,1,6),
					bg=COL,
					inches=FALSE,
					add=TRUE
				)
			}
		}
		dev.flush()
		Sys.sleep(max(tau/10,0.05))
	}
	return(invisible(x))
}


#' Represents a model in a visual way
#'
#' This draws an illustration that represents the interactions within
#' a model.
#'
#' @param x an sbmodel (list of data-frames)
#' @param rotate (logical) the compound labels will be angled if TRUE]
#' @param reorder an attempt is made to place the reaction nodes
#'     symmetrically with respect to the reactants
#' @param arrows direct arrows are drawn between reactantsm reaction
#'     node, and products. If FALSE, a curved line is drawn.
#' @param ... ignored
#' @return nothing (called for sideeffect)
#' @export
#' @examples
#' m <- uqsa_example("AKAP79") |> model_from_tsv()
#' plot(m)
plot.sbmodel <- function(x,rotate=TRUE,reorder=FALSE,arrows=TRUE,...){
	m <- x
	N <- stoichiometric_matrix(m)
	## the conservation law matrix, could be used to group species
	## Z <- pracma::rref(t(pracma::null(t(N))))
	R_ <- 0.45
	R1 <- 0.5
	R0 <- 1.0
	R2 <- 1.5
	R3 <- 1.65
	n <- NROW(m$Compound)
	dphi <- 2*pi/n
	phi <- seq(0,2*pi-dphi,length.out=n)
	C <- cos(phi)
	S <- sin(phi)
	r_ <- R_*cbind(C,S)
	r1 <- R1*cbind(C,S)
	r2 <- R2*cbind(C,S)
	r3 <- R3*cbind(C,S)
	d <- paste0("#",apply(m$Compound,1,digest::digest,algo="crc32"))
	frame()
	plot.window(xlim=c(-2,2),ylim=c(-2,2))
	title(main=comment(m))
	symbols(
		r2[,1],r2[,2],
		circles=rep(0.07,n),
		inches=FALSE,
		xlab="",
		ylab="",
		bg=d,
		add=TRUE
	)
	symbols(
		r1[,1],r1[,2],
		circles=rep(0.03,n),
		inches=FALSE,
		bg=d,
		add=TRUE)
	symbols(
		0,0,
		circles=(0.6),
		inches=FALSE,
		add=TRUE,
		fg=NA,
		bg="#40404010"
	)
	text(0,0,labels="reactants")
	if (as.logical(rotate)){
		for (k in seq(n)){
			text(
				r3[k,1],r3[k,2],
				labels=rownames(m$Compound)[k],
				srt = phi[k] * 180 / pi,
				adj = 0
			)
			text(
				r_[k,1],r_[k,2],
				labels=rownames(m$Compound)[k],
				srt = phi[k] * 180 / pi,
				adj = 1.1,
				col = "#80808080",
				cex=0.5
			)
		}
	} else {
		for (k in seq(n)){
			text(
				r3[k,1],r3[k,2]+0.07*S[k]*(abs(S[k])>0.99),
				labels=rownames(m$Compound)[k],
				adj=0.5*cbind((1-C[k]),(1-S[k]))
			)
		}
	}
	nr <- NROW(m$Reaction)
	dphi <- 2*pi/nr
	r_phi <- seq(0,2*pi-dphi,length.out=nr)
	if (isTRUE(reorder)){
		for (j in seq(nr)) {
			w <- which(N[, j] != 0)
			MPHI <- mean(phi[w])
			r_phi[j] <- MPHI
		}
	}
	C <- R0*cos(r_phi)
	S <- R0*sin(r_phi)
	r0 <- cbind(C,S)
	symbols(C,S,squares=rep(0.05,nr),inches=FALSE,add=TRUE)
	d <- paste0("#",substr(apply(m$Reaction,1,digest::digest,algo="crc32"),start=1,stop=6),"80")
	if (isTRUE(arrows)){
		for (j in seq(NCOL(N))){
			if (any(N[,j]<0)) {
				w <- which(N[,j]<0)
				for (i in w){
					arrows(r1[i,1],r1[i,2],r0[j,1],r0[j,2],lwd=4,col=d[j],length=0)
				}
			}
			if (any(N[,j]>0)) {
				w <- which(N[,j]>0)
				for (i in w){
					arrows(r0[j,1],r0[j,2],r2[i,1],r2[i,2],lwd=4,col=d[j])
				}
			}
		}
	} else {
		for (j in seq(NCOL(N))){
			if (any(N[,j]<0)) {
				w <- which(N[,j]<0)
				for (i in w){
					p1 <- r1[i,]
					p2 <- r0[j,]
					m1 <- 1.3*p1
					m2 <- 0.7*p2
					xspline(
						c(p1[1],m1[1],m2[1],p2[1]),
						c(p1[2],m1[2],m2[2],p2[2]),
						open=TRUE,
						shape=1,
						border=d[j],
						lwd=4
					)
				}
			}
			if (any(N[,j]>0)) {
				w <- which(N[,j]>0)
				for (i in w){
					p1 <- r0[j,]
					p2 <- r2[i,]
					m1 <- 1.3*p1
					m2 <- 0.7*p2
					xspline(
						c(p1[1],m1[1],m2[1],p2[1]),
						c(p1[2],m1[2],m2[2],p2[2]),
						open=TRUE,
						shape=1,
						lty=2,
						border=d[j],
						lwd=4
					)
				}
			}
		}
	}
}


replace_inf_with_na <- function(x){
	x[!is.finite(x)] <- NA
	if (inherits(x,"errors")){
		errors(x)[!is.finite(errors(x))] <- NA
		x[!is.finite(errors(x))] <- NA
	}
	return(x)
}

#' Plot function for experiments
#'
#' This function uses plot.errors and the base plot functions like
#' [matplot].
#' @keywords internal
#' @param x experiment setup (a list)
#' @param y simulation results (a list)
#' @param func.names restrict the shown output functions to these
#'     named outputs
#' @param verbose prints additional information about name matching
#'     between experiments and simulations
#' @param ... forwarded to the more specific plot function
#'     [plot.errors]
#' @return plot object
#' @export
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- write_and_compile(as_ode(m))
#' ex <- experiments(m,o)
#' s <- simulator.c(ex,o)
#' p0 <- values(m$Parameter)
#' y <- s(p0)
#' plot(ex,y)
plot.experiments <- function(x, y, func.names=NULL, verbose=getOption("uqsa.verbose",FALSE),...){
	if (missing(y)) {
		simulations <- NULL
	} else {
		stopifnot(length(x)==length(y))
		simulations <- y
		nf <- dim(simulations[[1]]$func)[1]
	}
	if (missing(func.names)){
		out <- seq(NROW(x[[1]]$data))
	} else if (is.character(func.names)) {
		out <- which(rownames(x[[1]]$data) %in% func.names)
	} else if (is.numeric(func.names)){
		out <- func.names
	}
	p <- list()
	for (i in seq_along(x)){
		time <- x[[i]]$outputTimes
		if (!is.null(simulations)){
			n <- dim(simulations[[i]]$func)
		} else {
			n <- c(dim(x[[i]]$data),1)
		}
		oNames <- rownames(x[[i]]$data)
		for (j in out){
			a <- ceiling(255*exp(-0.01*n[3]))
			## only plot if the output was actually measured
			if (verbose) cli::cli_alert_info("{oNames[j]} %in% {paste(colnames(x[[i]]$measurements),collapse=' ')}")
			if (oNames[j] %in% colnames(x[[i]]$measurements)){
				d <- replace_inf_with_na(x[[i]]$data[j,])
				p <- c(p,
					plot(
						errors::as.errors(time),
						d,
						main=names(x)[i],
						xlab="time",
						ylab=oNames[j],
						lty=1,
						...
					)
				)
				if (!is.null(simulations)){
					matplot(
						time,
						simulations[[i]]$func[j,,],
						col=rgb(0,0,255,a,maxColorValue=255),
						lwd=2,
						add=TRUE,
						type="l",
						lty=1
					)
				}
			}
		}
	}
	return(invisible(p))
}
