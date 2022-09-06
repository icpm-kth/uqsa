x<-c(3, 4, 3.5, 3.4, 4.2)

f_kde <- kde(x, xmin=min(x)-.1,xmax=max(x)+.1)
#hist(x, breaks=10, xlim=c(min(x)-.1,max(x)+.1))
hist(x, breaks=10, xlim=c(0,8))
lines(f_kde$eval.points,f_kde$estimate)

f_kde <- kde(x, xmin=0,xmax=8)
#hist(x, breaks=10, xlim=c(0,8))
lines(f_kde$eval.points,f_kde$estimate)

###

x <- c(3, 3, 3, 4, 3.5, 3.4, 4.2)

f_kde <- kde(x, xmin=min(x)-.1,xmax=max(x)+.1)
hist(x, breaks=10, xlim=c(0,8))
#hist(x, breaks=10, xlim=c(min(x)-.1,max(x)+.1))
lines(f_kde$eval.points,f_kde$estimate)

f_kde <- kde(x, xmin=0,xmax=8)
#hist(x, breaks=10, xlim=c(0,8))
lines(f_kde$eval.points,f_kde$estimate)
