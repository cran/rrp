"rrp.dist" <-
function(X, treated=NULL, msplit=10, Rep=250,cut.in=15, check.bal=FALSE, 
plot=FALSE){
	n <- dim(X)[1]
	n.var <- dim(X)[2]
	
	RRP <- rep(as.integer(Rep), n*(n-1)/2)
	class(RRP) <- "dist"
	attr(RRP,"Size") <- n
	attr(RRP,"Diag") <- FALSE
	attr(RRP,"Upper") <- FALSE
	attr(RRP,"method") <- "RRP"
	attr(RRP,"minsplit") <- msplit
	attr(RRP,"replications") <- Rep
	attr(RRP,"cov.cut") <- cut.in
	attr(RRP,"balanced") <- check.bal
	if(!is.null(treated))
		attr(RRP,"treated") <- as.logical(treated)
	attr(RRP,"call") <- match.call()
	
	Y <- X
	W <- X
	for(i in 1:dim(X)[2])
		W[,i] <- as.numeric(X[,i])

	if(cut.in>0){
		for(i in 1:n.var){
			if(is.numeric(X[,i]) | is.integer(X[,i])){
				if(length(unique(X[,i])) > cut.in)
				  Y[,i] <- ordered(cut(X[,i], seq(min(X[,i], na.rm=TRUE), 
					max(X[,i], na.rm=TRUE), length=cut.in),include.lowest=TRUE))				
				else		
					Y[,i] <- ordered(X[,i])			
			}
		}
	}

	x <- cbind(z=numeric(n), Y)

	for(K in 1:Rep){
	if(plot){
		if(is.null(treated))
			plot(X[,1:2])	
		else
			plot(X[,1:2],col=treated+2,pch=ifelse(treated,20,17),cex=1)
	}
	
	x[,1] <- runif(n) 

	rpart(z~., data=x,  method="anova",  minsplit=msplit, xval=0, cp=0) -> mod

	idx <- which(mod$frame$var=="<leaf>") # we get the terminal nodes
	n.leaves <- length(idx)  # we count the terminal nodes
	leaves <- as.integer(dimnames(mod$frame)[[1]]) # we get the labels of all the nodes
	LEAF.num <- leaves[idx] # we extract only the ones that are terminal nodes
	LEAF.pred <- leaves[mod$where] # we associate to each observation the corresponding node number

	group <- sapply(LEAF.num,function(x) which(LEAF.pred==x)) # indexes of matched per leaf
	if(!is.null(treated) & check.bal){
		idx.t <- lapply(group,function(x) x[which(treated[x]==1)])
		idx.c <- lapply(group,function(x) x[which(treated[x]==0)])
	}

	update.proxy <- function(g){
		similar <- unlist(group[[g]])
		residual <- numeric(0)
		if(check.bal){
			if(length(idx.t[[g]])>0 & length(idx.c[[g]])>0){
				apply(W[idx.t[[g]],], 2, function(x) min(x, na.rm=TRUE)) -> mins
				apply(W[idx.t[[g]],], 2, function(x) max(x, na.rm=TRUE)) -> maxs
				test <- sapply(idx.c[[g]], function(x) prod(mins <= W[x,])*prod(maxs >= W[x,]))
				test <- as.logical(test)
				similar <- c(idx.t[[g]], (idx.c[[g]])[test])
				residual <- (idx.c[[g]])[!test]
			}
		}
		RRP <<- addDist(RRP, similar, -1)
		if(length(residual)>0)
			RRP <<- addDist(RRP, residual, -1)
		if(plot){
			minX <- min(W[similar,1], na.rm=TRUE)
			maxX <- max(W[similar,1], na.rm=TRUE)
			minY <- min(W[similar,2], na.rm=TRUE)
			maxY <- max(W[similar,2], na.rm=TRUE)
			rect( minX,minY,maxX,maxY,lty=3)
			if(length(residual)>0){
				minX <- min(W[residual,1], na.rm=TRUE)
				maxX <- max(W[residual,1], na.rm=TRUE)
				minY <- min(W[residual,2], na.rm=TRUE)
				maxY <- max(W[residual,2], na.rm=TRUE)
				rect( minX,minY,maxX,maxY,lty=3)
			}
		}
	}

	sapply(1:n.leaves, update.proxy) # we fill up the proximity matrix
   
	if(K %/% 100 == K/100) {
		cat("0\n")
	} else {
		if(K %/% 10 == K/10)
		  cat("+")
		 else
 		  cat(".")
		}
	} # end for(rep in n.rep)
	RRP <- RRP/Rep
}

