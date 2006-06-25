"rrp.dist" <- function (X, treated = NULL, msplit = 10, Rep = 250, cut.in = 15, 
    check.bal = FALSE, plot = FALSE, asdist=FALSE) 
{
    n <- dim(X)[1]
    n.var <- dim(X)[2]
	
	RRP <- newXPtr(n, Rep);
	t.att <- NULL
	Y <- X

    if (cut.in > 0) {
        for (i in 1:n.var) {
            if (is.numeric(X[, i]) | is.integer(X[, i])) {
                if (length(unique(X[, i])) > cut.in) 
                  Y[, i] <- ordered(cut(X[, i], seq(min(X[, i], 
                    na.rm = TRUE), max(X[, i], na.rm = TRUE), 
                    length = cut.in), include.lowest = TRUE))
                else Y[, i] <- ordered(X[, i])
            }
        }
    }
    for (i in 1:dim(X)[2]) X[, i] <- as.numeric(X[, i])
    x <- cbind(z = numeric(n), Y)
	rm(Y)
    for (K in 1:Rep) {
    	if (plot) {
            if (is.null(treated)) 
                plot(X[, 1:2])
            else plot(X[, 1:2], col = treated + 2, pch = ifelse(treated, 
                20, 17), cex = 1)
        }
        x$z <- runif(n)
        mod <- rpart(z ~ ., data = x, method = "anova", minsplit = msplit, 
            xval = 0, cp = 0)
        idx <- which(mod$frame$var == "<leaf>")
        n.leaves <- length(idx)
        leaves <- as.integer(dimnames(mod$frame)[[1]])
        LEAF.num <- leaves[idx]
        LEAF.pred <- leaves[mod$where]
        group <- sapply(LEAF.num, function(x) which(LEAF.pred == 
            x))
        if (!is.null(treated) & check.bal) {
            idx.t <- lapply(group, function(x) x[which(treated[x] == 1)])
            idx.c <- lapply(group, function(x) x[which(treated[x] == 0)])
        }
        check.for.bal <- function(g) {
            similar <- unlist(group[[g]])
            residual <- numeric(0)
                if (length(idx.t[[g]]) > 0 & length(idx.c[[g]]) > 
                  0) {
                  mins <- apply(X[idx.t[[g]], ], 2, function(x) min(x, 
                    na.rm = TRUE))
                  maxs <- apply(X[idx.t[[g]], ], 2, function(x) max(x, 
                    na.rm = TRUE))
                  test <- sapply(idx.c[[g]], function(x) prod(mins <= 
                    X[x, ]) * prod(maxs >= X[x, ]))
                  test <- as.logical(test)
                  similar <- c(idx.t[[g]], (idx.c[[g]])[test])
                  residual <- (idx.c[[g]])[!test]
                }
             new.g[[2*g-1]] <<- similar
             new.g[[2*g]] <<- residual
            if (plot) {
                minX <- min(X[similar, 1], na.rm = TRUE)
                maxX <- max(X[similar, 1], na.rm = TRUE)
				minY <- min(X[similar, 2], na.rm = TRUE)
                maxY <- max(X[similar, 2], na.rm = TRUE)
                rect(minX, minY, maxX, maxY, lty = 3)
                if (length(residual) > 0) {
                  minX <- min(X[residual, 1], na.rm = TRUE)
                  maxX <- max(X[residual, 1], na.rm = TRUE)
                  minY <- min(X[residual, 2], na.rm = TRUE)
                  maxY <- max(X[residual, 2], na.rm = TRUE)
                  rect(minX, minY, maxX, maxY, lty = 3)
                }
            }
        }

		new.g <- vector(2*length(group),mode="list")
	    if(check.bal){
		 sapply(1:n.leaves, check.for.bal)
		 group <- new.g
		 }
        addXPtr(RRP, group, -1.0)
		rm(mod)
        if (K%/%100 == K/100) {
            cat("0\n")
        }
        else {
            if (K%/%10 == K/10) 
                cat("+")
            else cat(".")
        }
    }

    tmp <- mulXPtr(RRP, list(1:n), 1.0/Rep);
    t.att <- NULL
	if(!is.null(treated))
	 t.att <-  as.logical(treated)
	 RRPcl <- class(RRP) # XPtrToDist apparently removes class info

	if(asdist)
	 tmp <- XPtrToDist(RRP)

	attributes(tmp) <- list(Size = n, Diag = FALSE, Upper = FALSE, 
	method = "RRP", minsplit = msplit, replications = Rep, 
	cov.cut = cut.in, balanced = check.bal, treated = t.att,
	call = match.call()) 

	if(asdist)
	 class(tmp) <- "dist"
	else
	 class(tmp) <- RRPcl  

     return(invisible(tmp))
}


