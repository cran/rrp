"rrp.impute" <-
function(data, D = NULL, k=1,  msplit=10, Rep=250, cut.in=15){

 n <- dim(data)[1]
 all <- 1:n
 if(is.null(D) | (!any(class(D) == "XPtr")))
	d <- rrp.dist(data,  msplit=msplit, Rep=Rep,cut.in=cut.in, check.bal=FALSE, 
		plot=FALSE)

 miss.obs <- which(apply(data, 1, function(x) length(which(is.na(x)))>0 ) == TRUE)
 D <- XPtrToDist(d) 
 D[which(D==1)] <- NA
 
 y <- data
 comp.obs <- all[-miss.obs]

 f <- function(x) {tmp <- order(x); tmp[1:min(k,length(tmp))]}
 nn <- applyXPtr(d, miss.obs, comp.obs, f)

 g <- function(xx) {tmp <- order(xx); 
   tmp <- tmp[1:min(k,length(tmp))]; xx[tmp]}
 ww <- applyXPtr(d, miss.obs, comp.obs, g)

 for(i in 1:length(miss.obs)){
  tmp <- nn[[i]] 
  wt <- ww[[i]]

  v.idx <- which(is.na(data[miss.obs[i],]))	
  for(v in v.idx){
    if(!is.factor(data[[v]])){
	  y[miss.obs[i],v] <- weighted.mean(data[comp.obs[tmp],v], w=1-wt, na.rm=TRUE)
	} else {
	
	  if(length(tmp)>1){
       votes <- table(data[comp.obs[tmp],v])
       idx <- which(votes == max(votes))   
       if(length(idx)>1)
        idx <- sample(idx,1)
	     y[miss.obs[i],v] <- names(votes)[idx]
      } else {
       y[miss.obs[i],v] <- data[comp.obs[tmp],v]
      }  
   
	}
   }
  }
  return(list(new.data = y, dist = D))
}

