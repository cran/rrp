"rrp.impute" <-
function(data, k=1,  msplit=10, Rep=250, cut.in=15){

 n <- dim(data)[1]
 k <- dim(data)[2]
 all <- 1:n
 
 D <- rrp.dist(data,  msplit=msplit, Rep=Rep,cut.in=cut.in, check.bal=FALSE, 
 plot=FALSE)
 
 miss.obs <- which(apply(data, 1, function(x) length(which(is.na(x)))>0 ) == TRUE)
 
 D[which(D==1)] <- NA
 M <- as.matrix(D)
 
 y <- data
 comp.obs <- all[-miss.obs]

 for(i in 1:length(miss.obs)){
  tmp <- as.integer(which(M[miss.obs[i], comp.obs] == min(M[miss.obs[i], comp.obs], na.rm=TRUE)))
  tmp <- tmp[1:min(k,length(tmp))]

  v.idx <- which(is.na(data[miss.obs[i],]))	
  for(v in v.idx){
    if(!is.factor(data[[v]])){
	  y[miss.obs[i],v] <- weighted.mean(data[comp.obs[tmp],v], w=1-M[miss.obs[i],comp.obs[tmp]], na.rm=TRUE)
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

