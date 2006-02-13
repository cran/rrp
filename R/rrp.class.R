"rrp.class" <-
function(x, cl, train, test, k=1){

 if(class(x) == "dist"){
  if(attr(x,"method") != "RRP")
   stop("`x' must have attribute method = `RRP'")
 } else {
   stop("`x' must be an object of class `dist'")
 }

 if(attr(x, "Size") != (length(train)+length(test)))
  stop("Wrong dimensions")

 M <- as.matrix(x)
 if(!is.factor(cl))
  cl <- factor(cl)

 pred <- factor(rep(NA,length(test)), levels=levels(cl))
 
 for(i in 1:length(test)){
  tmp <- as.integer(which(M[test[i],train] == min(M[test[i],train], na.rm=TRUE)))
  tmp <- tmp[1:min(k,length(tmp))]

  if(length(tmp)>1){
   votes <- table(cl[tmp])
   idx <- which(votes == max(votes))   
   if(length(idx)>1)
     idx <- sample(idx,1)
	 
    pred.cl <- names(votes)[idx]
    pred[i] <- pred.cl
  } else {
    pred[i] <- cl[tmp]
  }  
 }
 return(pred)
}

