"rrp.predict" <-
function(x, y, train, test, k=1){

 if(class(x) == "dist"){
  if(attr(x,"method") != "RRP")
   stop("`x' must have attribute method = `RRP'")
 } else {
   stop("`x' must be an object of class `dist'")
 }

 if(attr(x, "Size") != (length(train)+length(test)))
  stop("Wrong dimensions")

 M <- as.matrix(x)
 if(!is.numeric(y))
  stop("`y' must be numeric")

 pred <- as.numeric(rep(NA,length(test)))
 
 for(i in 1:length(test)){
  tmp <- as.integer(order(M[test[i],train]))
  tmp <- tmp[1:min(k,length(tmp))]
  pred[i] <- weighted.mean(y[tmp], w=1-M[test[i],train[tmp]], na.rm=TRUE)
 }
 return(pred)
}

