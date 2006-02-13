"addDist" <-
function(d, x, k)
{
	if(class(d) == "dist")
		return(.Call("addDist", d, as.integer(sort(unique(x))), k, as.integer(attr(d,"Size")), PACKAGE="rrp"))
	else stop("`d' must be a of class `dist'")				
}

