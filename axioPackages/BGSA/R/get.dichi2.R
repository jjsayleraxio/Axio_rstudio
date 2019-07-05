get.dichi2 <-
function(x, df, scale){
	nu <- df/2
	return( (((nu)^(nu))/gamma(nu)) * (scale^nu) * (x^(-(nu + 1))) * exp(-nu * scale/x) )
	}
