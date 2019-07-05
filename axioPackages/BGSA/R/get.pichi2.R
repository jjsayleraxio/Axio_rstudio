get.pichi2 <-
function(x, df, scale){
	nu <- df/2
	return( gamma_inc(nu, scale*nu/x) /gamma(nu) )
	}
