get.richi2 <-
function(n, df, scale){	
	return((df * scale)/rchisq(n, df = df))
	}
