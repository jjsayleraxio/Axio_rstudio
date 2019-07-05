BGSA <-
    function(file, fout='outputFile.Rdata', nIter=1000, burnIn=200, get.beta=FALSE){
  
# .Rdata file created by the function read.GSE() and it includes c, y, and setInd
# c is the class indicator, i.e., treated/non treated 
# y are gene expression values 
# setInd is the list of genes (theri indices) assigned to each gene set
# fout is the name of the output file where the posterior samples are saved
# nIter is the number of MCMC samples
# burnIn is the number of burn-in samples
  
# output includes the following elements. 
# post.J: J is a vector of binary latent variables such that J=0 => H_0 is true, J=1 => H_1 is true. post.J is the posterior mean of J as a measure for P(H_1 | data)
# p.val: This is an alternative measure (comparable to p-value) that we can use to decide whether a gene set is significant
# post.beta: this includes the posterior samples for beta's if get.beta was set to TRUE. 
  
  load(file)
  
  c = labels
  y = scale(data)
  setInd = setInd
  
  n = dim(y)[1]
  p = dim(y)[2]
  
  n1 = sum(c==1)
  n2 = sum(c==2)
  
  nSet = length(setInd)
  
# Adding a column of 1's for the intercept
  x = cbind(rep(1, n), c-1); 
  
# d is the number of regression parameters including the intercept.
  d = dim(x)[2]
  
  X = t(x)%*%x;
  L.x = chol(X);
  invX = chol2inv(L.x);
  L.x.inv = chol(invX);
  
# This give the mean of the posterior distribution.
  BETA.hat = invX%*%t(x)%*%y;
  
  x.a = matrix(x[, 1:(d-1)], n, d-1, byrow=TRUE)
  x.b = matrix(x[, d], n, 1)
  
  X.a = t(x.a)%*%x.a;
  L.a = chol(X.a);
  invX.a = chol2inv(L.a);
  L.a.inv = chol(invX.a);
  
  
  X.b = t(x.b)%*%x.b;
  L.b = chol(X.b);
  invX.b = chol2inv(L.b);
  L.b.inv = chol(invX.b);
  
  
# The following two matrices hold the posterior samples for sigma2 and beta.
  alpha = matrix(BETA.hat[1:(d-1), ], nrow=d-1, ncol=p, byrow=TRUE)
  beta = matrix(BETA.hat[d, ], 1, p)
  sigma2 = colSums((y - x%*%BETA.hat)^2)/(n-1)
  
# Here, we creat a list of sets  
  sets = list()
  G = NULL
  post.beta = list()
  post.beta.p = list()
  count=0
  for (i in 1:nSet){ 	
    
    indSet = list()
    ind = setInd[[i]]
    indSet$ind = ind[!is.na(ind)]	
    
    indSet$p = length(indSet$ind)  
    indSet$y = y[, indSet$ind]
    indSet$n = length(c)
    indSet$alpha = alpha[, indSet$ind]
    indSet$beta = beta[, indSet$ind]
    indSet$sigma2 = sigma2[indSet$ind]
    indSet$tau2 = 1
    count = count+1 
    post.beta[[count]] = indSet$beta	
    post.beta.p[[count]] = rep(1, indSet$p)	
    sets[[count]] = indSet
    
  }
  
  nSet = length(sets)    
  
# These are for storing posterior samples 
  
# J is a vector of binary latent variables such that J=0 => H_0 is true, J=1 => H_1 is true. We can use the posterior mean of J as a measure for P(H_1 | data)
  post.J = matrix(NA, nrow=nIter, ncol = nSet)
  post.tau2 = matrix(NA, nrow=nIter, ncol = nSet)
  post.lambda = rep(NA, nIter)
  post.sc = matrix(NA, nIter, 2)
  post.nu = rep(NA, nIter)
  
# This is an alternative measure (comparable to p-value) that we can use to decide whether a gene set is significant
  p.val =  matrix(NA, nrow=nIter, ncol = nSet)
  
#Initial values for MCMC and fixed parameters
  lambda = 0.5	# lambda is P(H_1)
  J = rbinom(nSet, 1, lambda) 	
  
  mu0 = 0	
  
# lambda ~ Beta(aBeta, bBeta)			
  aBeta = 1
  bBeta = 1
  
# tau^2 ~ (1-lambda) x Inv-chi2(nu, sc0)	+ lambda x Inv-chi2(nu, sc0+sc1)
  nu0 = 5
  sc0 = 0.1
  sc1 = 0.2
  
  muNu = 1
  sigNu = 1
  muSc = 1
  sigSc = 1
  
  samp = list()
  
#MCMC samples
  for (iter in 1:nIter){
    print(iter)		
    tau2.0 = NULL
    tau2.1 = NULL
    
    
    for(i in 1:nSet){
      
      y = sets[[i]]$y
      alpha = sets[[i]]$alpha
      beta = sets[[i]]$beta
      sigma2 = sets[[i]]$sigma2
      tau2 = sets[[i]]$tau2
      p = sets[[i]]$p			
      
      y.a = y - x.b%*%beta
      
      alpha.hat = invX.a%*%t(x.a)%*%y.a
      
      sigma2.mat = matrix(sigma2, d-1, p, byrow=TRUE)
      
      u = matrix(rnorm((d-1)*p), p, d-1);
      alpha = (t(u%*%L.a.inv) + alpha.hat/sqrt(sigma2.mat))*sqrt(sigma2.mat)
      
      y.b = y - x.a%*%alpha
      y.b = y.b[c==2, ]
      
      V = 1/(1/(tau2*sigma2) + n2/sigma2)
      mu_n = V*(mu0/(tau2*sigma2) + colSums(y.b)/sigma2) 
      sigma_n = sqrt(V)
      beta = rnorm(p, mean = mu_n, sd = sigma_n)
      
      BETA = rbind(alpha, beta)
      
      eps = y - x%*%BETA
      eps.bar = colMeans(eps)
      
      nu_n = n-1
      sigma02_n = colSums( (eps - eps.bar)^2 ) /(n-1)
      z = rchisq(p, nu_n);
      sigma2 = nu_n*sigma02_n/z;
      
      beta.scaled = beta/sqrt(sigma2)
      nu_n = nu0 + p
      V = mean(beta.scaled^2)
      
      if(J[i]==0){
        sigma02_n = (nu0*sc0+p*V)/(nu0+p)
        z = rchisq(1, nu_n)
        tau2 = nu_n*sigma02_n/z
        tau2.0 = c(tau2.0, tau2)
      }else{
        sigma02_n = (nu0*(sc0+sc1)+p*V)/(nu0+p)
        z = rchisq(1, nu_n)
        tau2 = nu_n*sigma02_n/z
        tau2.1 = c(tau2.1, tau2)
      }	  							
      
      sets[[i]]$alpha = alpha
      sets[[i]]$beta = beta
      sets[[i]]$sigma2 = sigma2
      sets[[i]]$tau2 = tau2
      
      post.beta[[i]] = rbind(post.beta[[i]], sets[[i]]$beta)
      post.tau2[iter, i] = sets[[i]]$tau2
      
      
    }
    
    
    
    thisSc0 = NULL
    if(is.null(tau2.0)){
      thisSc0 = rgamma(20, muSc, sigSc)
      sc0 = thisSc0[1]
    }else{
      for(rep in 1:20){	
        sc0 = getSigTau2(c(tau2.0, tau2.1), nu0, sc0, c(rep(0, length(tau2.0)), rep(sc1, length(tau2.1))), muSc, sigSc)
        thisSc0[rep] = sc0
      }
    }
    
    thisSc1 = NULL
    if(is.null(tau2.1)){
      thisSc1 = rgamma(20, muSc, sigSc)
      sc1 = thisSc1[1]
    }else{
      for(rep in 1:20){	
        sc1 = getSigTau2(tau2.1, nu0, sc1, sc0, muSc, sigSc)
        thisSc1[rep] = sc1
      }
    }
    
    thisNu0 = NULL
    for(rep in 1:20){	
      nu0 = getNuTau2(c(tau2.0, tau2.1), nu0, c(rep(sc0, length(tau2.0)), rep(sc0+sc1, length(tau2.1))), muNu, sigNu)
      thisNu0[rep] = nu0
    }
    
    post.sc[iter, ] = c(sc0, sc1)
    post.nu[iter] = nu0
    
    beta.H0 = rnorm(1000, 0, sqrt(get.richi2(1000, nu0, sc0)))
    
    for(i in 1:nSet){
      q0 = get.dichi2.log(sets[[i]]$tau2, nu0, sc0)+log(1-lambda)
      q1 = get.dichi2.log(sets[[i]]$tau2, nu0, sc0+sc1)+log(lambda)
      m = apply(cbind(q0, q1), 1, max)
      A = m + log(exp(q0-m) + exp(q1-m))
      prob.1 = exp(q1 - A)
      picked = rbinom(1, 1, prob.1)
      if (picked == 0){
        J[i] = 0    	
      }else{
        J[i] = 1	
      }
    }
    
    post.J[iter, ] = J
    
    if (iter > 100){
      newNj1 = sum(J==1)
      lambda = rbeta(1, aBeta + newNj1, bBeta + nSet - newNj1)
    }
    post.lambda[iter] = lambda
    
    p.val[iter, ] = pchisq(nu0*sc0/post.tau2[iter, ], nu0)
    
    if(get.beta){ 
      for(s in 1:nSet){ 
        p.val.beta = colMeans(matrix(abs(beta.H0), 1000, sets[[s]]$p) > abs(matrix(sets[[s]]$beta/sqrt(sets[[s]]$sigma2), 1000, sets[[s]]$p, byrow=TRUE)) )
        post.beta.p[[s]] = rbind(post.beta.p[[s]], p.val.beta)
      }
    }
    
  }	
  
  
  if(get.beta){ 
    for(s in 1:nSet){ 
      post.beta.p[[s]] = colMeans(post.beta.p[[s]])
    }
  }
  
  
  return(list(post.J = colMeans(post.J[burnIn:nIter, ]), p.val = colMeans(p.val[burnIn:nIter, ]), post.beta = post.beta))
  
}
