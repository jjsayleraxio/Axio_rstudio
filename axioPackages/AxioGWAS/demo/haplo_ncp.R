haplo.power.qt.ncp<-
function (haplo, haplo.freq, base.index, haplo.beta, y.mu, y.var) 
{
    haplo <- as.matrix(haplo)
    n.loci <- ncol(haplo)
    hapPair.lst <- get.hapPair(haplo, haplo.freq, base.index)
    p.g <- hapPair.lst$p.g
    x.haplo <- hapPair.lst$x.haplo
    haplo.indx <- hapPair.lst$haplo.indx
    df <- ncol(x.haplo)
    x.mu <- apply(x.haplo * p.g, 2, sum)
    x.mu.mat <- matrix(rep(x.mu, nrow(x.haplo)), nrow = nrow(x.haplo), 
        byrow = TRUE)
    delta.x <- x.haplo - x.mu.mat
    t2 <- delta.x * p.g
    vx.complete <- t(t2) %*% delta.x
    geno.mat <- NULL
    for (i in 1:n.loci) {
        t1 <- haplo[haplo.indx[, 1], i]
        t2 <- haplo[haplo.indx[, 2], i]
        a1 <- ifelse(t1 < t2, t1, t2)
        a2 <- ifelse(t2 > t1, t2, t1)
        geno.mat <- cbind(geno.mat, a1, a2)
    }
    geno.hash <- haplo.hash(geno.mat)
    haplo.group <- geno.hash$hash
    ord <- order(haplo.group)
    p.g <- p.g[ord]
    haplo.group <- haplo.group[ord]
    x.haplo <- as.matrix(x.haplo[ord, ])
    p.haplo.group <- tapply(p.g, haplo.group, sum)
    nrep <- tapply(haplo.group, haplo.group, length)
    denom <- rep(p.haplo.group, nrep)
    post <- p.g/denom
    n.group <- length(unique(haplo.group))
    nx <- ncol(x.haplo)
    vx.incomplete <- matrix(numeric(nx^2), nrow = nx)
    for (i in 1:n.group) {
        zed <- (haplo.group == i)
        tmp.x.mu <- as.vector(apply(x.haplo[zed, , drop = FALSE] * 
            post[zed], 2, sum))
        tmp.delta <- (tmp.x.mu - x.mu)
        vx.incomplete <- vx.incomplete + ((tmp.delta %o% tmp.delta) * 
            p.haplo.group[i])
    }
    r2.phase.unknown <- (t(haplo.beta[-base.index]) %*% vx.incomplete %*% 
        haplo.beta[-base.index])/y.var
    r2.phase.known <- (t(haplo.beta[-base.index]) %*% vx.complete %*% 
        haplo.beta[-base.index])/y.var
    ncp.chi.phase.known <- r2.phase.known
    ncp.chi.phase.unknown <- r2.phase.unknown
    ncp.f.phase.known <- r2.phase.known/(1 - r2.phase.known)
    ncp.f.phase.unknown <- r2.phase.unknown/(1 - r2.phase.unknown)
    return(list(ncp.chi.phased.haplo = ncp.chi.phase.known, ncp.f.phased.haplo = ncp.f.phase.known, 
        ncp.chi.unphased.haplo = ncp.chi.phase.unknown, ncp.f.unphased.haplo = ncp.f.phase.unknown, 
        df = df))
}

haplo.power.qt.1 <-function (haplo, haplo.freq, base.index, haplo.beta, y.mu, y.var, alpha, sample.size = NULL, power = NULL) 
{
    if (is.null(power) & is.null(sample.size)) {
        stop("Must specify either power or sample.size")
    }
    if (!is.null(power) & !is.null(sample.size)) {
        stop("Must specify only one of power or sample.size")
    }
    ncp <- haplo.power.qt.ncp.1(haplo = haplo, haplo.freq = haplo.freq, 
        base.index = base.index, haplo.beta = haplo.beta, y.mu = y.mu, 
        y.var = y.var)
    if (is.null(sample.size) & !is.null(power)) {
        ss.phased.haplo <- f.sample.size(nc = ncp$ncp.f.phased.haplo, 
            df1 = ncp$df, alpha = alpha, power = power)
        ss.unphased.haplo <- f.sample.size(nc = ncp$ncp.f.unphased.haplo, 
            df1 = ncp$df, alpha = alpha, power = power)
        power.phased.haplo <- power
        power.unphased.haplo <- power
    }
    if (is.null(power) & !is.null(sample.size)) {
        power.phased.haplo <- f.power(n = sample.size, nc = ncp$ncp.f.phased.haplo, 
            df1 = ncp$df, alpha = alpha)
        power.unphased.haplo <- f.power(n = sample.size, nc = ncp$ncp.f.unphased.haplo, 
            df1 = ncp$df, alpha = alpha)
        ss.phased.haplo <- sample.size
        ss.unphased.haplo <- sample.size
    }
    return(list(ss.phased.haplo = ss.phased.haplo, ss.unphased.haplo = ss.unphased.haplo, 
        power.phased.haplo = power.phased.haplo, power.unphased.haplo = power.unphased.haplo))
}



load( "t:/Software_Projects/GeneticsPC/data/.RData")
library(haplo.stats)

haplo.freq.1 <- c( 0.7 , 0.3 )

tmp <- find.haplo.beta.qt(haplo.1,haplo.freq.1,base.index.1,haplo.risk.1, r2=0.5, y.mu=0, y.var=30)
haplo.beta.1 <- tmp$beta
haplo.power.qt(haplo.1, haplo.freq.1, base.index.1, haplo.beta.1, y.mu=0, y.var=30, alpha=.05, sample.size = 100)  

rsq<-seq(from=0, to=0.3, by= 0.005)
ss<-c(48, 75, 100)
#ss<-ss/2
pow<-matrix(0,61,3)
for(i in 1:61)
{
  tmp <- find.haplo.beta.qt(haplo.1,haplo.freq.1,base.index.1,haplo.risk.1, r2=rsq[i], y.mu=0, y.var=1)
  haplo.beta.1 <- tmp$beta  
  for(j in 1:3)
  {
    pow[i,j]<-haplo.power.qt(haplo.1, haplo.freq.1, base.index.1, haplo.beta.1, y.mu=0, y.var=30^2, alpha=.05, sample.size = ss[j])$power.phased.haplo
  }
}

png("t:/2008/Powercalcs/Xia/rsquare.png")
plot(rsq,pow[,3], xlab="R^2", ylab="Power", main="Power vs. R^2 for Sample Size = 48, 75, 100", type="l" , col = 1)
lines(rsq, pow[,2], lty=2 , col = 2)
lines(rsq, pow[,1], lty=3 , col = 3)
legend( "topleft" , c("100","75","48") , lty=c(1,2,3) , col = c(1,2,3) , bty = "n" )
dev.off()



haplo.power.qt.ncp.1(haplo.1, haplo.freq.1, base.index.1, haplo.beta.1, y.mu=0, y.var=30^2)