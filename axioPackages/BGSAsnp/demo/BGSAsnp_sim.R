# This code performs two simulations. First, it takes 2500 SNP related to 505 genes, which are assigned 50 pathways, and use GSEAsnp, ALIGATOR, and BGSAsnp to perform pathway analysis. For the second simulation, it removes 150 of the genes, assumes that their correct class can be predicted with probability 0.7, and apply GSEAsnp, ALIGATOR, BGSAsnp without prediction, and BGSAsnp with prediction. The final results (incuding AUC) are saved as res.Rdata in the current directory. auc1 is a matrix for the results of the first simulation, and auc2 is another matrix for the results of the second simulation. 

rm(list=ls())
library(caTools)


# These functions are obtained from SNPath package for running GSEAsnp and ALIGATOR. You need to obtain 'miscfun.r', 'gseaSnp.r', and 'aligator.r' from the SNPath package available online at http://linchen.fhcrc.org/grass.html, and save them in the current directory. 


wd = getwd()
source(paste(wd, 'miscfun.r', sep='/'))
source(paste(wd, 'gseaSnp.r', sep='/'))
source(paste(wd, 'aligator.r', sep='/'))
source(paste(wd, 'BGSAsnp.R', sep='/'))
source(paste(wd, 'BGSAsnpmiss.R', sep='/'))


# This is the function to get the SNP level statistics
get.snp.z <- function(snp.dat, y){
freq=rowMeans(snpDat/2)
r0=rowSums(snpDat[,y==1]==0); s0=rowSums(snpDat[,y==0]==0)
r1=rowSums(snpDat[,y==1]==1); s1=rowSums(snpDat[,y==0]==1)
r2=rowSums(snpDat[,y==1]==2); s2=rowSums(snpDat[,y==0]==2)
r=r0+r1+r2; s=s0+s1+s2
n0=r0+s0; n1=r1+s1; n2=r2+s2; n=n0+n1+n2
z.stat=((s*r1-r*s1)+2*(s*r2-s2*r))/sqrt(r*s/n)/sqrt(n*(n1+4*n2)-(n1+2*n2)^2)
return(z.stat)
}


# This is the function to get the gene level summary measures
get.gene.z <- function(snp.z, gene.snplist){
	nu0 <- 0
	sc0 <- 1
	
	nGene = length(gene.snplist)
	gene.z <- rep(NA, nGene)
	
	for(i in 1:nGene){
			
	beta = snp.z[ gene.snplist[[i]] ]
	p = length(beta)
	 V = mean(beta^2)
	 gene.z[i] = sign(beta[which.max(abs(beta))])*((nu0*sc0+V*p)/(nu0+p-2))
	}
	
	return(gene.z)
}



# Preparing the data for simulation
gene.info <- read.table(paste(wd, 'gene.info.txt', sep='/'), header=TRUE)
snp.info <- read.table(paste(wd, 'snp.info.txt', sep='/'), header=TRUE)
snpDat <- as.matrix(read.table(paste(wd, 'snpDatnew.txt', sep='/'), header=FALSE))

gene.snplist <- Snp2Gene.Adist(snp.info, gene.info)
snp.genelist <- gene2snp.Adist(snp.info, gene.info)

gene.names <- gene.info$gene.Name
nGene <- length(gene.names)


nSim <- 100
auc1 = matrix(NA, nSim, 3)
res1 = list()

auc2 = matrix(NA, nSim, 4)
res2 = list()

counter = 1

for(i in 1:nSim){
	
print(i)

set.seed(i)	

nSet <- 50
setInd <- sort(sample(1:50, 500, replace=TRUE))
sim.pathway=vector("list",nSet)
for(j in 1:nSet){
 	sim.pathway[[j]]= gene.names[setInd==j]
}
names(sim.pathway) <- paste('path', seq(1, 50), sep='')

pathway.genelist <- Gene2Set(gene.names, sim.pathway)

sig.genes <- NULL
for(jj in 1:5){
temp <- unlist(pathway.genelist[jj])	
sig.genes <- c(sig.genes, sample(temp, min(3, length(temp))))
}
ind <- unique(unlist(gene.snplist[sig.genes]))
samp.ind <- sample(ind, 10)
x <- snpDat[samp.ind, ]


beta <- rbeta(dim(x)[1], 2, 8)
beta <- (-1)^(rbinom(length(beta), 1, 0.5))*beta

eta <- as.vector(beta%*%x)
mu <- exp(eta)/(1+exp(eta))
y <- round(mu)


while(mean(y)<0.1 | mean(y)>0.9){
set.seed(i+1000*counter)
beta <- rbeta(dim(x)[1], 2, 8)
beta <- (-1)^(rbinom(length(beta), 1, 0.5))*beta
eta <- as.vector(beta%*%x)
mu <- exp(eta)/(1+exp(eta))
y <- round(mu)
counter = counter+1
}


snp.z <- get.snp.z(snpDat, y)	

z <- get.gene.z(snp.z, gene.snplist)

p.gsea <- gseaSnp(snp.dat=snpDat, snp.info=snp.info, gene.info=gene.info, gene.set=sim.pathway, y=y, snp.method="chiSq", gene.def="abs",dist=5, B=100)


pval <- calc.fun(snp.dat=snpDat, y=y, snp.method="chiSq")$pval
p.aligator <- aligator(snp.info=snp.info, gene.info=gene.info, gene.set=sim.pathway, snp.pval=pval, gene.def="rel", dist=5)


simRes = BGSA.snp(z, setInd=pathway.genelist, nIter = 5500, burnIn = 500)
p.bgsa = simRes$p.val


r = rep(0, nSet)
r[1:5] = 1

p = t(p.gsea) 	
auc1[i, 1] = colAUC(p, r)

p = t(p.aligator) 	
auc1[i, 2] = colAUC(p, r)

p = t(p.bgsa) 	
auc1[i, 3] = colAUC(p, r)

print(auc1[i, ])

temp <- list()

temp$p.gsea <- p.gsea
temp$p.aligator <- p.aligator
temp$p.bgsa <- p.bgsa
temp$auc <- auc1[i, ]
temp$y = y
temp$pathway <- pathway.genelist
temp$beta = beta

res1[[i]] <- temp


n.rm <- 150
unknown.genes <- matrix(NA, n.rm, 3)
names.unknown.genes <- rep('', n.rm)
for(j in 1: n.rm){
	temp <- sample(50, 1)
	while(length(sim.pathway[[temp]])<3){
			temp <- sample(50, 1)	
		}
		names.unknown.genes[j] <- as.character(sim.pathway[[temp]][1])
		sim.pathway[[temp]] <- sim.pathway[[temp]][-1]
		unknown.genes[j, 1] <- temp
		unknown.genes[j, 2:3] <- sample(seq(1, 50)[-temp], 2)
}
row.names(unknown.genes) <- names.unknown.genes
unknown.genes.prob <- matrix(c(0.7, .2, .1), n.rm, 3, byrow=TRUE)

unknown.genes.z = rep(NA, length(names.unknown.genes))
for(j in 1:length(names.unknown.genes)){
ind <- which(!is.na(match(gene.names, names.unknown.genes[j])))
if(length(ind)>0){
unknown.genes.z[j] <- z[ind] 	
}
}
ind = which(!is.na(unknown.genes.z))
unknown.genes.z <- unknown.genes.z[ind]
unknown.genes <- unknown.genes[ind, ]
unknown.genes.prob <- unknown.genes.prob[ind, ]


gene.names2 <- setdiff(gene.names, rownames(unknown.genes))

ind <- which(!is.na(match(gene.info$gene.Name, gene.names2)))

gene.info2 <- gene.info[ind, ]

gene.snplist2 <- gene.snplist[ind]

ind <- unlist(gene.snplist[ind])

snp.info2 <- snp.info[ind, ]

snpDat2 <- snpDat[ind, ]



pathway.genelist2 <- Gene2Set(gene.names2, sim.pathway)


snp.z2 <- get.snp.z(snpDat2, y)	

z2 <- get.gene.z(snp.z2, gene.snplist2)


p.gsea <- gseaSnp(snp.dat=snpDat2, snp.info=snp.info2, gene.info=gene.info2, gene.set=sim.pathway, y=y, snp.method="chiSq", gene.def="abs",dist=5, B=100)


pval <- calc.fun(snp.dat=snpDat2, y=y, snp.method="chiSq")$pval
p.aligator <- aligator(snp.info=snp.info2, gene.info=gene.info2, gene.set=sim.pathway, snp.pval=pval, gene.def="rel", dist=5)


simRes = BGSA.snp(z2, setInd=pathway.genelist2, nIter = 5500, burnIn = 500)
p.bgsa = simRes$p.val



simRes = BGSA.snpmiss(z2, setInd=pathway.genelist2, unknown.genes=unknown.genes, unknown.genes.prob = unknown.genes.prob, unknown.genes.z, nIter = 5500, burnIn = 500)
p.bgsa.miss = simRes$p.val



p = t(p.gsea) 	
auc2[i, 1] = colAUC(p, r)

p = t(p.aligator) 	
auc2[i, 2] = colAUC(p, r)

p = t(p.bgsa) 	
auc2[i, 3] = colAUC(p, r)


p = t(p.bgsa.miss) 	
auc2[i, 4] = colAUC(p, r)


print(auc2[i, ])

temp <- list()

temp$p.gsea <- p.gsea
temp$p.aligator <- p.aligator
temp$p.bgsa <- p.bgsa
temp$p.bgsa.miss <- p.bgsa.miss
temp$auc <- auc2[i, ]
temp$y = y
temp$pathway <- pathway.genelist
temp$beta = beta

res2[[i]] <- temp

save(res1, res2, auc1, auc2, file=paste(wd, 'res.Rdata', sep='/'))

}


