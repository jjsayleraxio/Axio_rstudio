read.GSE <-
function(expressionFile, labelsFile, outFile = 'gse.Rdata', BM.attributes = c("affy_hg_u133_plus_2","hgnc_symbol"), BM.filters = "affy_hg_u133_plus_2", gmt.file='c5mf.txt', lower=5, upper=200){
	
# expressionFile includes the gene expression values. Each row corresponds to one probe, and each column corresponds to one sample. The first column should have the probe ids.
# labelsFile should include one line providing the corresponding group (e.g., case/control) for each sample
# outfile is the name of the output file, which must have .Rdata extension. The output file will include the gene expression data (called y), the labels, the gene names, the set identifers (setInd), and the name of pathways. 
# BM.attributes and BM.filters specity the attributes and filters for the getBM function of the biomaRt package
# gmt.file, is a text file obtained from the GSEA database; it should include the pathways.
# lower and upper are the minimum and maximum number of genes assigned to a pathway. Pathways that include fewer genes or more genes will be excluded from the analysis. 

data <- read.table(expressionFile)
probe <- data[, 1]
data <- t(data[, -1])
labels <- as.numeric(as.factor(unlist(read.table(labelsFile))))


ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

map = getBM(attributes=BM.attributes,filters=BM.filters,values=probe, mart= ensembl)


m = match(probe, map[, 1])
genenames = toupper(as.character(map[m, 2]))


geneset.obj<- GSA.read.gmt(gmt.file)
pathways <- geneset.obj$geneset.names

setInd = lapply(geneset.obj$genesets, match, genenames)

allInd = unlist(setInd)
assignedInd = allInd[!is.na(allInd)]
u.assignedInd = unique(assignedInd)

data = data[, u.assignedInd]
genenames = genenames[u.assignedInd]

setInd = lapply(geneset.obj$genesets, match, genenames)

nSet = length(setInd)
l.set = NULL
for(i in 1:nSet){
	l.set[i] = sum(!is.na(setInd[[i]]))
}

setInd2 = list()
pathnames = NULL
counter=0
for(i in 1:nSet){
	if(l.set[i]>=lower & l.set[i] <= upper){ # only including pathways with more than 10 genes
	counter = counter + 1	
	setInd2[[counter]] <- setInd[[i]][!is.na(setInd[[i]])]
	pathnames[counter] <- pathways[i]
	}
}
setInd = setInd2


allInd = unlist(setInd)
assignedInd = allInd[!is.na(allInd)]
u.assignedInd = unique(assignedInd)

data = data[, u.assignedInd]
genenames = genenames[u.assignedInd]

setInd = lapply(geneset.obj$genesets, match, genenames)

nSet = length(setInd)
l.set = NULL
for(i in 1:nSet){
	l.set[i] = sum(!is.na(setInd[[i]]))
}

setInd2 = list()
pathnames = NULL
counter=0
for(i in 1:nSet){
	if(l.set[i]>=lower & l.set[i] <= upper){ # only including pathways with more than 10 genes
	counter = counter + 1	
	setInd2[[counter]] <- setInd[[i]][!is.na(setInd[[i]])]
	pathnames[counter] <- pathways[i]
	}
}
setInd = setInd2

nSet = length(setInd)
l.set = NULL
for(i in 1:nSet){
	l.set[i] = sum(!is.na(setInd[[i]]))
}

save(y = data, labels, genenames, setInd, pathnames, file = outFile)



}
