# updated to use new BiocManager instead of BiocInstaller

#Set so that an external arguments to this script are the list of views to use to identify bioconductor libraries to keep
tags <- commandArgs(trailingOnly = TRUE)
print(tags)

#install most recent version of Bioconductor

url <- "http://www.bioconductor.org/packages/release/bioc"

#always install these packages - installed in parent container
builtins <- c("Matrix", "KernSmooth", "mgcv")

t <- tempfile()
download.file(paste0(url,"/VIEWS"), t)
dcf <- as.data.frame(read.dcf(t), stringsAsFactors=FALSE)

BioC_Packages<-dcf$Package

FindPackageWithTag<-function(PackageName,PackageDF,Tag){
  PackageTags<-gsub("\n","",PackageDF[which(PackageDF$Package==PackageName),"biocViews"])
  if(!is.na(PackageTags)){
    tagsList<-strsplit(PackageTags,split = ", ")[[1]]
    if(any(tagsList%in%Tag)){
      return(PackageName)
    }
  }
}

bioCpackages<-c(do.call('rbind',lapply(BioC_Packages,FindPackageWithTag,dcf,tags)))

ap.db <- available.packages(contrib.url(BiocManager::repositories()))
ap <- rownames(ap.db)

#filter for only available packages and union with builtins
bioCpackages<-union(bioCpackages[which(bioCpackages%in%ap)],builtins)

#filter out already installed packages
if(any(bioCpackages%in%rownames(installed.packages()))){
  bioCpackages<-bioCpackages[-which(bioCpackages%in%rownames(installed.packages()))]
}

#number of bioconductor packages to install
length(bioCpackages)

#install packages - updated to use BiocManager
BiocManager::install(bioCpackages)

#Check that installed packages are consistent with R and Bioconductor version installed
suppressWarnings(BiocManager::valid(fix=TRUE, ask=FALSE))

