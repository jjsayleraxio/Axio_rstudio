
axioPackages<-c("AxioGWAS","AxioNanostring","AxioPC","AxioPlot","AxioQC","AxioSerializer","AxioqPCR","BGSA","BGSAsnp","svmrfe")
prerequisites<-c('haplo.stats','snpStats','car','mlogit','HH','snowfall','biomaRt','RSQLite','bit','NanoStringNorm','XML','Biobase','pwr','hexbin','XLConnect','RSQLite','blob','lz4','GSA','biomaRt','e1071','randomForest','ipred','cglNanoString','GeneticsPlot','cglQC')

install.packages(c("Unicode","fst","stringi","stringr","tidyverse"))
install.packages("flock", repos="http://R-Forge.R-project.org", type="source")
install.packages(prerequisites)

for(axioPackage in axioPackages){
  print(paste("installing",axioPackage))
  install.packages( file.path("axioPackages",axioPackage) , repos = NULL , type = "source" )
}
