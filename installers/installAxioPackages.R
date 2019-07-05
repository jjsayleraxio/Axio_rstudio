
axioPackages<-c("AxioGWAS","AxioNanostring","AxioPC","AxioPlot","AxioQC","AxioSerializer","AxioqPCR","BGSA","BGSAsnp","svmrfe")

install.packages(c("Unicode","fst","stringi","stringr","tidyverse"))
install.packages("flock", repos="http://R-Forge.R-project.org", type="source")

for(axioPackage in axioPackages){
  print(paste("installing",axioPackage))
  install.packages( file.path("axioPackages",axioPackage) , repos = NULL , type = "source" )
}
