FROM bioconductor/release_base2:R3.6.0_Bioc3.9

LABEL maintainer="josephs@axioresearch.com" \
      description="RStudio docker image used at Axio Research" \
      version="1.0"

RUN apt-get update -y && apt-get -y install git subversion python python3 python-pip python3-pip

WORKDIR /

RUN git clone https://github.com/jjsayleraxio/Axio_rstudio.git

WORKDIR /Axio_rstudio

RUN sh users/userconf.sh && R -f installers/BioC_installr.R --args Sequencing && R -f installers/installAxioPackages.R

ENV PASSWORD=axio
ENV HTTR_LOCALHOST=/etc/R/Renviron.site
ENV HTTR_PORT=/etc/R/Renviron.site

RUN chgrp -R users /usr/local/lib/R/site-library && chmod -R g+w /usr/local/lib/R/site-library 
