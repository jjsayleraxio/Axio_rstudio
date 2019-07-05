FROM bioconductor/release_base2:R3.6.0_Bioc3.9

LABEL maintainer="josephs@axioresearch.com" \
      description="RStudio docker image used at Axio Research" \
      version="1.0"

RUN apt-get update -y && apt-get -y install git subversion && mkdir /app/

COPY installers/ /tmp/

COPY axioPackages/ /app/

COPY users/userconf.sh /tmp/

RUN sh /tmp/userconf.sh && chgrp -R users /usr/local/lib/R/site-library && chmod -R g+w /usr/local/lib/R/site-library #&& R -f /tmp/BioC_installr.R --args Sequencing && R -f /tmp/installAxioPackages.R

ENV PASSWORD=axio
ENV HTTR_LOCALHOST=/etc/R/Renviron.site
ENV HTTR_PORT=/etc/R/Renviron.site
