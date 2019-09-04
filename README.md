<img align="right" src="https://raw.githubusercontent.com/jjsayleraxio/AxioShiny/master/images/axio-logo.png">
<br>

# Axio_rstudio
![GitHub release (latest by date)](https://img.shields.io/github/v/release/jjsayleraxio/Axio_rstudio?logo=github&style=flat)
[![GitHub issues](https://img.shields.io/github/issues/jjsayleraxio/Axio_rstudio?logo=github&style=flat)](https://github.com/jjsayleraxio/Axio_rstudio/issues)
![Docker Build Status](https://img.shields.io/docker/build/jjsaxio/axio_rstudio?logo=docker&style=flat)

Custom rstudio container used at Axio Research

__Bioconductor release_base2 image with R3.6.0 and BioC 3.9__
___
This is a custom rstudio container used by Axio Research. It contains ~Bioconductor packages as well as~ custom built Axio Research R packages.

To run:

```
docker run -dit --name [container name] -p [HOST PORT]:8787 -v [HOST DIR]:[CONTAINER DIR] jjsaxio/axio_rstudio:latest
```

The volume (-v) is optional.
