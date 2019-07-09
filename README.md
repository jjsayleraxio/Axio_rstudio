# Axio_rstudio
Custom rstudio container used at Axio Research

__Bioconductor release_base2 image with R3.5.0 and BioC 3.7__

___

This is a custom rstudio container used by Axio Research. It contains ~Bioconductor packages as well as~ custom built Axio Research R packages.

To run:

```
docker run -dit --name [container name] -p [HOST PORT]:8787 -v [HOST DIR]:[CONTAINER DIR] jjsaxio/axio_rstudio:latest
```

The volume (-v) is optional.
