# Axio_rstudio
Custom rstudio container used at Axio Research

This is a custom rstudio container used by Axio Research. It contains ~Bioconductor packages as well as~ custom built Axio Research R packages.

To run:

```
docker run -dit --name [container name] -p [HOST PORT]:8787 -v [HOST DIR]:[CONTAINER DIR] jjsaxio/axio_rstudio:latest
```

The volume (-v) is optional.
