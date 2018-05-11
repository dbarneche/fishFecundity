# Fish reproductive-energy output increases disproportionately with body size

This repository contains code and data needed to reproduce the article:

**Barneche DR, Robertson DR, White CR, Marshall DJ** (accepted) Fish reproductive-energy output increases disproportionately with body size. *Science*. [doi: 10.1126/science.aao6868](http://science.sciencemag.org/content/360/6389/642.full)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1213118.svg)](https://doi.org/10.5281/zenodo.1213118)

## Instructions

All analyses were done in `R`. To compile the paper, including figures and tables we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies = TRUE)
```
(run `install.packages("devtools")` to install devtools if needed.)

The `remake` package also depends on `storr`, install it like this:
```r
devtools::install_github("richfitz/storr", dependencies = TRUE)
```

Next you need to open an R session with working directory set to the root of the project.

We use a number of packages, missing packages can be easily installed by remake:

```r
remake::install_missing_packages()
```

Then, to generate all figures, analyses, and tables, simply run:

```r
remake::make()
```

All output will be automatically placed in a directory called `output` (it is going to be automatically created for you).

Also notice that all the combined Bayesian models in this paper will take a several days (up to a month) to run on a regular computer.

If you find remake confusing and prefer to run plain R, you can use remake to build a script `build.R` that produces a given output, e.g.

```r
remake::make_script(filename = 'build.R')
```

### This paper was produced using the following software and associated packages:
```
R version 3.4.3 (2017-11-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X El Capitan 10.11.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] tools     stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] png_0.1-7         Hmisc_4.1-1       Formula_1.2-2     survival_2.41-3   lattice_0.20-35   HDInterval_0.1.3  LoLinR_0.0.0.9000 brms_2.1.0        ggplot2_2.2.1     Rcpp_0.12.15     
[11] rfishbase_2.1.2.4 ape_5.0           rotl_3.0.3        plyr_1.8.4       

loaded via a namespace (and not attached):
 [1] Brobdingnag_1.2-4    httr_1.3.1           tidyr_0.8.0          jsonlite_1.5         splines_3.4.3        gtools_3.5.0         StanHeaders_2.17.2   threejs_0.3.1        shiny_1.0.5         
[10] assertthat_0.2.0     stats4_3.4.3         latticeExtra_0.6-28  progress_1.1.2       backports_1.1.2      pillar_1.1.0         glue_1.2.0           digest_0.6.15        checkmate_1.8.5     
[19] RColorBrewer_1.1-2   colorspace_1.3-2     htmltools_0.3.6      httpuv_1.3.5         Matrix_1.2-12        dygraphs_1.1.1.4     XML_3.98-1.9         pkgconfig_2.0.1      rncl_0.8.2          
[28] rstan_2.17.3         purrr_0.2.4          xtable_1.8-2         mvtnorm_1.0-7        scales_0.5.0         htmlTable_1.11.2     tibble_1.4.2         bayesplot_1.4.0      DT_0.4              
[37] shinyjs_1.0          nnet_7.3-12          lazyeval_0.2.1       magrittr_1.5         mime_0.5             nlme_3.1-131         foreign_0.8-69       xts_0.10-1           colourpicker_1.0    
[46] data.table_1.10.4-3  rsconnect_0.8.5      loo_1.1.0            prettyunits_1.0.2    matrixStats_0.53.0   stringr_1.2.0        munsell_0.4.3        cluster_2.0.6        bindrcpp_0.2        
[55] compiler_3.4.3       rlang_0.1.6          grid_3.4.3           rstudioapi_0.7       htmlwidgets_1.0      crosstalk_1.0.0      igraph_1.1.2         miniUI_0.1.1         base64enc_0.1-3     
[64] gtable_0.2.0         rentrez_1.1.0        inline_0.3.14        abind_1.4-5          markdown_0.8         reshape2_1.4.3       R6_2.2.2             gridExtra_2.3        rstantools_1.4.0    
[73] zoo_1.8-1            knitr_1.19           bridgesampling_0.4-0 dplyr_0.7.4          bindr_0.1            shinystan_2.4.0      shinythemes_1.1.1    stringi_1.1.6        parallel_3.4.3      
[82] rpart_4.1-11         acepack_1.4.1        coda_0.19-1          lmtest_0.9-35       
```

### How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the green button `clone or download` and then click on `Download ZIP`  

## Bug reporting
* Please [report any issues or bugs](https://github.com/dbarneche/fishFecundity/issues).
