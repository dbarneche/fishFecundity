###################################################################
# THIS IS THE CODE USED TO DOWNLOAD THE FISHBASE DATASET
# IT CAN BE USED TO UPDATE THE DATA IN CASE THERE IS A NEW VERSION
# OF THE DATASETS AVAILABLE ON FISHBASE, PARTICULARLY BECAUSE SOME
# LINKS MAY BECOME DEPRECATED OVER TIME. NOTICE THOUGH THAT BY
# DOWNLOADING UPDATED VERSIONS OF THESE DATASETS FROM FISHBASE,
# THE FINAL ANALYSES MAY LIKELY YIELD DIFFERENT RESULTS BECAUSE
# THE INPUT DATA WILL HAVE CHANGED. ALSO IT IS POSSIBLE THAT THE
# FISHBASE SYSTEM WILL GO DOWN FOR MAINTENANCE WHILE THE DOWNLOAD
# IS HAPPENING, SO KEEP AN EYE OUT BECAUSE YOU MAY HAVE TO RUN THIS
# A FEW TIMES TO GET IT TO WORK. FASTER INTERNET CONNECTIONS HELP.
# PLEASE CONTACT THE AUTHOR VIA THE ISSUES PAGE ON 
# https://github.com/dbarneche/fishFecundity/issues
# SHOULD YOU HAVE PROBLEMS USING THIS CODE.
###################################################################
giveTaxonTable  <-  function (pathToHtmlTable) {
    pagetree  <-  XML::htmlTreeParse(pathToHtmlTable, error = function (...){})
    urltable  <-  XML::xmlToList(pagetree$children$html$children$body$children$table$children[[2]])
    data.frame(species  =  unname(apply(urltable, 2, function (x)unlist(x[1])[1])),
               urls     =  unname(apply(urltable, 2, function (x)unlist(x[1])[2])),
               family   =  unname(apply(urltable, 2, function (x)unlist(x[3])[1])),
               stringsAsFactors = FALSE)
}

consistentTree  <-  function (...) {
    sptable  <-  readAndParseTree(...)
    while (inherits(sptable, 'try-error')) {
        sptable  <-  readAndParseTree(...)
    }
    sptable
}

readAndParseTree  <-  function (url, ...) {
    n  <-  1
    message('Attempt #', n, ' to open: ', url, '\n\n', sep = '')
    sptree  <-  openReadCloseConnection(url, ...)
    while (inherits(sptree, 'try-error')) {
        n  <-  n + 1
        message('Attempt #', n, ' to open: ', url, '\n\n', sep = '')
        sptree  <-  openReadCloseConnection(url, ...)
    }
    try(xml2::read_html(sptree, verbose = TRUE), silent = TRUE)
}

openReadCloseConnection  <-  function (url, verbose) {
    try(RCurl::getURL(url = url, followLocation = TRUE, .opts = list(timeout = 20, maxredirs = 2, verbose = verbose)), silent = TRUE)
}

downloadFishBaseLwData  <-  function (fishBaseHtmlFile, verbose = TRUE) {
    taxonInfo  <-  giveTaxonTable(fishBaseHtmlFile)
    lwData     <-  data.frame()
    for (i in 1:nrow(taxonInfo)) {
        lwData  <-  rbind(lwData, getLwData(taxonInfo[i, ], verbose = verbose))
    }
    cleanLwData(lwData)
}

getLwData  <-  function (data, verbose = TRUE) {
    if (verbose) {
        message('Extracting length-weight data for:', data$species, '\n')
    }
    sptable  <-  consistentTree(data$urls, verbose = verbose)
    htmlTab  <-  sptable %>% rvest::html_nodes('table')
    while (length(htmlTab) == 0) {
        message('Empty page: re-attempting to download length-weight data for:', data$species, '\n')
        sptable  <-  consistentTree(data$urls, verbose = verbose)
        htmlTab  <-  sptable %>% rvest::html_nodes('table')
    }
    tabPos        <-  grep('Score', htmlTab)
    filteredTab   <-  htmlTab %>% magrittr::extract2(tabPos)
    scoresName    <-  filteredTab %>% rvest::html_nodes('input') %>% rvest::html_attr(name = 'name')
    scoresVec     <-  filteredTab %>% rvest::html_nodes('input') %>% rvest::html_attr(name = 'value')
    scoresVec     <-  scoresVec[scoresName == 'score[]']
    gtab          <-  filteredTab %>% rvest::html_table(header = TRUE)
    names(gtab)   <-  c('score', 'a', 'b', 'doubtful', 'sex', 'length_cm', 'length_type', 'r2', 'sd_b', 'sd_log10_a', 'n', 'country', 'locality')
    gtab$score    <-  as.numeric(scoresVec)
    gtab$species  <-  data$species
    gtab$family   <-  data$family
    gtab
}

cleanLwData  <-  function (data) {
    for (j in c('score', 'a', 'b', 'sd_b', 'sd_log10_a', 'n')) {
        data[[j]]  <-  as.numeric(data[[j]])
    }
    re             <-  '([[:alpha:].]+)[[:space:]]'
    data$doubtful  <-  tolower(gsub(re, '\\1', data$doubtful))
    data$doubtful  <-  gsub('[[:space:]]', NA, data$doubtful)
    data$sex       <-  tolower(data$sex)
    
    data$length_cm      <-  fixASCII(data$length_cm)
    rangeList           <-  sapply(data$length_cm, strsplit, split = '[-][[:space:].]')
    data$length_cm_min  <-  sapply(rangeList, evaluateAndReturnRange, pos = 1)
    data$length_cm_max  <-  sapply(rangeList, evaluateAndReturnRange, pos = 2)
    
    re                <-  '([[:alpha:].]+)[[:space:]]'
    data$length_type  <-  toupper(gsub(re, '\\1', data$length_type))
    data$length_type[data$length_type == 'NA']  <-  NA
    data$length_type  <-  gsub('[[:space:]]', NA, data$length_type)
    
    data$r2  <-  as.numeric(gsub('&nbsp', '', data$r2))
    
    data
}

fixASCII  <-  function (charVector) {
    Encoding(charVector)  <-  'latin1'
    iconv(charVector, 'latin1', 'ASCII', sub = '')
}

evaluateAndReturnRange  <-  function (x, pos) {
    lenX  <-  length(x) == 2
    if (lenX) {
        x  <-  gsub(',', '', x)
        as.numeric(x[pos])
    } else {
        NA
    }
}

prepareFishBaseUrl  <-  function (species) {
    paste0('http://www.fishbase.de/summary/', sub(' ', '-', species), '.html')
}

getBayesianLwData  <-  function (data, verbose = TRUE) {
    sppUrl   <-  prepareFishBaseUrl(data$species)
    sptable  <-  consistentTree(sppUrl, verbose = verbose)
    htmlTab  <-  sptable %>% rvest::html_nodes('div')
    tabPos   <-  grep('Bayesian length-weight', htmlTab)
    if (length(tabPos) == 0) {
        message(data$species, ' does not contain Bayesian lw-data, consider adding it manually\n')
        data[, c('a', 'a_lower', 'a_upper','b', 'b_lower', 'b_upper', 'obs')]  <-  NA
    } else {
        lwInfo   <-  rvest::html_text(htmlTab[tabPos[length(tabPos)]], trim = TRUE)
        obs      <-  strsplit(lwInfo, '), ', fixed = TRUE)[[1]]
        numbers  <-  as.numeric(unlist(regmatches(obs[1:2], gregexpr("[[:digit:]]+\\.*[[:digit:]]*", obs[1:2]))))
        data[, c('a', 'a_lower', 'a_upper','b', 'b_lower', 'b_upper')]  <-  numbers
    }
    data
}

library(plyr)
library(XML)
library(xml2)
library(RCurl)
library(rvest)
library(dplyr)
library(magrittr)

# R version 3.4.3 (2017-11-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6

# Matrix products: default
# BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# locale:
# [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] magrittr_1.5    dplyr_0.7.4     rvest_0.3.2     RCurl_1.95-4.10 bitops_1.0-6    xml2_1.2.0      XML_3.98-1.9   
# [8] plyr_1.8.4     

# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.15     assertthat_0.2.0 R6_2.2.2         pillar_1.1.0     httr_1.3.1       rlang_0.1.6     
#  [7] bindrcpp_0.2     glue_1.2.0       compiler_3.4.3   pkgconfig_2.0.1  bindr_0.1        tibble_1.4.2    

lwData  <-  downloadFishBaseLwData('data/fishBaseData/lw_parameters.html')
write.csv(lwData, 'data/fishBaseData/lwData.csv', row.names = FALSE)

source('R/analyses.R')
lwMissingSpecies  <-  readFile('data/fishBaseData/lwMissingSpecies.csv')
lwMissingData     <-  plyr::ddply(lwMissingSpecies, .(species), getBayesianLwData)
write.csv(lwMissingData, 'data/fishBaseData/lwMissingData', row.names = FALSE, quote = FALSE)
