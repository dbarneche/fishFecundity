####################
# GENERAL STAN SPECS
####################
options(mc.cores = parallel::detectCores())

###################
# GENERAL FUNCTIONS
###################
readFile  <-  function (filePath, ...) {
    read.csv(filePath, header = TRUE, stringsAsFactors = FALSE, ...)
}

niceThousands  <-  function (number) {
    formatC(number, format = 'fg', big.mark = ',')
}

numbers2words  <-  function (x) {
    # https://github.com/ateucher/useful_code/blob/master/R/numbers2words.r
    # https://gist.githubusercontent.com/psychemedia/150cb9901529da58124a/raw/a12bfd600af11065255f39787c965b7a4be086d0/numbers2words.R
    # Function by John Fox found here: 
    # http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
    # Tweaks by AJH to add commas and "and"
    helper  <-  function (x) {
        digits   <-  rev(strsplit(as.character(x), '')[[1]])
        nDigits  <-  length(digits)
        if (nDigits == 1) {
            as.vector(ones[digits])
        } else if (nDigits == 2) {
            if (x <= 19) {
                as.vector(teens[digits[1]])
            } else {
                trim(paste(tens[digits[2]], Recall(as.numeric(digits[1]))))
            }
        } else if (nDigits == 3) {
            trim(paste(ones[digits[3]], 'hundred and', Recall(makeNumber(digits[2:1]))))
        } else {
            nSuffix  <-  ((nDigits + 2) %/% 3) - 1
            if (nSuffix > length(suffixes)) {
                stop(paste(x, 'is too large!'))
            }
            trim(paste(Recall(makeNumber(digits[nDigits:(3*nSuffix + 1)])), suffixes[nSuffix], ',', Recall(makeNumber(digits[(3*nSuffix):1]))))
        }
    }
    trim  <-  function (text) {
        # Tidy leading/trailing whitespace, space before comma
        text  <-  gsub('^\ ', '', gsub('\ *$', '', gsub('\ ,', ',', text)))
        # Clear any trailing ' and'
        text  <-  gsub(' and$','',text)
        # Clear any trailing comma
        gsub('\ *,$', '', text)
    }  
    makeNumber  <-  function(...) {
        as.numeric(paste(..., collapse = ''))
    }
    # Disable scientific notation
    opts  <-  options(scipen = 100) 
    on.exit(options(opts))
    ones  <-  c('', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine') 
    names(ones)  <-  0:9 
    teens  <-  c('ten', 'eleven', 'twelve', 'thirteen', 'fourteen', 'fifteen', 'sixteen', ' seventeen', 'eighteen', 'nineteen')
    names(teens)  <-  0:9 
    tens <- c('twenty', 'thirty', 'forty', 'fifty', 'sixty', 'seventy', 'eighty', 'ninety')
    names(tens)  <-  2:9
    x  <-  round(x)
    suffixes  <-  c('thousand', 'million', 'billion', 'trillion')
    if (length(x) > 1) {
        return(trim(sapply(x, helper)))
    }
    helper(x)
}

####################
# FECUNDITY DATABASE
# MANIPULATION
####################
readAndAddMass  <-  function (filePath, addLnVolume = FALSE, ...) {
    data  <-  readFile(filePath)
    if (addLnVolume) {
        data$lnVolume  <-  log((4 / 3) * pi * (exp(data$lnEggSize) / 2)^3)
    }
    addMassToData(data, ...)
}

addMassToData  <-  function (data, lwData, lwMissing, verboseMissing = TRUE) {
    ####################################
    # KEEP a & b VALUES WITH THE 
    # HIGHEST SCORE FOR EACH SPECIES
    ####################################
    subL   <-  lwData[, c('species', 'family', 'score', 'a', 'b', 'length_cm_min', 'length_cm_max')]
    subL   <-  plyr::ddply(subL, .(species), function (x)x[which.max(x$score), ])

    ##############################
    # MERGE WITH FECUNDITY DATA
    ##############################
    mgdDat     <-  merge(x = data, y = subL, by.x = 'Species', by.y = 'species', all.x = TRUE)

    if (verboseMissing) {
        message('The following species did not have original length-weight conversion data on FishBase:', paste0('\n', unique(mgdDat$Species[is.na(mgdDat$a)]), ';'), '\n\nGrabbing Bayesian closest LW...\n\n')
    }

    mgdDat[is.na(mgdDat$a), c('a', 'b')]  <-  lwMissing[match(mgdDat$Species[is.na(mgdDat$a)], lwMissing$species), c('a', 'b')]
    mgdDat$mass_g  <-  mgdDat$a * (mgdDat$FemaleSize_mm / 10) ^ mgdDat$b # size needs to be in cm
    mgdDat$mass_g[is.na(mgdDat$mass_g)]  <-  mgdDat$FemaleMass_g[is.na(mgdDat$mass_g)]
    mgdDat$lnMass  <-  log(mgdDat$mass_g)
    mgdDat
}

addTotalEggVolume  <-  function (data) {
    data$lnTotalVolume  <-  data$lnVolume + data$lnFecundity
    data
}

##############################
# EGG ENERGY DATA MANIPULATION
##############################
cleanEggDryWeightData  <-  function (path) {
    data                  <-  readFile(path)
    data$volume_mm3       <-  (4 / 3) * pi * (data$eggDiameter_mm / 2)^3
    data$energy_j_per_mg  <-  25 # from Kamler 2005, marine species: 24.9, 25.3 and 25.2 joules per mg so it's very stable -- reasonable to assume 25 J per g
    data$country          <-  NA
    data                  <-  data[, c('Species', 'volume_mm3', 'eggDryWeight_mg', 'energy_j_per_mg', 'country', 'ReferencePdf')]
    data$eggEnergy_j      <-  data$eggDryWeight_mg * data$energy_j_per_mg
    rownames(data)        <-  NULL
    data$lnVolume         <-  log(data$volume_mm3)
    data$lnTotalEnergy    <-  log(data$eggEnergy_j)
    data
}

cleanEggEnergyData  <-  function (path) {
    data                <-  readFile(path)
    data                <-  data[, c('Species', 'volume_mm3', 'eggDryWeight_mg', 'energy_j_per_mg', 'country')]
    data$ReferencePdf   <-  'robertson_dr_2018_unpublished.pdf'
    # outliers
    data$eggEnergy_j    <-  data$eggDryWeight_mg * data$energy_j_per_mg
    data$lnVolume       <-  log(data$volume_mm3)
    data$lnTotalEnergy  <-  log(data$eggEnergy_j)
    data
}

cleanAndBindEnergyAndWeightData  <-  function (pathEnergyData, pathWeightData) {
    energy     <-  cleanEggEnergyData(pathEnergyData)
    dryWeight  <-  cleanEggDryWeightData(pathWeightData)
    rbind(energy, dryWeight)
}

###############
# IDEAL STUDIES
###############
crunchBradfordData  <-  function () {
    fecundity  <-  readFile('data/trioData/bradford_rg_1992_c49_2045_fig4_spring.csv')
    dryWeight  <-  readFile('data/trioData/bradford_rg_1992_c49_2045_fig5a_spring.csv')
    
    fecundity$fecundity        <-  round(fecundity$fecundity_divided_1e3 * 1e3)
    dryWeight$eggDryWeight_mg  <-  dryWeight$eggDryWeight_mg_per_100eggs / 1e2

    # match data to closest size
    fecundity$matchedSize      <-  NA
    fecundity$eggDryWeight_mg  <-  NA
    for (j in 1:nrow(fecundity)) {
        matchedLine                   <-  which.min(abs(dryWeight$length_tl_mm - fecundity$length_tl_mm[j]))
        fecundity$matchedSize[j]      <-  dryWeight$length_tl_mm[matchedLine]
        fecundity$eggDryWeight_mg[j]  <-  dryWeight$eggDryWeight_mg[matchedLine]
    }
    fecundity$clutchWeight_mg  <-  fecundity$fecundity * fecundity$eggDryWeight_mg
    fecundity
}

####################
# PHYLOGENETIC TREES
####################
readTree  <-  function (treeTopology) {
    ape::read.tree(text = treeTopology)
}

rockFishTree   <-  function () {
    message('Acanthoclinus fuscus species reinserted right next to Moronidae\n')
    readTree('(Acanthoclinus_fuscus_ott3633802);')
}

snookTree   <-  function () {
    message('Centropomus undecimalis and Lates calcarifer species reinserted right next to Carangidae\n')
    readTree('(Centropomus_undecimalis_ott317368,Lates_calcarifer_ott6362446);')
}

cardinalTree   <-  function () {
    message('Apogon doederleini species reinserted into Apogonidae\n')
    readTree('(Apogon_doederleini_ott687634);')
}

damselTree  <-  function () {
    damselTopology  <-  '((((Microspathodon_bairdii_ott237630,(Microspathodon_chrysurus_ott847660,Microspathodon_dorsalis_ott205758)),Hypsypops_rubicundus_ott847666),((Stegastes_planifrons_ott665837,Stegastes_acapulcoensis_ott3635554),(Stegastes_partitus_ott345269,(Stegastes_flavilatus_ott3635541,((Stegastes_variabilis_ott323173,Stegastes_leucostictus_ott100830),(Stegastes_adustus_ott323181,(Stegastes_fuscus_ott3635543,Stegastes_diencaeus_ott729161))))))),((Chromis_atripectoralis_ott741423,(Chromis_multilineata_ott437016,Chromis_atrilobata_ott436999)),(((Abudefduf_septemfasciatus_ott129790,Abudefduf_sordidus_ott1053071),(Abudefduf_bengalensis_ott318931,(Abudefduf_vaigiensis_ott1053067,(Abudefduf_saxatilis_ott405751,Abudefduf_troschelii_ott961357)))),((Acanthochromis_polyacanthus_ott100410,Amphiprion_melanopus_ott45635),(Pomacentrus_amboinensis_ott880239,Pomacentrus_coelestis_ott622060)))));'
    message('Pomacentridae species reinserted right next to Labridae using topology published in Frédérich et al. 2013 Am Nat 181, 94--113\n')
    readTree(damselTopology)
}

########################################
# TREE EXTRACTION AND MATCHING WITH DATA
########################################
downloadAndProcessTree  <-  function () {
    actn  <-  rotl::tnrs_match_names(names = 'Actinopterygii', context_name = 'Animals')
    tree  <-  rotl::tol_subtree(ott_id = actn$ott_id)
    attr(tree, 'order')  <-  'cladewise'
    # add Dipnoi species as an outer group, necessary to root the tree
    out   <-  rotl::tol_subtree(ott_id = rotl::tnrs_match_names(names = 'Lepidosiren paradoxa', context_name = 'Animals')$ott_id)
    tree2 <-  readTree(paste0('(', out, ');'))
    tp    <-  ape::getMRCA(tree, tree$tip.label)
    tree  <-  ape::bind.tree(tree, tree2, where = tp)
    tree  <-  ape::root(tree, outgroup = out, resolve.root = TRUE)
    tree  <-  ape::compute.brlen(tree, method = 'Grafen')
    tree  <-  graftMissingTree(tree, treeToBeInsertedFun = list('damselTree', 'snookTree', 'rockFishTree', 'cardinalTree'), anchorFamily = list('Labridae', 'Carangidae', 'Moronidae', 'Apogonidae'))
    invisible(tree)
}

graftMissingTree  <-  function (tree, treeToBeInsertedFun, anchorFamily) {
    spp  <-  sapply(tree$tip.label, function (x)paste(strsplit(x, '_')[[1]][1:2], collapse = ' '))
    for (i in seq_along(treeToBeInsertedFun)) {
        newTree       <-  get(treeToBeInsertedFun[[i]])()
        fbSpp         <-  rfishbase::species_list(Family = anchorFamily[[i]])
        anchorPoint   <-  ape::getMRCA(tree, tree$tip.label[na.omit(match(fbSpp, spp))])
        tree          <-  ape::bind.tree(tree, newTree, where = anchorPoint)
    }
    # compute branch lengths on way out
    ape::compute.brlen(tree, method = 'Grafen')
}

cleanAndExtractPhyloTreeAndMatchedData  <-  function (data, tree, otlSubs) {
    all   <-  createSpeciesOtl(data, tree, otlSubs)
    tree  <-  replaceDuplicatedNodes(all$tree)
    tree  <-  treeCleanUpForAnalysis(all$data, tree)
    list('data' = all$data,
         'tree' = tree)
}

createSpeciesOtl  <-  function (data, tree, otlSubs) {
    # create a column in datafile to match OTL name
    # there are species from data that are found as subspecies in the
    # Actinopterygii tree from the Open Tree of Life. Some others are
    # not found in the tree: they are either misnomers or do not exist
    # at all, though have closely-related species from the same genera
    # so change all those to facilitate things. These are found in the
    # data.frame otlSubs
    spp        <-  unique(data$Species)
    ottFrames  <-  paste0(gsub(' ', '_', spp), '_ott')
    out        <-  data.frame()
    for (i in seq_along(ottFrames)) {
        x  <-  grep(ottFrames[i], tree$tip.label, fixed = TRUE)
        if (length(x) == 0) {
            x  <-  grep(otlSubs$ottFrame[otlSubs$Species == spp[i]], tree$tip.label, fixed = TRUE)
            tree$tip.label[x]  <-  modifyTreeTip(tree$tip.label[x], spp[i])
        }
        out  <-  rbind(out, data.frame(Species = spp[i], otlSpecies = tree$tip.label[x], stringsAsFactors = FALSE))
    }
    if (nrow(out) != nrow(out[complete.cases(out), ])) {
        stop('missing species in tree tips')
    }
    data$otlSpecies  <-  out$otlSpecies[match(data$Species, out$Species)]
    data$animal      <-  as.factor(data$otlSpecies)
    list('data' = data,
         'tree' = tree)
}

modifyTreeTip  <-  function (treeTip, species) {
    brokenName  <-  strsplit(treeTip, '_')[[1]]
    paste0(gsub(' ', '_', species), '_', grep('ott', brokenName, value = TRUE))
}

treeCleanUpForAnalysis  <-  function (data, tree) {
    # find species that are in the data, but not in the tree, and remove them from data
    missingFromData  <-  setdiff(tree$tip.label, unique(data$otlSpecies))
    ape::drop.tip(tree, missingFromData)
}

replaceDuplicatedNodes  <-  function (tree) {
    if (any(duplicated(tree$node.label))) {
        tree$node.label  <-  paste0('ott', seq_along(tree$node.label))
    }
    tree
}

####################
# STATISTICAL MODELS
####################
removeSingleSpecies  <-  function (data) {
    tree    <-  data$tree
    data    <-  data$data
    toKeep  <-  table(data$otlSpecies)
    data    <-  data[data$otlSpecies %in% names(toKeep)[toKeep > 1], ]
    tree    <-  ape::drop.tip(tree, setdiff(tree$tip.label, data$otlSpecies))
    list('data' = data,
         'tree' = tree)
}

removeSingleSpeciesAndFitModel  <-  function (modelFunction, ...) {
    data  <-  removeSingleSpecies(...)
    get(modelFunction)(data)
}

fitFecundityMassScalingModel  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnFecundity ~ lnMass + (1 | animal) + (1 + lnMass | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(1, 2), 'b'), prior(normal(3, 3), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 15e3, warmup = 7500, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitFecundityMassScalingModelWithSpawningInteraction  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnFecundity ~ lnMass + spawningMode * lnMass + (1 | animal) + (1 + lnMass | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(1, 2), 'b'), prior(normal(3, 3), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 15e3, warmup = 7500, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitFecundityMassScalingModelWithPhylogenySlope  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnFecundity ~ lnMass + (1 + lnMass | animal) + (1 + lnMass | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(1, 2), 'b'), prior(normal(3, 3), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 15e3, warmup = 7500)
}

fecundityModelComparison  <-  function (complexModel, simpleModel) {
    looCompTabIc  <-  brms::LOO(complexModel, simpleModel, pointwise = FALSE, cores = 4)
    looCompTab    <-  loo::compare(looCompTabIc[[1]], looCompTabIc[[2]])
    pVal          <-  2 * pnorm(-abs((looCompTab['elpd_diff'] - 0) / looCompTab['se']))
    list(looCompTabIc  =  looCompTabIc,
         looCompTab    =  looCompTab,
         pVal          =  pVal
        )
}

fitFecundityMassScalingModelNoFieldLab  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    data  <-  data[data$typeOfWork_Fecundity != 'Field/Lab', ]
    tree  <-  ape::drop.tip(tree, setdiff(tree$tip.label, data$otlSpecies))
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnFecundity ~ lnMass + (1 | animal) + (1 + lnMass | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(1, 2), 'b'), prior(normal(3, 3), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 15e3, warmup = 7500, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitTotalVolumeMassScalingModel  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnTotalVolume ~ lnMass + (1 | animal) + (1 + lnMass | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(1, 2), 'b'), prior(normal(3, 3), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 15e3, warmup = 7500, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitEggMassScalingModel  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnVolume ~ lnMass + (1 | animal) + (1 + lnMass | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitEggMassScalingModelWithSpawningInteraction  <-  function (data, spawningTable) {
    tree  <-  data$tree
    data  <-  data$data
    A     <-  ape::vcv(tree, corr = FALSE)
    data$matchSppName  <-  sapply(data$otlSpecies, function (x) {paste(strsplit(x, '_')[[1]][1:2], collapse = ' ')})
    data$spawningMode  <-  spawningTable$spawningModeDetailed[match(data$matchSppName, spawningTable$species)]
    brms::brm(lnVolume ~ lnMass + lnMass * spawningMode + (1 | animal) + (1 + lnMass | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitEggMassScalingModelNoFieldLab  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    data  <-  data[data$typeOfWork_eggSize != 'Field/Lab', ]
    tree  <-  ape::drop.tip(tree, setdiff(tree$tip.label, data$otlSpecies))
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnVolume ~ lnMass + (1 | animal) + (1 + lnMass | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitEnergyVolumeScalingModel  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnTotalEnergy ~ lnVolume + (1 | animal) + (1 + lnVolume | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitEnergyVolumeScalingModelWithSpawningInteraction  <-  function (data, spawningTable) {
    tree  <-  data$tree
    data  <-  data$data
    A     <-  ape::vcv(tree, corr = FALSE)
    data$matchSppName  <-  sapply(data$otlSpecies, function (x) {paste(strsplit(x, '_')[[1]][1:2], collapse = ' ')})
    data$spawningMode  <-  spawningTable$spawningModeDetailed[match(data$matchSppName, spawningTable$species)]
    brms::brm(lnTotalEnergy ~ lnVolume + lnVolume * spawningMode + (1 | animal) + (1 + lnVolume | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

fitEnergyVolumeScalingModelNoCorsica  <-  function (data) {
    tree  <-  data$tree
    data  <-  data$data
    data  <-  data[data$country != 'corsica', ]
    tree  <-  ape::drop.tip(tree, setdiff(tree$tip.label, data$otlSpecies))
    A     <-  ape::vcv(tree, corr = FALSE)
    brms::brm(lnTotalEnergy ~ lnVolume + (1 | animal) + (1 + lnVolume | otlSpecies), data = data, family = gaussian(), cov_ranef = list(animal = A), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 4, cores = 4, iter = 6000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 20))
}

####################
# MS DATA EXTRACTION
####################
heritabilityPosterior  <-  function (model, simplified = TRUE) {
    post  <-  brms::posterior_samples(model, pars = c('sd_animal__Intercept', 'sd_otlSpecies__Intercept', 'sigma'))
    x     <-  post$sd_animal__Intercept^2 / (post$sd_animal__Intercept^2 + post$sd_otlSpecies__Intercept^2 + post$sigma^2) * 100
    if (simplified) {    
        c('mean' = mean(x), quantile(x, probs = c(0.025, 0.5, 0.975)))
    } else {
        x
    }
}

sppSpecificMassFecundityScaling  <-  function (model) {
    # extract full posterior distributions for slopes
    posteriorSlopeFixed   <-  brms::posterior_samples(model, pars = 'b_lnMass')
    posteriorSlopeRandom  <-  brms::posterior_samples(model, pars = ',lnMass]')
    names(posteriorSlopeRandom)  <-  gsub(',lnMass]', '', gsub('r_otlSpecies[', '', names(posteriorSlopeRandom), fixed = TRUE))

    # get random slope coefficients for caterpillar plot
    postSlopeCoef  <-  data.frame()
    for (k in 1:ncol(posteriorSlopeRandom)) {
        posterior      <-  posteriorSlopeFixed + posteriorSlopeRandom[, names(posteriorSlopeRandom)[k]]
        cis95          <-  quantile(posterior[, 1], probs = c(0.025, 0.975))
        postSlopeCoef  <-  rbind(postSlopeCoef, data.frame(species = names(posteriorSlopeRandom)[k], mean = mean(posterior[, 1]), lower_95ci = cis95[1], upper_95ci = cis95[2], stringsAsFactors = FALSE, row.names = NULL))
    }
    postSlopeCoef
}

ovaryMassCoefs  <-  function (...) {
    ovaryEnergyMcmcMat       <-  totalOvaryEnergyWithCIs(...)
    ovaryEnergyMcmcMat[, 'intercept']  <-  exp(ovaryEnergyMcmcMat[, 'intercept'])
    ovaryEnergy              <-  apply(ovaryEnergyMcmcMat, 2, mean)
    qtlInt                   <-  quantile(ovaryEnergyMcmcMat[, 'intercept'], probs = c(0.025, 0.975))
    qtlSlp                   <-  quantile(ovaryEnergyMcmcMat[, 'slope'], probs = c(0.025, 0.975))
    data.frame('mean' = c(ovaryEnergy['intercept'], ovaryEnergy['slope']), 'lower95ci' = c(qtlInt[1], qtlSlp[1]), 'upper95ci' = c(qtlInt[2], qtlSlp[2]), stringsAsFactors = FALSE, row.names = c('intercept', 'slope'))
}

modelFixedStats  <-  function (model) {
    brms::fixef(model, estimate = 'mean')[, c('2.5%ile', 'Estimate', '97.5%ile')]
}

hpdiBrms  <-  function (brmsOutput, credMass = 0.9999) {
    mat  <-  brms::posterior_samples(brmsOutput, pars = 'b_')
    cis  <-  data.frame(sapply(mat, HDInterval::hdi, credMass = credMass))
    lgm  <-  matrix(TRUE, nrow(mat), ncol(mat))
    out  <-  as.list(mat)
    for (i in 1:ncol(lgm)) {
        lgm[, i]  <-  out[[i]] >= cis[[names(out)[i]]][1] & out[[i]] <= cis[[names(out)[i]]][2]
    }
    toKeep  <-  apply(lgm, 1, all)
    for (i in seq_along(out)) {
        out[[i]]  <-  out[[i]][toKeep]
    }
    out
}

totalOvaryEnergyWithCIs  <-  function (fecundityMassScalingModel, eggMassScalingModel, energyVolumeScalingModel, seed = 1) {
    set.seed(seed)
    fecMat   <-  hpdiBrms(fecundityMassScalingModel)
    eggMat   <-  hpdiBrms(eggMassScalingModel)
    engMat   <-  hpdiBrms(energyVolumeScalingModel)
    iter         <-  10000
    indexFecMat  <-  sample(1:length(fecMat[[1]]), iter, replace = TRUE)
    indexEggMat  <-  sample(1:length(eggMat[[1]]), iter, replace = TRUE)
    indexEngMat  <-  sample(1:length(engMat[[1]]), iter, replace = TRUE)

    iterDat  <-  data.frame()
    for (i in seq_len(iter)) {
        itcpt    <-  engMat[['b_Intercept']][indexEngMat[i]] + eggMat[['b_Intercept']][indexEggMat[i]] * engMat[['b_lnVolume']][indexEngMat[i]] + fecMat[['b_Intercept']][indexFecMat[i]]
        slope    <-  (engMat[['b_lnVolume']][indexEngMat[i]] * eggMat[['b_lnMass']][indexEggMat[i]]) + fecMat[['b_lnMass']][indexFecMat[i]]
        iterDat  <-  rbind(iterDat, data.frame(intercept = itcpt, slope = slope, stringsAsFactors = FALSE))
    }
    iterDat
}

createTaxonomyTab  <-  function (fecundityMassScalingModel, eggMassScalingModel, energyVolScalingModel, spawningTable) {
    allOtl  <-  sort(unique(c(fecundityMassScalingModel$data$otlSpecies, eggMassScalingModel$data$otlSpecies, energyVolScalingModel$data$otlSpecies)))

    # order and family-arrange phylogeny
    spp     <-  unname(sapply(allOtl, function (x) {paste(strsplit(x, '_')[[1]][1:2], collapse = ' ')}))
    orders  <-  spawningTable$order[match(spp, spawningTable$species)]
    family  <-  spawningTable$family[match(spp, spawningTable$species)]
    refTab  <-  data.frame(order = orders, family = family, otlSpecies = allOtl, stringsAsFactors = FALSE)
    refTab  <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
    rownames(refTab)  <-  NULL
    refTab
}

getSppSpecReproOutAll  <-  function (fecundityMassScalingModel, eggMassScalingModel, energyVolScalingModel, taxonomyTab) {
    taxonomyTab$avInt  <-  taxonomyTab$lowerCiInt  <-  taxonomyTab$upperCiInt  <-  taxonomyTab$avSlope  <-  taxonomyTab$lowerCiSlope  <-  taxonomyTab$upperCiSlope  <-  taxonomyTab$fecundityPars  <-  taxonomyTab$eggVolPars  <-  taxonomyTab$eggEngPars  <-  NA
    for (i in 1:nrow(taxonomyTab)) {
        x    <-  totOvEngCISppMix(fecundityMassScalingModel, eggMassScalingModel, energyVolScalingModel, taxonomyTab$otlSpecies[i])
        cis  <-  sapply(x$data, quantile, probs = c(0.025, 0.975))
        x$data$intercept  <-  exp(x$data$intercept)
        avs  <-  sapply(x$data, mean)
        taxonomyTab$avInt[i]          <-  avs['intercept']
        taxonomyTab$lowerCiInt[i]     <-  cis['2.5%', 'intercept']
        taxonomyTab$upperCiInt[i]     <-  cis['97.5%', 'intercept']
        taxonomyTab$avSlope[i]        <-  avs['slope']
        taxonomyTab$lowerCiSlope[i]   <-  cis['2.5%', 'slope']
        taxonomyTab$upperCiSlope[i]   <-  cis['97.5%', 'slope']
        taxonomyTab$fecundityPars[i]  <-  x$tags$fecTag
        taxonomyTab$eggVolPars[i]     <-  x$tags$eggTag
        taxonomyTab$eggEngPars[i]     <-  x$tags$engTag
    }
    taxonomyTab
}

