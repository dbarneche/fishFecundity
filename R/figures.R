######################
# AUXILLIARY FUNCTIONS
######################
generatePngs  <-  function (originFileName) {
    noExt  <-  tools::file_path_sans_ext(originFileName)
    system(paste0('sips -s formatOptions best -s format png ', originFileName, ' --out ', noExt, '.png'))
}

readAndTrimPng  <-  function (pngFileName, dirPic = 'figures') {
    pngFile  <-  png::readPNG(file.path(dirPic, paste0(pngFileName, '.png')))
    trimPngBlanks(pngFile)
}

trimPngBlanks  <-  function (pngMatrix) {
    pngMatrix[!apply(pngMatrix[, , 1], 1, function (x)all(x == 1)), !apply(pngMatrix[, , 1], 2, function (x)all(x == 1)), ]
}

toPdf <- function (expr, filename, ...) {
    toDev(expr, pdf, filename, ...)
}

toDev <- function (expr, dev, filename, ..., verbose = TRUE) {
    if (verbose) {
        cat(sprintf('Creating %s\n', filename))
    }
    dev(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}

proportionalPng  <-  function (logo, px, py, log = FALSE, middleY = FALSE, ...) {
    # utility function for embedding png images at specified fractional sizes in R plots
    # places the logo centred on a specified fraction of the the usr space, 
    # and sizes appropriately (respects aspect ratio)
    # works for any plot size, including logged axes
    # logo: a png object obtained with readPNG
    # px: is a vector with 2 values (start and end) specifying the relative x positions (i.e. between 0 and 1)
    # py: is the initial (i.e. bottom) relative y position (i.e. between 0 and 1)
    # log: should the x and/or y axes be logged? Default to FALSE; Alternative values are 'x', 'y' or 'xy'    
    if (!is.numeric(px) | !is.numeric(py) | length(px) != 2) {
    	stop('wrong position coordinates [0,1]')
    }
    usr  <-  par('usr')
    pin  <-  par('pin')
    # first get proportions of coordinates right
    pxRg    <-  px[2] - px[1]
    xProp   <-  usr[1] + px * (usr[2] - usr[1]) # x range from relative px 
    yProp   <-  usr[3] + py * (usr[4] - usr[3]) # minimum y from relative py
    # now get aspect ratio to calculate maximum y
    pinRatio  <-  pin[2] / pin[1] # aspect ratio of actual plot region, depends on device and plot size
    dims      <-  dim(logo)[1:2] # number of x-y pixels for the logo (aspect ratio)
    AR        <-  dims[1] / dims[2]
    yProp     <-  c(yProp, usr[3] + (py + pxRg * AR / pinRatio) * (usr[4] - usr[3])) # maximum y from relative py correcting for x and plot ratios
    if (middleY) {
        yProp     <-  yProp - (yProp[2] - yProp[1]) / 2 # place minimum py as middle point of png    
    }
    if (log == 'x') {
        xProp  <-  10^(xProp)
    }
    if (log == 'y') { 
        yProp  <-  10^(yProp)
    }
    if (log == 'xy') {
        xProp  <-  10^(xProp)
        yProp  <-  10^(yProp)
    }
    rasterImage(logo, xProp[1], yProp[1], xProp[2], yProp[2], interpolate = TRUE, ...)
}

changePngColour  <-  function (pngObject, col, ...) {
    # change colour
    # currently works for RGBa only
    rgbVals  <-  col2rgb(col, ...) / 255
    for(i in 1:3) {
        pngObject[, , i]  <-  rgbVals[i]
    }
    pngObject
}

lengthFromWeight  <-  function (massInGrams, a = 0.00610, b = 3.114) {
    # from female Gadus morhua
    (massInGrams / a)^(1 / b)
}

hpdiBrmsSpSpecific  <-  function (brmsOutput, credMass = 0.9999, selectedSpecies, sizeName) {
    prs  <-  c('b_Intercept', paste0('b_', sizeName), paste0('r_otlSpecies[', selectedSpecies, ',', sizeName, ']'), paste0('r_otlSpecies[', selectedSpecies, ',Intercept]'), paste0('r_animal[', selectedSpecies, ',Intercept]'))
    mat  <-  brms::posterior_samples(brmsOutput, pars = prs, exact_match = TRUE)
    
    dat  <-  vector(mode = 'list', length = length(selectedSpecies))
    for (k in seq_along(selectedSpecies)) {
        sdat  <-  data.frame(mat[, 1] + mat[, grep(paste0('r_otlSpecies[', selectedSpecies[k], ',Intercept]'), colnames(mat), fixed = TRUE)] + mat[, grep(paste0('r_animal[', selectedSpecies[k], ',Intercept]'), colnames(mat), fixed = TRUE)], mat[, 2] + mat[, grep(paste0('r_otlSpecies[', selectedSpecies[k], ',', sizeName, ']'), colnames(mat), fixed = TRUE)])
        names(sdat)  <-  c(paste0(selectedSpecies[k], '_Intercept'), paste0(selectedSpecies[k], '_', sizeName))
        dat[[k]]  <-  sdat
    }
    dat  <-  do.call(cbind.data.frame, dat)
    cis  <-  data.frame(sapply(dat, HDInterval::hdi, credMass = credMass), check.names = FALSE)
    lgm  <-  matrix(TRUE, nrow(dat), ncol(dat))
    out  <-  as.list(dat)
    for (i in 1:ncol(lgm)) {
        lgm[, i]  <-  out[[i]] >= cis[[names(out)[i]]][1] & out[[i]] <= cis[[names(out)[i]]][2]
    }
    toKeep  <-  apply(lgm, 1, all)
    for (i in seq_along(out)) {
        out[[i]]  <-  out[[i]][toKeep]
    }
    out
}

totOvEngCISpp  <-  function (fecundityMassScalingModel, eggMassScalingModel, energyVolumeScalingModel, selectedSpecies, seed = 1) {
    set.seed(seed)
    fecMat   <-  hpdiBrmsSpSpecific(fecundityMassScalingModel, selectedSpecies = selectedSpecies, sizeName = 'lnMass')
    eggMat   <-  hpdiBrmsSpSpecific(eggMassScalingModel, selectedSpecies = selectedSpecies, sizeName = 'lnMass')
    engMat   <-  hpdiBrmsSpSpecific(energyVolumeScalingModel, selectedSpecies = selectedSpecies, sizeName = 'lnVolume')
    iter         <-  10000
    indexFecMat  <-  sample(1:length(fecMat[[1]]), iter, replace = TRUE)
    indexEggMat  <-  sample(1:length(eggMat[[1]]), iter, replace = TRUE)
    indexEngMat  <-  sample(1:length(engMat[[1]]), iter, replace = TRUE)

    iterDat  <-  data.frame()
    for (i in seq_len(iter)) {
        subDat  <-  vector(mode = 'list', length = length(selectedSpecies))
        for (j in seq_along(selectedSpecies)) {
            itcpt    <-  engMat[[paste0(selectedSpecies[j], '_Intercept')]][indexEngMat[i]] + eggMat[[paste0(selectedSpecies[j], '_Intercept')]][indexEggMat[i]] * engMat[[paste0(selectedSpecies[j], '_lnVolume')]][indexEngMat[i]] + fecMat[[paste0(selectedSpecies[j], '_Intercept')]][indexFecMat[i]]
            slope    <-  (engMat[[paste0(selectedSpecies[j], '_lnVolume')]][indexEngMat[i]] * eggMat[[paste0(selectedSpecies[j], '_lnMass')]][indexEggMat[i]]) + fecMat[[paste0(selectedSpecies[j], '_lnMass')]][indexFecMat[i]]
            subDat[[j]]         <-  data.frame(itcpt, slope)
            names(subDat[[j]])  <-  c(paste0(selectedSpecies[j], '_Intercept'), paste0(selectedSpecies[j], '_Slope'))
        }
        iterDat  <-  rbind(iterDat, do.call(cbind.data.frame, subDat))
    }
    iterDat
}

ovaryMassCoefsSpSpecific  <-  function (...) {
    ovaryEnergyMcmcMat       <-  totOvEngCISpp(...)
    ovaryEnergyMcmcMat[, 1]  <-  exp(ovaryEnergyMcmcMat[, 1])
    ovaryEnergy              <-  apply(ovaryEnergyMcmcMat, 2, mean)
    qtlInt                   <-  quantile(ovaryEnergyMcmcMat[, 1], probs = c(0.025, 0.975))
    qtlSlp                   <-  quantile(ovaryEnergyMcmcMat[, 2], probs = c(0.025, 0.975))
    data.frame('mean' = c(ovaryEnergy[1], ovaryEnergy[2]), 'lower95ci' = c(qtlInt[1], qtlSlp[1]), 'upper95ci' = c(qtlInt[2], qtlSlp[2]), stringsAsFactors = FALSE, row.names = c('intercept', 'slope'))
}

totOvEngCISppMix  <-  function (fecundityMassScalingModel, eggMassScalingModel, energyVolScalingModel, selectedSpecies, seed = 1) {
    set.seed(seed)
    if (any(selectedSpecies %in% fecundityMassScalingModel$data$otlSpecies)) {
        fecMat   <-  hpdiBrmsSpSpecific(fecundityMassScalingModel, selectedSpecies = selectedSpecies, sizeName = 'lnMass')
        fecTag  <-  'species-specific'
    } else {
        fecMat   <-  hpdiBrms(fecundityMassScalingModel)
        fecTag  <-  'population-level'
    }
    
    if (any(selectedSpecies %in% eggMassScalingModel$data$otlSpecies)) {
        eggMat   <-  hpdiBrmsSpSpecific(eggMassScalingModel, selectedSpecies = selectedSpecies, sizeName = 'lnMass')
        eggTag  <-  'species-specific'
    } else {
        eggMat   <-  hpdiBrms(eggMassScalingModel)
        eggTag  <-  'population-level'
    }

    if (any(selectedSpecies %in% energyVolScalingModel$data$otlSpecies)) {
        engMat   <-  hpdiBrmsSpSpecific(energyVolScalingModel, selectedSpecies = selectedSpecies, sizeName = 'lnVolume')
        engTag  <-  'species-specific'
    } else {
        engMat   <-  hpdiBrms(energyVolScalingModel)
        engTag  <-  'population-level'
    }
    
    iter         <-  10000
    indexFecMat  <-  sample(1:length(fecMat[[1]]), iter, replace = TRUE)
    indexEggMat  <-  sample(1:length(eggMat[[1]]), iter, replace = TRUE)
    indexEngMat  <-  sample(1:length(engMat[[1]]), iter, replace = TRUE)

    iterDat  <-  data.frame()
    for (i in seq_len(iter)) {
        itcpt    <-  engMat[[1]][indexEngMat[i]] + eggMat[[1]][indexEggMat[i]] * engMat[[2]][indexEngMat[i]] + fecMat[[1]][indexFecMat[i]]
        slope    <-  (engMat[[2]][indexEngMat[i]] * eggMat[[2]][indexEggMat[i]]) + fecMat[[2]][indexFecMat[i]]
        iterDat  <-  rbind(iterDat, data.frame(intercept = itcpt, slope = slope, stringsAsFactors = FALSE))
    }
    list(data = iterDat,
        tags = data.frame(fecTag = fecTag, eggTag = eggTag, engTag = engTag, stringsAsFactors = FALSE))
}

#########
# FIGURES
#########
makeFigure1  <-  function (dest, ...) {
    toPdf(fig1(...), dest, width = 11, height = 5.5)
}

fig1  <-  function (...) {
    postOvMassCoef  <-  ovaryMassCoefsSpSpecific(...)
    fishLogo        <-  png::readPNG('figures/gadus_morhua_diane_peebles.png')
    modelIntercept  <-  round(postOvMassCoef['intercept', 'mean'], 2)
    modelSlope      <-  round(postOvMassCoef['slope', 'mean'], 2)
    redFishBig      <-  (modelIntercept * 30e3 ^ modelSlope) / 1e6 # in megajoules
    redFishSmall    <-  (modelIntercept * 2e3 ^ modelSlope) / 1e6 # in megajoules
    par(omi = c(0.5, 1, 0.2, 0.2), mgp = c(3, 0.5, 0), tck = -0.03)
    layout(matrix(c(1, 1, 2, 2, 2), nrow = 1, ncol = 5))
    par(mai = c(0.7, 0.3, 0.5, 0.3), cex = 1, cex.lab = 1.3)
    plot(NA, xlab = substitute('Female mass x 10'^3 * ', ' * italic('M'['i']) * ', (g)'), ylab = substitute('Reproductive output, ' * italic('R'['i']) * ', (MJ)'), xlim = c(-2.5, 40), ylim = c(-8, 80), axes = FALSE, xpd = NA, xaxs = 'i', yaxs = 'i')
    LoLinR::proportionalLabel(0, 1.05, 'A', adj = c(0, 0.5), cex = 1, xpd = NA, font = 2)
    axis(2, las = 1, at = 0, cex.axis = 1, mgp = c(3, 0.7, 0))
    for (k in seq(20, 80, length.out = 4)) {
        axis(2, las = 1, at = k, cex.axis = 1, mgp = c(3, 0.7, 0))
    }
    axis(1, at = 0)
    for (i in seq(10, 40, 10)) {
        axis(1, at = i, cex.axis = 1)
    }
    LoLinR::proportionalLabel(0.01, 1, substitute(italic('R'['i']) == b%*%italic('M'['i'])^beta[1], list(b = substitute(a%*%'10'^-6, list(a = modelIntercept)))), adj = c(0, 1), xpd = NA)
    lines(0:40, (modelIntercept * (0:40*1e3) ^ modelSlope) / 1e6, col = 'tomato3', lwd = 3)
    points(c(2, 30), (modelIntercept * c(2e3, 30e3) ^ modelSlope) / 1e6, col = 'tomato3', pch = 16, cex = 1.3)
    LoLinR::proportionalLabel(0.4, 0.58, 'Hyper-allometric\nscaling', adj = c(0.5, 0.5), font = 3, col = 'tomato3')
    LoLinR::proportionalLabel(0.4, 0.49, substitute(beta[1] == b, list(b = modelSlope)), adj = c(0.5, 0.5), font = 3, col = 'tomato3')
    lines(0:40, (modelIntercept * (0:40*1e3) ^ 1) / 1e6, col = 'dodgerblue3', lwd = 3, lty = 2)
    points(c(2, 30), (modelIntercept * c(2e3, 30e3) ^ 1) / 1e6, col = 'dodgerblue3', pch = 15, cex = 1.3)
    LoLinR::proportionalLabel(0.5, 0.24, 'Isometric\nscaling', adj = c(0.5, 0.5), font = 3, col = 'dodgerblue3')
    LoLinR::proportionalLabel(0.5, 0.16, substitute(beta[1] == 1), adj = c(0.5, 0.5), font = 3, col = 'dodgerblue3')
    box(bty = 'l')
    
    findPosX  <-  function (target) {
        usr  <-  par('usr')
        (target - usr[1]) / (usr[2] - usr[1])
    }

    LoLinR::proportionalLabel(findPosX(2), 0.02, '(2 kg)', adj = c(0.5, 0), cex = 1)
    LoLinR::proportionalLabel(findPosX(30), 0.18, '(30 kg)', adj = c(0.5, 0.5), cex = 1)
    
    sizePlot1X      <-  par('pin')[1]
    deltaXSmlFish1  <-  0.3
    deltaXBigFish1  <-  (lengthFromWeight(30e3) * deltaXSmlFish1) / lengthFromWeight(2e3)
    proportionalPng(fishLogo, c(0.01, 0.01 + deltaXSmlFish1), c(0.1))
    proportionalPng(fishLogo, c(0.01, 0.01 + deltaXBigFish1), c(0.62), xpd = NA)

    par(mai = c(0, 0.1, 0, 0))
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = '', ylab = '', axes = FALSE, xaxs = 'i', yaxs = 'i')
    
    sizePlot2X      <-  par('pin')[1]
    deltaXBigFish2  <-  deltaXBigFish1 * sizePlot1X / sizePlot2X
    deltaXSmlFish2  <-  deltaXSmlFish1 * sizePlot1X / sizePlot2X

    # one big fish
    proportionalPng(fishLogo, c(0.05, 0.05 + deltaXBigFish2), c(0.65), xpd = NA)
    lines(c(0.48, 0.48), c(0.68, 0.9), lty = 2, col = 'grey30')
    LoLinR::proportionalLabel(0.5, 0.84, substitute(italic('n') == 1), adj = c(0, 0.5), xpd = NA, cex = 1)
    LoLinR::proportionalLabel(0.5, 0.79, substitute(italic(sum(R['i'])) %~~% b*' MJ', list(b = round(redFishBig, 1))), adj = c(0, 0.5), xpd = NA, cex = 1)
    LoLinR::proportionalLabel(0.5, 0.74, substitute(italic(sum('M'['i'])) == z*' kg', list(z = 30e3/1e3)), adj = c(0, 0.5), xpd = NA, cex = 1)
    
    # shoal of small big fish
    codShl  <-  readFile('figures/codShoalCoords.csv')
    xShoal  <-  codShl$x
    yShoal  <-  codShl$y
    for (i in seq_along(xShoal)) {
        proportionalPng(fishLogo, c(xShoal[i], xShoal[i] + deltaXSmlFish2), yShoal[i], xpd = NA)
    }
    lines(c(0.78, 0.78), c(0.15, 0.65), lty = 2, col = 'grey30')
    LoLinR::proportionalLabel(0.8, 0.5, substitute(italic('n') == a, list(a = round(redFishBig / redFishSmall))), adj = c(0, 0.5), xpd = NA, cex = 1)
    LoLinR::proportionalLabel(0.8, 0.44, substitute(italic(sum(R['i'])) %~~% b*' MJ', list(b = round(redFishSmall * (redFishBig / redFishSmall), 1))), adj = c(0, 0.5), xpd = NA, cex = 1)
    LoLinR::proportionalLabel(0.8, 0.38, substitute(italic(sum('M'['i'])) == z*' kg', list(z = (2e3 * round(redFishBig / redFishSmall)) / 1e3)), adj = c(0, 0.5), xpd = NA, cex = 1)

    LoLinR::proportionalLabel(0.05, 0.932, 'B', adj = c(0, 0.5), cex = 1, xpd = NA, font = 2)
}

makeFigure2  <-  function (dest, ...) {
    toPdf(fig2(...), dest, width = 13.5, height = 4.5)
}

fig2  <-  function (data1, fecundityMassScalingModel, data2, eggMassScalingModel, data3, energyVolScalingModel) {
    # fecundity
    data       <-  data1$data
    # extract fixed and random effects and plot
    fixefs     <-  brms::fixef(fecundityMassScalingModel, estimate = c('mean'))
    mcmcSamp   <-  brms::posterior_samples(fecundityMassScalingModel, pars = 'b_lnMass')
    cis        <-  data.frame(cis = HDInterval::hdi(mcmcSamp[, 'b_lnMass']),  credMass = 0.95)
    ranefs     <-  brms::ranef(fecundityMassScalingModel, estimate = c('mean'), var = FALSE)
    randomPhy  <-  data.frame(Intercept = ranefs$animal[, 'Estimate', 'Intercept'])
    randomSpp  <-  data.frame(Intercept = ranefs$otlSpecies[, 'Estimate', 'Intercept'], lnMass = ranefs$otlSpecies[, 'Estimate', 'lnMass'])

    data$correctedFecundity  <-  data$lnFecundity - randomPhy$Intercept[match(data$otlSpecies, rownames(randomPhy))] - randomSpp$Intercept[match(data$otlSpecies, rownames(randomSpp))] - randomSpp$lnMass[match(data$otlSpecies, rownames(randomSpp))] * data$lnMass

    par(mfrow = c(1, 3), mai = c(1.02, 1.12, 0.82, 0.3),  omi = rep(0, 4), cex = 1, cex.lab = 1.3, cex.axis = 1.2, mgp = c(3, 0.5, 0), tck = -0.03)
    plot(NA, xlab = 'Female mass (g)', ylab = 'Fecundity (# eggs/female)', xlim = range(data$lnMass), ylim = log(c(1e-1, 1e7)), axes = FALSE, xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1, at = log(c(1e-2, 1, 1e2, 9.5e3, 9e5)), labels = c(expression(paste(1, ' x ', 10^-2)), '1', expression(paste(1, ' x ', 10^2)), expression(paste(9.5, ' x ', 10^3)), expression(paste(9, ' x ', 10^5))), cex.axis = 1)
    for (k in seq(-1, 7, 2)) {
        axis(2, las = 1, at = log(10^k), labels = substitute(10^a, list(a = k)), cex.axis = 1)
    }
    points(data$lnMass, data$correctedFecundity, col = LoLinR::transparentColor('tomato3', 0.4), pch = 16, cex = 1.5)

    LoLinR::proportionalLabel(0.05, 0.9, substitute(italic(y) == b%*%italic(x)^a, list(b = LoLinR::rounded(exp(fixefs['Intercept', 'Estimate']), 2), a = LoLinR::rounded(fixefs['lnMass', 'Estimate'], 2))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.8, substitute(z :~a - b, list(z = '95% C.I.', a = LoLinR::rounded(cis['lower', 'cis'], 2), b = LoLinR::rounded(cis['upper', 'cis'], 2))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.7, substitute(italic('n') == a, list(a = niceThousands(nrow(data)))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.6, paste0(length(unique(data$otlSpecies)), ' spp.'), adj = c(0, 0.5), cex = 0.9)
    lines(range(data$lnMass), fixefs['Intercept', 'Estimate'] + range(data$lnMass) * fixefs['lnMass', 'Estimate'], lty = 2, lwd = 1.5)
    LoLinR::proportionalLabel(0, 1.05, 'A', adj = c(0, 0.5), cex = 1, xpd = NA, font = 2)

    # egg volume
    data  <-  data2$data
    # extract fixed and random effects and plot
    fixefs     <-  brms::fixef(eggMassScalingModel, estimate = c('mean'))
    mcmcSamp   <-  brms::posterior_samples(eggMassScalingModel, pars = 'b_lnMass')
    cis        <-  data.frame(cis = HDInterval::hdi(mcmcSamp[, 'b_lnMass']),  credMass = 0.95)
    ranefs     <-  brms::ranef(eggMassScalingModel, estimate = c('mean'), var = FALSE)
    randomPhy  <-  data.frame(Intercept = ranefs$animal[, 'Estimate', 'Intercept'])
    randomSpp  <-  data.frame(Intercept = ranefs$otlSpecies[, 'Estimate', 'Intercept'], lnMass = ranefs$otlSpecies[, 'Estimate', 'lnMass'])

    data$correctedlnEggVol  <-  data$lnVolume - randomPhy$Intercept[match(data$otlSpecies, rownames(randomPhy))] - randomSpp$Intercept[match(data$otlSpecies, rownames(randomSpp))] - randomSpp$lnMass[match(data$otlSpecies, rownames(randomSpp))] * data$lnMass

    plot(NA, xlab = 'Female mass (g)', ylab = substitute('Egg volume (mm'^3*')'), xlim = range(data$lnMass), ylim = range(data$correctedlnEggVol), axes = FALSE, xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1, at = log(c(1e-2, 8.8e-1, 68, 5.2e3, 4e5)), labels = c(expression(paste(1, ' x ', 10^-2)), expression(paste(8.8, ' x ', 10^-1)), '68', expression(paste(5.2, ' x ', 10^3)), expression(paste(4, ' x ', 10^5))), cex.axis = 1)    
    axis(2, at = seq(min(data$correctedlnEggVol), max(data$correctedlnEggVol), length.out = 5), labels = round(exp(seq(min(data$correctedlnEggVol), max(data$correctedlnEggVol), length.out = 5)), 2), las = 1, cex.axis = 1)

    points(data$lnMass, data$correctedlnEggVol, col = LoLinR::transparentColor('seagreen4', 0.4), pch = 17, cex = 1.5)

    LoLinR::proportionalLabel(0.05, 0.9, substitute(italic(y) == b%*%italic(x)^a, list(b = round(exp(fixefs['Intercept', 'Estimate']), 2), a = round(fixefs['lnMass', 'Estimate'], 2))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.8, substitute(z :~a - b, list(z = '95% C.I.', a = LoLinR::rounded(cis['lower', 'cis'], 2), b = LoLinR::rounded(cis['upper', 'cis'], 2))), adj = c(0, 0.5), cex = 0.9)    
    LoLinR::proportionalLabel(0.05, 0.7, substitute(italic('n') == a, list(a = nrow(data))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.6, paste0(length(unique(data$otlSpecies)), ' spp.'), adj = c(0, 0.5), cex = 0.9)
    lines(range(data$lnMass), fixefs['Intercept', 'Estimate'] + range(data$lnMass) * fixefs['lnMass', 'Estimate'], lty = 2, lwd = 1.5)
    LoLinR::proportionalLabel(0, 1.05, 'B', adj = c(0, 0.5), cex = 1, xpd = NA, font = 2)
    
    # egg energy
    data  <-  data3$data
    # extract fixed and random effects and plot
    fixefs     <-  brms::fixef(energyVolScalingModel, estimate = c('mean'))
    mcmcSamp   <-  brms::posterior_samples(energyVolScalingModel, pars = 'b_lnVolume')
    cis        <-  data.frame(cis = HDInterval::hdi(mcmcSamp[, 'b_lnVolume']),  credMass = 0.95)
    ranefs     <-  brms::ranef(energyVolScalingModel, estimate = c('mean'), var = FALSE)
    randomPhy  <-  data.frame(Intercept = ranefs$animal[, 'Estimate', 'Intercept'])
    randomSpp  <-  data.frame(Intercept = ranefs$otlSpecies[, 'Estimate', 'Intercept'], lnVolume = ranefs$otlSpecies[, 'Estimate', 'lnVolume'])

    data$correctedEnergy  <-  data$lnTotalEnergy - randomPhy$Intercept[match(data$otlSpecies, rownames(randomPhy))] - randomSpp$Intercept[match(data$otlSpecies, rownames(randomSpp))] - randomSpp$lnVolume[match(data$otlSpecies, rownames(randomSpp))] * data$lnVolume
    
    plot(NA, xlab = substitute('Egg volume (mm'^3*')'), ylab = 'Egg energy (J)', xlim = range(data$lnVolume), ylim = range(data$correctedEnergy), axes = FALSE, xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1, at = log(c(0.07, 0.27, 1.14, 4.8, 20.14)), labels = c(0.07, 0.27, 1.14, 4.8, 20.14), cex.axis = 1)
    axis(2, at = seq(min(data$correctedEnergy), max(data$correctedEnergy), length.out = 5), labels = round(exp(seq(min(data$correctedEnergy), max(data$correctedEnergy), length.out = 5)), 1), las = 1, cex.axis = 1)

    points(data$lnVolume, data$correctedEnergy, col = LoLinR::transparentColor('seagreen4', 0.4), pch = 17, cex = 1.5)
    LoLinR::proportionalLabel(0.05, 0.9, substitute(italic(y) == b%*%italic(x)^a, list(b = round(exp(fixefs['Intercept', 'Estimate']), 2), a = round(fixefs['lnVolume', 'Estimate'], 2))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.8, substitute(z :~a - b, list(z = '95% C.I.', a = LoLinR::rounded(cis['lower', 'cis'], 2), b = LoLinR::rounded(cis['upper', 'cis'], 2))), adj = c(0, 0.5), cex = 0.9)        
    LoLinR::proportionalLabel(0.05, 0.7, substitute(italic('n') == a, list(a = niceThousands(nrow(data)))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.6, paste0(length(unique(data$otlSpecies)), ' spp.'), adj = c(0, 0.5), cex = 0.9)
    lines(range(data$lnVolume), fixefs['Intercept', 'Estimate'] + range(data$lnVolume) * fixefs['lnVolume', 'Estimate'], lty = 2, lwd = 1.5)
    LoLinR::proportionalLabel(0, 1.05, 'C', adj = c(0, 0.5), cex = 1, xpd = NA, font = 2)
}

makeFigure3  <-  function (dest, ...) {
    toPdf(fig3(...), dest, width = 12, height = 7.5)
}

fig3  <-  function (sppSpecReproOutAll, fecundityMassScalingModel, eggMassScalingModel, energyVolScalingModel, otlsPngPicTab, changeColor = TRUE) {
    sppSpecReproOutAll  <-  sppSpecReproOutAll[with(sppSpecReproOutAll, order(order, family, otlSpecies)), ]
    postOvMassCoef      <-  ovaryMassCoefs(fecundityMassScalingModel, eggMassScalingModel, energyVolScalingModel)
    modelSlope          <-  round(postOvMassCoef['slope', 'mean'], 2)

    par(mfrow = c(1, 3), cex = 1, cex.lab = 1.3, cex.axis = 1)
    par(mai = c(0.9, 1.4, 0.3, 0.05))

    plotSeq        <-  c(seq(1, nrow(sppSpecReproOutAll), nrow(sppSpecReproOutAll) / 3), nrow(sppSpecReproOutAll) + 1)
    n  <-  1
    for (j in seq_along(plotSeq)[-4]) {
        xMin  <-  plotSeq[j]
        xMax  <-  plotSeq[j + 1] - 1
        plot(NA, xlab = '', ylab = '', xlim = c(min(sppSpecReproOutAll$lowerCiSlope), max(sppSpecReproOutAll$upperCiSlope)), ylim = rev(c((xMin - 1), (xMax + 1))), axes = FALSE, xpd = NA, yaxs = 'i')
        if (j == 2) {
            LoLinR::proportionalLabel(0.5, -0.1, 'Species-specific exponent', cex = 1.3, xpd  = NA, adj = c(0.5, 0.5))
        }
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
        LoLinR::whiteGrid()
        box()
        axis(1)
        lines(rep(1, 2),  par('usr')[3:4], lty = 2, lwd = 1.5)
        lines(rep(modelSlope, 2),  par('usr')[3:4], lty = 3, lwd = 1, col = 'grey30')
        rowNumbers  <-  rep(NA, nrow(sppSpecReproOutAll) / 3)
        rowNumbers[seq(1, nrow(sppSpecReproOutAll) / 3, 5)] <-  seq(xMin, xMax, 5)

        subSeq  <-  xMin:xMax
        for (k in seq_along(subSeq)) {
            lines(sppSpecReproOutAll[subSeq[k], c('lowerCiSlope', 'upperCiSlope')], rep(subSeq[k], 2), lwd = 1.2, col = 'grey30')
            if (sppSpecReproOutAll[subSeq[k], 'avSlope'] > 1) {
                ptCol  <-  'tomato3'
                ptShp  <-  16
            } else if (sppSpecReproOutAll[subSeq[k], 'avSlope'] < 1) {
                ptCol  <-  'seagreen4'
                ptShp  <-  17
            } else {
                ptCol  <-  'dodgerblue3'
                ptShp  <-  15
            }
            points(sppSpecReproOutAll[subSeq[k], 'avSlope'], subSeq[k], col = ptCol, pch = ptShp, cex = 0.6)
            axis(2, las = 1, at = subSeq[k], labels = rowNumbers[k], xpd = NA, cex.axis = 0.8, mgp = c(3, 0.7, 0))
            if (sppSpecReproOutAll$otlSpecies[subSeq[k]] == otlsPngPicTab$otlSpecies[n]) {
                yPos      <-  1 - abs(subSeq[k]  - par('usr')[4]) / abs(par('usr')[3]  - par('usr')[4])
                pngObjct  <-  readAndTrimPng(otlsPngPicTab[n, 2])
                if (changeColor) {
                    pngObjct  <-  changePngColour(pngObjct, 'grey30')
                }
                proportionalPng(pngObjct, c(-0.2 - otlsPngPicTab[n, 3], -0.2), c(yPos), xpd = NA, middleY = TRUE)
                if (n != nrow(otlsPngPicTab)) {
                    n  <-  n + 1    
                }
            }
            
        }
    }
}

makeFigureS1  <-  function (dest, ...) {
    toPdf(figS1(...), dest, width = 7, height = 7)
}

figS1  <-  function (engVolPhyData, energyVolScalingModel) {
    # egg energy
    data  <-  engVolPhyData$data
    data  <-  data[data$ReferencePdf == 'robertson_dr_2018_unpublished.pdf', ]
    data$lnDryWeight  <-  log(data$eggDryWeight_mg)
    par(omi = c(0.5, 0.8, 0.5, 0.5), cex = 1, cex.lab = 1.3, cex.axis = 1.2)        
    plot(NA, xlab = substitute('Egg dry weight (mg)'), ylab = 'Egg energy (J)', xlim = range(data$lnDryWeight), ylim = range(data$lnTotalEnergy), axes = FALSE, xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1, at = seq(min(data$lnDryWeight), max(data$lnDryWeight), length.out = 5), labels = round(exp(seq(min(data$lnDryWeight), max(data$lnDryWeight), length.out = 5)), 2), cex.axis = 1)
    axis(2, at = seq(min(data$lnTotalEnergy), max(data$lnTotalEnergy), length.out = 5), labels = round(exp(seq(min(data$lnTotalEnergy), max(data$lnTotalEnergy), length.out = 5)), 1), las = 1, cex.axis = 1)
    points(data$lnDryWeight, data$lnTotalEnergy, col = LoLinR::transparentColor('dodgerblue3', 0.4), pch = 15, cex = 1.5)
    mod  <-  coef(summary(lm(lnTotalEnergy ~ lnDryWeight, data = data)))
    LoLinR::proportionalLabel(0.05, 0.9, substitute(italic(y) == b%*%italic(x)^a, list(b = round(exp(mod['(Intercept)', 'Estimate']), 2), a = round(mod['lnDryWeight', 'Estimate'], 2))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.8, substitute(italic('n') == a, list(a = niceThousands(nrow(data)))), adj = c(0, 0.5), cex = 0.9)
    LoLinR::proportionalLabel(0.05, 0.7, paste0(length(unique(data$otlSpecies)), ' spp.'), adj = c(0, 0.5), cex = 0.9)
    lines(range(data$lnDryWeight), mod['(Intercept)', 'Estimate'] + range(data$lnDryWeight) * mod['lnDryWeight', 'Estimate'], lty = 2, lwd = 1.5)
    mean(data$lnTotalEnergy / (data$lnDryWeight / 1e3))
}
