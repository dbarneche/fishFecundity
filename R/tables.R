makeTableS3  <-  function (dest, data, fecundityMassScalingModel, refTab, spawningTable) {
    data  <-  data$data

    # order and family-arrange phylogeny
	refTab       <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
	refTab$f3    <-  seq_len(nrow(refTab))
	refTab       <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
    refTab       <-  refTab[refTab$otlSpecies %in% data$otlSpecies, ]
    rownames(refTab)  <-  NULL
	
	# extract full posterior distributions for slopes 
    fixefs                <-  brms::fixef(fecundityMassScalingModel, estimate = c('mean'))
	posteriorSlopeFixed   <-  brms::posterior_samples(fecundityMassScalingModel, pars = 'b_lnMass')
	posteriorSlopeRandom  <-  brms::posterior_samples(fecundityMassScalingModel, pars = ',lnMass]')
	names(posteriorSlopeRandom)  <-  gsub(',lnMass]', '', gsub('r_otlSpecies[', '', names(posteriorSlopeRandom), fixed = TRUE))
    posteriorSlopeRandom  <-  posteriorSlopeRandom[, match(refTab$otlSpecies, names(posteriorSlopeRandom))]

	tableS3          <-  data.frame()
	for (k in 1:ncol(posteriorSlopeRandom)) {
		su  <-  refTab[refTab$otlSpecies == names(posteriorSlopeRandom)[k], ]
		xs  <-  data$lnMass[data$otlSpecies == names(posteriorSlopeRandom)[k]]
		if (length(xs) > 1) {
			obs  <-  sum(data$otlSpecies == names(posteriorSlopeRandom)[k])
			fcR  <-  paste0(niceThousands(range(data$Fecundity_nOfEggs_per_female[data$otlSpecies == names(posteriorSlopeRandom)[k]])), collapse = ' – ')
			msR  <-  paste0(LoLinR::rounded(range(data$mass_g[data$otlSpecies == names(posteriorSlopeRandom)[k]]), 2), collapse = ' – ')
		} else {
			obs  <-  1
			fcR  <-  niceThousands(data$Fecundity_nOfEggs_per_female[data$otlSpecies == names(posteriorSlopeRandom)[k]])
			msR  <-  LoLinR::rounded(data$mass_g[data$otlSpecies == names(posteriorSlopeRandom)[k]], 2)
		}
		posterior  <-  posteriorSlopeFixed + posteriorSlopeRandom[, names(posteriorSlopeRandom)[k]]
		cis95      <-  quantile(posterior[, 1], probs = c(0.025, 0.975))
    	tableS3    <-  rbind(tableS3, data.frame('otlSpecies' = names(posteriorSlopeRandom)[k], 'F3' = su$f3, '$\\beta_1$' = LoLinR::rounded(mean(posterior[, 1]), 2), '2.5% C.I.' = LoLinR::rounded(cis95[1], 2), '97.5% C.I.' = LoLinR::rounded(cis95[2], 2), 'n' = obs, 'Fecundity range' = fcR, 'Mass range (g)' = msR,  'ReferencesPdfs' = paste0(unique(data$ReferencePdf[data$otlSpecies == names(posteriorSlopeRandom)[k]]), collapse = ' ; '), 'Family' = su$family, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL))
	}

	splitSpp              <-  lapply(tableS3$otlSpecies, function(x)strsplit(x, '_')[[1]])
	ottPos                <-  sapply(splitSpp, function (x)grep('ott', x, fixed = TRUE))
	tableS3[['Species']]  <-  sapply(splitSpp, function(x)paste(x[1:2], collapse = ' '))
	tableS3[['OTL']]      <-  mapply(function(x, y)x[y], x = splitSpp, y = ottPos)
	tableS3[['SM']]       <-  spawningTable$spawningModeCode[match(tableS3$Species, spawningTable$species)]
	tableS3               <-  tableS3[with(tableS3, order(Family, Species)), ]
	write.csv(tableS3[, c(10:13, 2:8)], dest, row.names = FALSE, na = '-')
}

makeTableS4  <-  function (dest, data, eggMassScalingModel, refTab, spawningTable) {
    data  <-  data$data

    # order and family-arrange phylogeny
	refTab       <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
	refTab$f3    <-  seq_len(nrow(refTab))
	refTab       <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
    refTab       <-  refTab[refTab$otlSpecies %in% data$otlSpecies, ]
    rownames(refTab)  <-  NULL

	# extract full posterior distributions for slopes 
    fixefs                <-  brms::fixef(eggMassScalingModel, estimate = c('mean'))
	posteriorSlopeFixed   <-  brms::posterior_samples(eggMassScalingModel, pars = 'b_lnMass')
	posteriorSlopeRandom  <-  brms::posterior_samples(eggMassScalingModel, pars = ',lnMass]')
	names(posteriorSlopeRandom)  <-  gsub(',lnMass]', '', gsub('r_otlSpecies[', '', names(posteriorSlopeRandom), fixed = TRUE))

	tableS4  <-  data.frame()
	for (k in 1:ncol(posteriorSlopeRandom)) {
		su  <-  refTab[refTab$otlSpecies == names(posteriorSlopeRandom)[k], ]
		xs  <-  data$lnMass[data$otlSpecies == names(posteriorSlopeRandom)[k]]
		if (length(xs) > 1) {
			obs      <-  sum(data$otlSpecies == names(posteriorSlopeRandom)[k])
			eggVolR  <-  paste0(LoLinR::rounded(range(exp(data$lnVolume[data$otlSpecies == names(posteriorSlopeRandom)[k]])), 3), collapse = ' – ')
			msR      <-  paste0(LoLinR::rounded(range(data$mass_g[data$otlSpecies == names(posteriorSlopeRandom)[k]]), 2), collapse = ' – ')
		} else {
			obs      <-  1
			eggVolR  <-  LoLinR::rounded(exp(data$lnVolume[data$otlSpecies == names(posteriorSlopeRandom)[k]]), 3)
			msR      <-  LoLinR::rounded(data$mass_g[data$otlSpecies == names(posteriorSlopeRandom)[k]], 2)
		}
		posterior  <-  posteriorSlopeFixed + posteriorSlopeRandom[, names(posteriorSlopeRandom)[k]]
		cis95      <-  quantile(posterior[, 1], probs = c(0.025, 0.975))
    	tableS4    <-  rbind(tableS4, data.frame('otlSpecies' = names(posteriorSlopeRandom)[k], 'F3' = su$f3, '$\\beta_1$' = LoLinR::rounded(mean(posterior[, 1]), 2), '2.5% C.I.' = LoLinR::rounded(cis95[1], 2), '97.5% C.I.' = LoLinR::rounded(cis95[2], 2), 'n' = obs, 'Egg-volume range (mm\\textsuperscript{3})' = eggVolR, 'Mass range (g)' = msR,  'ReferencesPdfs' = paste0(unique(data$ReferencePdf[data$otlSpecies == names(posteriorSlopeRandom)[k]]), collapse = ' ; '), 'Family' = su$family, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL))
	}

	splitSpp              <-  lapply(tableS4$otlSpecies, function(x)strsplit(x, '_')[[1]])
	ottPos                <-  sapply(splitSpp, function (x)grep('ott', x, fixed = TRUE))
	tableS4[['Species']]  <-  sapply(splitSpp, function(x)paste(x[1:2], collapse = ' '))
	tableS4[['OTL']]      <-  mapply(function(x, y)x[y], x = splitSpp, y = ottPos)
	tableS4[['SM']]       <-  spawningTable$spawningModeCode[match(tableS4$Species, spawningTable$species)]
	tableS4               <-  tableS4[with(tableS4, order(Family, Species)), ]
	write.csv(tableS4[, c(10:13, 2:8)], dest, row.names = FALSE, na = '-')
}

makeTableS5  <-  function (dest, data, energyVolScalingModel, refTab, spawningTable) {
    data  <-  data$data

    # order and family-arrange phylogeny
	refTab       <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
	refTab$f3    <-  seq_len(nrow(refTab))
	refTab       <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
    refTab       <-  refTab[refTab$otlSpecies %in% data$otlSpecies, ]
    rownames(refTab)  <-  NULL

	# extract full posterior distributions for slopes 
    fixefs                <-  brms::fixef(energyVolScalingModel, estimate = c('mean'))
	posteriorSlopeFixed   <-  brms::posterior_samples(energyVolScalingModel, pars = 'b_lnVolume')
	posteriorSlopeRandom  <-  brms::posterior_samples(energyVolScalingModel, pars = ',lnVolume]')
	names(posteriorSlopeRandom)  <-  gsub(',lnVolume]', '', gsub('r_otlSpecies[', '', names(posteriorSlopeRandom), fixed = TRUE))

	tableS5          <-  data.frame()
	for (k in 1:ncol(posteriorSlopeRandom)) {
		su  <-  refTab[refTab$otlSpecies == names(posteriorSlopeRandom)[k], ]
		xs  <-  data$lnVolume[data$otlSpecies == names(posteriorSlopeRandom)[k]]
		if (length(xs) > 1) {
			obs      <-  sum(data$otlSpecies == names(posteriorSlopeRandom)[k])
			enR      <-  paste0(LoLinR::rounded(range(exp(data$lnTotalEnergy[data$otlSpecies == names(posteriorSlopeRandom)[k]])), 2), collapse = ' – ')
			eggVolR  <-  paste0(LoLinR::rounded(range(exp(data$lnVolume[data$otlSpecies == names(posteriorSlopeRandom)[k]])), 3), collapse = ' – ')
			
		} else {
			obs      <-  1
			enR      <-  LoLinR::rounded(exp(data$lnTotalEnergy[data$otlSpecies == names(posteriorSlopeRandom)[k]]), 2)
			eggVolR  <-  LoLinR::rounded(exp(data$lnVolume[data$otlSpecies == names(posteriorSlopeRandom)[k]]), 3)
		}
		posterior  <-  posteriorSlopeFixed + posteriorSlopeRandom[, names(posteriorSlopeRandom)[k]]
		cis95      <-  quantile(posterior[, 1], probs = c(0.025, 0.975))
    	tableS5    <-  rbind(tableS5, data.frame('otlSpecies' = names(posteriorSlopeRandom)[k], 'F3' = su$f3, '$\\beta_1$' = LoLinR::rounded(mean(posterior[, 1]), 2), '2.5% C.I.' = LoLinR::rounded(cis95[1], 2), '97.5% C.I.' = LoLinR::rounded(cis95[2], 2), 'n' = obs, 'Egg-energy range (J)' = enR,  'Egg-volume range (mm\\textsuperscript{3})' = eggVolR, 'ReferencesPdfs' = paste0(unique(data$ReferencePdf[data$otlSpecies == names(posteriorSlopeRandom)[k]]), collapse = ' ; '), 'Family' = su$family, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL))
	}

	splitSpp              <-  lapply(tableS5$otlSpecies, function(x)strsplit(x, '_')[[1]])
	ottPos                <-  sapply(splitSpp, function (x)grep('ott', x, fixed = TRUE))
	tableS5[['Species']]  <-  sapply(splitSpp, function(x)paste(x[1:2], collapse = ' '))
	tableS5[['OTL']]      <-  mapply(function(x, y)x[y], x = splitSpp, y = ottPos)
	tableS5[['SM']]       <-  spawningTable$spawningModeCode[match(tableS5$Species, spawningTable$species)]
	tableS5               <-  tableS5[with(tableS5, order(Family, Species)), ]
	write.csv(tableS5[, c(10:13, 2:8)], dest, row.names = FALSE, na = '-')
}

makeTableS6  <-  function (dest, fecFemSizeVolPhyData, totalVolumeMassScalingNoSingleObsModel, refTab, spawningTable) {
    data  <-  removeSingleSpecies(fecFemSizeVolPhyData)$data

    # order and family-arrange phylogeny
	refTab       <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
	refTab$f3    <-  seq_len(nrow(refTab))
	refTab       <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
    refTab       <-  refTab[refTab$otlSpecies %in% data$otlSpecies, ]
    rownames(refTab)  <-  NULL
	
	# extract full posterior distributions for slopes 
    fixefs                <-  brms::fixef(totalVolumeMassScalingNoSingleObsModel, estimate = c('mean'))
	posteriorSlopeFixed   <-  brms::posterior_samples(totalVolumeMassScalingNoSingleObsModel, pars = 'b_lnMass')
	posteriorSlopeRandom  <-  brms::posterior_samples(totalVolumeMassScalingNoSingleObsModel, pars = ',lnMass]')
	names(posteriorSlopeRandom)  <-  gsub(',lnMass]', '', gsub('r_otlSpecies[', '', names(posteriorSlopeRandom), fixed = TRUE))
    posteriorSlopeRandom  <-  posteriorSlopeRandom[, match(refTab$otlSpecies, names(posteriorSlopeRandom))]

	tableS6  <-  data.frame()
	for (k in 1:ncol(posteriorSlopeRandom)) {
		su  <-  refTab[refTab$otlSpecies == names(posteriorSlopeRandom)[k], ]
		xs  <-  data$lnMass[data$otlSpecies == names(posteriorSlopeRandom)[k]]
		if (length(xs) > 1) {
			obs  <-  sum(data$otlSpecies == names(posteriorSlopeRandom)[k])
			fcR  <-  paste0(niceThousands(range(exp(data$lnTotalVolume)[data$otlSpecies == names(posteriorSlopeRandom)[k]])), collapse = ' – ')
			msR  <-  paste0(LoLinR::rounded(range(data$mass_g[data$otlSpecies == names(posteriorSlopeRandom)[k]]), 2), collapse = ' – ')
		} else {
			obs  <-  1
			fcR  <-  niceThousands(exp(data$lnTotalVolume)[data$otlSpecies == names(posteriorSlopeRandom)[k]])
			msR  <-  LoLinR::rounded(data$mass_g[data$otlSpecies == names(posteriorSlopeRandom)[k]], 2)
		}
		posterior  <-  posteriorSlopeFixed + posteriorSlopeRandom[, names(posteriorSlopeRandom)[k]]
		cis95      <-  quantile(posterior[, 1], probs = c(0.025, 0.975))
    	tableS6    <-  rbind(tableS6, data.frame('otlSpecies' = names(posteriorSlopeRandom)[k], 'F3' = su$f3, '$\\beta_1$' = LoLinR::rounded(mean(posterior[, 1]), 2), '2.5% C.I.' = LoLinR::rounded(cis95[1], 2), '97.5% C.I.' = LoLinR::rounded(cis95[2], 2), 'n' = obs, 'Total egg-volume range (mm\\textsuperscript{3})' = fcR, 'Mass range (g)' = msR,  'ReferencesPdfs' = paste0(unique(data$ReferencePdf[data$otlSpecies == names(posteriorSlopeRandom)[k]]), collapse = ' ; '), 'Family' = su$family, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL))
	}

	splitSpp              <-  lapply(tableS6$otlSpecies, function(x)strsplit(x, '_')[[1]])
	ottPos                <-  sapply(splitSpp, function (x)grep('ott', x, fixed = TRUE))
	tableS6[['Species']]  <-  sapply(splitSpp, function(x)paste(x[1:2], collapse = ' '))
	tableS6[['OTL']]      <-  mapply(function(x, y)x[y], x = splitSpp, y = ottPos)
	tableS6[['SM']]       <-  spawningTable$spawningModeCode[match(tableS6$Species, spawningTable$species)]
	tableS6               <-  tableS6[with(tableS6, order(Family, Species)), ]
	write.csv(tableS6[, c(10:13, 2:8)], dest, row.names = FALSE, na = '-')
}

makeTableS7  <-  function (dest, refTab, spawningTable) {
	data     <-  refTab[with(refTab, order(order, family, otlSpecies)), ]
	data$f3  <-  seq_len(nrow(data))
	data     <-  data[with(data, order(order, family, otlSpecies)), ]
	tableS7  <-  data.frame('Family' = data$family, 'otlSpecies' = data$otlSpecies, 'F3' = data$f3, '$\\beta_1$' = LoLinR::rounded(data$avSlope, 2), '2.5% C.I.' = LoLinR::rounded(data$lowerCiSlope, 2), '97.5% C.I.' = LoLinR::rounded(data$upperCiSlope, 2), 'D1' = ifelse(data$fecundityPars == 'population-level', 'No', 'Yes'), 'D2' = ifelse(data$eggVolPars == 'population-level', 'No', 'Yes'), 'D3' = ifelse(data$eggEngPars == 'population-level', 'No', 'Yes'), stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
	splitSpp              <-  lapply(tableS7$otlSpecies, function(x)strsplit(x, '_')[[1]])
	ottPos                <-  sapply(splitSpp, function (x)grep('ott', x, fixed = TRUE))
	tableS7[['Species']]  <-  sapply(splitSpp, function(x)paste(x[1:2], collapse = ' '))
	tableS7[['OTL']]      <-  mapply(function(x, y)x[y], x = splitSpp, y = ottPos)
	tableS7[['SM']]       <-  spawningTable$spawningModeCode[match(tableS7$Species, spawningTable$species)]
	allDatasets           <-  apply(tableS7[, paste0('D', 1:3)] == 'Yes', 1, all)
	tableS7[['Species']][allDatasets]  <-  paste0(tableS7[['Species']][allDatasets], '\\textbf{*}')
	tableS7               <-  tableS7[with(tableS7, order(Family, Species)), ]
	write.csv(tableS7[, c(1, 10:12, 3:9)], dest, row.names = FALSE, na = '-')
}

makeTableS8  <-  function (dest, buckleyTrioData, kokitaTrioData, bradfordTrioData) {    
    # buckleyTrioData
    buckleyTrioData$eggDryWeight_g     <-  buckleyTrioData$eggDryWeight_ug / 1e6
    buckleyTrioData$clucthDryWeight_g  <-  buckleyTrioData$eggDryWeight_g * buckleyTrioData$fecundity
    bucMod  <-  coef(summary(lm(log(clucthDryWeight_g) ~ log(wetWeight_g), data = buckleyTrioData)))

    # kokitaTrioData
	kokitaTrioData$mass_g  <-  0.037 * ((kokitaTrioData$lengthSl_mm / 0.79) / 10) ^ 2.63 # from FishBase
    kokMod  <-  coef(summary(lm(log(clutchDryWeight_g) ~ log(mass_g), data = kokitaTrioData)))

    # bradfordTrioData
	head(bradfordTrioData)
	bradfordTrioData$mass_g          <-  0.00420 * (bradfordTrioData$length_tl_mm / 10) ^ 3.178 # from FishBase
	bradfordTrioData$clutchWeight_g  <-  bradfordTrioData$clutchWeight_mg / 1e3
    brdMod  <-  coef(summary(lm(log(clutchWeight_g) ~ log(mass_g), data = bradfordTrioData)))

    tableS8  <-  data.frame('Species'        = c('Pseudopleuronectes americanus', 'Clupea harengus', 'Pomacentrus coelestis'),
    						'$\\beta_0$'       =  LoLinR::rounded(exp(c(bucMod[1, 1], kokMod[1, 1], brdMod[1, 1])), 2),
    						'$\\beta_1$'       =  LoLinR::rounded(c(bucMod[2, 1], kokMod[2, 1], brdMod[2, 1]), 2),
    						'$\\beta_1$ S.E.'  =  LoLinR::rounded(c(bucMod[2, 2], kokMod[2, 2], brdMod[2, 2]), 2),
    						'n'              =  c(nrow(buckleyTrioData), nrow(kokitaTrioData), nrow(bradfordTrioData)),
    						'ref'            =  c('buckley_lj_1991_m74_125.pdf', 'bradford_rg_1992_c49_2045.pdf', 'kokita_t_2003_m143_593.pdf'),
    						stringsAsFactors = FALSE, check.names = FALSE)
	write.csv(tableS8, dest, row.names = FALSE)
}
