geth2 <- function(model) {
  vcov <- as.data.frame(VarCorr(model))[,c("grp", "vcov")]
  rownames(vcov) <- vcov$grp
  vg <- vcov["name2", "vcov"]
  if(is.na(vcov["loc:name2",  "vcov"])) {
    vge <- 0
  } else {
    vge <- vcov["loc:name2",  "vcov"]
  }
  ve <- vcov["Residual","vcov"]
  h2 <- vg/(vg + (vge/2) + (ve/(2*nrep)))
  return(h2)
}

varg <- function(model) {
  vcov <- as.data.frame(VarCorr(model))[,c("grp", "vcov")]
  rownames(vcov) <- vcov$grp
  vg <- vcov["name2", "vcov"]
  return(vg)
}

coh2 <- function(wave_1, wave_2, trait, flip) {
  wave_1 <- paste0('wave_', wave_1)
  wave_2 <- paste0('wave_', wave_2)
  # ---------------------data---------------------
  wave_data <- phenotypes[, c("name2", "loc", "set", "block", "taxa", "range", "row", trait)]

  if (flip){
    wave_data['freq_ratio'] <- phenotypes[wave_2] / phenotypes[wave_1] #* -1.0
  } else{
    wave_data['freq_ratio'] <- phenotypes[wave_1] / phenotypes[wave_2] #* -1.0
  }
  #wave_data['freq_ratio'] <- scale(wave_data['freq_ratio'])

  #for testing: coheritability with a different phenotype
  # wave_data['freq_ratio'] <- scale(phenotypes[trait])
  wave_data['sum'] <- wave_data['freq_ratio'] + wave_data[trait]


tryCatch({
    # ---------------------models---------------------
    modelW <- lmer(freq_ratio ~ loc + set + loc:set + (1 | loc:block) + (1 | name2) + (1 | loc:name2),
                   data = wave_data)
    messagew <- ifelse(is.null(modelW@optinfo$conv$lme4$messages), "...", modelW@optinfo$conv$lme4$messages)

    if(isSingular(modelW, 10e-8) | grepl('failed to converge|boundary', messagew)){
      xw <- as.data.frame(VarCorr(modelW))
      xw[xw$grp == "loc:name2","vcov"] < 10e-8
      modelW <- lmer(freq_ratio ~ loc + set + loc:set + (1 | loc:block) + (1 | name2),
                     data = wave_data)
    }

    modelS <- lmer(sum ~ loc + set + loc:set + (1 | loc:block) + (1 | name2) + (1 | loc:name2),
                   data = wave_data)
    messages <- ifelse(is.null(modelS@optinfo$conv$lme4$messages), "...", modelS@optinfo$conv$lme4$messages)

    if(isSingular(modelS, 10e-8) | grepl('failed to converge|boundary', messages)){
      xs <- as.data.frame(VarCorr(modelS))
      xs[xs$grp == "loc:name2","vcov"] < 10e-8
      modelS <- lmer(sum ~ loc + set + loc:set + (1 | loc:block) + (1 | name2),
                     data = wave_data)
    }

    if(isSingular(modelW, 10e-8) | isSingular(modelS, 10e-8)){
      #zz <- data.frame(trait = trait, ratio = paste0(wave_1,"/", wave_2), coh2= NA)
      output2 <- data.frame(
        wave_1 = wave_1,
        wave_2= wave_2,
        trait = trait,
        coh2= "NA",
        h2 = "NA", #heritability of trait
        hw2 = "NA", #heritability of wavelength ratios
        corg = "NA", #correlation between two models,
        corgblup = "NA",#correlation between blups
        covs = "NA", #genetic covariance,
        varw = "NA", #wavelength model variance
        vars = "NA", #sum model variance,
        vartrait = "NA" #trait variance
      )
      cat(wave_1, wave_2, "did not converge...\n")
    } else {
      # ---------------------h2 and var---------------------
      hw <- geth2(modelW)
      varw <- varg(modelW)
      vars <- varg(modelS)

      # ---------------------calculation genetic covariance---------------------#
      #varS = varg + varW + 2covS
      covs <- (vars - vartrait - varw) / 2

      #genetic correlation
      corg <- covs / (sqrt(vartrait) * sqrt(varw))

      #genetic correlation blups
      blups <- blups |>
        left_join(ranef(modelW)$name2 |> rownames_to_column("name2"), by = join_by(name2))

      corgblup <- cor(blups[,-1])[1,2]/(sqrt(hw)*sqrt(h2))

      # ---------------------co-heritability---------------------  #
      coh2 <- sqrt(h2) * sqrt(hw) * corg
      output2 <- data.frame(
        wave_1 = wave_1,
        wave_2= wave_2,
        trait = trait,
        coh2= coh2,
        h2 = h2, #heritability of trait
        hw2 = hw, #heritability of wavelength ratios
        corg = corg, #correlation between two models,
        corgblup = corgblup,#correlation between blups
        covs = covs, #genetic covariance
        varw = varw, #wavelength model variance
        vars = vars, #sum model variance,
        vartrait = vartrait #trait variance
      )
      #corg <- 10
      #if(corg > 1 | corg < 1){
      #  cat(wave_1, wave_2, "correlation not [-1,1]\n" )
      #  print(output2)
      #}
    }
    return(output2)
  },
  error = function(cnd) {
  message(cnd)
    output2 <- data.frame(
      wave_1 = wave_1,
      wave_2= wave_2,
      trait = trait,
      coh2= "NA",
      h2 = "NA", #heritability of trait
      hw2 = "NA", #heritability of wavelength ratios
      corg = "NA", #correlation between two models,
      corgblup = "NA",#correlation between blups
      covs = "NA", #genetic covariance
      varw = "NA", #wavelength model variance
      vars = "NA", #sum model variance,
      vartrait = "NA" #trait variance
    )
    cat(wave_1, wave_2, "did not converge...\n")
    return(output2)
  })
}
