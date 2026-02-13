#---- Data Processing ----
library(tidyverse)
library(data.table)
library(here)
library(asreml)
library(asremlPlus)
library(ASRgenomics)
library(fs)
library(cvTools)
library(janitor)
library(AGHmatrix)
source("./function/aux_function.R")

asreml.options(workspace = "800mb",
               pworkspace = "300mb",
               maxit = 300)
# load kinship
kin <- fread("./data/kinship_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)
# load phenotypes
wave <- fread("./data/phenotypes_whole.csv", data.table = FALSE) %>% clean_names()
# common IDs
common_ids <- intersect(rownames(kin), unique(wave$taxa))
# filter both
kin  <- kin[common_ids, common_ids, drop = FALSE]
wave <- wave[wave$taxa %in% common_ids, ]
# save filtered phenotypes (cleanly)
write.csv(wave, "./data/phenotypes_whole_filtered.csv", row.names = FALSE)
# tune + invert filtered kinship
Gb   <- G.tuneup(G = kin, bend = TRUE, eig.tol = 1e-06)$Gb
GINV <- G.inverse(G = Gb, sparseform = TRUE)
# saveRDS(kin,  "./data/kin_common_ids.rds")
# saveRDS(Gb,   "./data/Gb_common_ids.rds")
# saveRDS(GINV, "./data/GINV_common_ids.rds")
# write.table(common_ids, "./data/common_ids.txt", row.names=FALSE, col.names="taxa", quote=FALSE)

Nwave_pheno <- read.csv("./output/Nwave_pheno.csv")
# ---------------------processing of data---------------------
Nwave_pheno<-
  Nwave_pheno|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for each wavelength for location MW
variables <- colnames(Nwave_pheno)[11:13]
models <- vector("list",length(variables))
Nratio_bluesMW <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~corgh(loc):ar1(range):ar1(row),
      data = Nwave_pheno,subset= loc== "MW", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      Nratio_bluesMW <- bind_rows(Nratio_bluesMW, temp)
      cat('\n')

    }  else {
      Nratio_bluesMW <- bind_rows(Nratio_bluesMW, data.frame(name2 = NA,
                                                             predicted.value = NA,
                                                             wave = variables[i]))
      cat('\n')
    }
  },

  error = function(err) {
    message("An error occured")
    print(err)
  })
}

Nratio_bluesMW<- Nratio_bluesMW %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesMW, "./output/Nratio_bluesMW.csv", row.names = FALSE)
plot(model)
#calculating blues for each wavelength
variables <- colnames(Nwave_pheno)[11:13]
models <- vector("list",length(variables))
Nratio_bluesEF <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = Nwave_pheno,subset= loc== "EF", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      Nratio_bluesEF <- bind_rows(Nratio_bluesEF, temp)
      cat('\n')

    }  else {
      Nratio_bluesEF <- bind_rows(Nratio_bluesEF, data.frame(name2 = NA,
                                                             predicted.value = NA,
                                                             wave = variables[i]))
      cat('\n')
    }
  },

  error = function(err) {
    message("An error occured")
    print(err)
  })
}
Nratio_bluesEF<- Nratio_bluesEF %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesEF, "./output/Nratio_bluesEF.csv", row.names = FALSE)
Nratio_bluesEF<- read.csv("./output/Nratio_bluesEF.csv")

# ---------------------NareaEF blues fitting model---------------------
nareabluesEF <- asreml(
  fixed = narea ~ name2 +set,
  random =  ~  block,
  data = Nwave_pheno, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
nareabluesEF<- update.asreml(nareabluesEF)
# ---------------------storing prediction---------------------
nareabluesEF <- data.frame(
  name2= nareabluesEF$predictions$pvals$name2,
  narea = round(nareabluesEF$predictions$pvals$predicted.value, 3))
NbluesEF<-left_join(Nratio_bluesEF, nareabluesEF, by= "name2")
fwrite( NbluesEF, "./output/NbluesEF.csv", row.names = FALSE)
# ---------------------NareaMW blues fitting model---------------------
nareabluesMW <- asreml(
  fixed = narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = Nwave_pheno, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
nareabluesMW<- update.asreml(nareabluesMW)
# ---------------------storing prediction---------------------
nareabluesMW <- data.frame(
  name2= nareabluesMW$predictions$pvals$name2,
  narea = round(nareabluesMW$predictions$pvals$predicted.value, 3))
NbluesMW<-left_join(Nratio_bluesMW, nareabluesMW, by= "name2")
fwrite( NbluesMW, "./output/NbluesMW.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
NbluesEF<- read.csv("./output/NbluesEF.csv")
NbluesMW<- read.csv("./output/NbluesMW.csv")

Nblues <- bind_rows(
  "EF" = NbluesEF,
  "MW" = NbluesMW,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Nblues <-Nblues|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
N_blues<- subset(Nblues, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
N_blues<- droplevels(N_blues[N_blues$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(N_blues,"./output/N_blues.csv", row.names = F)

#narea heritability for Ef
nareaEF <- asreml(
  fixed = narea ~ set,
  random = ~name2 + block,
  data = Nwave_pheno, subset= loc== "EF",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
nareaEF<- update.asreml(nareaEF)

# ---------------------calculating and storing heritability---------------------
h2<- (1 - ((nareaEF$predictions$avsed["mean"] ^ 2) /(2 * summary(nareaEF)$varcomp["name2", "component"]))) # narea h2= 0.342

#h2 of narea in MW loc
nareaMW <- asreml(
  fixed = narea ~ set,
  random = ~name2 + block,
  data = Nwave_pheno, subset= loc== "MW",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
nareaMW<- update.asreml(nareaMW)
# ---------------------calculating and storing heritability---------------------
h2<- (1 - ((nareaMW$predictions$avsed["mean"] ^ 2) /(2 * summary(nareaMW)$varcomp["name2", "component"]))) # narea h2=0.42 for MW
#SLA blues
Swave_pheno <- read.csv("./output/Swave_pheno.csv")
# ---------------------processing of data---------------------
Swave_pheno<-
  Swave_pheno|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for wave ratio for loc MW
variables <- colnames(Swave_pheno)[11:13]
models <- vector("list",length(variables))
Sratio_bluesMW <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = Swave_pheno,subset= loc== "MW", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      Sratio_bluesMW <- bind_rows(Sratio_bluesMW, temp)
      cat('\n')

    }  else {
      Sratio_bluesMW <- bind_rows(Sratio_bluesMW, data.frame(name2 = NA,
                                                             predicted.value = NA,
                                                             wave = variables[i]))
      cat('\n')
    }
  },

  error = function(err) {
    message("An error occured")
    print(err)
  })
}
Sratio_bluesMW %>%
  group_by(name2, wave) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1L) %>%
  print()
Sratio_bluesMW<- Sratio_bluesMW %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesMW, "./output/Sratio_bluesMW.csv", row.names = FALSE)

#calculating blues for each waveratio for loc EF for sla
variables <- colnames(Swave_pheno)[11:13]
models <- vector("list",length(variables))
Sratio_bluesEF <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = Swave_pheno,subset= loc== "EF", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      Sratio_bluesEF <- bind_rows(Sratio_bluesEF, temp)
      cat('\n')

    }  else {
      Sratio_bluesEF <- bind_rows(Sratio_bluesEF, data.frame(name2 = NA,
                                                             predicted.value = NA,
                                                             wave = variables[i]))
      cat('\n')
    }
  },

  error = function(err) {
    message("An error occured")
    print(err)
  })
}

Sratio_bluesEF<- Sratio_bluesEF %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesEF, "./output/Sratio_bluesEF.csv", row.names = FALSE)


# ---------------------SareaEF blues fitting model---------------------
slabluesEF <- asreml(
  fixed = sla ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = Swave_pheno, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
slabluesEF<- update.asreml(slabluesEF)
# ---------------------storing prediction---------------------
slabluesEF <- data.frame(
  name2= slabluesEF$predictions$pvals$name2,
  sla = round(slabluesEF$predictions$pvals$predicted.value, 3))
SbluesEF <-left_join(Sratio_bluesEF, slabluesEF, by= "name2")
fwrite( SbluesEF, "./output/SbluesEF.csv", row.names = FALSE)
# ---------------------SlaMW blues fitting model---------------------
slabluesMW <- asreml(
  fixed = sla ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = Swave_pheno, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
slabluesMW<- update.asreml(slabluesMW)
# ---------------------storing prediction---------------------
slabluesMW <- data.frame(
  name2= slabluesMW$predictions$pvals$name2,
  sla = round(slabluesMW$predictions$pvals$predicted.value, 3))
SbluesMW<-left_join(Sratio_bluesMW, slabluesMW, by= "name2")
fwrite( SbluesMW, "./output/SbluesMW.csv", row.names = FALSE)
SbluesEF<- read.csv("./output/SbluesEF.csv")
SbluesMW<- read.csv("./output/SbluesMW.csv")

Sblues <- bind_rows(
  "EF" = SbluesEF,
  "MW" = SbluesMW,
  .id = "env")

# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Sblues <-Sblues|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
S_blues<- subset(Sblues, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
S_blues<- droplevels(S_blues[S_blues$taxa %in%rownames(kin), ]) #filtering indv that are present both in kin and phenotypic data
write.csv(S_blues, "./output/S_blues.csv", row.names = F)

#sla heritability
slaEF <- asreml(
  fixed = sla ~ set,
  random = ~name2 + block,
  data = Swave_pheno, subset= loc== "EF",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
slaEF<- update.asreml(slaEF)
# ---------------------calculating and storing heritability---------------------
h2<- (1 - ((slaEF$predictions$avsed["mean"] ^ 2) /(2 * summary(slaEF)$varcomp["name2", "component"]))) # sla h2= 0.332 for EF location

#sla h2 for MW
#sla heritability
slaMW <- asreml(
  fixed = sla ~ set,
  random = ~name2 + block,
  data = Swave_pheno, subset= loc== "MW",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
slaMW<- update.asreml(slaMW)
# ---------------------calculating and storing heritability---------------------
h2<- (1 - ((slaMW$predictions$avsed["mean"] ^ 2) /(2 * summary(slaMW)$varcomp["name2", "component"]))) # sla h2= 0.312 for MW location

#plsr-narea(pn)
pnwave_pheno <- read.csv("./output/pnwave_pheno.csv")
# ---------------------processing of data---------------------
pnwave_pheno<-
  pnwave_pheno|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for wave ratio for loc MW
variables <- colnames(pnwave_pheno)[11:13]
models <- vector("list",length(variables))
pnratio_bluesMW <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = pnwave_pheno,subset= loc== "MW", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      pnratio_bluesMW <- bind_rows(pnratio_bluesMW, temp)
      cat('\n')

    }  else {
      pnratio_bluesMW <- bind_rows(pnratio_bluesMW, data.frame(name2 = NA,
                                                               predicted.value = NA,
                                                               wave = variables[i]))
      cat('\n')
    }
  },

  error = function(err) {
    message("An error occured")
    print(err)
  })
}
pnratio_bluesMW<- pnratio_bluesMW %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesMW, "./output/pnratio_bluesMW.csv", row.names = FALSE)

#calculating blues for each waveratio for loc EF fornarea
variables <- colnames(pnwave_pheno)[11:13]
models <- vector("list",length(variables))
pnratio_bluesEF <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = pnwave_pheno,subset= loc== "EF", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      pnratio_bluesEF <- bind_rows(pnratio_bluesEF, temp)
      cat('\n')

    }  else {
      pnratio_bluesEF <- bind_rows(pnratio_bluesEF, data.frame(name2 = NA,
                                                               predicted.value = NA,
                                                               wave = variables[i]))
      cat('\n')
    }
  },

  error = function(err) {
    message("An error occured")
    print(err)
  })
}

pnratio_bluesEF<- pnratio_bluesEF %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesEF, "./output/pnratio_bluesEF.csv", row.names = FALSE)
pnratio_bluesEF<- read.csv("./output/pnratio_bluesEF.csv")
pnratio_bluesMW<- read.csv("./output/pnratio_bluesMW.csv")

# ---------------------pnEF blues fitting model---------------------
pnbluesEF <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_pheno, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesEF<- update.asreml(pnbluesEF)
# ---------------------storing prediction---------------------
pnbluesEF <- data.frame(
  name2= pnbluesEF$predictions$pvals$name2,
  plsr_narea = round(pnbluesEF$predictions$pvals$predicted.value, 3))
pnbluesEF <-left_join(pnratio_bluesEF, pnbluesEF, by= "name2")
fwrite( pnbluesEF, "./output/pnbluesEF.csv", row.names = FALSE)
# ---------------------pnMW blues fitting model---------------------
pnbluesMW <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_pheno, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesMW<- update.asreml(pnbluesMW)
# ---------------------storing prediction---------------------
pnbluesMW <- data.frame(
  name2= pnbluesMW$predictions$pvals$name2,
  plsr_narea = round(pnbluesMW$predictions$pvals$predicted.value, 3))
pnbluesMW<-left_join(pnratio_bluesMW, pnbluesMW, by= "name2")
fwrite( pnbluesMW, "./output/pnbluesMW.csv", row.names = FALSE)
pnbluesEF<- read.csv("./output/pnbluesEF.csv")
pnbluesMW<- read.csv("./output/pnbluesMW.csv")

pnblues <- bind_rows(
  "EF" = pnbluesEF,
  "MW" = pnbluesMW,
  .id = "env")

# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues <-pnblues|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues<- subset(pnblues, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
pn_blues<- droplevels(pn_blues[pn_blues$taxa %in%rownames(kin), ]) #filtering indv that are present both in kin and phenotypic data
write.csv(pn_blues, "./output/pn_blues.csv", row.names = F)
#h2 for pn narea for both location
pnMW <- asreml(
  fixed = fs_plsr_narea ~ set,
  random = ~name2 + block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_pheno, subset= loc== "MW",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
pnMW<- update.asreml(pnMW)
# ---------------------calculating and storing heritability---------------------
h2<- (1 - ((pnMW$predictions$avsed["mean"] ^ 2) /(2 * summary(pnMW)$varcomp["name2", "component"]))) # pn h2= 0.42 for MW location

pnEF <- asreml(
  fixed = fs_plsr_narea ~ set,
  random = ~name2 + block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_pheno, subset= loc== "EF",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
pnEF<- update.asreml(pnEF)
# ---------------------calculating and storing heritability---------------------
h2<- (1 - ((pnEF$predictions$avsed["mean"] ^ 2) /(2 * summary(pnEF)$varcomp["name2", "component"]))) # pn h2= 0.383 for EF location

#plsr-sla(ps)
pswave_pheno <- read.csv("./output/pswave_pheno.csv")
# ---------------------processing of data---------------------
pswave_pheno<-
  pswave_pheno|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for wave ratio for loc MW
variables <- colnames(pswave_pheno)[11:13]
models <- vector("list",length(variables))
psratio_bluesMW <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = pswave_pheno,subset= loc== "MW", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      psratio_bluesMW <- bind_rows(psratio_bluesMW, temp)
      cat('\n')

    }  else {
      psratio_bluesMW <- bind_rows(psratio_bluesMW, data.frame(name2 = NA,
                                                               predicted.value = NA,
                                                               wave = variables[i]))
      cat('\n')
    }
  },

  error = function(err) {
    message("An error occured")
    print(err)
  })
}
psratio_bluesMW<- psratio_bluesMW %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesMW, "./output/psratio_bluesMW.csv", row.names = FALSE)

#calculating blues for each waveratio for loc EF forsla
variables <- colnames(pswave_pheno)[11:13]
models <- vector("list",length(variables))
psratio_bluesEF <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = pswave_pheno,subset= loc== "EF", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      psratio_bluesEF <- bind_rows(psratio_bluesEF, temp)
      cat('\n')

    }  else {
      psratio_bluesEF <- bind_rows(psratio_bluesEF, data.frame(name2 = NA,
                                                               predicted.value = NA,
                                                               wave = variables[i]))
      cat('\n')
    }
  },

  error = function(err) {
    message("An error occured")
    print(err)
  })
}

psratio_bluesEF<- psratio_bluesEF %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesEF, "./output/psratio_bluesEF.csv", row.names = FALSE)


# ---------------------psEF blues fitting model---------------------
psbluesEF <- asreml(
  fixed = plsr_sla_sorghum ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pswave_pheno, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
psbluesEF<- update.asreml(psbluesEF)
# ---------------------storing prediction---------------------
psbluesEF <- data.frame(
  name2= psbluesEF$predictions$pvals$name2,
  plsr_sla = round(psbluesEF$predictions$pvals$predicted.value, 3))
psbluesEF <-left_join(psratio_bluesEF, psbluesEF, by= "name2")
fwrite( psbluesEF, "./output/psbluesEF.csv", row.names = FALSE)
# ---------------------psMW blues fitting model---------------------
psbluesMW <- asreml(
  fixed = plsr_sla_sorghum ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pswave_pheno, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
psbluesMW<- update.asreml(psbluesMW)
# ---------------------storing prediction---------------------
psbluesMW <- data.frame(
  name2= psbluesMW$predictions$pvals$name2,
  plsr_sla = round(psbluesMW$predictions$pvals$predicted.value, 3))
psbluesMW<-left_join(psratio_bluesMW, psbluesMW, by= "name2")
fwrite( psbluesMW, "./output/psbluesMW.csv", row.names = FALSE)
psbluesEF<- read.csv("./output/psbluesEF.csv")
psbluesMW<- read.csv("./output/psbluesMW.csv")

psblues <- bind_rows(
  "EF" = psbluesEF,
  "MW" = psbluesMW,
  .id = "env")

# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
psblues <-psblues|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
ps_blues<- subset(psblues, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
ps_blues<- droplevels(ps_blues[ps_blues$taxa %in%rownames(kin), ]) #filtering indv that are present both in kin and phenotypic data
write.csv(ps_blues, "./output/ps_blues.csv", row.names = F)


#h2 for plsr sla for both locations
#h2 for pn narea for both location
psMW <- asreml(
  fixed = plsr_sla_sorghum ~ set,
  random = ~name2 + block,
  residual =  ~ar1(range):ar1(row),
  data = pswave_pheno, subset= loc== "MW",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
psMW<- update.asreml(psMW)
# ---------------------calculating and storing heritability---------------------
h2<- (1 - ((psMW$predictions$avsed["mean"] ^ 2) /(2 * summary(psMW)$varcomp["name2", "component"]))) # ps h2= 0.24 for MW location

psEF <- asreml(
  fixed = plsr_sla_sorghum ~ set,
  random = ~name2 + block,
  residual =  ~ar1(range):ar1(row),
  data = pswave_pheno, subset= loc== "EF",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
psEF<- update.asreml(psEF)
# ---------------------calculating and storing heritability---------------------
h2<- (1 - ((psEF$predictions$avsed["mean"] ^ 2) /(2 * summary(psEF)$varcomp["name2", "component"]))) # ps h2= 0.287 for EF location
