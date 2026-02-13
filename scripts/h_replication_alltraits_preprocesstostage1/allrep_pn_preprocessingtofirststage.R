# --------------------------------------------------------
# ---------------------Processing Coh2---------------------
# --------------------------------------------------------
library(data.table)
library(tidyverse)
library(janitor)
library(Hmisc)
library(corrplot)
library(factoextra)
library(here)
library(asreml)
library(fs)
library(knitr)
library(kableExtra)
library(dplyr)
library(sparkline)
library(htmlwidgets)
library(gt)
source("./function/aux_function.R")

pnratio_transform_rep1<-read.csv("./output/pnratio_transform_rep1.csv")
pnratio_transform_rep2<-read.csv("./output/pnratio_transform_rep2.csv")
pnratio_transform_rep3<-read.csv("./output/pnratio_transform_rep3.csv")
pnratio_transform_rep4<-read.csv("./output/pnratio_transform_rep4.csv")
pnratio_transform_rep5<-read.csv("./output/pnratio_transform_rep5.csv")
pnratio_combined<- rbind(pnratio_transform_rep1,pnratio_transform_rep2,pnratio_transform_rep3,pnratio_transform_rep4,pnratio_transform_rep5 )
write.csv(pnratio_combined,"./output/pnratio_combined.csv", row.names = F)
pnratio_combined<- read.csv("./output/pnratio_combined.csv")


#for five replications  first stage analysis
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
phenotypes<- phenotypes_whole %>% clean_names()
pnratio_transform<- read.csv("./output/pnratio_combined.csv")

# Compute and append ratios for each wave pair
for (i in 1:nrow(pnratio_transform)) {
  wave1 <- pnratio_transform$wave_1[i]
  wave2 <- pnratio_transform$wave_2[i]
  if (wave1 %in% colnames(phenotypes) && wave2 %in% colnames(phenotypes)) {
    # Calculate the ratio
    ratio_column_name <- paste(wave1, wave2, sep = "_")
    phenotypes[[ratio_column_name]] <- phenotypes[[wave1]] / phenotypes[[wave2]]
  } else {
    warning(paste("Missing columns:", wave1, "or", wave2, "in phenotypes."))
  }
}

phenotypes<- phenotypes %>% select(-c(14:2166),-narea,-n_perc,-sla)
write.csv(phenotypes,"./data/phenotypes_allReppn.csv", row.names = F)
phenotypes_allRep<- read.csv("./data/phenotypes_allReppn.csv")
#obtaining blues for replication pn
pnwave_pheno <- read.csv("./data/phenotypes_allReppn.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(pnwave_pheno)[11:25]
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
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
pnratiobluesMW_rep1<-pnratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(1:4)
pnratiobluesMW_rep2<-pnratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2, 5:7)
pnratiobluesMW_rep3<-pnratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2, 8:10)
pnratiobluesMW_rep4<-pnratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2, 11:13)
pnratiobluesMW_rep5<-pnratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2, 14:16)
fwrite(pnratiobluesMW_rep5, "./output/pnratiobluesMW_rep5.csv", row.names = FALSE)
fwrite(pnratiobluesMW_rep4, "./output/pnratiobluesMW_rep4.csv", row.names = FALSE)
fwrite(pnratiobluesMW_rep3, "./output/pnratiobluesMW_rep3.csv", row.names = FALSE)
fwrite(pnratiobluesMW_rep2, "./output/pnratiobluesMW_rep2.csv", row.names = FALSE)
fwrite(pnratiobluesMW_rep1, "./output/pnratiobluesMW_rep1.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pnwave_pheno)[11:25]
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
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
pnratiobluesEF_rep1<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(name2,2:4)
pnratiobluesEF_rep2<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2,5:7)
pnratiobluesEF_rep3<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2,8:10)
pnratiobluesEF_rep4<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2,11:13)
pnratiobluesEF_rep5<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2,14:16)

fwrite( pnratiobluesEF_rep5, "./output/pnratiobluesEF_rep5.csv", row.names = FALSE)
fwrite( pnratiobluesEF_rep4, "./output/pnratiobluesEF_rep4.csv", row.names = FALSE)
fwrite( pnratiobluesEF_rep3, "./output/pnratiobluesEF_rep3.csv", row.names = FALSE)
fwrite( pnratiobluesEF_rep2, "./output/pnratiobluesEF_rep2.csv", row.names = FALSE)
fwrite( pnratiobluesEF_rep1, "./output/pnratiobluesEF_rep1.csv", row.names = FALSE)

# ---------------------pnEF blues fitting model---------------------
pnwave_pheno <- read.csv("./data/phenotypes_whole.csv") %>% clean_names()
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
pnbluesEF <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  data = pnwave_pheno, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesEF<- update.asreml(pnbluesEF)
# ---------------------storing prediction---------------------
pnbluesEF <- data.frame(
  name2= pnbluesEF$predictions$pvals$name2,
  pn = round(pnbluesEF$predictions$pvals$predicted.value, 3))
pnbluesEF_rep1<-left_join(pnratiobluesEF_rep1, pnbluesEF, by= "name2")
pnbluesEF_rep2<-left_join(pnratiobluesEF_rep2, pnbluesEF, by= "name2")
pnbluesEF_rep3<-left_join(pnratiobluesEF_rep3, pnbluesEF, by= "name2")
pnbluesEF_rep4<-left_join(pnratiobluesEF_rep4, pnbluesEF, by= "name2")
pnbluesEF_rep5<-left_join(pnratiobluesEF_rep5, pnbluesEF, by= "name2")
fwrite(pnbluesEF_rep4, "./output/pnbluesEF_rep4.csv", row.names = FALSE)
fwrite(pnbluesEF_rep5, "./output/pnbluesEF_rep5.csv", row.names = FALSE)
fwrite(pnbluesEF_rep3, "./output/pnbluesEF_rep3.csv", row.names = FALSE)
fwrite(pnbluesEF_rep2, "./output/pnbluesEF_rep2.csv", row.names = FALSE)
fwrite(pnbluesEF_rep1, "./output/pnbluesEF_rep1.csv", row.names = FALSE)

# ---------------------NareaMW blues fitting model---------------------
pnbluesMW <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  residual =  ~ ar1(range):ar1(row),
  data = pnwave_pheno, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesMW<- update.asreml(pnbluesMW)
# ---------------------storing prediction---------------------
pnbluesMW <- data.frame(
  name2= pnbluesMW$predictions$pvals$name2,
  pn = round(pnbluesMW$predictions$pvals$predicted.value, 3))
pnbluesMW_rep1<-left_join(pnratiobluesMW_rep1, pnbluesMW, by= "name2")
pnbluesMW_rep2<-left_join(pnratiobluesMW_rep2, pnbluesMW, by= "name2")
pnbluesMW_rep3<-left_join(pnratiobluesMW_rep3, pnbluesMW, by= "name2")
pnbluesMW_rep4<-left_join(pnratiobluesMW_rep4, pnbluesMW, by= "name2")
pnbluesMW_rep5<-left_join(pnratiobluesMW_rep5, pnbluesMW, by= "name2")
fwrite( pnbluesMW_rep4, "./output/pnbluesMW_rep4.csv", row.names = FALSE)
fwrite( pnbluesMW_rep5, "./output/pnbluesMW_rep5.csv", row.names = FALSE)
fwrite( pnbluesMW_rep3, "./output/pnbluesMW_rep3.csv", row.names = FALSE)
fwrite( pnbluesMW_rep2, "./output/pnbluesMW_rep2.csv", row.names = FALSE)
fwrite( pnbluesMW_rep1, "./output/pnbluesMW_rep1.csv", row.names = FALSE)

#combining blues for trait and waveratio for EF and MW location for rep1 to 5

Names_WEST <- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST) == "Name2"] <- "name2"
# load kinship
kin <- fread("./data/kinship_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

for (i in 1:5) {
  process_pn_blues(rep_id = i, kin = kin, names_west = Names_WEST)
}




# pn- lowest coheritable syntehtic trait selection for rep 1-5
get_coh2_0_syntrait(
  trait_col = "fs_plsr_narea",
  breakdown_prefix = "pn",
  output_prefix = "pnwave",
  keep_cols = 1:9
)


#calculating blues for the selected low coheritable synthetic trait
pnwave_lowcoh2<- read.csv("./output/pnwave_rep1_lowcoh2.csv")
# ---------------------processing of data---------------------
pnwave_lowcoh2<-
  pnwave_lowcoh2|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for each wavelength for location MW
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesMW1 <- bind_rows(pnratio_bluesMW1, temp)
pnratio_bluesMW1<- pnratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesMW1, "./output/pnratio_bluesMW1_rep1.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesEF1 <- bind_rows(pnratio_bluesEF1, temp)

pnratio_bluesEF1<- pnratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesEF1, "./output/pnratio_bluesEF1_rep1.csv", row.names = FALSE)
pnratio_bluesEF1<- read.csv("./output/pnratio_bluesEF1_rep1.csv")

# ---------------------pnEF blues fitting model---------------------
pnbluesEF1 <- asreml(
  fixed =  fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesEF1<- update.asreml(pnbluesEF1)
# ---------------------storing prediction---------------------
pnbluesEF1 <- data.frame(
  name2= pnbluesEF1$predictions$pvals$name2,
  pn = round(pnbluesEF1$predictions$pvals$predicted.value, 3))
pnbluesEF1<-left_join(pnratio_bluesEF1, pnbluesEF1, by= "name2")
fwrite( pnbluesEF1, "./output/pnbluesEF1_rep1.csv", row.names = FALSE)
# ---------------------pnMW blues fitting model---------------------
pnbluesMW1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesMW1<- update.asreml(pnbluesMW1)
# ---------------------storing prediction---------------------
pnbluesMW1 <- data.frame(
  name2= pnbluesMW1$predictions$pvals$name2,
  pn = round(pnbluesMW1$predictions$pvals$predicted.value, 3))
pnbluesMW1<-left_join(pnratio_bluesMW1, pnbluesMW1, by= "name2")
fwrite( pnbluesMW1, "./output/pnbluesMW1_rep1.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
pnbluesEF1<- read.csv("./output/pnbluesEF1_rep1.csv")
pnbluesMW1<- read.csv("./output/pnbluesMW1_rep1.csv")

pnblues1 <- bind_rows(
  "EF" = pnbluesEF1,
  "MW" = pnbluesMW1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues1 <-pnblues1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues1<- subset(pnblues1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
pn_blues1<- droplevels(pn_blues1[pn_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(pn_blues1,"./output/pn_blues1_rep1.csv", row.names = F)


##calculating blues for the selected low coheritable synthetic trait
pnwave_lowcoh2<- read.csv("./output/pnwave_rep2_lowcoh2.csv")
# ---------------------processing of data---------------------
pnwave_lowcoh2<-
  pnwave_lowcoh2|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for each wavelength for location MW
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesMW1 <- bind_rows(pnratio_bluesMW1, temp)
pnratio_bluesMW1<- pnratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesMW1, "./output/pnratio_bluesMW1_rep2.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesEF1 <- bind_rows(pnratio_bluesEF1, temp)

pnratio_bluesEF1<- pnratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesEF1, "./output/pnratio_bluesEF1_rep2.csv", row.names = FALSE)
pnratio_bluesEF1<- read.csv("./output/pnratio_bluesEF1_rep2.csv")

# ---------------------pnEF blues fitting model---------------------
pnbluesEF1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesEF1<- update.asreml(pnbluesEF1)
# ---------------------storing prediction---------------------
pnbluesEF1 <- data.frame(
  name2= pnbluesEF1$predictions$pvals$name2,
  pn = round(pnbluesEF1$predictions$pvals$predicted.value, 3))
pnbluesEF1<-left_join(pnratio_bluesEF1, pnbluesEF1, by= "name2")
fwrite( pnbluesEF1, "./output/pnbluesEF1_rep2.csv", row.names = FALSE)
# ---------------------pnMW blues fitting model---------------------
pnbluesMW1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesMW1<- update.asreml(pnbluesMW1)
# ---------------------storing prediction---------------------
pnbluesMW1 <- data.frame(
  name2= pnbluesMW1$predictions$pvals$name2,
  pn = round(pnbluesMW1$predictions$pvals$predicted.value, 3))
pnbluesMW1<-left_join(pnratio_bluesMW1, pnbluesMW1, by= "name2")
fwrite( pnbluesMW1, "./output/pnbluesMW1_rep2.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
pnbluesEF1<- read.csv("./output/pnbluesEF1_rep2.csv")
pnbluesMW1<- read.csv("./output/pnbluesMW1_rep2.csv")

pnblues1 <- bind_rows(
  "EF" = pnbluesEF1,
  "MW" = pnbluesMW1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues1 <-pnblues1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues1<- subset(pnblues1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
pn_blues1<- droplevels(pn_blues1[pn_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(pn_blues1,"./output/pn_blues1_rep2.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
pnwave_lowcoh2<- read.csv("./output/pnwave_rep3_lowcoh2.csv")
# ---------------------processing of data---------------------
pnwave_lowcoh2<-
  pnwave_lowcoh2|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for each wavelength for location MW
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesMW1 <- bind_rows(pnratio_bluesMW1, temp)
pnratio_bluesMW1<- pnratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesMW1, "./output/pnratio_bluesMW1_rep3.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesEF1 <- bind_rows(pnratio_bluesEF1, temp)

pnratio_bluesEF1<- pnratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesEF1, "./output/pnratio_bluesEF1_rep3.csv", row.names = FALSE)
pnratio_bluesEF1<- read.csv("./output/pnratio_bluesEF1_rep3.csv")

# ---------------------pnEF blues fitting model---------------------
pnbluesEF1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesEF1<- update.asreml(pnbluesEF1)
# ---------------------storing prediction---------------------
pnbluesEF1 <- data.frame(
  name2= pnbluesEF1$predictions$pvals$name2,
  pn = round(pnbluesEF1$predictions$pvals$predicted.value, 3))
pnbluesEF1<-left_join(pnratio_bluesEF1, pnbluesEF1, by= "name2")
fwrite( pnbluesEF1, "./output/pnbluesEF1_rep3.csv", row.names = FALSE)
# ---------------------pnMW blues fitting model---------------------
pnbluesMW1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesMW1<- update.asreml(pnbluesMW1)
# ---------------------storing prediction---------------------
pnbluesMW1 <- data.frame(
  name2= pnbluesMW1$predictions$pvals$name2,
  pn = round(pnbluesMW1$predictions$pvals$predicted.value, 3))
pnbluesMW1<-left_join(pnratio_bluesMW1, pnbluesMW1, by= "name2")
fwrite( pnbluesMW1, "./output/pnbluesMW1_rep3.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
pnbluesEF1<- read.csv("./output/pnbluesEF1_rep3.csv")
pnbluesMW1<- read.csv("./output/pnbluesMW1_rep3.csv")

pnblues1 <- bind_rows(
  "EF" = pnbluesEF1,
  "MW" = pnbluesMW1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues1 <-pnblues1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues1<- subset(pnblues1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
pn_blues1<- droplevels(pn_blues1[pn_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(pn_blues1,"./output/pn_blues1_rep3.csv", row.names = F)



###calculating blues for the selected low coheritable synthetic trait
pnwave_lowcoh2<- read.csv("./output/pnwave_rep4_lowcoh2.csv")
# ---------------------processing of data---------------------
pnwave_lowcoh2<-
  pnwave_lowcoh2|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for each wavelength for location MW
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesMW1 <- bind_rows(pnratio_bluesMW1, temp)
pnratio_bluesMW1<- pnratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesMW1, "./output/pnratio_bluesMW1_rep4.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesEF1 <- bind_rows(pnratio_bluesEF1, temp)

pnratio_bluesEF1<- pnratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesEF1, "./output/pnratio_bluesEF1_rep4.csv", row.names = FALSE)
pnratio_bluesEF1<- read.csv("./output/pnratio_bluesEF1_rep4.csv")

# ---------------------pnEF blues fitting model---------------------
pnbluesEF1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesEF1<- update.asreml(pnbluesEF1)
# ---------------------storing prediction---------------------
pnbluesEF1 <- data.frame(
  name2= pnbluesEF1$predictions$pvals$name2,
  pn = round(pnbluesEF1$predictions$pvals$predicted.value, 3))
pnbluesEF1<-left_join(pnratio_bluesEF1, pnbluesEF1, by= "name2")
fwrite( pnbluesEF1, "./output/pnbluesEF1_rep4.csv", row.names = FALSE)
# ---------------------pnMW blues fitting model---------------------
pnbluesMW1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesMW1<- update.asreml(pnbluesMW1)
# ---------------------storing prediction---------------------
pnbluesMW1 <- data.frame(
  name2= pnbluesMW1$predictions$pvals$name2,
  pn = round(pnbluesMW1$predictions$pvals$predicted.value, 3))
pnbluesMW1<-left_join(pnratio_bluesMW1, pnbluesMW1, by= "name2")
fwrite( pnbluesMW1, "./output/pnbluesMW1_rep4.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
pnbluesEF1<- read.csv("./output/pnbluesEF1_rep4.csv")
pnbluesMW1<- read.csv("./output/pnbluesMW1_rep4.csv")

pnblues1 <- bind_rows(
  "EF" = pnbluesEF1,
  "MW" = pnbluesMW1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues1 <-pnblues1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues1<- subset(pnblues1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
pn_blues1<- droplevels(pn_blues1[pn_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(pn_blues1,"./output/pn_blues1_rep4.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
pnwave_lowcoh2<- read.csv("./output/pnwave_rep5_lowcoh2.csv")
# ---------------------processing of data---------------------
pnwave_lowcoh2<-
  pnwave_lowcoh2|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for each wavelength for location MW
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesMW1 <- bind_rows(pnratio_bluesMW1, temp)
pnratio_bluesMW1<- pnratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesMW1, "./output/pnratio_bluesMW1_rep5.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
pnratio_bluesEF1 <- bind_rows(pnratio_bluesEF1, temp)

pnratio_bluesEF1<- pnratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesEF1, "./output/pnratio_bluesEF1_rep5.csv", row.names = FALSE)
pnratio_bluesEF1<- read.csv("./output/pnratio_bluesEF1_rep5.csv")

# ---------------------pnEF blues fitting model---------------------
pnbluesEF1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesEF1<- update.asreml(pnbluesEF1)
# ---------------------storing prediction---------------------
pnbluesEF1 <- data.frame(
  name2= pnbluesEF1$predictions$pvals$name2,
  pn = round(pnbluesEF1$predictions$pvals$predicted.value, 3))
pnbluesEF1<-left_join(pnratio_bluesEF1, pnbluesEF1, by= "name2")
fwrite( pnbluesEF1, "./output/pnbluesEF1_rep5.csv", row.names = FALSE)
# ---------------------pnMW blues fitting model---------------------
pnbluesMW1 <- asreml(
  fixed = fs_plsr_narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pnwave_lowcoh2, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
pnbluesMW1<- update.asreml(pnbluesMW1)
# ---------------------storing prediction---------------------
pnbluesMW1 <- data.frame(
  name2= pnbluesMW1$predictions$pvals$name2,
  pn = round(pnbluesMW1$predictions$pvals$predicted.value, 3))
pnbluesMW1<-left_join(pnratio_bluesMW1, pnbluesMW1, by= "name2")
fwrite( pnbluesMW1, "./output/pnbluesMW1_rep5.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
pnbluesEF1<- read.csv("./output/pnbluesEF1_rep5.csv")
pnbluesMW1<- read.csv("./output/pnbluesMW1_rep5.csv")

pnblues1 <- bind_rows(
  "EF" = pnbluesEF1,
  "MW" = pnbluesMW1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues1 <-pnblues1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues1<- subset(pnblues1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
pn_blues1<- droplevels(pn_blues1[pn_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(pn_blues1,"./output/pn_blues1_rep5.csv", row.names = F)
