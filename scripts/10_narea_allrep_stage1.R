#loading library
library(data.table)
library(tidyverse)
library(janitor)
library(Hmisc)
library(corrplot)
library(factoextra)
library(here)
library(fs)
library(knitr)
library(kableExtra)
library(dplyr)
source("./function/aux_function.R")

#selecting the synthetic trait for all five reps for sla and combining with respective phenotypic file

# nareapreprocessing_allrep(rep_id,  # while calling script change your rep_id=1 to 5 serially to get all rep pre-processing output for narea
#                                       breakdown_dir = "./output",
#                                       phenotypes_path = "./data/phenotypes_whole.csv",
#                                       sample_path = "./data/sample_per_rep.csv",
#                                       out_dir = "./output",
#                                       k_groups = 3) #if you want to change your path and file name for input & output dir, change likewise

#load all the transformed dataset for each reps
Nratio_transform_rep1<- read.csv("./output/Nratio_transform_rep1.csv")
Nratio_transform_rep2<- read.csv("./output/Nratio_transform_rep2.csv")
Nratio_transform_rep3<- read.csv("./output/Nratio_transform_rep3.csv")
Nratio_transform_rep4<- read.csv("./output/Nratio_transform_rep4.csv")
Nratio_transform_rep5<- read.csv("./output/Nratio_transform_rep5.csv")

Nratio_combined<- rbind(Nratio_transform_rep1,Nratio_transform_rep2,Nratio_transform_rep3,Nratio_transform_rep4,Nratio_transform_rep5 )
write.csv(Nratio_combined,"./output/Nratio_combined.csv", row.names = F)
Nratio_combined<- read.csv("./output/Nratio_combined.csv")

#for five replications  first stage analysis
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
phenotypes<- phenotypes_whole %>% clean_names()
Nratio_transform<- read.csv("./output/Nratio_combined.csv")
Nratio_transform<- read.csv("./output/Nratio_transform_rep1.csv")
Nratio_transform<- read.csv("./output/Nratio_transform_rep2.csv")
Nratio_transform<- read.csv("./output/Nratio_transform_rep3.csv")
Nratio_transform<- read.csv("./output/Nratio_transform_rep4.csv")
Nratio_transform<- read.csv("./output/Nratio_transform_rep5.csv")

# Compute and append ratios for each wave pair
for (i in 1:nrow(Nratio_transform)) {
  wave1 <- Nratio_transform$wave_1[i]
  wave2 <- Nratio_transform$wave_2[i]
  if (wave1 %in% colnames(phenotypes) && wave2 %in% colnames(phenotypes)) {
    # Calculate the ratio
    ratio_column_name <- paste(wave1, wave2, sep = "_")
    phenotypes[[ratio_column_name]] <- phenotypes[[wave1]] / phenotypes[[wave2]]
  } else {
    warning(paste("Missing columns:", wave1, "or", wave2, "in phenotypes."))
  }
}

phenotypes<- phenotypes %>% select(-c(11:2166))
write.csv(phenotypes,"./data/phenotypes_allRep.csv", row.names = F)
phenotypes_allRep<- read.csv("./data/phenotypes_allRep.csv")
#obtaining blues for replication 1 narea
Nwave_pheno <- read.csv("./data/phenotypes_allRep.csv")
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
variables <- colnames(Nwave_pheno)[11:25]
models <- vector("list",length(variables))
Nratio_bluesMW <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
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
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
NratiobluesMW_rep1<-Nratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(1:4)
NratiobluesMW_rep2<-Nratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2, 5:7)
NratiobluesMW_rep3<-Nratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2, 8:10)
NratiobluesMW_rep4<-Nratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2, 11:13)
NratiobluesMW_rep5<-Nratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2, 14:16)
fwrite(NratiobluesMW_rep5, "./output/NratiobluesMW_rep5.csv", row.names = FALSE)
fwrite(NratiobluesMW_rep4, "./output/NratiobluesMW_rep4.csv", row.names = FALSE)
fwrite(NratiobluesMW_rep3, "./output/NratiobluesMW_rep3.csv", row.names = FALSE)
fwrite(NratiobluesMW_rep2, "./output/NratiobluesMW_rep2.csv", row.names = FALSE)
fwrite(NratiobluesMW_rep1, "./output/NratiobluesMW_rep1.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Nwave_pheno)[11:25]
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
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
NratiobluesEF_rep1<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(name2,2:4)
NratiobluesEF_rep2<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2,5:7)
NratiobluesEF_rep3<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2,8:10)
NratiobluesEF_rep4<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2,11:13)
NratiobluesEF_rep5<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2,14:16)

fwrite( NratiobluesEF_rep5, "./output/NratiobluesEF_rep5.csv", row.names = FALSE)
fwrite( NratiobluesEF_rep4, "./output/NratiobluesEF_rep4.csv", row.names = FALSE)
fwrite( NratiobluesEF_rep3, "./output/NratiobluesEF_rep3.csv", row.names = FALSE)
fwrite( NratiobluesEF_rep2, "./output/NratiobluesEF_rep2.csv", row.names = FALSE)
fwrite( NratiobluesEF_rep1, "./output/NratiobluesEF_rep1.csv", row.names = FALSE)

# ---------------------NareaEF blues fitting model---------------------
Nwave_pheno <- read.csv("./data/phenotypes_whole.csv") %>% clean_names()
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
NbluesEF_rep1<-left_join(NratiobluesEF_rep1, nareabluesEF, by= "name2")
NbluesEF_rep2<-left_join(NratiobluesEF_rep2, nareabluesEF, by= "name2")
NbluesEF_rep3<-left_join(NratiobluesEF_rep3, nareabluesEF, by= "name2")
NbluesEF_rep4<-left_join(NratiobluesEF_rep4, nareabluesEF, by= "name2")
NbluesEF_rep5<-left_join(NratiobluesEF_rep5, nareabluesEF, by= "name2")
fwrite(NbluesEF_rep5, "./output/NbluesEF_rep5.csv", row.names = FALSE)
fwrite(NbluesEF_rep4, "./output/NbluesEF_rep4.csv", row.names = FALSE)
fwrite(NbluesEF_rep3, "./output/NbluesEF_rep3.csv", row.names = FALSE)
fwrite(NbluesEF_rep2, "./output/NbluesEF_rep2.csv", row.names = FALSE)
fwrite(NbluesEF_rep1, "./output/NbluesEF_rep1.csv", row.names = FALSE)

# ---------------------NareaMW blues fitting model---------------------
nareabluesMW <- asreml(
  fixed = narea ~ name2 +set,
  random =  ~  block,
  residual =  ~ ar1(range):ar1(row),
  data = Nwave_pheno, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
nareabluesMW<- update.asreml(nareabluesMW)
# ---------------------storing prediction---------------------
nareabluesMW <- data.frame(
  name2= nareabluesMW$predictions$pvals$name2,
  narea = round(nareabluesMW$predictions$pvals$predicted.value, 3))
NbluesMW_rep1<-left_join(NratiobluesMW_rep1, nareabluesMW, by= "name2")
NbluesMW_rep2<-left_join(NratiobluesMW_rep2, nareabluesMW, by= "name2")
NbluesMW_rep3<-left_join(NratiobluesMW_rep3, nareabluesMW, by= "name2")
NbluesMW_rep4<-left_join(NratiobluesMW_rep4, nareabluesMW, by= "name2")
NbluesMW_rep5<-left_join(NratiobluesMW_rep5, nareabluesMW, by= "name2")
fwrite( NbluesMW_rep2, "./output/NbluesMW_rep2.csv", row.names = FALSE)
fwrite( NbluesMW_rep1, "./output/NbluesMW_rep1.csv", row.names = FALSE)
fwrite( NbluesMW_rep3, "./output/NbluesMW_rep3.csv", row.names = FALSE)
fwrite( NbluesMW_rep4, "./output/NbluesMW_rep4.csv", row.names = FALSE)
fwrite( NbluesMW_rep5, "./output/NbluesMW_rep5.csv", row.names = FALSE)

#combining blues for trait and waveratio for EF and MW location for rep1-5
Names_WEST <- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST) == "Name2"] <- "name2"
# load kinship
kin <- fread("./data/kinship_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)
# Processing all rep blues files
for (i in 1:5) {
  process_nareablues(
    rep_id = i,
    kin = kin,
    names_west = Names_WEST
  )
}


#selecting synthetic trait that has lowest coh2 i.e 0
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = TRUE, sep = ",")
indv_sample <- read.csv("./data/sample_per_rep.csv")
selected_list <- vector("list", 5)


##calling this function to get narea- lowest coheritable syntehtic trait selection for rep 1-5
get_coh2_0_syntrait(
  trait_col = "narea",
  breakdown_prefix = "narea",
  output_prefix = "Nwave",
  keep_cols = 1:10
)


#calculating blues for the selected low coheritable synthetic trait for narea
# narea (inputs: Nwave_repX_lowcoh2.csv; outputs: Nratio_*, Nblues*)
run_blues_all_reps_oldnames("narea", "Nwave", "N") #gives the blues/stage1 data for all rep for narea


#filtering the blues based on names west dataset
#load the kinship matrix
kin <- fread("./data/kinship_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)
clean_with_names_west_and_kin_all_reps(out_prefix = "N", kin = kin)


