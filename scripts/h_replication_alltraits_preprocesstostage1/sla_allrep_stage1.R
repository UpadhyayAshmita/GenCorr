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

#selecting the synthetic trait for all five reps for sla and combining with respective phenotypic file
slapreprocessing_allrep<- function(rep_id,
                                   breakdown_dir = "./output",
                                   phenotypes_path = "./data/phenotypes_whole.csv",
                                   sample_path = "./data/sample_per_rep.csv",
                                   out_dir = "./output",
                                   k_groups = 3)


Sratio_transform_rep1<-read.csv("./output/Sratio_transform_rep1.csv")
Sratio_transform_rep2<-read.csv("./output/Sratio_transform_rep2.csv")
Sratio_transform_rep3<-read.csv("./output/Sratio_transform_rep3.csv")
Sratio_transform_rep4<-read.csv("./output/Sratio_transform_rep4.csv")
Sratio_transform_rep5<-read.csv("./output/Sratio_transform_rep5.csv")
Sratio_combined<- rbind(Sratio_transform_rep1,Sratio_transform_rep2,Sratio_transform_rep3,Sratio_transform_rep4,Sratio_transform_rep5 )
write.csv(Sratio_combined,"./output/Sratio_combined.csv", row.names = F)
Sratio_combined<- read.csv("./output/Sratio_combined.csv")

#for five replications  first stage analysis
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
phenotypes<- phenotypes_whole %>% clean_names()
Sratio_transform<- read.csv("./output/Sratio_combined.csv")

# Compute and append ratios for each wave pair
for (i in 1:nrow(Sratio_transform)) {
  wave1 <- Sratio_transform$wave_1[i]
  wave2 <- Sratio_transform$wave_2[i]
  if (wave1 %in% colnames(phenotypes) && wave2 %in% colnames(phenotypes)) {
    # Calculate the ratio
    ratio_column_name <- paste(wave1, wave2, sep = "_")
    phenotypes[[ratio_column_name]] <- phenotypes[[wave1]] / phenotypes[[wave2]]
  } else {
    warning(paste("Missing columns:", wave1, "or", wave2, "in phenotypes."))
  }
}

phenotypes<- phenotypes %>% select(-c(13:2166),-narea,-n_perc)
write.csv(phenotypes,"./data/phenotypes_allRepsla.csv", row.names = F)
phenotypes_allRep<- read.csv("./data/phenotypes_allRepsla.csv")
#obtaining blues for replication  sla
Swave_pheno <- read.csv("./data/phenotypes_allRepsla.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(Swave_pheno)[11:25]
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

Sratio_bluesMW<- Sratio_bluesMW %>% pivot_wider(names_from = wave, values_from= predicted.value)
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
SratiobluesMW_rep1<-Sratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(1:4)
SratiobluesMW_rep2<-Sratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2, 5:7)
SratiobluesMW_rep3<-Sratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2, 8:10)
SratiobluesMW_rep4<-Sratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2, 11:13)
SratiobluesMW_rep5<-Sratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2, 14:16)
fwrite(SratiobluesMW_rep5, "./output/SratiobluesMW_rep5.csv", row.names = FALSE)
fwrite(SratiobluesMW_rep4, "./output/SratiobluesMW_rep4.csv", row.names = FALSE)
fwrite(SratiobluesMW_rep3, "./output/SratiobluesMW_rep3.csv", row.names = FALSE)
fwrite(SratiobluesMW_rep2, "./output/SratiobluesMW_rep2.csv", row.names = FALSE)
fwrite(SratiobluesMW_rep1, "./output/SratiobluesMW_rep1.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Swave_pheno)[11:25]
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
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/sample_per_rep.csv")
SratiobluesEF_rep1<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(name2,2:4)
SratiobluesEF_rep2<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2,5:7)
SratiobluesEF_rep3<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2,8:10)
SratiobluesEF_rep4<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2,11:13)
SratiobluesEF_rep5<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2,14:16)

fwrite( SratiobluesEF_rep4, "./output/SratiobluesEF_rep4.csv", row.names = FALSE)
fwrite( SratiobluesEF_rep5, "./output/SratiobluesEF_rep5.csv", row.names = FALSE)
fwrite( SratiobluesEF_rep3, "./output/SratiobluesEF_rep3.csv", row.names = FALSE)
fwrite( SratiobluesEF_rep2, "./output/SratiobluesEF_rep2.csv", row.names = FALSE)
fwrite( SratiobluesEF_rep1, "./output/SratiobluesEF_rep1.csv", row.names = FALSE)

# ---------------------slaEF blues fitting model---------------------
Swave_pheno <- read.csv("./data/phenotypes_whole.csv") %>% clean_names()
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
slabluesEF <- asreml(
  fixed = sla ~ name2 +set,
  random =  ~  block,
  data = Swave_pheno, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
slabluesEF<- update.asreml(slabluesEF)
# ---------------------storing prediction---------------------
slabluesEF <- data.frame(
  name2= slabluesEF$predictions$pvals$name2,
  sla = round(slabluesEF$predictions$pvals$predicted.value, 3))
SbluesEF_rep1<-left_join(SratiobluesEF_rep1, slabluesEF, by= "name2")
SbluesEF_rep2<-left_join(SratiobluesEF_rep2, slabluesEF, by= "name2")
SbluesEF_rep3<-left_join(SratiobluesEF_rep3, slabluesEF, by= "name2")
SbluesEF_rep4<-left_join(SratiobluesEF_rep4, slabluesEF, by= "name2")
SbluesEF_rep5<-left_join(SratiobluesEF_rep5, slabluesEF, by= "name2")
fwrite(SbluesEF_rep5, "./output/SbluesEF_rep5.csv", row.names = FALSE)
fwrite(SbluesEF_rep4, "./output/SbluesEF_rep4.csv", row.names = FALSE)
fwrite(SbluesEF_rep3, "./output/SbluesEF_rep3.csv", row.names = FALSE)
fwrite(SbluesEF_rep2, "./output/SbluesEF_rep2.csv", row.names = FALSE)
fwrite(SbluesEF_rep1, "./output/SbluesEF_rep1.csv", row.names = FALSE)

# ---------------------NareaMW blues fitting model---------------------
slabluesMW <- asreml(
  fixed = sla ~ name2 +set,
  random =  ~  block,
  residual =  ~ ar1(range):ar1(row),
  data = Swave_pheno, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
slabluesMW<- update.asreml(slabluesMW)
# ---------------------storing prediction---------------------
slabluesMW <- data.frame(
  name2= slabluesMW$predictions$pvals$name2,
  sla = round(slabluesMW$predictions$pvals$predicted.value, 3))
SbluesMW_rep1<-left_join(SratiobluesMW_rep1, slabluesMW, by= "name2")
SbluesMW_rep2<-left_join(SratiobluesMW_rep2, slabluesMW, by= "name2")
SbluesMW_rep3<-left_join(SratiobluesMW_rep3, slabluesMW, by= "name2")
SbluesMW_rep4<-left_join(SratiobluesMW_rep4, slabluesMW, by= "name2")
SbluesMW_rep5<-left_join(SratiobluesMW_rep5, slabluesMW, by= "name2")
fwrite( SbluesMW_rep5, "./output/SbluesMW_rep5.csv", row.names = FALSE)
fwrite( SbluesMW_rep4, "./output/SbluesMW_rep4.csv", row.names = FALSE)
fwrite( SbluesMW_rep3, "./output/SbluesMW_rep3.csv", row.names = FALSE)
fwrite( SbluesMW_rep2, "./output/SbluesMW_rep2.csv", row.names = FALSE)
fwrite( SbluesMW_rep1, "./output/SbluesMW_rep1.csv", row.names = FALSE)

#combining blues for trait and waveratio for EF and MW location for rep1-5
Names_WEST <- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST) == "Name2"] <- "name2"
# load kinship
kin <- fread("./data/kinship_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

for (i in 1:5) {
  process_sla_blues(rep_id = i, kin = kin, names_west = Names_WEST)
}


# getting lowest coheritability syntehtic trait for rep 1-5 for sla
get_coh2_0_syntrait(
  trait_col = "sla",
  breakdown_prefix = "sla",
  output_prefix = "Swave",
  keep_cols = 1:9
)

#calculating blues for the selected low coheritable synthetic trait
# sla (inputs: Swave_repX_lowcoh2.csv; outputs: Sratio_*, Sblues*)
run_blues_all_reps_oldnames("sla", "Swave", "S")


#filtering the blues based on names west dataset
#load the kinship matrix
kin <- fread("./data/kinship_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)
clean_with_names_west_and_kin_all_reps(out_prefix = "S", kin = kin)

