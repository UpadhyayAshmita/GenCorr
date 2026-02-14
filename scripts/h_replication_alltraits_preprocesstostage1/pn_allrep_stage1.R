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

pnpreprocessing_allrep(rep_id= 1, breakdown_dir = "./output", # change rep id = 1 to 5 to get all 5 reps pn output
                               phenotypes_path = "./data/phenotypes_whole.csv",
                               sample_path = "./data/sample_per_rep.csv",
                               out_dir = "./output",
                               k_groups = 3,
                               pn_trait_col = "fs_plsr_narea")

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

#selecting synthetic trait that has lowest coh2 i.e 0
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = TRUE, sep = ",")
indv_sample <- read.csv("./data/sample_per_rep.csv")
selected_list <- vector("list", 5)
# calling this function to get pn- lowest coheritable syntehtic trait selection for rep 1-5
get_coh2_0_syntrait(
  trait_col = "fs_plsr_narea",
  breakdown_prefix = "pn",
  output_prefix = "pnwave",
  keep_cols = 1:9
)


#calculating blues for the selected low coheritable synthetic trait
# pn (inputs: pnwave_repX_lowcoh2.csv; outputs: pnratio_*, pnblues*)
run_blues_all_reps_oldnames("fs_plsr_narea", "pnwave", "pn")
