library(asreml)
library(dplyr)
# Select the one wave ratio with the lowest coh2
selected <- N |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Nwave_lowcoh2<- wave %>% select(1:10, narea,wave_1812,wave_1574)
wave_1812_wave_1574 <- paste("wave_1812", "wave_1574", sep = "_")
Nwave_lowcoh2[[wave_1812_wave_1574]] <- Nwave_lowcoh2[["wave_1812"]] / Nwave_lowcoh2[["wave_1574"]]
Nwave_lowcoh2<- Nwave_lowcoh2 %>% select(-wave_1812,-wave_1574)
write.csv(Nwave_lowcoh2, "./output/Nwave_lowcoh2.csv", row.names = F)

#for sla
# Select the one wave ratio with the lowest coh2
selected <- S |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Swave_lowcoh2<- wave %>% select(1:9, sla,wave_1812,wave_1799)
wave_1812_wave_1799 <- paste("wave_1812", "wave_1799", sep = "_")
Swave_lowcoh2[[wave_1812_wave_1799]] <- Swave_lowcoh2[["wave_1812"]] / Swave_lowcoh2[["wave_1799"]]
Swave_lowcoh2<- Swave_lowcoh2 %>% select(-wave_1812,-wave_1799)
write.csv(Swave_lowcoh2, "./output/Swave_lowcoh2.csv", row.names = F)

#for plsr-narea
# Select the one wave ratio with the lowest coh2
selected <- pn |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pnwave_lowcoh2<- wave %>% select(1:9, fs_plsr_narea,wave_2285,wave_1459)
wave_2285_wave_1459 <- paste("wave_2285", "wave_1459", sep = "_")
pnwave_lowcoh2[[wave_2285_wave_1459]] <- pnwave_lowcoh2[["wave_2285"]] / pnwave_lowcoh2[["wave_1459"]]
pnwave_lowcoh2<- pnwave_lowcoh2 %>% select(-wave_2285,-wave_1459)
write.csv(pnwave_lowcoh2, "./output/pnwave_lowcoh2.csv", row.names = F)

#for plsr-sla
# Select the one wave ratio with the lowest coh2
selected <- ps |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pswave_lowcoh2<- wave %>% select(1:9, plsr_sla_sorghum,wave_1641,wave_638)
wave_1641_wave_638<- paste("wave_1641", "wave_638", sep = "_")

pswave_lowcoh2[[wave_1641_wave_638]] <- pswave_lowcoh2[["wave_1641"]] / pswave_lowcoh2[["wave_638"]]
pswave_lowcoh2<- pswave_lowcoh2 %>% select(-wave_1641,-wave_638)
write.csv(pswave_lowcoh2, "./output/pswave_lowcoh2.csv", row.names = F)

Nwave_lowcoh2<- read.csv("./output/Nwave_lowcoh2.csv")
# ---------------------processing of data---------------------
Nwave_lowcoh2<-
  Nwave_lowcoh2|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for each wavelength for location MW
variables <- colnames(Nwave_lowcoh2)[11]
models <- vector("list",length(variables))
Nratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
    model <- asreml(
      fixed = get(variables) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = Nwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
      model <- update.asreml(model)
      #---------------------storing prediction---------------------
      temp <- model$predictions$pvals[, 1:2]
      temp$wave <- variables
      Nratio_bluesMW1 <- bind_rows(Nratio_bluesMW1, temp)
Nratio_bluesMW1<- Nratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesMW1, "./output/Nratio_bluesMW_lowcoh2.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Nwave_lowcoh2)[11]
models <- vector("list",length(variables))
Nratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
      fixed = get(variables) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = Nwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
      model <- update.asreml(model)
  #---------------------storing prediction---------------------
      temp <- model$predictions$pvals[, 1:2]
      temp$wave <- variables
      Nratio_bluesEF1 <- bind_rows(Nratio_bluesEF1, temp)

Nratio_bluesEF1<- Nratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesEF1, "./output/Nratio_bluesEF_lowcoh2.csv", row.names = FALSE)
Nratio_bluesEF<- read.csv("./output/Nratio_bluesEF_lowcoh2.csv")

# ---------------------NareaEF blues fitting model---------------------
nareabluesEF1 <- asreml(
  fixed = narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = Nwave_lowcoh2, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
nareabluesEF1<- update.asreml(nareabluesEF1)
# ---------------------storing prediction---------------------
nareabluesEF1 <- data.frame(
  name2= nareabluesEF1$predictions$pvals$name2,
  narea = round(nareabluesEF1$predictions$pvals$predicted.value, 3))
NbluesEF1<-left_join(Nratio_bluesEF1, nareabluesEF1, by= "name2")
fwrite( NbluesEF1, "./output/NbluesEF_lowcoh2.csv", row.names = FALSE)
# ---------------------NareaMW blues fitting model---------------------
nareabluesMW1 <- asreml(
  fixed = narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = Nwave_lowcoh2, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
nareabluesMW1<- update.asreml(nareabluesMW1)
# ---------------------storing prediction---------------------
nareabluesMW1 <- data.frame(
  name2= nareabluesMW1$predictions$pvals$name2,
  narea = round(nareabluesMW1$predictions$pvals$predicted.value, 3))
NbluesMW1<-left_join(Nratio_bluesMW1, nareabluesMW1, by= "name2")
fwrite( NbluesMW1, "./output/NbluesMW.csv_lowcoh2", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
NbluesEF1<- read.csv("./output/NbluesEF.csv")
NbluesMW1<- read.csv("./output/NbluesMW.csv")

Nblues1 <- bind_rows(
  "EF" = NbluesEF1,
  "MW" = NbluesMW1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Nblues1 <-Nblues1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
N_blues1<- subset(Nblues1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
N_blues1<- droplevels(N_blues1[N_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(N_blues1,"./output/N_blues_lowcoh2.csv", row.names = F)
#SLA blues
Swave_lowcoh2 <- read.csv("./output/Swave_lowcoh2.csv")
# ---------------------processing of data---------------------
Swave_lowcoh2<-
  Swave_lowcoh2|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

variable <- colnames(Swave_lowcoh2)[11]
cat(variable, '\n')
Sratio_bluesMW1 <-data.frame()
#fit the model
model <- asreml(
  fixed = get(variable) ~ set + name2,
  random = ~block,
  residual = ~ar1(range):ar1(row),
  data = Swave_lowcoh2,
  subset = loc == "MW",
  na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

model <- update.asreml(model)
model<- update.asreml(model)
# Process predictions
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variable
Sratio_bluesMW1 <- bind_rows(Sratio_bluesMW1, temp)
Sratio_bluesMW1<- Sratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesMW1, "./output/Sratio_bluesMW_lowcoh2.csv", row.names = FALSE)
Sratio_bluesMW1<- read.csv("./output/Sratio_bluesMW_lowcoh2.csv")
#calculating blues for each waveratio for loc EF for sla
variable <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesEF1 <-data.frame()
#fit the model
model <- asreml(
  fixed = get(variable) ~ set + name2,
  random = ~block,
  residual = ~ar1(range):ar1(row),
  data = Swave_lowcoh2,
  subset = loc == "EF",
  na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

model <- update.asreml(model)
model<- update.asreml(model)
# Process predictions
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variable
Sratio_bluesEF1 <- bind_rows(Sratio_bluesEF1, temp)

Sratio_bluesEF1<- Sratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesEF1, "./output/Sratio_bluesEF_lowcoh2.csv", row.names = FALSE)


# ---------------------SareaEF blues fitting model---------------------
slabluesEF1 <- asreml(
  fixed = sla ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = Swave_lowcoh2, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
slabluesEF1<- update.asreml(slabluesEF1)
# ---------------------storing prediction---------------------
slabluesEF1 <- data.frame(
  name2= slabluesEF1$predictions$pvals$name2,
  sla = round(slabluesEF1$predictions$pvals$predicted.value, 3))
SbluesEF1 <-left_join(Sratio_bluesEF1, slabluesEF1, by= "name2")
fwrite( SbluesEF1, "./output/SbluesEF_lowcoh2.csv", row.names = FALSE)
# ---------------------SlaMW blues fitting model---------------------
slabluesMW1 <- asreml(
  fixed = sla ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = Swave_lowcoh2, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
slabluesMW1<- update.asreml(slabluesMW1)
# ---------------------storing prediction---------------------
slabluesMW1 <- data.frame(
  name2= slabluesMW1$predictions$pvals$name2,
  sla = round(slabluesMW1$predictions$pvals$predicted.value, 3))
SbluesMW1<-left_join(Sratio_bluesMW1, slabluesMW1, by= "name2")
fwrite( SbluesMW1, "./output/SbluesMW_lowcoh2.csv", row.names = FALSE)
SbluesEF1<- read.csv("./output/SbluesEF_lowcoh2.csv")
SbluesMW1<- read.csv("./output/SbluesMW_lowcoh2.csv")

Sblues1 <- bind_rows(
  "EF" = SbluesEF1,
  "MW" = SbluesMW1,
  .id = "env")

# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Sblues1 <-Sblues1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
S_blues1<- subset(Sblues1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
S_blues1<- droplevels(S_blues1[S_blues1$taxa %in%rownames(kin), ]) #filtering indv that are present both in kin and phenotypic data
write.csv(S_blues1, "./output/S_blues_lowcoh2.csv", row.names = F)

#plsr-narea(pn)
pnwave_lowcoh2 <- read.csv("./output/pnwave_lowcoh2.csv")
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

#calculating blues for wave ratio for loc MW
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
    model <- asreml(
      fixed = get(variables) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = pnwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
  model <- update.asreml(model)

  # Process predictions
  temp <- model$predictions$pvals[, 1:2]
  temp$wave <- variables
      #---------------------storing prediction---------------------
  pnratio_bluesMW1 <- bind_rows(pnratio_bluesMW1, temp)

pnratio_bluesMW1<- pnratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesMW1, "./output/pnratio_bluesMW_lowcoh2.csv", row.names = FALSE)

#calculating blues for each waveratio for loc EF fornarea
variables <- colnames(pnwave_lowcoh2)[11]
models <- vector("list",length(variables))
pnratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  residual =  ~ar1(range):ar1(row),
  data = pnwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)

# Process predictions
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
#---------------------storing prediction---------------------
pnratio_bluesEF1 <- bind_rows(pnratio_bluesEF1, temp)
pnratio_bluesEF1<- pnratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( pnratio_bluesEF1, "./output/pnratio_bluesEF_lowcoh2.csv", row.names = FALSE)
pnratio_bluesEF1<- read.csv("./output/pnratio_bluesEF_lowcoh2.csv")
pnratio_bluesMW1<- read.csv("./output/pnratio_bluesMW_lowcoh2.csv")

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
  plsr_narea = round(pnbluesEF1$predictions$pvals$predicted.value, 3))
pnbluesEF1 <-left_join(pnratio_bluesEF1, pnbluesEF1, by= "name2")
fwrite( pnbluesEF1, "./output/pnbluesEF_lowcoh2.csv", row.names = FALSE)
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
  plsr_narea = round(pnbluesMW1$predictions$pvals$predicted.value, 3))
pnbluesMW1<-left_join(pnratio_bluesMW1, pnbluesMW1, by= "name2")
fwrite( pnbluesMW1, "./output/pnbluesMW_lowcoh2.csv", row.names = FALSE)
pnbluesEF1<- read.csv("./output/pnbluesEF_lowcoh2.csv")
pnbluesMW1<- read.csv("./output/pnbluesMW_lowcoh2.csv")

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
pn_blues1<- droplevels(pn_blues1[pn_blues1$taxa %in%rownames(kin), ]) #filtering indv that are present both in kin and phenotypic data
write.csv(pn_blues1, "./output/pn_blues_lowcoh2.csv", row.names = F)


#plsr-sla(ps)
pswave_lowcoh2 <- read.csv("./output/pswave_lowcoh2.csv")
# ---------------------processing of data---------------------
pswave_lowcoh2<-
  pswave_lowcoh2|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for wave ratio for loc MW
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesMW1<-data.frame()
# ---------------------fitting model---------------------
    model <- asreml(
      fixed = get(variables) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = pswave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    model <- update.asreml(model)
      #---------------------storing prediction---------------------
      temp <- model$predictions$pvals[, 1:2]
      temp$wave <- variables
      psratio_bluesMW1 <- bind_rows(psratio_bluesMW1, temp)
      cat('\n')

psratio_bluesMW1<- psratio_bluesMW1%>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesMW1, "./output/psratio_bluesMW_lowcoh2.csv", row.names = FALSE)

#calculating blues for each waveratio for loc EF forsla
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
    model <- asreml(
      fixed = get(variables) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = pswave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
      model <- update.asreml(model)
     #---------------------storing prediction---------------------
      temp <- model$predictions$pvals[, 1:2]
      temp$wave <- variables
      psratio_bluesEF1 <- bind_rows(psratio_bluesEF1, temp)
psratio_bluesEF1<- psratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesEF1, "./output/psratio_bluesEF_lowcoh2.csv", row.names = FALSE)


# ---------------------psEF blues fitting model---------------------
psbluesEF1 <- asreml(
  fixed = plsr_sla_sorghum ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pswave_lowcoh2, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
psbluesEF1<- update.asreml(psbluesEF1)
# ---------------------storing prediction---------------------
psbluesEF1 <- data.frame(
  name2= psbluesEF1$predictions$pvals$name2,
  plsr_sla = round(psbluesEF1$predictions$pvals$predicted.value, 3))
psbluesEF1 <-left_join(psratio_bluesEF1, psbluesEF1, by= "name2")
fwrite( psbluesEF1, "./output/psbluesEF_lowcoh2.csv", row.names = FALSE)
# ---------------------psMW blues fitting model---------------------
psbluesMW1 <- asreml(
  fixed = plsr_sla_sorghum ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = pswave_lowcoh2, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
psbluesMW1<- update.asreml(psbluesMW1)
# ---------------------storing prediction---------------------
psbluesMW1 <- data.frame(
  name2= psbluesMW1$predictions$pvals$name2,
  plsr_sla = round(psbluesMW1$predictions$pvals$predicted.value, 3))
psbluesMW1<-left_join(psratio_bluesMW1, psbluesMW1, by= "name2")
fwrite( psbluesMW1, "./output/psbluesMW_lowcoh2.csv", row.names = FALSE)
psbluesEF1<- read.csv("./output/psbluesEF_lowcoh2.csv")
psbluesMW1<- read.csv("./output/psbluesMW_lowcoh2.csv")

psblues1 <- bind_rows(
  "EF" = psbluesEF1,
  "MW" = psbluesMW1,
  .id = "env")

# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
psblues1 <-psblues1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
ps_blues1<- subset(psblues1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
ps_blues1<- droplevels(ps_blues1[ps_blues1$taxa %in%rownames(kin), ]) #filtering indv that are present both in kin and phenotypic data
write.csv(ps_blues1, "./output/ps_blues_lowcoh2.csv", row.names = F)


