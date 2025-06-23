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
library(fs)
library(knitr)
library(kableExtra)
library(dplyr)
library(sparkline)
library(htmlwidgets)
source("./scripts/aux_function.R")


#replication 1
N1 <- fread("./data/narea_breakdown1.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_5)
wave <-phenotypes |>
  clean_names()

sum(is.na(N$corg))
sum(N$corg< -1 | N$corg > 1, na.rm = T)
sum(N$corgblup< -1 | N$corgblup > 1, na.rm = T)

N <- N |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

# sum(abs(N$corg - N$corgblup)>0.2, na.rm = T)
# 176644
sum(is.na(N$corg))
ntraits <- ceiling(sum(!is.na(N$coh2))*0.01)

# Calculate the wave ratio for the single row in 'sel'
#wave_ratio_name <- paste0("wave$", lowest_ntraits$wave_1, "_", lowest_ntraits$wave_2)
#assign(wave_ratio_name, wave[[lowest_ntraits$wave_1]] / wave[[lowest_ntraits$wave_2]])

#lowest_ntraits$rowname <- paste0(lowest_ntraits$wave_1, "_", lowest_ntraits$wave_2)
#Nratio_low <- wave |> select(wave_524_wave_681)

sel <- N |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Nratios <- wave |> select(wave_1693_wave_1696:wave_862_wave_376)
Nratios <- t(Nratios)
colnames(Nratios) <- wave$plot_id
#distance
dist <- get_dist(Nratios, method = "pearson")
# fviz_dist(dist)
hcN <- hclust(dist, method = "average")
sub_grp <- cutree(hcN, k = 3)
table(sub_grp)
Nratios <- Nratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Nratios)

Nratios<- as.data.frame(Nratios)
write.csv(Nratios,"./data/Nratios_rep1.csv", row.names=F)
Nratios<- read.csv("./data/Nratios_rep1.csv")

# Compute the absolute difference between wave_1 and wave_2
Nratios <- Nratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Nratio_sel <- Nratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Nratio_transform <- Nratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Nratio_transform, "./data/Nratio_transform_rep1.csv", row.names=F)
Nratio_transform<- read.csv("./data/Nratio_transform_rep1.csv")
com_col <- colnames(wave)[colnames(wave) %in% Nratio_transform$wave_sel]
Nwave_pheno<- wave %>% select(1:10, com_col)
write.csv(Nwave_pheno,"./data/Nwave_pheno_rep1.csv", row.names = F)
Nwave<- read.csv("./data/Nwave_pheno_rep1.csv")
#plotting

l <- Nratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()

# gg <- N |>
#   ggplot(aes(x = wave_1, y = wave_2, fill = coh2)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red",na.value="black") +
#   labs(x = "Wave2", y = "Wave1", title = "Co-h2") +
#   geom_label(data = l, aes(label = label))

heatmap_rep1 <-N |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                  wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Coh2-Replictaion1") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_rep1.jpeg")


#replication 2
N <- fread("./data/narea_breakdown2.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_2)
wave <-phenotypes |>
  clean_names()

sum(is.na(N$corg))
sum(N$corg< -1 | N$corg > 1, na.rm = T)
sum(N$corgblup< -1 | N$corgblup > 1, na.rm = T)

N <- N |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(N$corg))
ntraits <- ceiling(sum(!is.na(N$coh2))*0.01)
sel <- N |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Nratios <- wave |> select(wave_1659_wave_1716:wave_1855_wave_2231)
Nratios <- t(Nratios)
colnames(Nratios) <- wave$plot_id
#distance
dist <- get_dist(Nratios, method = "pearson")
# fviz_dist(dist)
hcN <- hclust(dist, method = "average")
sub_grp <- cutree(hcN, k = 3)
table(sub_grp)
Nratios <- Nratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Nratios)

Nratios<- as.data.frame(Nratios)
write.csv(Nratios,"./data/Nratios_rep2.csv", row.names=F)
Nratios<- read.csv("./data/Nratios_rep2.csv")

# Compute the absolute difference between wave_1 and wave_2
Nratios <- Nratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Nratio_sel <- Nratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Nratio_transform <- Nratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Nratio_transform, "./data/Nratio_transform_rep2.csv", row.names=F)
Nratio_transform<- read.csv("./data/Nratio_transform_rep2.csv")
com_col <- colnames(wave)[colnames(wave) %in% Nratio_transform$wave_sel]
Nwave_pheno<- wave %>% select(1:10, com_col)
write.csv(Nwave_pheno,"./data/Nwave_pheno_rep2.csv", row.names = F)
Nwave<- read.csv("./data/Nwave_pheno_rep2.csv")
#plotting

l <- Nratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()

heatmap_rep2 <-N |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Coh2-Replictaion2") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep2,"./figure/Heatmap_rep2.jpeg")

#replication 3
N <- fread("./data/narea_breakdown3.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_3)
wave <-phenotypes |>
  clean_names()

sum(is.na(N$corg))
sum(N$corg< -1 | N$corg > 1, na.rm = T)
sum(N$corgblup< -1 | N$corgblup > 1, na.rm = T)

N <- N |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(N$corg))
ntraits <- ceiling(sum(!is.na(N$coh2))*0.01)
sel <- N |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Nratios <- wave |> select(wave_1693_wave_1696:wave_862_wave_376)
Nratios <- t(Nratios)
colnames(Nratios) <- wave$plot_id
#distance
dist <- get_dist(Nratios, method = "pearson")
# fviz_dist(dist)
hcN <- hclust(dist, method = "average")
sub_grp <- cutree(hcN, k = 3)
table(sub_grp)
Nratios <- Nratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Nratios)

Nratios<- as.data.frame(Nratios)
write.csv(Nratios,"./data/Nratios_rep3.csv", row.names=F)
Nratios<- read.csv("./data/Nratios_rep3.csv")

# Compute the absolute difference between wave_1 and wave_2
Nratios <- Nratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Nratio_sel <- Nratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Nratio_transform <- Nratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Nratio_transform, "./data/Nratio_transform_rep3.csv", row.names=F)
Nratio_transform<- read.csv("./data/Nratio_transform_rep3.csv")
com_col <- colnames(wave)[colnames(wave) %in% Nratio_transform$wave_sel]
Nwave_pheno<- wave %>% select(1:10, com_col)
write.csv(Nwave_pheno,"./data/Nwave_pheno_rep3.csv", row.names = F)
Nwave<- read.csv("./data/Nwave_pheno_rep3.csv")
#plotting

l <- Nratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()

heatmap_rep3 <-N |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Coh2-Replictaion3") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep3,"./figure/Heatmap_rep3.jpeg")


#replication-4
N <- fread("./data/narea_breakdown4.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_4)
wave <-phenotypes |>
  clean_names()

sum(is.na(N$corg))
sum(N$corg< -1 | N$corg > 1, na.rm = T)
sum(N$corgblup< -1 | N$corgblup > 1, na.rm = T)

N <- N |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(N$corg))
ntraits <- ceiling(sum(!is.na(N$coh2))*0.01)
sel <- N |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Nratios <- wave |> select(wave_1692_wave_1705:wave_1691_wave_1718)
Nratios <- t(Nratios)
colnames(Nratios) <- wave$plot_id
#distance
dist <- get_dist(Nratios, method = "pearson")
# fviz_dist(dist)
hcN <- hclust(dist, method = "average")
sub_grp <- cutree(hcN, k = 3)
table(sub_grp)
Nratios <- Nratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Nratios)

Nratios<- as.data.frame(Nratios)
write.csv(Nratios,"./data/Nratios_rep4.csv", row.names=F)
Nratios<- read.csv("./data/Nratios_rep4.csv")

# Compute the absolute difference between wave_1 and wave_2
Nratios <- Nratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Nratio_sel <- Nratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Nratio_transform <- Nratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Nratio_transform, "./data/Nratio_transform_rep4.csv", row.names=F)
Nratio_transform<- read.csv("./data/Nratio_transform_rep4.csv")
com_col <- colnames(wave)[colnames(wave) %in% Nratio_transform$wave_sel]
Nwave_pheno<- wave %>% select(1:10, com_col)
write.csv(Nwave_pheno,"./data/Nwave_pheno_rep4.csv", row.names = F)
Nwave<- read.csv("./data/Nwave_pheno_rep4.csv")
#plotting

l <- Nratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()

heatmap_rep4 <-N |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Coh2-Replictaion4") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep4,"./figure/Heatmap_rep4.jpeg")


#replication 5
N <- fread("./data/narea_breakdown5.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_5)
wave <-phenotypes |>
  clean_names()

sum(is.na(N$corg))
sum(N$corg< -1 | N$corg > 1, na.rm = T)
sum(N$corgblup< -1 | N$corgblup > 1, na.rm = T)

N <- N |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(N$corg))
ntraits <- ceiling(sum(!is.na(N$coh2))*0.01)
sel <- N |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Nratios <- wave |> select(wave_1039_wave_390:wave_840_wave_402)
Nratios <- t(Nratios)
colnames(Nratios) <- wave$plot_id
#distance
dist <- get_dist(Nratios, method = "pearson")
# fviz_dist(dist)
hcN <- hclust(dist, method = "average")
sub_grp <- cutree(hcN, k = 3)
table(sub_grp)
Nratios <- Nratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Nratios)

Nratios<- as.data.frame(Nratios)
write.csv(Nratios,"./data/Nratios_rep5.csv", row.names=F)
Nratios<- read.csv("./data/Nratios_rep5.csv")

# Compute the absolute difference between wave_1 and wave_2
Nratios <- Nratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Nratio_sel <- Nratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Nratio_transform <- Nratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Nratio_transform, "./data/Nratio_transform_rep5.csv", row.names=F)
Nratio_transform<- read.csv("./data/Nratio_transform_rep5.csv")
com_col <- colnames(wave)[colnames(wave) %in% Nratio_transform$wave_sel]
Nwave_pheno<- wave %>% select(1:10, com_col)
write.csv(Nwave_pheno,"./data/Nwave_pheno_rep5.csv", row.names = F)
Nwave<- read.csv("./data/Nwave_pheno_rep5.csv")
#plotting

l <- Nratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()

heatmap_rep5 <-N |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Coh2-Replictaion5") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()
ggsave(plot=heatmap_rep5,"./figure/Heatmap_rep5.jpeg")


Nratio_combined<- read.csv("./data/Nratio_combined.csv")
#calculating avg coh2 from 5 reps and plot the heatmaps
N1 <- fread("./data/narea_breakdown1.csv", data.table = F)
N2 <- fread("./data/narea_breakdown2.csv", data.table = F)
N3<- fread("./data/narea_breakdown3.csv", data.table = F)
N4 <- fread("./data/narea_breakdown4.csv", data.table = F)
N5 <- fread("./data/narea_breakdown5.csv", data.table = F)
N1_sub <- N1[c("wave_1", "wave_2", "coh2")]
N2_sub <- N2[c("wave_1", "wave_2", "coh2")]
N3_sub <- N3[c("wave_1", "wave_2", "coh2")]
N4_sub <- N4[c("wave_1", "wave_2", "coh2")]
N5_sub <- N5[c("wave_1", "wave_2", "coh2")]
names(N1_sub)[3] <- "coh2_rep1"
names(N2_sub)[3] <- "coh2_rep2"
names(N3_sub)[3] <- "coh2_rep3"
names(N4_sub)[3] <- "coh2_rep4"
names(N5_sub)[3] <- "coh2_rep5"
reps<- list(N1_sub, N2_sub,N3_sub,N4_sub,N5_sub)
final_reps <- Reduce(function(x, y) merge(x, y, by = c("wave_1", "wave_2")), reps)

# Calculate the average of coh2 values from different reps
final_reps <- final_reps %>%
  mutate(coh2_avg = rowMeans(select(., matches("^coh2_rep")), na.rm = TRUE))

# Define the symbols you want to use
symbols <- c("*", "+", "o", "x", "s")
Nratio_combined <- Nratio_combined %>%
  mutate(label = rep(symbols,each=3, length.out = nrow(Nratio_combined)))
# Process Nratio_combined to prepare it for merging
l <- Nratio_combined %>%
  mutate(wave_1 = as.numeric(gsub("_.*", "", gsub("wave_", "", wave_1))),
         wave_2 = as.numeric(gsub(".*_", "", wave_2))) %>%
  select(wave_1, wave_2, label) %>%
  as.data.frame()

# Process final_reps for plotting
heatmap_rep1 <- final_reps %>%
  mutate(wave_1 = as.numeric(gsub("wave_", "", wave_1)),
         wave_2 = as.numeric(gsub("wave_", "", wave_2))) %>%
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2_avg)) +
  scale_fill_gradient(low = "white", high = "red", na.value = "black") +
  labs(x = "Wave1", y = "Wave2", title = "Avg_coh2_Combined_reps") +
  geom_label(data = l, aes(label = label),size=4) +
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_avgreps.jpeg")




Nratio_combined<- rbind(Nratio_transform_rep1,Nratio_transform_rep2,Nratio_transform_rep3,Nratio_transform_rep4,Nratio_transform_rep5 )
write.csv(Nratio_combined,"./data/Nratio_combined.csv", row.names = F)

#for five replications  first stage analysis
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% clean_names()
Nratio_transform<- read.csv("./data/Nratio_combined.csv")
Nratio_transform<- read.csv("./data/Nratio_transform_rep1.csv")
Nratio_transform<- read.csv("./data/Nratio_transform_rep2.csv")
Nratio_transform<- read.csv("./data/Nratio_transform_rep3.csv")
Nratio_transform<- read.csv("./data/Nratio_transform_rep4.csv")
Nratio_transform<- read.csv("./data/Nratio_transform_rep5.csv")

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
indv_sample<- read.csv("./data/indv_sample.csv")
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
indv_sample<- read.csv("./data/indv_sample.csv")
NratiobluesEF_rep1<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(name2,2:4)
NratiobluesEF_rep2<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2,5:7)
NratiobluesEF_rep3<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2,8:10)
NratiobluesEF_rep4<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2,11:13)
NratiobluesEF_rep5<-Nratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2,14:16)

fwrite( NratiobluesEF_rep5, "./output/NratiobluesEF_rep5.csv", row.names = FALSE)

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

#combining blues for trait and waveratio for EF and MW location for rep1
NbluesEF_rep1<- read.csv("./output/NbluesEF_rep1.csv")
NbluesMW_rep1<- read.csv("./output/NbluesMW_rep1.csv")

Nblues_rep1 <- bind_rows(
  "EF" = NbluesEF_rep1,
  "MW" = NbluesMW_rep1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Nblues_rep1 <-Nblues_rep1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
N_blues_rep1<- subset(Nblues_rep1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
N_blues_rep1<- droplevels(N_blues_rep1[N_blues_rep1$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(N_blues_rep1,"./output/N_blues_rep1.csv", row.names = F)

#combining blues for trait and waveratio for EF and MW location for rep2
NbluesEF_rep2<- read.csv("./output/NbluesEF_rep2.csv")
NbluesMW_rep2<- read.csv("./output/NbluesMW_rep2.csv")

Nblues_rep2 <- bind_rows(
  "EF" = NbluesEF_rep2,
  "MW" = NbluesMW_rep2,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Nblues_rep2 <-Nblues_rep2|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
N_blues_rep2<- subset(Nblues_rep2, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
N_blues_rep2<- droplevels(N_blues_rep2[N_blues_rep2$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(N_blues_rep2,"./output/N_blues_rep2.csv", row.names = F)
# #combining blues for trait and waveratio for EF and MW location for rep2
NbluesEF_rep3<- read.csv("./output/NbluesEF_rep3.csv")
NbluesMW_rep3<- read.csv("./output/NbluesMW_rep3.csv")
Nblues_rep3 <- bind_rows(
   "EF" = NbluesEF_rep3,
   "MW" = NbluesMW_rep3,
   .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
 Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Nblues_rep3 <-Nblues_rep3|>
   left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
   mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
N_blues_rep3<- subset(Nblues_rep3, !is.na(taxa)) %>% select(-name2,-Corrected_names)
N_blues_rep3<- droplevels(N_blues_rep3[N_blues_rep3$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(N_blues_rep3,"./output/N_blues_rep3.csv", row.names = F)

#combining blues for trait and waveratio for EF and MW location for rep4
NbluesEF_rep4<- read.csv("./output/NbluesEF_rep4.csv")
NbluesMW_rep4<- read.csv("./output/NbluesMW_rep4.csv")
Nblues_rep4 <- bind_rows(
  "EF" = NbluesEF_rep4,
  "MW" = NbluesMW_rep4,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Nblues_rep4 <-Nblues_rep4|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
N_blues_rep4<- subset(Nblues_rep4, !is.na(taxa)) %>% select(-name2,-Corrected_names)
N_blues_rep4<- droplevels(N_blues_rep4[N_blues_rep4$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(N_blues_rep4,"./output/N_blues_rep4.csv", row.names = F)
#combining blues for trait and waveratio for EF and MW location for rep5
NbluesEF_rep5<- read.csv("./output/NbluesEF_rep5.csv")
NbluesMW_rep5<- read.csv("./output/NbluesMW_rep5.csv")
Nblues_rep5 <- bind_rows(
  "EF" = NbluesEF_rep5,
  "MW" = NbluesMW_rep5,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Nblues_rep5 <-Nblues_rep5|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
N_blues_rep5<- subset(Nblues_rep5, !is.na(taxa)) %>% select(-name2,-Corrected_names)
N_blues_rep5<- droplevels(N_blues_rep5[N_blues_rep5$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(N_blues_rep5,"./output/N_blues_rep5.csv", row.names = F)


#----------------2nd step --------------------
N_blues<- read.csv("N_blues_rep1.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))
sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
kin<- fread('kin_additive.txt', data.table= F)
rownames(kin) <- colnames(kin)
kin<- as.matrix(kin)

GINV <- readRDS(file = "./data/GINV.rds")
#----- ST Narea for rep1 --------
# ---------------------narea-running cross validation for EFMW loc---------------------
result_N_rep1 <-
  crossv(
    sort = sort,
    train = N_bluesEF,
    validation = N_bluesMW,
    mytrait = "narea",
    kin=kin,
    scheme = "EFMW")

corr_N_EFMW<- data.frame(result_N_rep1$ac)
fwrite(corr_N_EFMW, "./output/corr_N_EFMWrep1.csv")

# ---------------------narea-running cross validation for MWEF loc---------------------
result_N_rep1 <-
  crossv(
    sort = sort,
    train = N_bluesMW,
    validation = N_bluesEF,
    mytrait = "narea",
    kin=kin,
    scheme = "MWEF_rep1")

corr_N_MWEF<- data.frame(result_N_rep1$ac)
fwrite(corr_N_MWEF, "./output/corr_N_MWEFrep1.csv")

#----- ST Narea for rep2 --------
N_blues<- read.csv("N_blues_rep2.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))
sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)

# ---------------------narea-running cross validation for EFMW loc---------------------
result_N_rep2 <-
  crossv(
    sort = sort,
    train = N_bluesEF,
    validation = N_bluesMW,
    mytrait = "narea",
    kin=kin,
    scheme = "EFMW_rep2")

corr_N_EFMW<- data.frame(result_N_rep2$ac)
fwrite(corr_N_EFMW, "./output/corr_N_EFMWrep2.csv")

# ---------------------narea-running cross validation for MWEF loc---------------------
result_N_rep2 <-
  crossv(
    sort = sort,
    train = N_bluesMW,
    validation = N_bluesEF,
    mytrait = "narea",
    kin = kin,
    scheme = "MWEF_rep2")

corr_N_MWEF<- data.frame(result_N_rep2$ac)
fwrite(corr_N_MWEF, "./output/corr_N_MWEFrep2.csv")

#----- ST Narea for rep3 --------
N_blues<- read.csv("N_blues_rep3.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))
sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)

# ---------------------narea-running cross validation for EFMW loc---------------------
result_N_rep3 <-
  crossv(
    sort = sort,
    train = N_bluesEF,
    validation = N_bluesMW,
    mytrait = "narea",
    kin=kin,
    scheme = "EFMW_rep3")

corr_N_EFMW<- data.frame(result_N_rep3$ac)
fwrite(corr_N_EFMW, "./output/corr_N_EFMWrep3.csv")

# ---------------------narea-running cross validation for MWEF loc---------------------
result_N_rep3 <-
  crossv(
    sort = sort,
    train = N_bluesMW,
    validation = N_bluesEF,
    mytrait = "narea",
    kin = kin,
    scheme = "MWEF_rep3")

corr_N_MWEF<- data.frame(result_N_rep3$ac)
fwrite(corr_N_MWEF, "./output/narea/corr_N_MWEFrep3.csv")

#----- ST Narea for rep4 --------
N_blues<- read.csv("N_blues_rep4.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))
sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)

# ---------------------narea-running cross validation for EFMW loc---------------------
result_N_rep4 <-
  crossv(
    sort = sort,
    train = N_bluesEF,
    validation = N_bluesMW,
    mytrait = "narea",
    kin=kin,
    scheme = "EFMW_rep4")

corr_N_EFMW<- data.frame(result_N_rep4$ac)
fwrite(corr_N_EFMW, "./output/corr_N_EFMWrep4.csv")

# ---------------------narea-running cross validation for MWEF loc---------------------
result_N_rep4 <-
  crossv(
    sort = sort,
    train = N_bluesMW,
    validation = N_bluesEF,
    mytrait = "narea",
    kin = kin,
    scheme = "MWEF_rep4")

corr_N_MWEF<- data.frame(result_N_rep4$ac)
fwrite(corr_N_MWEF, "./output/narea/corr_N_MWEFrep4.csv")

#----- ST Narea for rep5 --------
N_blues<- read.csv("N_blues_rep5.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))
sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)

# ---------------------narea-running cross validation for EFMW loc---------------------
result_N_rep5 <-
  crossv(
    sort = sort,
    train = N_bluesEF,
    validation = N_bluesMW,
    mytrait = "narea",
    kin=kin,
    scheme = "EFMW_rep5")

corr_N_EFMW<- data.frame(result_N_rep5$ac)
fwrite(corr_N_EFMW, "./output/corr_N_EFMWrep5.csv")

# ---------------------narea-running cross validation for MWEF loc---------------------
result_N_rep5 <-
  crossv(
    sort = sort,
    train = N_bluesMW,
    validation = N_bluesEF,
    mytrait = "narea",
    kin = kin,
    scheme = "MWEF_rep5")

corr_N_MWEF<- data.frame(result_N_rep5$ac)
fwrite(corr_N_MWEF, "./output/narea/corr_N_MWEFrep5.csv")

#obtaining synthetic trait table and value in each rep
Nratio_combined <- read.csv("./data/Nratio_combined.csv", header = TRUE) %>% select(-3)
Nratio_combined$ST <- rep(c("ST1", "ST2", "ST3"), times = 5)
Nratio_combined$Rep <- rep(1:5, each = 3)
Nratio_combined<- Nratio_combined %>% select(-wave_1, -wave_2)
Nratio_combined$coh2_Rep1<- c(0.63,0.63,0.63,0.63,0.43,0.54,0.62,0.48,0.56,0.63,0.61,NA,0.32,0.54,0.62)
Nratio_combined$coh2_Rep2<- c(0.55,NA,NA,0.66,0.58,0.59,NA,NA,NA,NA,0.52,NA,0.40,NA,NA)
Nratio_combined$coh2_Rep3<- c(0.46,0.56,0.44,0.64,0.47,0.59,0.52,0.56,0.51,0.55,0.47,0.32,0.49,0.43,0.55)
Nratio_combined$coh2_Rep4<- c(0.46,0.51,0.50,0.55,0.35,0.32,0.48,0.32,0.36,0.50,0.49,0.52,0.28,0.22,0.53)
Nratio_combined$coh2_Rep5<- c(0.34,0.46,0.42,0.55,0.40,0.48,0.43,0.42,0.54,0.45,0.39,0.31,0.48,0.48,0.49)

gt_table <- gt(Nratio_combined)
# Apply the yellow highlighting
for (i in 1:5) {
  column_name <- paste0("coh2_Rep", i)
  gt_table <- gt_table %>%
    tab_style(
      style = cell_fill(color = "yellow"),
      locations = cells_body(
        columns = column_name,
        rows = Rep == i
      )
    )
}
htmltools::save_html(gt_table, file = "coh2_all_Narea_rep.html")
webshot("coh2_all_Narea_rep.html", "./figure/coh2_all_Narea_rep.png")


#selecting synthetic trait that has lowest coh2 i.e 0 and using it to fit the model

#Rep1
# Select the one wave ratio with the lowest coh2
N1 <- fread("./data/narea_breakdown1.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_1))
wave <-phenotypes |>
  clean_names()
selected <- N1 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Nwave_lowcoh2<- wave %>% select(1:10, narea,wave_1933,wave_1434)
wave_1933_wave_1434 <- paste("wave_1933", "wave_1434", sep = "_")
Nwave_lowcoh2[[wave_1933_wave_1434]] <- Nwave_lowcoh2[["wave_1933"]] / Nwave_lowcoh2[["wave_1434"]]
Nwave_lowcoh2<- Nwave_lowcoh2 %>% select(-wave_1933,-wave_1434)
write.csv(Nwave_lowcoh2, "./output/Nwave_rep1_lowcoh2.csv", row.names = F)

#Rep2
N2 <- fread("./data/narea_breakdown2.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_2))
wave <-phenotypes |>
  clean_names()

selected <- N2 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Nwave_lowcoh2<- wave %>% select(1:10, narea,wave_1878,wave_2161)
wave_1878_wave_2161 <- paste("wave_1878", "wave_2161", sep = "_")
Nwave_lowcoh2[[wave_1878_wave_2161]] <- Nwave_lowcoh2[["wave_1878"]] / Nwave_lowcoh2[["wave_2161"]]
Nwave_lowcoh2<- Nwave_lowcoh2 %>% select(-wave_1878,-wave_2161)
write.csv(Nwave_lowcoh2, "./output/Nwave_rep2_lowcoh2.csv", row.names = F)

#Rep3
N3 <- fread("./data/narea_breakdown3.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_3))
wave <-phenotypes |>
  clean_names()

selected <- N3 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Nwave_lowcoh2<- wave %>% select(1:10, narea,wave_1968,wave_2380)
wave_1968_wave_2380 <- paste("wave_1968", "wave_2380", sep = "_")
Nwave_lowcoh2[[wave_1968_wave_2380]] <- Nwave_lowcoh2[["wave_1968"]] / Nwave_lowcoh2[["wave_2380"]]
Nwave_lowcoh2<- Nwave_lowcoh2 %>% select(-wave_1968,-wave_2380)
write.csv(Nwave_lowcoh2, "./output/Nwave_rep3_lowcoh2.csv", row.names = F)


#Rep4
N4 <- fread("./data/narea_breakdown4.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_4))
wave <-phenotypes |>
  clean_names()

selected <- N4 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Nwave_lowcoh2<- wave %>% select(1:10, narea,wave_1377,wave_571)
wave_1377_wave_571 <- paste("wave_1377", "wave_571", sep = "_")
Nwave_lowcoh2[[wave_1377_wave_571]] <- Nwave_lowcoh2[["wave_1377"]] / Nwave_lowcoh2[["wave_571"]]
Nwave_lowcoh2<- Nwave_lowcoh2 %>% select(-wave_1377,-wave_571)
write.csv(Nwave_lowcoh2, "./output/Nwave_rep4_lowcoh2.csv", row.names = F)

#Rep5
N5 <- fread("./data/narea_breakdown5.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_5))
wave <-phenotypes |>
  clean_names()

selected <- N5 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Nwave_lowcoh2<- wave %>% select(1:10, narea,wave_2109,wave_537)
wave_2109_wave_537 <- paste("wave_2109", "wave_537", sep = "_")
Nwave_lowcoh2[[wave_2109_wave_537]] <- Nwave_lowcoh2[["wave_2109"]] / Nwave_lowcoh2[["wave_537"]]
Nwave_lowcoh2<- Nwave_lowcoh2 %>% select(-wave_2109,-wave_537)
write.csv(Nwave_lowcoh2, "./output/Nwave_rep5_lowcoh2.csv", row.names = F)



#calculating blues for the selected low coheritable synthetic trait
Nwave_lowcoh2<- read.csv("./output/Nwave_rep1_lowcoh2.csv")
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
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesMW1 <- bind_rows(Nratio_bluesMW1, temp)
Nratio_bluesMW1<- Nratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesMW1, "./output/Nratio_bluesMW1_rep1.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Nwave_lowcoh2)[11]
models <- vector("list",length(variables))
Nratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesEF1 <- bind_rows(Nratio_bluesEF1, temp)

Nratio_bluesEF1<- Nratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesEF1, "./output/Nratio_bluesEF1_rep1.csv", row.names = FALSE)
Nratio_bluesEF1<- read.csv("./output/Nratio_bluesEF1_rep1.csv")

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
fwrite( NbluesEF1, "./output/NbluesEF1_rep1.csv", row.names = FALSE)
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
fwrite( NbluesMW1, "./output/NbluesMW1_rep1.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
NbluesEF1<- read.csv("./output/NbluesEF1_rep1.csv")
NbluesMW1<- read.csv("./output/NbluesMW1_rep1.csv")

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
write.csv(N_blues1,"./output/N_blues1_rep1.csv", row.names = F)


##calculating blues for the selected low coheritable synthetic trait
Nwave_lowcoh2<- read.csv("./output/Nwave_rep2_lowcoh2.csv")
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
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesMW1 <- bind_rows(Nratio_bluesMW1, temp)
Nratio_bluesMW1<- Nratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesMW1, "./output/Nratio_bluesMW1_rep2.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Nwave_lowcoh2)[11]
models <- vector("list",length(variables))
Nratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesEF1 <- bind_rows(Nratio_bluesEF1, temp)

Nratio_bluesEF1<- Nratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesEF1, "./output/Nratio_bluesEF1_rep2.csv", row.names = FALSE)
Nratio_bluesEF1<- read.csv("./output/Nratio_bluesEF1_rep2.csv")

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
fwrite( NbluesEF1, "./output/NbluesEF1_rep2.csv", row.names = FALSE)
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
fwrite( NbluesMW1, "./output/NbluesMW1_rep2.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
NbluesEF1<- read.csv("./output/NbluesEF1_rep2.csv")
NbluesMW1<- read.csv("./output/NbluesMW1_rep2.csv")

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
write.csv(N_blues1,"./output/N_blues1_rep2.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
Nwave_lowcoh2<- read.csv("./output/Nwave_rep3_lowcoh2.csv")
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
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesMW1 <- bind_rows(Nratio_bluesMW1, temp)
Nratio_bluesMW1<- Nratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesMW1, "./output/Nratio_bluesMW1_rep3.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Nwave_lowcoh2)[11]
models <- vector("list",length(variables))
Nratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesEF1 <- bind_rows(Nratio_bluesEF1, temp)

Nratio_bluesEF1<- Nratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesEF1, "./output/Nratio_bluesEF1_rep3.csv", row.names = FALSE)
Nratio_bluesEF1<- read.csv("./output/Nratio_bluesEF1_rep3.csv")

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
fwrite( NbluesEF1, "./output/NbluesEF1_rep3.csv", row.names = FALSE)
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
fwrite( NbluesMW1, "./output/NbluesMW1_rep3.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
NbluesEF1<- read.csv("./output/NbluesEF1_rep3.csv")
NbluesMW1<- read.csv("./output/NbluesMW1_rep3.csv")

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
write.csv(N_blues1,"./output/N_blues1_rep3.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
Nwave_lowcoh2<- read.csv("./output/Nwave_rep4_lowcoh2.csv")
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
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesMW1 <- bind_rows(Nratio_bluesMW1, temp)
Nratio_bluesMW1<- Nratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesMW1, "./output/Nratio_bluesMW1_rep4.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Nwave_lowcoh2)[11]
models <- vector("list",length(variables))
Nratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesEF1 <- bind_rows(Nratio_bluesEF1, temp)

Nratio_bluesEF1<- Nratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesEF1, "./output/Nratio_bluesEF1_rep4.csv", row.names = FALSE)
Nratio_bluesEF1<- read.csv("./output/Nratio_bluesEF1_rep4.csv")

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
fwrite( NbluesEF1, "./output/NbluesEF1_rep4.csv", row.names = FALSE)
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
fwrite( NbluesMW1, "./output/NbluesMW1_rep4.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
NbluesEF1<- read.csv("./output/NbluesEF1_rep4.csv")
NbluesMW1<- read.csv("./output/NbluesMW1_rep4.csv")

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
write.csv(N_blues1,"./output/N_blues1_rep4.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
Nwave_lowcoh2<- read.csv("./output/Nwave_rep5_lowcoh2.csv")
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
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesMW1 <- bind_rows(Nratio_bluesMW1, temp)
Nratio_bluesMW1<- Nratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesMW1, "./output/Nratio_bluesMW1_rep5.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Nwave_lowcoh2)[11]
models <- vector("list",length(variables))
Nratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Nwave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Nratio_bluesEF1 <- bind_rows(Nratio_bluesEF1, temp)

Nratio_bluesEF1<- Nratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Nratio_bluesEF1, "./output/Nratio_bluesEF1_rep5.csv", row.names = FALSE)
Nratio_bluesEF1<- read.csv("./output/Nratio_bluesEF1_rep5.csv")

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
fwrite( NbluesEF1, "./output/NbluesEF1_rep5.csv", row.names = FALSE)
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
fwrite( NbluesMW1, "./output/NbluesMW1_rep5.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
NbluesEF1<- read.csv("./output/NbluesEF1_rep5.csv")
NbluesMW1<- read.csv("./output/NbluesMW1_rep5.csv")

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
write.csv(N_blues1,"./output/N_blues1_rep5.csv", row.names = F)
