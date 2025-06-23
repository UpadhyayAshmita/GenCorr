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


#replication 1
S <- fread("./data/sla_breakdown1.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_1)

wave <-phenotypes |>
  clean_names()

sum(is.na(S$corg))
sum(S$corg< -1 | S$corg > 1, na.rm = T)
sum(S$corgblup< -1 | S$corgblup > 1, na.rm = T)

S <- S |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(S$corg))
ntraits <- ceiling(sum(!is.na(S$coh2))*0.01)
sel <- S |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Sratios <- wave |> select(wave_1720_wave_1362:wave_1656_wave_1304)
Sratios <- t(Sratios)
colnames(Sratios) <- wave$plot_id
#distance
dist <- get_dist(Sratios, method = "pearson")
# fviz_dist(dist)
hcS <- hclust(dist, method = "average")
sub_grp <- cutree(hcS, k = 3)
table(sub_grp)
Sratios <- Sratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Sratios)

Sratios<- as.data.frame(Sratios)
write.csv(Sratios,"./data/Sratios_rep1.csv", row.names=F)
Sratios<- read.csv("./data/Sratios_rep1.csv")

# Compute the absolute difference between wave_1 and wave_2
Sratios <- Sratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Sratio_sel <- Sratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Sratio_transform <- Sratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Sratio_transform, "./data/Sratio_transform_rep1.csv", row.names=F)
Sratio_transform<- read.csv("./data/Sratio_transform_rep1.csv")
com_col <- colnames(wave)[colnames(wave) %in% Sratio_transform$wave_sel]
Swave_pheno<- wave %>% select(1:9, sla,com_col)
write.csv(Swave_pheno,"./data/Swave_pheno_rep1.csv", row.names = F)
Swave<- read.csv("./data/Swave_pheno_rep1.csv")
#plotting

l <- Sratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep1 <-S |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Sla_Coh2-Replictaion1") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_rep1_sla.jpeg")

S1<- S
#replication 2
S <- fread("./data/sla_breakdown2.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_2)

wave <-phenotypes |>
  clean_names()

sum(is.na(S$corg))
sum(S$corg< -1 | S$corg > 1, na.rm = T)
sum(S$corgblup< -1 | S$corgblup > 1, na.rm = T)

S <- S |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(S$corg))
ntraits <- ceiling(sum(!is.na(S$coh2))*0.01)
sel <- S |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Sratios <- wave |> select(wave_1717_wave_1682:wave_1736_wave_1191)
Sratios <- t(Sratios)
colnames(Sratios) <- wave$plot_id
#distance
dist <- get_dist(Sratios, method = "pearson")
# fviz_dist(dist)
hcS <- hclust(dist, method = "average")
sub_grp <- cutree(hcS, k = 3)
table(sub_grp)
Sratios <- Sratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Sratios)

Sratios<- as.data.frame(Sratios)
write.csv(Sratios,"./data/Sratios_rep2.csv", row.names=F)
Sratios<- read.csv("./data/Sratios_rep2.csv")

# Compute the absolute difference between wave_1 and wave_2
Sratios <- Sratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Sratio_sel <- Sratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Sratio_transform <- Sratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Sratio_transform, "./data/Sratio_transform_rep2.csv", row.names=F)
Sratio_transform<- read.csv("./data/Sratio_transform_rep2.csv")
com_col <- colnames(wave)[colnames(wave) %in% Sratio_transform$wave_sel]
Swave_pheno<- wave %>% select(1:9, sla,com_col)
write.csv(Swave_pheno,"./data/Swave_pheno_rep2.csv", row.names = F)
Swave<- read.csv("./data/Swave_pheno_rep2.csv")
#plotting

l <- Sratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep1 <-S |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Sla_Coh2-Replictaion2") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_rep2_sla.jpeg")
Sratio_combined<- read.csv("./data/Sratio_combined.csv")

#replication 3
S <- fread("./data/sla_breakdown3.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_3)

wave <-phenotypes |>
  clean_names()

sum(is.na(S$corg))
sum(S$corg< -1 | S$corg > 1, na.rm = T)
sum(S$corgblup< -1 | S$corgblup > 1, na.rm = T)

S <- S |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(S$corg))
ntraits <- ceiling(sum(!is.na(S$coh2))*0.01)
sel <- S |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Sratios <- wave |> select(wave_751_wave_805:wave_404_wave_1068)
Sratios <- t(Sratios)
colnames(Sratios) <- wave$plot_id
#distance
dist <- get_dist(Sratios, method = "pearson")
# fviz_dist(dist)
hcS <- hclust(dist, method = "average")
sub_grp <- cutree(hcS, k = 3)
table(sub_grp)
Sratios <- Sratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Sratios)

Sratios<- as.data.frame(Sratios)
write.csv(Sratios,"./data/Sratios_rep3.csv", row.names=F)
Sratios<- read.csv("./data/Sratios_rep3.csv")

# Compute the absolute difference between wave_1 and wave_2
Sratios <- Sratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Sratio_sel <- Sratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Sratio_transform <- Sratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Sratio_transform, "./data/Sratio_transform_rep3.csv", row.names=F)
Sratio_transform<- read.csv("./data/Sratio_transform_rep3.csv")
com_col <- colnames(wave)[colnames(wave) %in% Sratio_transform$wave_sel]
Swave_pheno<- wave %>% select(1:9, sla,com_col)
write.csv(Swave_pheno,"./data/Swave_pheno_rep3.csv", row.names = F)
Swave<- read.csv("./data/Swave_pheno_rep3.csv")
#plotting

l <- Sratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep1 <-S |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Sla_Coh2-Replictaion3") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_rep3_sla.jpeg")


#replication4
S <- fread("./data/sla_breakdown4.csv", data.table = F)

phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_4)

wave <-phenotypes |>
  clean_names()

sum(is.na(S$corg))
sum(S$corg< -1 | S$corg > 1, na.rm = T)
sum(S$corgblup< -1 | S$corgblup > 1, na.rm = T)

S <- S |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(S$corg))
ntraits <- ceiling(sum(!is.na(S$coh2))*0.01)
sel <- S |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Sratios <- wave |> select(wave_752_wave_754:wave_380_wave_874)
Sratios <- t(Sratios)
colnames(Sratios) <- wave$plot_id
#distance
dist <- get_dist(Sratios, method = "pearson")
# fviz_dist(dist)
hcS <- hclust(dist, method = "average")
sub_grp <- cutree(hcS, k = 3)
table(sub_grp)
Sratios <- Sratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Sratios)

Sratios<- as.data.frame(Sratios)
write.csv(Sratios,"./data/Sratios_rep4.csv", row.names=F)
Sratios<- read.csv("./data/Sratios_rep4.csv")

# Compute the absolute difference between wave_1 and wave_2
Sratios <- Sratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Sratio_sel <- Sratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Sratio_transform <- Sratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Sratio_transform, "./data/Sratio_transform_rep4.csv", row.names=F)
Sratio_transform<- read.csv("./data/Sratio_transform_rep4.csv")
com_col <- colnames(wave)[colnames(wave) %in% Sratio_transform$wave_sel]
Swave_pheno<- wave %>% select(1:9, sla,com_col)
write.csv(Swave_pheno,"./data/Swave_pheno_rep4.csv", row.names = F)
Swave<- read.csv("./data/Swave_pheno_rep4.csv")
#plotting

l <- Sratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep1 <-S |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Sla_Coh2-Replictaion4") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_rep4_sla.jpeg")

#replication 5
S <- fread("./data/sla_breakdown5.csv", data.table = F)

phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_5)

wave <-phenotypes |>
  clean_names()

sum(is.na(S$corg))
sum(S$corg< -1 | S$corg > 1, na.rm = T)
sum(S$corgblup< -1 | S$corgblup > 1, na.rm = T)

S <- S |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(S$corg))
ntraits <- ceiling(sum(!is.na(S$coh2))*0.01)
sel <- S |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Sratios <- wave |> select(wave_752_wave_754:wave_380_wave_874)
Sratios <- t(Sratios)
colnames(Sratios) <- wave$plot_id
#distance
dist <- get_dist(Sratios, method = "pearson")
# fviz_dist(dist)
hcS <- hclust(dist, method = "average")
sub_grp <- cutree(hcS, k = 3)
table(sub_grp)
Sratios <- Sratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(Sratios)

Sratios<- as.data.frame(Sratios)
write.csv(Sratios,"./data/Sratios_rep5.csv", row.names=F)
Sratios<- read.csv("./data/Sratios_rep5.csv")

# Compute the absolute difference between wave_1 and wave_2
Sratios <- Sratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Sratio_sel <- Sratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Sratio_transform <- Sratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Sratio_transform, "./data/Sratio_transform_rep5.csv", row.names=F)
Sratio_transform<- read.csv("./data/Sratio_transform_rep5.csv")
com_col <- colnames(wave)[colnames(wave) %in% Sratio_transform$wave_sel]
Swave_pheno<- wave %>% select(1:9, sla,com_col)
write.csv(Swave_pheno,"./data/Swave_pheno_rep5.csv", row.names = F)
Swave<- read.csv("./data/Swave_pheno_rep5.csv")
#plotting

l <- Sratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep1 <-S |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Sla_Coh2-Replictaion5") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_rep5_sla.jpeg")

Sratio_transform_rep1<-read.csv("./data/Sratio_transform_rep1.csv")
Sratio_transform_rep2<-read.csv("./data/Sratio_transform_rep2.csv")
Sratio_transform_rep3<-read.csv("./data/Sratio_transform_rep3.csv")
Sratio_transform_rep4<-read.csv("./data/Sratio_transform_rep4.csv")
Sratio_transform_rep5<-read.csv("./data/Sratio_transform_rep5.csv")
Sratio_combined<- rbind(Sratio_transform_rep1,Sratio_transform_rep2,Sratio_transform_rep3,Sratio_transform_rep4,Sratio_transform_rep5 )
write.csv(Sratio_combined,"./data/Sratio_combined.csv", row.names = F)
Sratio_combined<- read.csv("./data/Sratio_combined.csv")

#calculating avg coh2 from 5 reps and plot the heatmaps
S1 <- fread("./data/sla_breakdown1.csv", data.table = F)
S2 <- fread("./data/sla_breakdown2.csv", data.table = F)
S3<- fread("./data/sla_breakdown3.csv", data.table = F)
S4 <- fread("./data/sla_breakdown4.csv", data.table = F)
S5 <- fread("./data/sla_breakdown5.csv", data.table = F)
S1_sub <- S1[c("wave_1", "wave_2", "coh2")]
S2_sub <- S2[c("wave_1", "wave_2", "coh2")]
S3_sub <- S3[c("wave_1", "wave_2", "coh2")]
S4_sub <- S4[c("wave_1", "wave_2", "coh2")]
S5_sub <- S5[c("wave_1", "wave_2", "coh2")]
names(S1_sub)[3] <- "coh2_rep1"
names(S2_sub)[3] <- "coh2_rep2"
names(S3_sub)[3] <- "coh2_rep3"
names(S4_sub)[3] <- "coh2_rep4"
names(S5_sub)[3] <- "coh2_rep5"
reps<- list(S1_sub, S2_sub,S3_sub,S4_sub,S5_sub)
final_reps <- Reduce(function(x, y) merge(x, y, by = c("wave_1", "wave_2")), reps)

# Calculate the average of coh2 values from different reps
final_reps <- final_reps %>%
  mutate(coh2_avg = rowMeans(select(., matches("^coh2_rep")), na.rm = TRUE))

# Define the symbols you want to use
symbols <- c("*", "+", "o", "x", "s")
Sratio_combined <- Sratio_combined %>%
  mutate(label = rep(symbols,each=3, length.out = nrow(Sratio_combined)))
# Process Sratio_combined to prepare it for merging
l <- Sratio_combined %>%
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
  labs(x = "Wave1", y = "Wave2", title = "Avg_coh2_Combinedsla_reps") +
  geom_label(data = l, aes(label = label),size=4) +
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmapsla_avgreps.jpeg")
#for five replications  first stage analysis
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% clean_names()
Sratio_transform<- read.csv("./data/Sratio_combined.csv")

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
indv_sample<- read.csv("./data/indv_sample.csv")
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
indv_sample<- read.csv("./data/indv_sample.csv")
SratiobluesEF_rep1<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(name2,2:4)
SratiobluesEF_rep2<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2,5:7)
SratiobluesEF_rep3<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2,8:10)
SratiobluesEF_rep4<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2,11:13)
SratiobluesEF_rep5<-Sratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2,14:16)

fwrite( SratiobluesEF_rep4, "./output/SratiobluesEF_rep4.csv", row.names = FALSE)

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

#combining blues for trait and waveratio for EF and MW location for rep1
SbluesEF_rep1<- read.csv("./output/SbluesEF_rep1.csv")
SbluesMW_rep1<- read.csv("./output/SbluesMW_rep1.csv")

Sblues_rep1 <- bind_rows(
  "EF" = SbluesEF_rep1,
  "MW" = SbluesMW_rep1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Sblues_rep1 <-Sblues_rep1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
S_blues_rep1<- subset(Sblues_rep1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
S_blues_rep1<- droplevels(S_blues_rep1[S_blues_rep1$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(S_blues_rep1,"./output/S_blues_rep1.csv", row.names = F)

#combining blues for trait and waveratio for EF and MW location for rep2
SbluesEF_rep2<- read.csv("./output/SbluesEF_rep2.csv")
SbluesMW_rep2<- read.csv("./output/SbluesMW_rep2.csv")

Sblues_rep2 <- bind_rows(
  "EF" = SbluesEF_rep2,
  "MW" = SbluesMW_rep2,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Sblues_rep2 <-Sblues_rep2|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
S_blues_rep2<- subset(Sblues_rep2, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
S_blues_rep2<- droplevels(S_blues_rep2[S_blues_rep2$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(S_blues_rep2,"./output/S_blues_rep2.csv", row.names = F)
# #combining blues for trait and waveratio for EF and MW location for rep2
SbluesEF_rep3<- read.csv("./output/SbluesEF_rep3.csv")
SbluesMW_rep3<- read.csv("./output/SbluesMW_rep3.csv")
Sblues_rep3 <- bind_rows(
  "EF" = SbluesEF_rep3,
  "MW" = SbluesMW_rep3,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Sblues_rep3 <-Sblues_rep3|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
S_blues_rep3<- subset(Sblues_rep3, !is.na(taxa)) %>% select(-name2,-Corrected_names)
S_blues_rep3<- droplevels(S_blues_rep3[S_blues_rep3$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(S_blues_rep3,"./output/S_blues_rep3.csv", row.names = F)

#combining blues for trait and waveratio for EF and MW location for rep4
SbluesEF_rep4<- read.csv("./output/SbluesEF_rep4.csv")
SbluesMW_rep4<- read.csv("./output/SbluesMW_rep4.csv")
Sblues_rep4 <- bind_rows(
  "EF" = SbluesEF_rep4,
  "MW" = SbluesMW_rep4,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Sblues_rep4 <-Sblues_rep4|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
S_blues_rep4<- subset(Sblues_rep4, !is.na(taxa)) %>% select(-name2,-Corrected_names)
S_blues_rep4<- droplevels(S_blues_rep4[S_blues_rep4$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(S_blues_rep4,"./output/S_blues_rep4.csv", row.names = F)
#combining blues for trait and waveratio for EF and MW location for rep5
SbluesEF_rep5<- read.csv("./output/SbluesEF_rep5.csv")
SbluesMW_rep5<- read.csv("./output/SbluesMW_rep5.csv")
Sblues_rep5 <- bind_rows(
  "EF" = SbluesEF_rep5,
  "MW" = SbluesMW_rep5,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
Sblues_rep5 <-Sblues_rep5|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
S_blues_rep5<- subset(Sblues_rep5, !is.na(taxa)) %>% select(-name2,-Corrected_names)
S_blues_rep5<- droplevels(S_blues_rep5[S_blues_rep5$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(S_blues_rep5,"./output/S_blues_rep5.csv", row.names = F)


#----------------2nd step --------------------
S_blues<- read.csv("S_blues_rep1.csv")
S_bluesEF<- S_blues %>% filter(env== "EF")
S_bluesMW<- S_blues %>% filter(env== "MW")
S_bluesEF<- S_bluesEF%>% mutate(taxa= factor(taxa))
S_bluesMW<- S_bluesMW%>% mutate(taxa= factor(taxa))
N_blues<- read.csv("N_blues_rep1.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
kin<- fread('kin_additive.txt', data.table= F)
rownames(kin) <- colnames(kin)
kin<- as.matrix(kin)

GINV <- readRDS(file = "./data/GINV.rds")
#----- ST Narea for rep1 --------
# ---------------------sla-running cross validation for EFMW loc---------------------
result_S_rep1 <-
  crossv(
    sort = sort,
    train = S_bluesEF,
    validation = S_bluesMW,
    mytrait = "sla",
    kin=kin,
    scheme = "EFMW")

corr_S_EFMW<- data.frame(result_S_rep1$ac)
fwrite(corr_S_EFMW, "./output/corr_S_EFMWrep1.csv")

# ---------------------sla-running cross validation for MWEF loc---------------------
result_S_rep1 <-
  crossv(
    sort = sort,
    train = S_bluesMW,
    validation = S_bluesEF,
    mytrait = "sla",
    kin=kin,
    scheme = "MWEF_rep1")

corr_S_MWEF<- data.frame(result_S_rep1$ac)
fwrite(corr_S_MWEF, "./output/corr_S_MWEFrep1.csv")

#----- ST Narea for rep2 --------
S_blues<- read.csv("S_blues_rep2.csv")
S_bluesEF<- S_blues %>% filter(env== "EF")
S_bluesMW<- S_blues %>% filter(env== "MW")
S_bluesEF<- S_bluesEF%>% mutate(taxa= factor(taxa))
S_bluesMW<- S_bluesMW%>% mutate(taxa= factor(taxa))
N_blues<- read.csv("N_blues_rep2.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)

# ---------------------sla-running cross validation for EFMW loc---------------------
result_S_rep2 <-
  crossv(
    sort = sort,
    train = S_bluesEF,
    validation = S_bluesMW,
    mytrait = "sla",
    kin=kin,
    scheme = "EFMW_rep2")

corr_S_EFMW<- data.frame(result_S_rep2$ac)
fwrite(corr_S_EFMW, "./output/corr_S_EFMWrep2.csv")

# ---------------------sla-running cross validation for MWEF loc---------------------
result_S_rep2 <-
  crossv(
    sort = sort,
    train = S_bluesMW,
    validation = S_bluesEF,
    mytrait = "sla",
    kin = kin,
    scheme = "MWEF_rep2")

corr_S_MWEF<- data.frame(result_S_rep2$ac)
fwrite(corr_S_MWEF, "./output/corr_S_MWEFrep2.csv")

#----- ST Narea for rep3 --------
S_blues<- read.csv("S_blues_rep3.csv")
S_bluesEF<- S_blues %>% filter(env== "EF")
S_bluesMW<- S_blues %>% filter(env== "MW")
S_bluesEF<- S_bluesEF%>% mutate(taxa= factor(taxa))
S_bluesMW<- S_bluesMW%>% mutate(taxa= factor(taxa))

N_blues<- read.csv("N_blues_rep3.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------sla-running cross validation for EFMW loc---------------------
result_S_rep3 <-
  crossv(
    sort = sort,
    train = S_bluesEF,
    validation = S_bluesMW,
    mytrait = "sla",
    kin=kin,
    scheme = "EFMW_rep3")

corr_S_EFMW<- data.frame(result_S_rep3$ac)
fwrite(corr_S_EFMW, "./output/corr_S_EFMWrep3.csv")

# ---------------------sla-running cross validation for MWEF loc---------------------
result_S_rep3 <-
  crossv(
    sort = sort,
    train = S_bluesMW,
    validation = S_bluesEF,
    mytrait = "sla",
    kin = kin,
    scheme = "MWEF_rep3")

corr_S_MWEF<- data.frame(result_S_rep3$ac)
fwrite(corr_S_MWEF, "./output/sla/corr_S_MWEFrep3.csv")

#----- ST Narea for rep4 --------
S_blues<- read.csv("S_blues_rep4.csv")
S_bluesEF<- S_blues %>% filter(env== "EF")
S_bluesMW<- S_blues %>% filter(env== "MW")
S_bluesEF<- S_bluesEF%>% mutate(taxa= factor(taxa))
S_bluesMW<- S_bluesMW%>% mutate(taxa= factor(taxa))
N_blues<- read.csv("N_blues_rep4.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------sla-running cross validation for EFMW loc---------------------
result_S_rep4 <-
  crossv(
    sort = sort,
    train = S_bluesEF,
    validation = S_bluesMW,
    mytrait = "sla",
    kin=kin,
    scheme = "EFMW_rep4")

corr_S_EFMW<- data.frame(result_S_rep4$ac)
fwrite(corr_S_EFMW, "./output/corr_S_EFMWrep4.csv")

# ---------------------sla-running cross validation for MWEF loc---------------------
result_S_rep4 <-
  crossv(
    sort = sort,
    train = S_bluesMW,
    validation = S_bluesEF,
    mytrait = "sla",
    kin = kin,
    scheme = "MWEF_rep4")

corr_S_MWEF<- data.frame(result_S_rep4$ac)
fwrite(corr_S_MWEF, "./output/sla/corr_S_MWEFrep4.csv")

#----- ST Narea for rep5 --------
S_blues<- read.csv("S_blues_rep5.csv")
S_bluesEF<- S_blues %>% filter(env== "EF")
S_bluesMW<- S_blues %>% filter(env== "MW")
S_bluesEF<- S_bluesEF%>% mutate(taxa= factor(taxa))
S_bluesMW<- S_bluesMW%>% mutate(taxa= factor(taxa))

N_blues<- read.csv("N_blues_rep5.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------sla-running cross validation for EFMW loc---------------------
result_S_rep5 <-
  crossv(
    sort = sort,
    train = S_bluesEF,
    validation = S_bluesMW,
    mytrait = "sla",
    kin=kin,
    scheme = "EFMW_rep5")

corr_S_EFMW<- data.frame(result_S_rep5$ac)
fwrite(corr_S_EFMW, "./output/corr_S_EFMWrep5.csv")

# ---------------------sla-running cross validation for MWEF loc---------------------
result_S_rep5 <-
  crossv(
    sort = sort,
    train = S_bluesMW,
    validation = S_bluesEF,
    mytrait = "sla",
    kin = kin,
    scheme = "MWEF_rep5")

corr_S_MWEF<- data.frame(result_S_rep5$ac)
fwrite(corr_S_MWEF, "./output/sla/corr_S_MWEFrep5.csv")

#obtaining synthetic trait table and value in each rep
Sratio_combined <- read.csv("./data/Sratio_combined.csv", header = TRUE) %>% select(-3)
Sratio_combined$ST <- rep(c("ST1", "ST2", "ST3"), times = 5)
Sratio_combined$Rep <- rep(1:5, each = 3)
Sratio_combined<- Sratio_combined %>% select(-wave_1, -wave_2)
Sratio_combined$coh2_Rep1<- c(0.63,0.65,0.63,0.61,0.55,0.49,0.55,0.62, 0.52,0.63,0.64,0.60,0.55,0.62,0.33)
Sratio_combined$coh2_Rep2<- c(0.65,0.64,0.57,0.64,0.63,0.60,0.50,0.61,0.54,0.60,NA,NA,0.56,NA,0.21)
Sratio_combined$coh2_Rep3<-c(0.43,0.45,0.40,0.42,0.53,0.46,0.51,0.57,0.48,0.43,0.39,0.50,0.52,0.58,0.26)
Sratio_combined$coh2_Rep4<- c(0.49,0.53,0.44,0.49,0.39,0.37,0.38,0.49,0.36,0.49,0.49,0.49,0.38,0.53,NA)
Sratio_combined$coh2_Rep5<- c(0.49,0.40,0.76,0.46,NA,0.48,0.45,0.34,0.45,0.34,NA,NA,0.59,0.53,0.52)
library(gtable)
gt_table <- gt(Sratio_combined)
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
htmltools::save_html(gt_table, file = "coh2_all_sla_rep.html")
webshot("coh2_all_sla_rep.html", "./figure/coh2_all_sla_rep.png")


#selecting synthetic trait that has 0 coh2 value from each replication
#Rep1
# Select the one wave ratio with the lowest coh2
S1 <- fread("./data/sla_breakdown1.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_1))
wave <-phenotypes |>
  clean_names()
selected <- S1 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Swave_lowcoh2<- wave %>% select(1:9, sla,wave_367,wave_676)
wave_367_wave_676 <- paste("wave_367", "wave_676", sep = "_")
Swave_lowcoh2[[wave_367_wave_676]] <- Swave_lowcoh2[["wave_367"]] / Swave_lowcoh2[["wave_676"]]
Swave_lowcoh2<- Swave_lowcoh2 %>% select(-wave_367,-wave_676)
write.csv(Swave_lowcoh2, "./output/Swave_rep1_lowcoh2.csv", row.names = F)

#Rep2
S2 <- fread("./data/sla_breakdown2.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_2))
wave <-phenotypes |>
  clean_names()

selected <- S2 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Swave_lowcoh2<- wave %>% select(1:9, sla,wave_618,wave_1713)
wave_618_wave_1713 <- paste("wave_618", "wave_1713", sep = "_")
Swave_lowcoh2[[wave_618_wave_1713]] <- Swave_lowcoh2[["wave_618"]] / Swave_lowcoh2[["wave_1713"]]
Swave_lowcoh2<- Swave_lowcoh2 %>% select(-wave_618,-wave_1713)
write.csv(Swave_lowcoh2, "./output/Swave_rep2_lowcoh2.csv", row.names = F)

#Rep3
S3 <- fread("./data/sla_breakdown3.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_3))
wave <-phenotypes |>
  clean_names()

selected <- S3 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Swave_lowcoh2<- wave %>% select(1:9, sla,wave_717,wave_490)
wave_717_wave_490 <- paste("wave_717", "wave_490", sep = "_")
Swave_lowcoh2[[wave_717_wave_490]] <- Swave_lowcoh2[["wave_717"]] / Swave_lowcoh2[["wave_490"]]
Swave_lowcoh2<- Swave_lowcoh2 %>% select(-wave_717,-wave_490)
write.csv(Swave_lowcoh2, "./output/Swave_rep3_lowcoh2.csv", row.names = F)


#Rep4
S4 <- fread("./data/sla_breakdown4.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_4))
wave <-phenotypes |>
  clean_names()

selected <- S4 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Swave_lowcoh2<- wave %>% select(1:9, sla,wave_1954,wave_2292)
wave_1954_wave_2292 <- paste("wave_1954", "wave_2292", sep = "_")
Swave_lowcoh2[[wave_1954_wave_2292]] <- Swave_lowcoh2[["wave_1954"]] / Swave_lowcoh2[["wave_2292"]]
Swave_lowcoh2<- Swave_lowcoh2 %>% select(-wave_1954,-wave_2292)
write.csv(Swave_lowcoh2, "./output/Swave_rep4_lowcoh2.csv", row.names = F)

#Rep5
S5 <- fread("./data/sla_breakdown5.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_5))
wave <-phenotypes |>
  clean_names()

selected <- S5 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
Swave_lowcoh2<- wave %>% select(1:10, sla,wave_915,wave_1137)
wave_915_wave_1137 <- paste("wave_915", "wave_1137", sep = "_")
Swave_lowcoh2[[wave_915_wave_1137]] <- Swave_lowcoh2[["wave_915"]] / Swave_lowcoh2[["wave_1137"]]
Swave_lowcoh2<- Swave_lowcoh2 %>% select(-wave_915,-wave_1137)
write.csv(Swave_lowcoh2, "./output/Swave_rep5_lowcoh2.csv", row.names = F)



#calculating blues for the selected low coheritable synthetic trait
Swave_lowcoh2<- read.csv("./output/Swave_rep1_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesMW1 <- bind_rows(Sratio_bluesMW1, temp)
Sratio_bluesMW1<- Sratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesMW1, "./output/Sratio_bluesMW1_rep1.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesEF1 <- bind_rows(Sratio_bluesEF1, temp)

Sratio_bluesEF1<- Sratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesEF1, "./output/Sratio_bluesEF1_rep1.csv", row.names = FALSE)
Sratio_bluesEF1<- read.csv("./output/Sratio_bluesEF1_rep1.csv")

# ---------------------slaEF blues fitting model---------------------
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
SbluesEF1<-left_join(Sratio_bluesEF1, slabluesEF1, by= "name2")
fwrite( SbluesEF1, "./output/SbluesEF1_rep1.csv", row.names = FALSE)
# ---------------------slaMW blues fitting model---------------------
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
fwrite( SbluesMW1, "./output/SbluesMW1_rep1.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
SbluesEF1<- read.csv("./output/SbluesEF1_rep1.csv")
SbluesMW1<- read.csv("./output/SbluesMW1_rep1.csv")

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
S_blues1<- droplevels(S_blues1[S_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(S_blues1,"./output/S_blues1_rep1.csv", row.names = F)


##calculating blues for the selected low coheritable synthetic trait
Swave_lowcoh2<- read.csv("./output/Swave_rep2_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesMW1 <- bind_rows(Sratio_bluesMW1, temp)
Sratio_bluesMW1<- Sratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesMW1, "./output/Sratio_bluesMW1_rep2.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesEF1 <- bind_rows(Sratio_bluesEF1, temp)

Sratio_bluesEF1<- Sratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesEF1, "./output/Sratio_bluesEF1_rep2.csv", row.names = FALSE)
Sratio_bluesEF1<- read.csv("./output/Sratio_bluesEF1_rep2.csv")

# ---------------------slaEF blues fitting model---------------------
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
SbluesEF1<-left_join(Sratio_bluesEF1, slabluesEF1, by= "name2")
fwrite( SbluesEF1, "./output/SbluesEF1_rep2.csv", row.names = FALSE)
# ---------------------slaMW blues fitting model---------------------
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
fwrite( SbluesMW1, "./output/SbluesMW1_rep2.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
SbluesEF1<- read.csv("./output/SbluesEF1_rep2.csv")
SbluesMW1<- read.csv("./output/SbluesMW1_rep2.csv")

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
S_blues1<- droplevels(S_blues1[S_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(S_blues1,"./output/S_blues1_rep2.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
Swave_lowcoh2<- read.csv("./output/Swave_rep3_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesMW1 <- bind_rows(Sratio_bluesMW1, temp)
Sratio_bluesMW1<- Sratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesMW1, "./output/Sratio_bluesMW1_rep3.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesEF1 <- bind_rows(Sratio_bluesEF1, temp)

Sratio_bluesEF1<- Sratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesEF1, "./output/Sratio_bluesEF1_rep3.csv", row.names = FALSE)
Sratio_bluesEF1<- read.csv("./output/Sratio_bluesEF1_rep3.csv")

# ---------------------slaEF blues fitting model---------------------
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
SbluesEF1<-left_join(Sratio_bluesEF1, slabluesEF1, by= "name2")
fwrite( SbluesEF1, "./output/SbluesEF1_rep3.csv", row.names = FALSE)
# ---------------------slaMW blues fitting model---------------------
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
fwrite( SbluesMW1, "./output/SbluesMW1_rep3.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
SbluesEF1<- read.csv("./output/SbluesEF1_rep3.csv")
SbluesMW1<- read.csv("./output/SbluesMW1_rep3.csv")

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
S_blues1<- droplevels(S_blues1[S_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(S_blues1,"./output/S_blues1_rep3.csv", row.names = F)



###calculating blues for the selected low coheritable synthetic trait
Swave_lowcoh2<- read.csv("./output/Swave_rep4_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesMW1 <- bind_rows(Sratio_bluesMW1, temp)
Sratio_bluesMW1<- Sratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesMW1, "./output/Sratio_bluesMW1_rep4.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Swave_lowcoh2)[11]
models <- vector("list",length(variables))
Sratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesEF1 <- bind_rows(Sratio_bluesEF1, temp)

Sratio_bluesEF1<- Sratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesEF1, "./output/Sratio_bluesEF1_rep4.csv", row.names = FALSE)
Sratio_bluesEF1<- read.csv("./output/Sratio_bluesEF1_rep4.csv")

# ---------------------slaEF blues fitting model---------------------
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
SbluesEF1<-left_join(Sratio_bluesEF1, slabluesEF1, by= "name2")
fwrite( SbluesEF1, "./output/SbluesEF1_rep4.csv", row.names = FALSE)
# ---------------------slaMW blues fitting model---------------------
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
fwrite( SbluesMW1, "./output/SbluesMW1_rep4.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
SbluesEF1<- read.csv("./output/SbluesEF1_rep4.csv")
SbluesMW1<- read.csv("./output/SbluesMW1_rep4.csv")

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
S_blues1<- droplevels(S_blues1[S_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(S_blues1,"./output/S_blues1_rep4.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
Swave_lowcoh2<- read.csv("./output/Swave_rep5_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(Swave_lowcoh2)[12]
models <- vector("list",length(variables))
Sratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesMW1 <- bind_rows(Sratio_bluesMW1, temp)
Sratio_bluesMW1<- Sratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesMW1, "./output/Sratio_bluesMW1_rep5.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(Swave_lowcoh2)[12]
models <- vector("list",length(variables))
Sratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = Swave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
Sratio_bluesEF1 <- bind_rows(Sratio_bluesEF1, temp)

Sratio_bluesEF1<- Sratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( Sratio_bluesEF1, "./output/Sratio_bluesEF1_rep5.csv", row.names = FALSE)
Sratio_bluesEF1<- read.csv("./output/Sratio_bluesEF1_rep5.csv")

# ---------------------slaEF blues fitting model---------------------
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
SbluesEF1<-left_join(Sratio_bluesEF1, slabluesEF1, by= "name2")
fwrite( SbluesEF1, "./output/SbluesEF1_rep5.csv", row.names = FALSE)
# ---------------------slaMW blues fitting model---------------------
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
fwrite( SbluesMW1, "./output/SbluesMW1_rep5.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
SbluesEF1<- read.csv("./output/SbluesEF1_rep5.csv")
SbluesMW1<- read.csv("./output/SbluesMW1_rep5.csv")

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
S_blues1<- droplevels(S_blues1[S_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(S_blues1,"S_blues1_rep5.csv", row.names = F)

