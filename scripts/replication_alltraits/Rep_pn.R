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
source("./scripts/aux_function.R")


#replication 1
pn <- fread("./data/pn_breakdown1.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_1)

wave <-phenotypes |>
  clean_names()

sum(is.na(pn$corg))
sum(pn$corg< -1 | pn$corg > 1, na.rm = T)
sum(pn$corgblup< -1 | pn$corgblup > 1, na.rm = T)

pn <- pn |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(pn$corg))
ntraits <- ceiling(sum(!is.na(pn$coh2))*0.01)
sel <- pn |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
pnratios <- wave |> select(wave_822_wave_726:wave_1035_wave_1757)
pnratios <- t(pnratios)
colnames(pnratios) <- wave$plot_id
#distance
dist <- get_dist(pnratios, method = "pearson")
# fviz_dist(dist)
hcpn <- hclust(dist, method = "average")
sub_grp <- cutree(hcpn, k = 3)
table(sub_grp)
pnratios <- pnratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(pnratios)

pnratios<- as.data.frame(pnratios)
write.csv(pnratios,"./data/pnratios_rep1.csv", row.names=F)
pnratios<- read.csv("./data/pnratios_rep1.csv")

# Compute the absolute difference between wave_1 and wave_2
pnratios <- pnratios %>%
  mutate(diff = abs(wave_1 - wave_2))
pnratio_sel <- pnratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

pnratio_transform <- pnratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(pnratio_transform, "./data/pnratio_transform_rep1.csv", row.names=F)
pnratio_transform<- read.csv("./data/pnratio_transform_rep1.csv")
com_col <- colnames(wave)[colnames(wave) %in% pnratio_transform$wave_sel]
pnwave_pheno<- wave %>% select(1:9,fs_plsr_narea,com_col)
write.csv(pnwave_pheno,"./data/pnwave_pheno_rep1.csv", row.names = F)
pnwave<- read.csv("./data/pnwave_pheno_rep1.csv")
#plotting

l <- pnratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep1 <-pn |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "pn_Coh2-Replictaion1") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_rep1_pn.jpeg")


#replication 2
pn <- fread("./data/pn_breakdown2.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_2)

wave <-phenotypes |>
  clean_names()

sum(is.na(pn$corg))
sum(pn$corg< -1 | pn$corg > 1, na.rm = T)
sum(pn$corgblup< -1 | pn$corgblup > 1, na.rm = T)

pn <- pn |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(pn$corg))
ntraits <- ceiling(sum(!is.na(pn$coh2))*0.01)
sel <- pn |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
pnratios <- wave |> select(wave_1093_wave_2395:wave_1321_wave_2391)
pnratios <- t(pnratios)
colnames(pnratios) <- wave$plot_id
#distance
dist <- get_dist(pnratios, method = "pearson")
# fviz_dist(dist)
hcpn <- hclust(dist, method = "average")
sub_grp <- cutree(hcpn, k = 3)
table(sub_grp)
pnratios <- pnratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(pnratios)

pnratios<- as.data.frame(pnratios)
write.csv(pnratios,"./data/pnratios_rep2.csv", row.names=F)
pnratios<- read.csv("./data/pnratios_rep2.csv")

# Compute the absolute difference between wave_1 and wave_2
pnratios <- pnratios %>%
  mutate(diff = abs(wave_1 - wave_2))
pnratio_sel <- pnratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

pnratio_transform <- pnratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(pnratio_transform, "./data/pnratio_transform_rep2.csv", row.names=F)
pnratio_transform<- read.csv("./data/pnratio_transform_rep2.csv")
com_col <- colnames(wave)[colnames(wave) %in% pnratio_transform$wave_sel]
pnwave_pheno<- wave %>% select(1:9,fs_plsr_narea,com_col)
write.csv(pnwave_pheno,"./data/pnwave_pheno_rep2.csv", row.names = F)
pnwave<- read.csv("./data/pnwave_pheno_rep2.csv")
#plotting
l <- pnratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep2 <-pn |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "pn_Coh2-Replictaion2") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep2,"./figure/Heatmap_rep2_pn.jpeg")

#replication 3
pn <- fread("./data/pn_breakdown3.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_3)

wave <-phenotypes |>
  clean_names()

sum(is.na(pn$corg))
sum(pn$corg< -1 | pn$corg > 1, na.rm = T)
sum(pn$corgblup< -1 | pn$corgblup > 1, na.rm = T)

pn <- pn |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(pn$corg))
ntraits <- ceiling(sum(!is.na(pn$coh2))*0.01)
sel <- pn |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
pnratios <- wave |> select(wave_800_wave_742:wave_1137_wave_741)
pnratios <- t(pnratios)
colnames(pnratios) <- wave$plot_id
#distance
dist <- get_dist(pnratios, method = "pearson")
# fviz_dist(dist)
hcpn <- hclust(dist, method = "average")
sub_grp <- cutree(hcpn, k = 3)
table(sub_grp)
pnratios <- pnratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(pnratios)

pnratios<- as.data.frame(pnratios)
write.csv(pnratios,"./data/pnratios_rep3.csv", row.names=F)
pnratios<- read.csv("./data/pnratios_rep3.csv")

# Compute the absolute difference between wave_1 and wave_2
pnratios <- pnratios %>%
  mutate(diff = abs(wave_1 - wave_2))
pnratio_sel <- pnratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

pnratio_transform <- pnratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(pnratio_transform, "./data/pnratio_transform_rep3.csv", row.names=F)
pnratio_transform<- read.csv("./data/pnratio_transform_rep3.csv")
com_col <- colnames(wave)[colnames(wave) %in% pnratio_transform$wave_sel]
pnwave_pheno<- wave %>% select(1:9,fs_plsr_narea,com_col)
write.csv(pnwave_pheno,"./data/pnwave_pheno_rep3.csv", row.names = F)
pnwave<- read.csv("./data/pnwave_pheno_rep3.csv")
#plotting

l <- pnratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep3 <-pn |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "pn_Coh2-Replictaion3") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep3,"./figure/Heatmap_rep3_pn.jpeg")


#replication4
pn <- fread("./data/pn_breakdown4.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_4)

wave <-phenotypes |>
  clean_names()

sum(is.na(pn$corg))
sum(pn$corg< -1 | pn$corg > 1, na.rm = T)
sum(pn$corgblup< -1 | pn$corgblup > 1, na.rm = T)

pn <- pn |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(pn$corg))
ntraits <- ceiling(sum(!is.na(pn$coh2))*0.01)
sel <- pn |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
pnratios <- wave |> select(wave_1697_wave_1711:wave_1021_wave_1746)
pnratios <- t(pnratios)
colnames(pnratios) <- wave$plot_id
#distance
dist <- get_dist(pnratios, method = "pearson")
# fviz_dist(dist)
hcpn <- hclust(dist, method = "average")
sub_grp <- cutree(hcpn, k = 3)
table(sub_grp)
pnratios <- pnratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(pnratios)

pnratios<- as.data.frame(pnratios)
write.csv(pnratios,"./data/pnratios_rep4.csv", row.names=F)
pnratios<- read.csv("./data/pnratios_rep4.csv")

# Compute the absolute difference between wave_1 and wave_2
pnratios <- pnratios %>%
  mutate(diff = abs(wave_1 - wave_2))
pnratio_sel <- pnratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

pnratio_transform <- pnratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(pnratio_transform, "./data/pnratio_transform_rep4.csv", row.names=F)
pnratio_transform<- read.csv("./data/pnratio_transform_rep4.csv")
com_col <- colnames(wave)[colnames(wave) %in% pnratio_transform$wave_sel]
pnwave_pheno<- wave %>% select(1:9,fs_plsr_narea,com_col)
write.csv(pnwave_pheno,"./data/pnwave_pheno_rep4.csv", row.names = F)
pnwave<- read.csv("./data/pnwave_pheno_rep4.csv")
#plotting

l <- pnratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep4 <-pn |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "pn_Coh2-Replictaion4") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep4,"./figure/Heatmap_rep4_pn.jpeg")
#replication 5
pn <- fread("./data/pn_breakdown5.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_5)

wave <-phenotypes |>
  clean_names()

sum(is.na(pn$corg))
sum(pn$corg< -1 | pn$corg > 1, na.rm = T)
sum(pn$corgblup< -1 | pn$corgblup > 1, na.rm = T)

pn <- pn |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(pn$corg))
ntraits <- ceiling(sum(!is.na(pn$coh2))*0.01)
sel <- pn |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
pnratios <- wave |> select(wave_1001_wave_736:wave_1026_wave_721)
pnratios <- t(pnratios)
colnames(pnratios) <- wave$plot_id
#distance
dist <- get_dist(pnratios, method = "pearson")
# fviz_dist(dist)
hcpn <- hclust(dist, method = "average")
sub_grp <- cutree(hcpn, k = 3)
table(sub_grp)
pnratios <- pnratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(pnratios)

pnratios<- as.data.frame(pnratios)
write.csv(pnratios,"./data/pnratios_rep5.csv", row.names=F)
pnratios<- read.csv("./data/pnratios_rep5.csv")

# Compute the absolute difference between wave_1 and wave_2
pnratios <- pnratios %>%
  mutate(diff = abs(wave_1 - wave_2))
pnratio_sel <- pnratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

pnratio_transform <- pnratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(pnratio_transform, "./data/pnratio_transform_rep5.csv", row.names=F)
pnratio_transform<- read.csv("./data/pnratio_transform_rep5.csv")
com_col <- colnames(wave)[colnames(wave) %in% pnratio_transform$wave_sel]
pnwave_pheno<- wave %>% select(1:9,fs_plsr_narea,com_col)
write.csv(pnwave_pheno,"./data/pnwave_pheno_rep5.csv", row.names = F)
pnwave<- read.csv("./data/pnwave_pheno_rep5.csv")
#plotting

l <- pnratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep5 <-pn |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "pn_Coh2-Replictaion5") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep5,"./figure/Heatmap_rep5_pn.jpeg")



pnratio_transform_rep1<-read.csv("./data/pnratio_transform_rep1.csv")
pnratio_transform_rep2<-read.csv("./data/pnratio_transform_rep2.csv")
pnratio_transform_rep3<-read.csv("./data/pnratio_transform_rep3.csv")
pnratio_transform_rep4<-read.csv("./data/pnratio_transform_rep4.csv")
pnratio_transform_rep5<-read.csv("./data/pnratio_transform_rep5.csv")
pnratio_combined<- rbind(pnratio_transform_rep1,pnratio_transform_rep2,pnratio_transform_rep3,pnratio_transform_rep4,pnratio_transform_rep5 )
write.csv(pnratio_combined,"./data/pnratio_combined.csv", row.names = F)
pnratio_combined<- read.csv("./data/pnratio_combined.csv")

#calculating avg coh2 from 5 reps and plot the heatmaps
pn1 <- fread("./data/pn_breakdown1.csv", data.table = F)
pn2 <- fread("./data/pn_breakdown2.csv", data.table = F)
pn3<- fread("./data/pn_breakdown3.csv", data.table = F)
pn4 <- fread("./data/pn_breakdown4.csv", data.table = F)
pn5 <- fread("./data/pn_breakdown5.csv", data.table = F)
pn1_sub <- pn1[c("wave_1", "wave_2", "coh2")]
pn2_sub <- pn2[c("wave_1", "wave_2", "coh2")]
pn3_sub <- pn3[c("wave_1", "wave_2", "coh2")]
pn4_sub <- pn4[c("wave_1", "wave_2", "coh2")]
pn5_sub <- pn5[c("wave_1", "wave_2", "coh2")]
names(pn1_sub)[3] <- "coh2_rep1"
names(pn2_sub)[3] <- "coh2_rep2"
names(pn3_sub)[3] <- "coh2_rep3"
names(pn4_sub)[3] <- "coh2_rep4"
names(pn5_sub)[3] <- "coh2_rep5"
reps<- list(pn1_sub, pn2_sub,pn3_sub,pn4_sub,pn5_sub)
final_reps <- Reduce(function(x, y) merge(x, y, by = c("wave_1", "wave_2")), reps)

# Calculate the average of coh2 values from different reps
final_reps <- final_reps %>%
  mutate(coh2_avg = rowMeans(select(., matches("^coh2_rep")), na.rm = TRUE))

# Define the symbols you want to use
symbols <- c("*", "+", "o", "x", "s")
pnratio_combined <- pnratio_combined %>%
  mutate(label = rep(symbols,each=3, length.out = nrow(pnratio_combined)))
# Process Sratio_combined to prepare it for merging
l <- pnratio_combined %>%
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
  labs(x = "Wave1", y = "Wave2", title = "Avg_coh2_Combinedpn_reps") +
  geom_label(data = l, aes(label = label),size=4) +
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmappn_avgreps.jpeg")

#for five replications  first stage analysis
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% clean_names()
pnratio_transform<- read.csv("./data/pnratio_combined.csv")

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
indv_sample<- read.csv("./data/indv_sample.csv")
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
indv_sample<- read.csv("./data/indv_sample.csv")
pnratiobluesEF_rep1<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(name2,2:4)
pnratiobluesEF_rep2<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2,5:7)
pnratiobluesEF_rep3<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2,8:10)
pnratiobluesEF_rep4<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2,11:13)
pnratiobluesEF_rep5<-pnratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2,14:16)

fwrite( pnratiobluesEF_rep5, "./output/pnratiobluesEF_rep5.csv", row.names = FALSE)

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

#combining blues for trait and waveratio for EF and MW location for rep1
pnbluesEF_rep1<- read.csv("./output/pnbluesEF_rep1.csv")
pnbluesMW_rep1<- read.csv("./output/pnbluesMW_rep1.csv")

pnblues_rep1 <- bind_rows(
  "EF" = pnbluesEF_rep1,
  "MW" = pnbluesMW_rep1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues_rep1 <-pnblues_rep1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues_rep1<- subset(pnblues_rep1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
pn_blues_rep1<- droplevels(pn_blues_rep1[pn_blues_rep1$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(pn_blues_rep1,"./output/pn_blues_rep1.csv", row.names = F)

#combining blues for trait and waveratio for EF and MW location for rep2
pnbluesEF_rep2<- read.csv("./output/pnbluesEF_rep2.csv")
pnbluesMW_rep2<- read.csv("./output/pnbluesMW_rep2.csv")

pnblues_rep2 <- bind_rows(
  "EF" = pnbluesEF_rep2,
  "MW" = pnbluesMW_rep2,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues_rep2 <-pnblues_rep2|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues_rep2<- subset(pnblues_rep2, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
pn_blues_rep2<- droplevels(pn_blues_rep2[pn_blues_rep2$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(pn_blues_rep2,"./output/pn_blues_rep2.csv", row.names = F)
#combining blues for trait and waveratio for EF and MW location for rep2
pnbluesEF_rep3<- read.csv("./output/pnbluesEF_rep3.csv")
pnbluesMW_rep3<- read.csv("./output/pnbluesMW_rep3.csv")
pnblues_rep3 <- bind_rows(
  "EF" = pnbluesEF_rep3,
  "MW" = pnbluesMW_rep3,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues_rep3 <-pnblues_rep3|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues_rep3<- subset(pnblues_rep3, !is.na(taxa)) %>% select(-name2,-Corrected_names)
pn_blues_rep3<- droplevels(pn_blues_rep3[pn_blues_rep3$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(pn_blues_rep3,"./output/pn_blues_rep3.csv", row.names = F)

#combining blues for trait and waveratio for EF and MW location for rep4
pnbluesEF_rep4<- read.csv("./output/pnbluesEF_rep4.csv")
pnbluesMW_rep4<- read.csv("./output/pnbluesMW_rep4.csv")
pnblues_rep4 <- bind_rows(
  "EF" = pnbluesEF_rep4,
  "MW" = pnbluesMW_rep4,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues_rep4 <-pnblues_rep4|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues_rep4<- subset(pnblues_rep4, !is.na(taxa)) %>% select(-name2,-Corrected_names)
pn_blues_rep4 <- droplevels(pn_blues_rep4[pn_blues_rep4$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(pn_blues_rep4,"./output/pn_blues_rep4.csv", row.names = F)
#combining blues for trait and waveratio for EF and MW location for rep5
pnbluesEF_rep5<- read.csv("./output/pnbluesEF_rep5.csv")
pnbluesMW_rep5<- read.csv("./output/pnbluesMW_rep5.csv")
pnblues_rep5 <- bind_rows(
  "EF" = pnbluesEF_rep5,
  "MW" = pnbluesMW_rep5,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
pnblues_rep5 <-pnblues_rep5|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
pn_blues_rep5<- subset(pnblues_rep5, !is.na(taxa)) %>% select(-name2,-Corrected_names)
pn_blues_rep5<- droplevels(pn_blues_rep5[pn_blues_rep5$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(pn_blues_rep5,"./output/pn_blues_rep5.csv", row.names = F)


#----------------2nd step --------------------
pn_blues<- read.csv("pn_blues_rep1.csv")
pn_bluesEF<- pn_blues %>% filter(env== "EF")
pn_bluesMW<- pn_blues %>% filter(env== "MW")
pn_bluesEF<- pn_bluesEF%>% mutate(taxa= factor(taxa))
pn_bluesMW<- pn_bluesMW%>% mutate(taxa= factor(taxa))
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
#----- ST plar Narea for rep1 --------
# ---------------------sla-running cross validation for EFMW loc---------------------
result_pn_rep1 <-
  crossv(
    sort = sort,
    train = pn_bluesEF,
    validation = pn_bluesMW,
    mytrait = "pn",
    kin=kin,
    scheme = "EFMW")

corr_pn_EFMW<- data.frame(result_pn_rep1$ac)
fwrite(corr_pn_EFMW, "./output/corr_pn_EFMWrep1.csv")

# ---------------------plsr narea-running cross validation for MWEF loc---------------------
result_pn_rep1 <-
  crossv(
    sort = sort,
    train = pn_bluesMW,
    validation = pn_bluesEF,
    mytrait = "pn",
    kin=kin,
    scheme = "MWEF_rep1")

corr_pn_MWEF<- data.frame(result_pn_rep1$ac)
fwrite(corr_pn_MWEF, "./output/corr_pn_MWEFrep1.csv")

#----- ST Narea for rep2 --------
pn_blues<- read.csv("pn_blues_rep2.csv")
pn_bluesEF<- pn_blues %>% filter(env== "EF")
pn_bluesMW<- pn_blues %>% filter(env== "MW")
pn_bluesEF<- pn_bluesEF%>% mutate(taxa= factor(taxa))
pn_bluesMW<- pn_bluesMW%>% mutate(taxa= factor(taxa))
N_blues<- read.csv("N_blues_rep2.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)

# ---------------------pn-running cross validation for EFMW loc---------------------
result_pn_rep2 <-
  crossv(
    sort = sort,
    train = pn_bluesEF,
    validation = pn_bluesMW,
    mytrait = "pn",
    kin=kin,
    scheme = "EFMW_rep2")

corr_pn_EFMW<- data.frame(result_pn_rep2$ac)
fwrite(corr_pn_EFMW, "./output/corr_pn_EFMWrep2.csv")

# ---------------------pn-running cross validation for MWEF loc---------------------
result_pn_rep2 <-
  crossv(
    sort = sort,
    train = pn_bluesMW,
    validation = pn_bluesEF,
    mytrait = "pn",
    kin = kin,
    scheme = "MWEF_rep2")

corr_pn_MWEF<- data.frame(result_pn_rep2$ac)
fwrite(corr_pn_MWEF, "./output/corr_pn_MWEFrep2.csv")

#----- ST Narea for rep3 --------
pn_blues<- read.csv("pn_blues_rep3.csv")
pn_bluesEF<- pn_blues %>% filter(env== "EF")
pn_bluesMW<- pn_blues %>% filter(env== "MW")
pn_bluesEF<- pn_bluesEF%>% mutate(taxa= factor(taxa))
pn_bluesMW<- pn_bluesMW%>% mutate(taxa= factor(taxa))

N_blues<- read.csv("N_blues_rep3.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------pn-running cross validation for EFMW loc---------------------
result_pn_rep3 <-
  crossv(
    sort = sort,
    train = pn_bluesEF,
    validation = pn_bluesMW,
    mytrait = "pn",
    kin=kin,
    scheme = "EFMW_rep3")

corr_pn_EFMW<- data.frame(result_pn_rep3$ac)
fwrite(corr_pn_EFMW, "./output/corr_pn_EFMWrep3.csv")

# ---------------------pn-running cross validation for MWEF loc---------------------
result_pn_rep3 <-
  crossv(
    sort = sort,
    train = pn_bluesMW,
    validation = pn_bluesEF,
    mytrait = "pn",
    kin = kin,
    scheme = "MWEF_rep3")

corr_pn_MWEF<- data.frame(result_pn_rep3$ac)
fwrite(corr_pn_MWEF, "./output/corr_pn_MWEFrep3.csv")

#----- ST Narea for rep4 --------
pn_blues<- read.csv("pn_blues_rep4.csv")
pn_bluesEF<- pn_blues %>% filter(env== "EF")
pn_bluesMW<- pn_blues %>% filter(env== "MW")
pn_bluesEF<- pn_bluesEF%>% mutate(taxa= factor(taxa))
pn_bluesMW<- pn_bluesMW%>% mutate(taxa= factor(taxa))
N_blues<- read.csv("N_blues_rep4.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------pn-running cross validation for EFMW loc---------------------
result_pn_rep4 <-
  crossv(
    sort = sort,
    train = pn_bluesEF,
    validation = pn_bluesMW,
    mytrait = "pn",
    kin=kin,
    scheme = "EFMW_rep4")

corr_pn_EFMW<- data.frame(result_pn_rep4$ac)
fwrite(corr_pn_EFMW, "./output/corr_pn_EFMWrep4.csv")

# ---------------------pn-running cross validation for MWEF loc---------------------
result_pn_rep4 <-
  crossv(
    sort = sort,
    train = pn_bluesMW,
    validation = pn_bluesEF,
    mytrait = "pn",
    kin = kin,
    scheme = "MWEF_rep4")

corr_pn_MWEF<- data.frame(result_pn_rep4$ac)
fwrite(corr_pn_MWEF, "./output/corr_pn_MWEFrep4.csv")

#----- ST Narea for rep5 --------
pn_blues<- read.csv("pn_blues_rep5.csv")
pn_bluesEF<- pn_blues %>% filter(env== "EF")
pn_bluesMW<- pn_blues %>% filter(env== "MW")
pn_bluesEF<- pn_bluesEF%>% mutate(taxa= factor(taxa))
pn_bluesMW<- pn_bluesMW%>% mutate(taxa= factor(taxa))

N_blues<- read.csv("N_blues_rep5.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------pn-running cross validation for EFMW loc---------------------
result_pn_rep5 <-
  crossv(
    sort = sort,
    train = pn_bluesEF,
    validation = pn_bluesMW,
    mytrait = "pn",
    kin=kin,
    scheme = "EFMW_rep5")

corr_pn_EFMW<- data.frame(result_pn_rep5$ac)
fwrite(corr_pn_EFMW, "./output/corr_pn_EFMWrep5.csv")

# ---------------------pn-running cross validation for MWEF loc---------------------
result_pn_rep5 <-
  crossv(
    sort = sort,
    train = pn_bluesMW,
    validation = pn_bluesEF,
    mytrait = "pn",
    kin = kin,
    scheme = "MWEF_rep5")

corr_pn_MWEF<- data.frame(result_pn_rep5$ac)
fwrite(corr_pn_MWEF, "./output/corr_pn_MWEFrep5.csv")

# #obtaining synthetic trait table and value in each rep
# pnratio_combined <- read.csv("./data/pnratio_combined.csv", header = TRUE) %>% select(-3)
# pnratio_combined$ST <- rep(c("ST1", "ST2", "ST3"), times = 5)
# pnratio_combined$Rep <- rep(1:5, each = 3)
# pnratio_combined<- pnratio_combined %>% select(-wave_1, -wave_2)
# pnratio_combined$coh2_Rep1<- c(0.63,0.65,0.63,0.61,0.55,0.49,0.55,0.62, 0.52,0.63,0.64,0.60,0.55,0.62,0.33)
# pnratio_combined$coh2_Rep2<- c(0.65,0.64,0.57,0.64,0.63,0.60,0.50,0.61,0.54,0.60,NA,NA,0.56,NA,0.21)
# pnratio_combined$coh2_Rep3<-c(0.43,0.45,0.40,0.42,0.53,0.46,0.51,0.57,0.48,0.43,0.39,0.50,0.52,0.58,0.26)
# pnratio_combined$coh2_Rep4<- c(0.49,0.53,0.44,0.49,0.39,0.37,0.38,0.49,0.36,0.49,0.49,0.49,0.38,0.53,NA)
# pnratio_combined$coh2_Rep5<- c(0.49,0.40,0.76,0.46,NA,0.48,0.45,0.34,0.45,0.34,NA,NA,0.59,0.53,0.52)
# library(gtable)
# gt_table <- gt(pnratio_combined)
# # Apply the yellow highlighting
# for (i in 1:5) {
#   column_name <- paste0("coh2_Rep", i)
#   gt_table <- gt_table %>%
#     tab_style(
#       style = cell_fill(color = "yellow"),
#       locations = cells_body(
#         columns = column_name,
#         rows = Rep == i
#       )
#     )
# }
# htmltools::save_html(gt_table, file = "coh2_all_pn_rep.html")
# webshot("coh2_all_pn_rep.html", "./figure/coh2_all_pn_rep.png")
#
#


#selecting synthetic trait that has 0 coh2 value from each replication
#Rep1
# Select the one wave ratio with the lowest coh2
pn1 <- fread("./data/pn_breakdown1.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_1))
wave <-phenotypes |>
  clean_names()
selected <- pn1 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pnwave_lowcoh2<- wave %>% select(1:9, fs_plsr_narea,wave_351,wave_1389)
wave_351_wave_1389 <- paste("wave_351", "wave_1389", sep = "_")
pnwave_lowcoh2[[wave_351_wave_1389]] <- pnwave_lowcoh2[["wave_351"]] / pnwave_lowcoh2[["wave_1389"]]
pnwave_lowcoh2<- pnwave_lowcoh2 %>% select(-wave_351,-wave_1389)
write.csv(pnwave_lowcoh2, "./output/pnwave_rep1_lowcoh2.csv", row.names = F)

#Rep2
pn2 <- fread("./data/pn_breakdown2.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_2))
wave <-phenotypes |>
  clean_names()

selected <- pn2 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pnwave_lowcoh2<- wave %>% select(1:9, fs_plsr_narea,wave_2326,wave_1956)
wave_2326_wave_1956 <- paste("wave_2326", "wave_1956", sep = "_")
pnwave_lowcoh2[[wave_2326_wave_1956]] <- pnwave_lowcoh2[["wave_2326"]] / pnwave_lowcoh2[["wave_1956"]]
pnwave_lowcoh2<- pnwave_lowcoh2 %>% select(-wave_2326,-wave_1956)
write.csv(pnwave_lowcoh2, "./output/pnwave_rep2_lowcoh2.csv", row.names = F)

#Rep3
pn3 <- fread("./data/pn_breakdown3.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_3))
wave <-phenotypes |>
  clean_names()

selected <- pn3 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pnwave_lowcoh2<- wave %>% select(1:9, fs_plsr_narea,wave_1289,wave_759)
wave_1289_wave_759 <- paste("wave_1289", "wave_759", sep = "_")
pnwave_lowcoh2[[wave_1289_wave_759]] <- pnwave_lowcoh2[["wave_1289"]] / pnwave_lowcoh2[["wave_759"]]
pnwave_lowcoh2<- pnwave_lowcoh2 %>% select(-wave_1289,-wave_759)
write.csv(pnwave_lowcoh2, "./output/pnwave_rep3_lowcoh2.csv", row.names = F)


#Rep4
pn4 <- fread("./data/pn_breakdown4.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_4))
wave <-phenotypes |>
  clean_names()

selected <- pn4 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pnwave_lowcoh2<- wave %>% select(1:9, fs_plsr_narea,wave_667,wave_573)
wave_667_wave_573 <- paste("wave_667", "wave_573", sep = "_")
pnwave_lowcoh2[[wave_667_wave_573]] <- pnwave_lowcoh2[["wave_667"]] / pnwave_lowcoh2[["wave_573"]]
pnwave_lowcoh2<- pnwave_lowcoh2 %>% select(-wave_667,-wave_573)
write.csv(pnwave_lowcoh2, "./output/pnwave_rep4_lowcoh2.csv", row.names = F)

#Rep5
pn5 <- fread("./data/pn_breakdown5.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_5))
wave <-phenotypes |>
  clean_names()

selected <- pn5 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pnwave_lowcoh2<- wave %>% select(1:9, fs_plsr_narea,wave_2278,wave_2111)
wave_2278_wave_2111 <- paste("wave_2278", "wave_2111", sep = "_")
pnwave_lowcoh2[[wave_2278_wave_2111]] <- pnwave_lowcoh2[["wave_2278"]] / pnwave_lowcoh2[["wave_2111"]]
pnwave_lowcoh2<- pnwave_lowcoh2 %>% select(-wave_2278,-wave_2111)
write.csv(pnwave_lowcoh2, "./output/pnwave_rep5_lowcoh2.csv", row.names = F)



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
