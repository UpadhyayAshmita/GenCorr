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
ps <- fread("./data/ps_breakdown1.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_1)

wave <-phenotypes |>
  clean_names()

sum(is.na(ps$corg))
sum(ps$corg< -1 | ps$corg > 1, na.rm = T)
sum(ps$corgblup< -1 | ps$corgblup > 1, na.rm = T)

ps <- ps |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(ps$corg))
ntraits <- ceiling(sum(!is.na(ps$coh2))*0.01)
sel <- ps |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
psratios <- wave |> select(wave_1645_wave_1362:wave_2295_wave_1369)
psratios <- t(psratios)
colnames(psratios) <- wave$plot_id
#distance
dist <- get_dist(psratios, method = "pearson")
# fviz_dist(dist)
hcps <- hclust(dist, method = "average")
sub_grp <- cutree(hcps, k = 3)
table(sub_grp)
psratios <- psratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(psratios)

psratios<- as.data.frame(psratios)
write.csv(psratios,"./data/psratios_rep1.csv", row.names=F)
psratios<- read.csv("./data/psratios_rep1.csv")

# Compute the absolute difference between wave_1 and wave_2
psratios <- psratios %>%
  mutate(diff = abs(wave_1 - wave_2))
psratio_sel <- psratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

psratio_transform <- psratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(psratio_transform, "./data/psratio_transform_rep1.csv", row.names=F)
psratio_transform<- read.csv("./data/psratio_transform_rep1.csv")
com_col <- colnames(wave)[colnames(wave) %in% psratio_transform$wave_sel]
pswave_pheno<- wave %>% select(1:9,plsr_sla_sorghum,com_col)
write.csv(pswave_pheno,"./data/pswave_pheno_rep1.csv", row.names = F)
pswave<- read.csv("./data/pswave_pheno_rep1.csv")
#plotting

l <- psratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep1 <-ps |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "ps_Coh2-Replictaion1") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/Heatmap_rep1_ps.jpeg")


#replication 2
ps <- fread("./data/ps_breakdown2.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_2)

wave <-phenotypes |>
  clean_names()

sum(is.na(ps$corg))
sum(ps$corg< -1 | ps$corg > 1, na.rm = T)
sum(ps$corgblup< -1 | ps$corgblup > 1, na.rm = T)

ps <- ps |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(ps$corg))
ntraits <- ceiling(sum(!is.na(ps$coh2))*0.01)
sel <- ps |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
psratios <- wave |> select(wave_1709_wave_1681:wave_389_wave_922)
psratios <- t(psratios)
colnames(psratios) <- wave$plot_id
#distance
dist <- get_dist(psratios, method = "pearson")
# fviz_dist(dist)
hcps <- hclust(dist, method = "average")
sub_grp <- cutree(hcps, k = 3)
table(sub_grp)
psratios <- psratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(psratios)

psratios<- as.data.frame(psratios)
write.csv(psratios,"./data/psratios_rep2.csv", row.names=F)
psratios<- read.csv("./data/psratios_rep2.csv")

# Compute the absolute difference between wave_1 and wave_2
psratios <- psratios %>%
  mutate(diff = abs(wave_1 - wave_2))
psratio_sel <- psratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

psratio_transform <- psratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(psratio_transform, "./data/psratio_transform_rep2.csv", row.names=F)
psratio_transform<- read.csv("./data/psratio_transform_rep2.csv")
com_col <- colnames(wave)[colnames(wave) %in% psratio_transform$wave_sel]
pswave_pheno<- wave %>% select(1:9,plsr_sla_sorghum,com_col)
write.csv(pswave_pheno,"./data/pswave_pheno_rep2.csv", row.names = F)
pswave<- read.csv("./data/pswave_pheno_rep2.csv")
#plotting
l <- psratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep2 <-ps |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "ps_Coh2-Replictaion2") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep2,"./figure/Heatmap_rep2_ps.jpeg")

#replication 3
ps <- fread("./data/ps_breakdown3.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_3)

wave <-phenotypes |>
  clean_names()

sum(is.na(ps$corg))
sum(ps$corg< -1 | ps$corg > 1, na.rm = T)
sum(ps$corgblup< -1 | ps$corgblup > 1, na.rm = T)

ps <- ps |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(ps$corg))
ntraits <- ceiling(sum(!is.na(ps$coh2))*0.01)
sel <- ps |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
psratios <- wave |> select(wave_2147_wave_1521:wave_2294_wave_1857)
psratios <- t(psratios)
colnames(psratios) <- wave$plot_id
#distance
dist <- get_dist(psratios, method = "pearson")
# fviz_dist(dist)
hcps <- hclust(dist, method = "average")
sub_grp <- cutree(hcps, k = 3)
table(sub_grp)
psratios <- psratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(psratios)

psratios<- as.data.frame(psratios)
write.csv(psratios,"./data/psratios_rep3.csv", row.names=F)
psratios<- read.csv("./data/psratios_rep3.csv")

# Compute the absolute difference between wave_1 and wave_2
psratios <- psratios %>%
  mutate(diff = abs(wave_1 - wave_2))
psratio_sel <- psratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

psratio_transform <- psratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(psratio_transform, "./data/psratio_transform_rep3.csv", row.names=F)
psratio_transform<- read.csv("./data/psratio_transform_rep3.csv")
com_col <- colnames(wave)[colnames(wave) %in% psratio_transform$wave_sel]
pswave_pheno<- wave %>% select(1:9,plsr_sla_sorghum,com_col)
write.csv(pswave_pheno,"./data/pswave_pheno_rep3.csv", row.names = F)
pswave<- read.csv("./data/pswave_pheno_rep3.csv")
#plotting

l <- psratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep3 <-ps |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "ps_Coh2-Replictaion3") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep3,"./figure/Heatmap_rep3_ps.jpeg")


#replication4
ps <- fread("./data/ps_breakdown4.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_4)

wave <-phenotypes |>
  clean_names()

sum(is.na(ps$corg))
sum(ps$corg< -1 | ps$corg > 1, na.rm = T)
sum(ps$corgblup< -1 | ps$corgblup > 1, na.rm = T)

ps <- ps |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(ps$corg))
ntraits <- ceiling(sum(!is.na(ps$coh2))*0.01)
sel <- ps |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
psratios <- wave |> select(wave_1208_wave_1276:wave_1725_wave_1135)
psratios <- t(psratios)
colnames(psratios) <- wave$plot_id
#distance
dist <- get_dist(psratios, method = "pearson")
# fviz_dist(dist)
hcps <- hclust(dist, method = "average")
sub_grp <- cutree(hcps, k = 3)
table(sub_grp)
psratios <- psratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(psratios)

psratios<- as.data.frame(psratios)
write.csv(psratios,"./data/psratios_rep4.csv", row.names=F)
psratios<- read.csv("./data/psratios_rep4.csv")

# Compute the absolute difference between wave_1 and wave_2
psratios <- psratios %>%
  mutate(diff = abs(wave_1 - wave_2))
psratio_sel <- psratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

psratio_transform <- psratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(psratio_transform, "./data/psratio_transform_rep4.csv", row.names=F)
psratio_transform<- read.csv("./data/psratio_transform_rep4.csv")
com_col <- colnames(wave)[colnames(wave) %in% psratio_transform$wave_sel]
pswave_pheno<- wave %>% select(1:9,plsr_sla_sorghum,com_col)
write.csv(pswave_pheno,"./data/pswave_pheno_rep4.csv", row.names = F)
pswave<- read.csv("./data/pswave_pheno_rep4.csv")
#plotting

l <- psratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep4 <-ps |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "ps_Coh2-Replictaion4") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep4,"./figure/Heatmap_rep4_ps.jpeg")
#replication 5
ps <- fread("./data/ps_breakdown5.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_5)

wave <-phenotypes |>
  clean_names()

sum(is.na(ps$corg))
sum(ps$corg< -1 | ps$corg > 1, na.rm = T)
sum(ps$corgblup< -1 | ps$corgblup > 1, na.rm = T)

ps <- ps |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(is.na(ps$corg))
ntraits <- ceiling(sum(!is.na(ps$coh2))*0.01)
sel <- ps |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
psratios <- wave |> select(wave_759_wave_760:wave_423_wave_478)
psratios <- t(psratios)
colnames(psratios) <- wave$plot_id
#distance
dist <- get_dist(psratios, method = "pearson")
# fviz_dist(dist)
hcps <- hclust(dist, method = "average")
sub_grp <- cutree(hcps, k = 3)
table(sub_grp)
psratios <- psratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=1) |> select(coh2, wave_1, wave_2)
print(psratios)

psratios<- as.data.frame(psratios)
write.csv(psratios,"./data/psratios_rep5.csv", row.names=F)
psratios<- read.csv("./data/psratios_rep5.csv")

# Compute the absolute difference between wave_1 and wave_2
psratios <- psratios %>%
  mutate(diff = abs(wave_1 - wave_2))
psratio_sel <- psratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

psratio_transform <- psratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(psratio_transform, "./data/psratio_transform_rep5.csv", row.names=F)
psratio_transform<- read.csv("./data/psratio_transform_rep5.csv")
com_col <- colnames(wave)[colnames(wave) %in% psratio_transform$wave_sel]
pswave_pheno<- wave %>% select(1:9,plsr_sla_sorghum,com_col)
write.csv(pswave_pheno,"./data/pswave_pheno_rep5.csv", row.names = F)
pswave<- read.csv("./data/pswave_pheno_rep5.csv")
#plotting

l <- psratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap_rep5 <-ps |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                            wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "ps_Coh2-Replictaion5") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap_rep5,"./figure/Heatmap_rep5_ps.jpeg")



psratio_transform_rep1<-read.csv("./data/psratio_transform_rep1.csv")
psratio_transform_rep2<-read.csv("./data/psratio_transform_rep2.csv")
psratio_transform_rep3<-read.csv("./data/psratio_transform_rep3.csv")
psratio_transform_rep4<-read.csv("./data/psratio_transform_rep4.csv")
psratio_transform_rep5<-read.csv("./data/psratio_transform_rep5.csv")
psratio_combined<- rbind(psratio_transform_rep1,psratio_transform_rep2,psratio_transform_rep3,psratio_transform_rep4,psratio_transform_rep5 )
write.csv(psratio_combined,"./data/psratio_combined.csv", row.names = F)
psratio_combined<- read.csv("./data/psratio_combined.csv")

#calculating avg coh2 from 5 reps and plot the heatmaps
ps1 <- fread("./data/ps_breakdown1.csv", data.table = F)
ps2 <- fread("./data/ps_breakdown2.csv", data.table = F)
ps3<- fread("./data/ps_breakdown3.csv", data.table = F)
ps4 <- fread("./data/ps_breakdown4.csv", data.table = F)
ps5 <- fread("./data/ps_breakdown5.csv", data.table = F)
ps1_sub <- ps1[c("wave_1", "wave_2", "coh2")]
ps2_sub <- ps2[c("wave_1", "wave_2", "coh2")]
ps3_sub <- ps3[c("wave_1", "wave_2", "coh2")]
ps4_sub <- ps4[c("wave_1", "wave_2", "coh2")]
ps5_sub <- ps5[c("wave_1", "wave_2", "coh2")]
names(ps1_sub)[3] <- "coh2_rep1"
names(ps2_sub)[3] <- "coh2_rep2"
names(ps3_sub)[3] <- "coh2_rep3"
names(ps4_sub)[3] <- "coh2_rep4"
names(ps5_sub)[3] <- "coh2_rep5"
reps<- list(ps1_sub, ps2_sub,ps3_sub,ps4_sub,ps5_sub)
final_reps <- Reduce(function(x, y) merge(x, y, by = c("wave_1", "wave_2")), reps)

# Calculate the average of coh2 values from different reps
final_reps <- final_reps %>%
  mutate(coh2_avg = rowMeans(select(., matches("^coh2_rep")), na.rm = TRUE))

# Define the symbols you want to use
symbols <- c("*", "+", "o", "x", "s")
psratio_combined <- psratio_combined %>%
  mutate(label = rep(symbols,each=3, length.out = nrow(psratio_combined)))
# Process Sratio_combined to prepare it for merging
l <- psratio_combined %>%
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
  labs(x = "Wave1", y = "Wave2", title = "Avg_coh2_Combinedps_reps") +
  geom_label(data = l, aes(label = label),size=4) +
  theme_minimal()

ggsave(plot=heatmap_rep1,"./figure/HeatmapPS_avgreps.jpeg")

#for five replications  first stage analysis
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% clean_names()
psratio_transform<- read.csv("./data/psratio_combined.csv")

# Compute and append ratios for each wave pair
for (i in 1:nrow(psratio_transform)) {
  wave1 <- psratio_transform$wave_1[i]
  wave2 <- psratio_transform$wave_2[i]
  if (wave1 %in% colnames(phenotypes) && wave2 %in% colnames(phenotypes)) {
    # Calculate the ratio
    ratio_column_name <- paste(wave1, wave2, sep = "_")
    phenotypes[[ratio_column_name]] <- phenotypes[[wave1]] / phenotypes[[wave2]]
  } else {
    warning(paste("Missing columns:", wave1, "or", wave2, "in phenotypes."))
  }
}

phenotypes<- phenotypes %>% select(-c(16:2166),-narea,-n_perc,-sla,-fs_plsr_narea,-fs_plsr_nmass)
write.csv(phenotypes,"./data/phenotypes_allRepps.csv", row.names = F)
phenotypes_allRep<- read.csv("./data/phenotypes_allRepps.csv")
#obtaining blues for replication ps
pswave_pheno <- read.csv("./data/phenotypes_allRepps.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(pswave_pheno)[11:25]
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
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
psratiobluesMW_rep1<-psratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(1:4)
psratiobluesMW_rep2<-psratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2, 5:7)
psratiobluesMW_rep3<-psratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2, 8:10)
psratiobluesMW_rep4<-psratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2, 11:13)
psratiobluesMW_rep5<-psratio_bluesMW %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2, 14:16)
fwrite(psratiobluesMW_rep5, "./output/psratiobluesMW_rep5.csv", row.names = FALSE)
fwrite(psratiobluesMW_rep4, "./output/psratiobluesMW_rep4.csv", row.names = FALSE)
fwrite(psratiobluesMW_rep3, "./output/psratiobluesMW_rep3.csv", row.names = FALSE)
fwrite(psratiobluesMW_rep2, "./output/psratiobluesMW_rep2.csv", row.names = FALSE)
fwrite(psratiobluesMW_rep1, "./output/psratiobluesMW_rep1.csv", row.names = FALSE)
psratiobluesMW_rep5<- read.csv("./output/psratiobluesMW_rep5.csv")
#calculating blues for each wavelength
variables <- colnames(pswave_pheno)[11:25]
models <- vector("list",length(variables))
psratio_bluesEF_a <-data.frame()
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
      psratio_bluesEF_a <- bind_rows(psratio_bluesEF_a, temp)
      cat('\n')

    }  else {
      psratio_bluesEF_a <- bind_rows(psratio_bluesEF_a, data.frame(name2 = NA,
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

variables <- colnames(pswave_pheno)[21]
models <- vector("list",length(variables))
psratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      #residual =  ~ar1(range):ar1(row),
      data = pswave_pheno,subset= loc== "EF", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp1 <- models[[i]]$predictions$pvals[, 1:2]
      temp1$wave <- variables[i]
      psratio_bluesEF1 <- bind_rows(psratio_bluesEF1, temp1)
      cat('\n')
    }}

psratio_bluesEF<- rbind(psratio_bluesEF_a, psratio_bluesEF1)
psratio_bluesEF<- psratio_bluesEF %>% pivot_wider(names_from = wave, values_from= predicted.value)
psratio_bluesEF<- psratio_bluesEF[-c(855), ]
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
psratiobluesEF_rep1<-psratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_1) %>% select(name2,2:4)
psratiobluesEF_rep2<-psratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_2) %>% select(name2,5:7)
psratiobluesEF_rep3<-psratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_3) %>% select(name2,8:10)
psratiobluesEF_rep4<-psratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_4) %>% select(name2,11:13)
psratiobluesEF_rep5<-psratio_bluesEF %>% filter(!Name2 %in% indv_sample$Rep_5) %>% select(name2,14:16)

fwrite( psratiobluesEF_rep5, "./output/psratiobluesEF_rep5.csv", row.names = FALSE)

# ---------------------psEF blues fitting model---------------------
pswave_pheno <- read.csv("./data/phenotypes_whole.csv") %>% clean_names()
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
psbluesEF <- asreml(
  fixed = plsr_sla_sorghum ~ name2 +set,
  random =  ~  block,
  data = pswave_pheno, subset = loc == "EF",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
psbluesEF<- update.asreml(psbluesEF)
# ---------------------storing prediction---------------------
psbluesEF <- data.frame(
  name2= psbluesEF$predictions$pvals$name2,
  ps = round(psbluesEF$predictions$pvals$predicted.value, 3))
psbluesEF_rep1<-left_join(psratiobluesEF_rep1, psbluesEF, by= "name2")
psbluesEF_rep2<-left_join(psratiobluesEF_rep2, psbluesEF, by= "name2")
psbluesEF_rep3<-left_join(psratiobluesEF_rep3, psbluesEF, by= "name2")
psbluesEF_rep4<-left_join(psratiobluesEF_rep4, psbluesEF, by= "name2")
psbluesEF_rep5<-left_join(psratiobluesEF_rep5, psbluesEF, by= "name2")
fwrite(psbluesEF_rep5, "./output/psbluesEF_rep5.csv", row.names = FALSE)

# ---------------------NareaMW blues fitting model---------------------
psbluesMW <- asreml(
  fixed = plsr_sla_sorghum ~ name2 +set,
  random =  ~  block,
  residual =  ~ ar1(range):ar1(row),
  data = pswave_pheno, subset = loc == "MW",
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2"))
psbluesMW<- update.asreml(psbluesMW)
# ---------------------storing prediction---------------------
psbluesMW <- data.frame(
  name2= psbluesMW$predictions$pvals$name2,
  ps = round(psbluesMW$predictions$pvals$predicted.value, 3))
psbluesMW_rep1<-left_join(psratiobluesMW_rep1, psbluesMW, by= "name2")
psbluesMW_rep2<-left_join(psratiobluesMW_rep2, psbluesMW, by= "name2")
psbluesMW_rep3<-left_join(psratiobluesMW_rep3, psbluesMW, by= "name2")
psbluesMW_rep4<-left_join(psratiobluesMW_rep4, psbluesMW, by= "name2")
psbluesMW_rep5<-left_join(psratiobluesMW_rep5, psbluesMW, by= "name2")
fwrite( psbluesMW_rep5, "./output/psbluesMW_rep5.csv", row.names = FALSE)

#combining blues for trait and waveratio for EF and MW location for rep1
psbluesEF_rep1<- read.csv("./output/psbluesEF_rep1.csv")
psbluesMW_rep1<- read.csv("./output/psbluesMW_rep1.csv")

psblues_rep1 <- bind_rows(
  "EF" = psbluesEF_rep1,
  "MW" = psbluesMW_rep1,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
psblues_rep1 <-psblues_rep1|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
ps_blues_rep1<- subset(psblues_rep1, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
ps_blues_rep1<- droplevels(ps_blues_rep1[ps_blues_rep1$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(ps_blues_rep1,"./output/ps_blues_rep1.csv", row.names = F)

#combining blues for trait and waveratio for EF and MW location for rep2
psbluesEF_rep2<- read.csv("./output/psbluesEF_rep2.csv")
psbluesMW_rep2<- read.csv("./output/psbluesMW_rep2.csv")

psblues_rep2 <- bind_rows(
  "EF" = psbluesEF_rep2,
  "MW" = psbluesMW_rep2,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
psblues_rep2 <-psblues_rep2|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
ps_blues_rep2<- subset(psblues_rep2, !is.na(taxa))%>%
  select(-name2,-Corrected_names)
ps_blues_rep2<- droplevels(ps_blues_rep2[ps_blues_rep2$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(ps_blues_rep2,"./output/ps_blues_rep2.csv", row.names = F)
#combining blues for trait and waveratio for EF and MW location for rep2
psbluesEF_rep3<- read.csv("./output/psbluesEF_rep3.csv")
psbluesMW_rep3<- read.csv("./output/psbluesMW_rep3.csv")
psblues_rep3 <- bind_rows(
  "EF" = psbluesEF_rep3,
  "MW" = psbluesMW_rep3,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
psblues_rep3 <-psblues_rep3|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
ps_blues_rep3<- subset(psblues_rep3, !is.na(taxa)) %>% select(-name2,-Corrected_names)
ps_blues_rep3<- droplevels(ps_blues_rep3[ps_blues_rep3$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(ps_blues_rep3,"./output/ps_blues_rep3.csv", row.names = F)

#combining blues for trait and waveratio for EF and MW location for rep4
psbluesEF_rep4<- read.csv("./output/psbluesEF_rep4.csv")
psbluesMW_rep4<- read.csv("./output/psbluesMW_rep4.csv")
psblues_rep4 <- bind_rows(
  "EF" = psbluesEF_rep4,
  "MW" = psbluesMW_rep4,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
psblues_rep4 <-psblues_rep4|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
ps_blues_rep4<- subset(psblues_rep4, !is.na(taxa)) %>% select(-name2,-Corrected_names)
ps_blues_rep4 <- droplevels(ps_blues_rep4[ps_blues_rep4$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(ps_blues_rep4,"./output/ps_blues_rep4.csv", row.names = F)
#combining blues for trait and waveratio for EF and MW location for rep5
psbluesEF_rep5<- read.csv("./output/psbluesEF_rep5.csv")
psbluesMW_rep5<- read.csv("./output/psbluesMW_rep5.csv")
psblues_rep5 <- bind_rows(
  "EF" = psbluesEF_rep5,
  "MW" = psbluesMW_rep5,
  .id = "env")
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
psblues_rep5 <-psblues_rep5|>
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) &
                         grepl("PI", Corrected_names), NA, Corrected_names))
ps_blues_rep5<- subset(psblues_rep5, !is.na(taxa)) %>% select(-name2,-Corrected_names)
ps_blues_rep5<- droplevels(ps_blues_rep5[ps_blues_rep5$taxa %in%rownames(kin_int), ])#filtering indv that are present both in kinship and phenotypic data
write.csv(ps_blues_rep5,"./output/ps_blues_rep5.csv", row.names = F)


#----------------2nd step --------------------
ps_blues<- read.csv("./output/ps_blues_rep1.csv")
ps_bluesEF<- ps_blues %>% filter(env== "EF")
ps_bluesMW<- ps_blues %>% filter(env== "MW")
ps_bluesEF<- ps_bluesEF%>% mutate(taxa= factor(taxa))
ps_bluesMW<- ps_bluesMW%>% mutate(taxa= factor(taxa))
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
result_ps_rep1 <-
  crossv(
    sort = sort,
    train = ps_bluesEF,
    validation = ps_bluesMW,
    mytrait = "ps",
    kin=kin,
    scheme = "EFMW")

corr_ps_EFMW<- data.frame(result_ps_rep1$ac)
fwrite(corr_ps_EFMW, "./output/corr_ps_EFMWrep1.csv")

# ---------------------plsr narea-running cross validation for MWEF loc---------------------
result_ps_rep1 <-
  crossv(
    sort = sort,
    train = ps_bluesMW,
    validation = ps_bluesEF,
    mytrait = "ps",
    kin=kin,
    scheme = "MWEF_rep1")

corr_ps_MWEF<- data.frame(result_ps_rep1$ac)
fwrite(corr_ps_MWEF, "./output/corr_ps_MWEFrep1.csv")

#----- ST Narea for rep2 --------
ps_blues<- read.csv("ps_blues_rep2.csv")
ps_bluesEF<- ps_blues %>% filter(env== "EF")
ps_bluesMW<- ps_blues %>% filter(env== "MW")
ps_bluesEF<- ps_bluesEF%>% mutate(taxa= factor(taxa))
ps_bluesMW<- ps_bluesMW%>% mutate(taxa= factor(taxa))
N_blues<- read.csv("N_blues_rep2.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)

# ---------------------ps-running cross validation for EFMW loc---------------------
result_ps_rep2 <-
  crossv(
    sort = sort,
    train = ps_bluesEF,
    validation = ps_bluesMW,
    mytrait = "ps",
    kin=kin,
    scheme = "EFMW_rep2")

corr_ps_EFMW<- data.frame(result_ps_rep2$ac)
fwrite(corr_ps_EFMW, "./output/corr_ps_EFMWrep2.csv")

# ---------------------ps-running cross validation for MWEF loc---------------------
result_ps_rep2 <-
  crossv(
    sort = sort,
    train = ps_bluesMW,
    validation = ps_bluesEF,
    mytrait = "ps",
    kin = kin,
    scheme = "MWEF_rep2")

corr_ps_MWEF<- data.frame(result_ps_rep2$ac)
fwrite(corr_ps_MWEF, "./output/corr_ps_MWEFrep2.csv")

#----- ST Narea for rep3 --------
ps_blues<- read.csv("ps_blues_rep3.csv")
ps_bluesEF<- ps_blues %>% filter(env== "EF")
ps_bluesMW<- ps_blues %>% filter(env== "MW")
ps_bluesEF<- ps_bluesEF%>% mutate(taxa= factor(taxa))
ps_bluesMW<- ps_bluesMW%>% mutate(taxa= factor(taxa))

N_blues<- read.csv("N_blues_rep3.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------ps-running cross validation for EFMW loc---------------------
result_ps_rep3 <-
  crossv(
    sort = sort,
    train = ps_bluesEF,
    validation = ps_bluesMW,
    mytrait = "ps",
    kin=kin,
    scheme = "EFMW_rep3")

corr_ps_EFMW<- data.frame(result_ps_rep3$ac)
fwrite(corr_ps_EFMW, "./output/corr_ps_EFMWrep3.csv")

# ---------------------ps-running cross validation for MWEF loc---------------------
result_ps_rep3 <-
  crossv(
    sort = sort,
    train = ps_bluesMW,
    validation = ps_bluesEF,
    mytrait = "ps",
    kin = kin,
    scheme = "MWEF_rep3")

corr_ps_MWEF<- data.frame(result_ps_rep3$ac)
fwrite(corr_ps_MWEF, "./output/corr_ps_MWEFrep3.csv")

#----- ST Narea for rep4 --------
ps_blues<- read.csv("ps_blues_rep4.csv")
ps_bluesEF<- ps_blues %>% filter(env== "EF")
ps_bluesMW<- ps_blues %>% filter(env== "MW")
ps_bluesEF<- ps_bluesEF%>% mutate(taxa= factor(taxa))
ps_bluesMW<- ps_bluesMW%>% mutate(taxa= factor(taxa))
N_blues<- read.csv("N_blues_rep4.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------ps-running cross validation for EFMW loc---------------------
result_ps_rep4 <-
  crossv(
    sort = sort,
    train = ps_bluesEF,
    validation = ps_bluesMW,
    mytrait = "ps",
    kin=kin,
    scheme = "EFMW_rep4")

corr_ps_EFMW<- data.frame(result_ps_rep4$ac)
fwrite(corr_ps_EFMW, "./output/corr_ps_EFMWrep4.csv")

# ---------------------ps-running cross validation for MWEF loc---------------------
result_ps_rep4 <-
  crossv(
    sort = sort,
    train = ps_bluesMW,
    validation = ps_bluesEF,
    mytrait = "ps",
    kin = kin,
    scheme = "MWEF_rep4")

corr_ps_MWEF<- data.frame(result_ps_rep4$ac)
fwrite(corr_ps_MWEF, "./output/corr_ps_MWEFrep4.csv")

#----- ST Narea for rep5 --------
ps_blues<- read.csv("ps_blues_rep5.csv")
ps_bluesEF<- ps_blues %>% filter(env== "EF")
ps_bluesMW<- ps_blues %>% filter(env== "MW")
ps_bluesEF<- ps_bluesEF%>% mutate(taxa= factor(taxa))
ps_bluesMW<- ps_bluesMW%>% mutate(taxa= factor(taxa))

N_blues<- read.csv("N_blues_rep5.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))


sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)
# ---------------------ps-running cross validation for EFMW loc---------------------
result_ps_rep5 <-
  crossv(
    sort = sort,
    train = ps_bluesEF,
    validation = ps_bluesMW,
    mytrait = "ps",
    kin=kin,
    scheme = "EFMW_rep5")

corr_ps_EFMW<- data.frame(result_ps_rep5$ac)
fwrite(corr_ps_EFMW, "./output/corr_ps_EFMWrep5.csv")

# ---------------------ps-running cross validation for MWEF loc---------------------
result_ps_rep5 <-
  crossv(
    sort = sort,
    train = ps_bluesMW,
    validation = ps_bluesEF,
    mytrait = "ps",
    kin = kin,
    scheme = "MWEF_rep5")

corr_ps_MWEF<- data.frame(result_ps_rep5$ac)
fwrite(corr_ps_MWEF, "./output/corr_ps_MWEFrep5.csv")

# # obtaining synthetic trait table and value in each rep
# ps  ratio_combined <- read.csv("./data/psratio_combined.csv", header = TRUE) %>% select(-3)
# psra tio_combined$ST <- rep(c("ST1", "ST2", "ST3"), times = 5)
# psrat  io_combined$Rep <- rep(1:5, each = 3)
# psratio_combined<- psratio_combined %>% select(-wave_1, -wave_2)
# psratio_combined$coh2_Rep1<- c(0.63,0.65,0.63,0.61,0.55,0.49,0.55,0.62, 0.52,0.63,0.64,0.60,0.55,0.62,0.33)
   # psratio_combined$coh2_Rep3<-c(0.43,0.45,0.40,0.42,0.53,0.46,0.51,0.57,0.48,0.43,0.39,0.50,0.52,0.58,0.26)
# psratio_combined$coh2_Rep4<- c(0.49,0.53,0.44,0.49,0.39,0.37,0.38,0.49,0.36,0.49,0.49,0.49,0.38,0.53,NA)
# psratio_combined$coh2_Rep5<- c(0.49,0.40,0.76,0.46,NA,0.48,0.45,0.34,0.45,0.34,NA,NA,0.59,0.53,0.52)
# library(gtable)
# gt_table <- gt(psratio_combined)
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
# htmltools::save_html(gt_table, file = "coh2_all_ps_rep.html")
# webshot("coh2_all_ps_rep.html", "./figure/coh2_all_ps_rep.psg")
#
#

#selecting synthetic trait that has 0 coh2 value from each replication
#Rep1
# Select the one wave ratio with the lowest coh2
ps1 <- fread("./data/ps_breakdown1.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_1))
wave <-phenotypes |>
  clean_names()
selected <- ps1 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pswave_lowcoh2<- wave %>% select(1:9, plsr_sla_sorghum,wave_1954,wave_2311)
wave_1954_wave_2311 <- paste("wave_1954", "wave_2311", sep = "_")
pswave_lowcoh2[[wave_1954_wave_2311]] <- pswave_lowcoh2[["wave_1954"]] / pswave_lowcoh2[["wave_2311"]]
pswave_lowcoh2<- pswave_lowcoh2 %>% select(-wave_1954,-wave_2311)
write.csv(pswave_lowcoh2, "./output/pswave_rep1_lowcoh2.csv", row.names = F)

#Rep2
ps2 <- fread("./data/ps_breakdown2.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_2))
wave <-phenotypes |>
  clean_names()

selected <- ps2 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pswave_lowcoh2<- wave %>% select(1:9, plsr_sla_sorghum,wave_722,wave_1782)
wave_722_wave_1782 <- paste("wave_722", "wave_1782", sep = "_")
pswave_lowcoh2[[wave_722_wave_1782]] <- pswave_lowcoh2[["wave_722"]] / pswave_lowcoh2[["wave_1782"]]
pswave_lowcoh2<- pswave_lowcoh2 %>% select(-wave_722,-wave_1782)
write.csv(pswave_lowcoh2, "./output/pswave_rep2_lowcoh2.csv", row.names = F)

#Rep3
ps3 <- fread("./data/ps_breakdown3.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_3))
wave <-phenotypes |>
  clean_names()

selected <- ps3 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pswave_lowcoh2<- wave %>% select(1:9,plsr_sla_sorghum,wave_467,wave_2178)
wave_467_wave_2178 <- paste("wave_467", "wave_2178", sep = "_")
pswave_lowcoh2[[wave_467_wave_2178]] <- pswave_lowcoh2[["wave_467"]] / pswave_lowcoh2[["wave_2178"]]
pswave_lowcoh2<- pswave_lowcoh2 %>% select(-wave_467,-wave_2178)
write.csv(pswave_lowcoh2, "./output/pswave_rep3_lowcoh2.csv", row.names = F)


#Rep4
ps4 <- fread("./data/ps_breakdown4.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_4))
wave <-phenotypes |>
  clean_names()

selected <- ps4 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pswave_lowcoh2<- wave %>% select(1:9, plsr_sla_sorghum,wave_705,wave_1825)
wave_705_wave_1825 <- paste("wave_705", "wave_1825", sep = "_")
pswave_lowcoh2[[wave_705_wave_1825]] <- pswave_lowcoh2[["wave_705"]] / pswave_lowcoh2[["wave_1825"]]
pswave_lowcoh2<- pswave_lowcoh2 %>% select(-wave_705,-wave_1825)
write.csv(pswave_lowcoh2, "./output/pswave_rep4_lowcoh2.csv", row.names = F)

#Rep5
ps5 <- fread("./data/ps_breakdown5.csv", data.table = F)
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(!(Name2 %in% indv_sample$Rep_5))
wave <-phenotypes |>
  clean_names()

selected <- ps5 |>
  filter(coh2 > 0) %>%
  slice_min(coh2, n = 1) |>
  mutate(across(where(is.double), \(x) round(x, 4)))
pswave_lowcoh2<- wave %>% select(1:9, plsr_sla_sorghum,wave_654,wave_1088)
wave_654_wave_1088 <- paste("wave_654", "wave_1088", sep = "_")
pswave_lowcoh2[[wave_654_wave_1088]] <- pswave_lowcoh2[["wave_654"]] / pswave_lowcoh2[["wave_1088"]]
pswave_lowcoh2<- pswave_lowcoh2 %>% select(-wave_654,-wave_1088)
write.csv(pswave_lowcoh2, "./output/pswave_rep5_lowcoh2.csv", row.names = F)



#calculating blues for the selected low coheritable synthetic trait
pswave_lowcoh2<- read.csv("./output/pswave_rep1_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesMW1 <- bind_rows(psratio_bluesMW1, temp)
psratio_bluesMW1<- psratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesMW1, "./output/psratio_bluesMW1_rep1.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesEF1 <- bind_rows(psratio_bluesEF1, temp)

psratio_bluesEF1<- psratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesEF1, "./output/psratio_bluesEF1_rep1.csv", row.names = FALSE)
psratio_bluesEF1<- read.csv("./output/psratio_bluesEF1_rep1.csv")

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
  ps = round(psbluesEF1$predictions$pvals$predicted.value, 3))
psbluesEF1<-left_join(psratio_bluesEF1, psbluesEF1, by= "name2")
fwrite( psbluesEF1, "./output/psbluesEF1_rep1.csv", row.names = FALSE)
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
  ps = round(psbluesMW1$predictions$pvals$predicted.value, 3))
psbluesMW1<-left_join(psratio_bluesMW1, psbluesMW1, by= "name2")
fwrite( psbluesMW1, "./output/psbluesMW1_rep1.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
psbluesEF1<- read.csv("./output/psbluesEF1_rep1.csv")
psbluesMW1<- read.csv("./output/psbluesMW1_rep1.csv")

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
ps_blues1<- droplevels(ps_blues1[ps_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(ps_blues1,"./output/ps_blues1_rep1.csv", row.names = F)


##calculating blues for the selected low coheritable synthetic trait
pswave_lowcoh2<- read.csv("./output/pswave_rep2_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesMW1 <- bind_rows(psratio_bluesMW1, temp)
psratio_bluesMW1<- psratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesMW1, "./output/psratio_bluesMW1_rep2.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesEF1 <- bind_rows(psratio_bluesEF1, temp)

psratio_bluesEF1<- psratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesEF1, "./output/psratio_bluesEF1_rep2.csv", row.names = FALSE)
psratio_bluesEF1<- read.csv("./output/psratio_bluesEF1_rep2.csv")

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
  ps = round(psbluesEF1$predictions$pvals$predicted.value, 3))
psbluesEF1<-left_join(psratio_bluesEF1, psbluesEF1, by= "name2")
fwrite( psbluesEF1, "./output/psbluesEF1_rep2.csv", row.names = FALSE)
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
  ps = round(psbluesMW1$predictions$pvals$predicted.value, 3))
psbluesMW1<-left_join(psratio_bluesMW1, psbluesMW1, by= "name2")
fwrite( psbluesMW1, "./output/psbluesMW1_rep2.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
psbluesEF1<- read.csv("./output/psbluesEF1_rep2.csv")
psbluesMW1<- read.csv("./output/psbluesMW1_rep2.csv")

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
ps_blues1<- droplevels(ps_blues1[ps_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(ps_blues1,"./output/ps_blues1_rep2.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
pswave_lowcoh2<- read.csv("./output/pswave_rep3_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesMW1 <- bind_rows(psratio_bluesMW1, temp)
psratio_bluesMW1<- psratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesMW1, "./output/psratio_bluesMW1_rep3.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesEF1 <- bind_rows(psratio_bluesEF1, temp)

psratio_bluesEF1<- psratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesEF1, "./output/psratio_bluesEF1_rep3.csv", row.names = FALSE)
psratio_bluesEF1<- read.csv("./output/psratio_bluesEF1_rep3.csv")

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
  ps = round(psbluesEF1$predictions$pvals$predicted.value, 3))
psbluesEF1<-left_join(psratio_bluesEF1, psbluesEF1, by= "name2")
fwrite( psbluesEF1, "./output/psbluesEF1_rep3.csv", row.names = FALSE)
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
  ps = round(psbluesMW1$predictions$pvals$predicted.value, 3))
psbluesMW1<-left_join(psratio_bluesMW1, psbluesMW1, by= "name2")
fwrite( psbluesMW1, "./output/psbluesMW1_rep3.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
psbluesEF1<- read.csv("./output/psbluesEF1_rep3.csv")
psbluesMW1<- read.csv("./output/psbluesMW1_rep3.csv")

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
ps_blues1<- droplevels(ps_blues1[ps_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(ps_blues1,"./output/ps_blues1_rep3.csv", row.names = F)



###calculating blues for the selected low coheritable synthetic trait
pswave_lowcoh2<- read.csv("./output/pswave_rep4_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesMW1 <- bind_rows(psratio_bluesMW1, temp)
psratio_bluesMW1<- psratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesMW1, "./output/psratio_bluesMW1_rep4.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesEF1 <- bind_rows(psratio_bluesEF1, temp)

psratio_bluesEF1<- psratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesEF1, "./output/psratio_bluesEF1_rep4.csv", row.names = FALSE)
psratio_bluesEF1<- read.csv("./output/psratio_bluesEF1_rep4.csv")

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
  ps = round(psbluesEF1$predictions$pvals$predicted.value, 3))
psbluesEF1<-left_join(psratio_bluesEF1, psbluesEF1, by= "name2")
fwrite( psbluesEF1, "./output/psbluesEF1_rep4.csv", row.names = FALSE)
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
  ps = round(psbluesMW1$predictions$pvals$predicted.value, 3))
psbluesMW1<-left_join(psratio_bluesMW1, psbluesMW1, by= "name2")
fwrite( psbluesMW1, "./output/psbluesMW1_rep4.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
psbluesEF1<- read.csv("./output/psbluesEF1_rep4.csv")
psbluesMW1<- read.csv("./output/psbluesMW1_rep4.csv")

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
ps_blues1<- droplevels(ps_blues1[ps_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(ps_blues1,"./output/ps_blues1_rep4.csv", row.names = F)


###calculating blues for the selected low coheritable synthetic trait
pswave_lowcoh2<- read.csv("./output/pswave_rep5_lowcoh2.csv")
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

#calculating blues for each wavelength for location MW
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesMW1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "MW", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesMW1 <- bind_rows(psratio_bluesMW1, temp)
psratio_bluesMW1<- psratio_bluesMW1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesMW1, "./output/psratio_bluesMW1_rep5.csv", row.names = FALSE)

#calculating blues for each wavelength
variables <- colnames(pswave_lowcoh2)[11]
models <- vector("list",length(variables))
psratio_bluesEF1 <-data.frame()
# ---------------------fitting model---------------------
model <- asreml(
  fixed = get(variables) ~set + name2,
  random = ~block,
  #residual =  ~ar1(range):ar1(row),
  data = pswave_lowcoh2,subset= loc== "EF", na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
model <- update.asreml(model)
#---------------------storing prediction---------------------
temp <- model$predictions$pvals[, 1:2]
temp$wave <- variables
psratio_bluesEF1 <- bind_rows(psratio_bluesEF1, temp)

psratio_bluesEF1<- psratio_bluesEF1 %>% pivot_wider(names_from = wave, values_from= predicted.value)
fwrite( psratio_bluesEF1, "./output/psratio_bluesEF1_rep5.csv", row.names = FALSE)
psratio_bluesEF1<- read.csv("./output/psratio_bluesEF1_rep5.csv")

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
  ps = round(psbluesEF1$predictions$pvals$predicted.value, 3))
psbluesEF1<-left_join(psratio_bluesEF1, psbluesEF1, by= "name2")
fwrite( psbluesEF1, "./output/psbluesEF1_rep5.csv", row.names = FALSE)
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
  ps = round(psbluesMW1$predictions$pvals$predicted.value, 3))
psbluesMW1<-left_join(psratio_bluesMW1, psbluesMW1, by= "name2")
fwrite( psbluesMW1, "./output/psbluesMW1_rep5.csv", row.names = FALSE)
#combining blues for trait and waveratio for Ef and MW location
psbluesEF1<- read.csv("./output/psbluesEF1_rep5.csv")
psbluesMW1<- read.csv("./output/psbluesMW1_rep5.csv")

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
ps_blues1<- droplevels(ps_blues1[ps_blues1$taxa %in%rownames(kin), ])#filtering indv that are present both in kin and phenotypic data
write.csv(ps_blues1,"./output/ps_blues1_rep5.csv", row.names = F)
