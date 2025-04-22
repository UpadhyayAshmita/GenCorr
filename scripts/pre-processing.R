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

N <- fread("./data/narea_breakdown.csv", data.table = F)
S<- fread("./data/sla_breakdown.csv",data.table=F)
pn<- fread("./data/pn_breakdown.csv",data.table=F)
ps<- fread("./data/ps_breakdown.csv", data.table = F)

wave <- fread("./data/phenotypes.csv", data.table = F) |>
  clean_names()

sum(is.na(N$corg))#396
sum(N$corg< -1 | N$corg > 1, na.rm = T) #634
sum(N$corgblup< -1 | N$corgblup > 1, na.rm = T)  #17362
#1800   1827 narea  0.022461370 0.6179698 0.0008999314  0.95246255 -0.2692682
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

# lowest_ntraits <- N %>%
#   filter(coh2 >= 0) %>%
#   slice_min(order_by = coh2, n = 1)

# # Calculate the wave ratio for the single row in 'sel'
# wave_ratio_name <- paste0("wave$", lowest_ntraits$wave_1, "_", lowest_ntraits$wave_2)
# assign(wave_ratio_name, wave[[lowest_ntraits$wave_1]] / wave[[lowest_ntraits$wave_2]])
#
# lowest_ntraits$rowname <- paste0(lowest_ntraits$wave_1, "_", lowest_ntraits$wave_2)
# Nratio_low <- wave |> select(wave_524_wave_681)

sel <- N |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))


for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Nratios <- wave |> select(wave_1675_wave_1723:wave_1841_wave_2221)
Nratios <- t(Nratios)
colnames(Nratios) <- wave$plot_id
#distance
dist <- get_dist(Nratios, method = "pearson")
# fviz_dist(dist)
hcN <- hclust(dist, method = "average")
# fviz_nbclust(Nratios, FUN = hcut, method = "silhouette")
#
# fviz_dend(hcN, k = 2,                 # Cut in 3 groups
#           cex = 0.5,                 # label size
#           main = "Nitrogen Area",
#           color_labels_by_k = TRUE,  # color labels by groups
#           ggtheme = theme_gray()     # Change theme
# )
#
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
  slice_max(coh2, n=3) |> select(coh2, wave_1, wave_2)
print(Nratios)

Nratios<- as.data.frame(Nratios)
Nratios<- read.csv("./data/Nratios.csv")

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
write.csv(Nratio_transform, "./data/Nratio_transform.csv", row.names=F)
Nratio_transform<- read.csv("./data/Nratio_transform.csv")
com_col <- colnames(wave)[colnames(wave) %in% Nratio_transform$wave_sel]
Nwave_pheno<- wave %>% select(1:10, com_col)
write.csv(Nwave_pheno,"./data/Nwave_pheno.csv", row.names = F)
#plotting
l <- Nratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap <-N |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                           wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Heatmap-Narea") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap,"Heatmap_narea.jpeg")
#sla
sum(is.na(S$corg))#396
sum(S$corg< -1 | S$corg > 1, na.rm = T) #13141
sum(S$corgblup< -1 | S$corgblup > 1, na.rm = T)  #35278
S <- S |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(abs(S$corg - S$corgblup)>0.2, na.rm = T)
# 176644
sum(is.na(S$corg))

ntraits <- ceiling(sum(!is.na(S$coh2))*0.01)

sel <- S |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))


for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
Sratios <- wave |> select(wave_1719_wave_1680:wave_1673_wave_1370)
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
  slice_max(coh2, n=3) |> select(coh2, wave_1, wave_2)
print(Sratios)
Sratios<- as.data.frame(Sratios)
Sratios<- read.csv("./data/Sratios.csv")
#absolute difference between wave_1 and wave_2
Sratios <- Sratios %>%
  mutate(diff = abs(wave_1 - wave_2))
Sratio_sel <- Sratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

Sratio_transform <- Sratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
write.csv(Sratio_transform, "./data/Sratio_transform.csv", row.names = F)
com_col <- colnames(wave)[colnames(wave) %in% Sratio_transform$wave_sel]
Swave_pheno<- wave %>% select(1:9, sla, com_col)
write.csv(Swave_pheno,"./data/Swave_pheno.csv", row.names = F)

#plotting
l <- Sratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap <-S |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                      wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Heatmap-SLA") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap,"Heatmap_sla.jpeg")

#plsr-narea
sum(is.na(pn$corg))#396
sum(pn$corg< -1 | pn$corg > 1, na.rm = T) #13141
sum(pn$corgblup< -1 | pn$corgblup > 1, na.rm = T)  #35278
pn <- pn |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(abs(pn$corg - pn$corgblup)>0.2, na.rm = T)
# 176644
sum(is.na(pn$corg))

ntraits <- ceiling(sum(!is.na(pn$coh2))*0.01)

sel <- pn |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))


for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
pnratios <- wave |> select(wave_1076_wave_389:wave_1286_wave_2304)
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
  slice_max(coh2, n=3) |> select(coh2, wave_1, wave_2)
print(pnratios)
pnratios<- as.data.frame(pnratios)
write.csv(pnratios, "./data/pnratios.csv", row.names= F)
#absolute difference between wave_1 and wave_2
pnratios <- pnratios %>%
  mutate(diff = abs(wave_1 - wave_2))
pnratio_sel <- pnratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()
pnratio_transform <- pnratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
#write.csv(pnratio_transform, "./data/pnratio_transform.csv", row.names=F)
#pnratio_transform<- read.csv("./data/pnratio_transform.csv")
com_col <- colnames(wave)[colnames(wave) %in% pnratio_transform$wave_sel]
pnwave_pheno<- wave %>% select(1:9,fs_plsr_narea,all_of(com_col))
write.csv(pnwave_pheno,"./data/pnwave_pheno.csv", row.names = F)

#plotting
l <- pnratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap <-pn |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                      wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "Heatmap-PLSR Narea") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap,"Heatmap_plsrnarea.jpeg")
#plsr-sla
sum(is.na(ps$corg))#396
sum(ps$corg< -1 | ps$corg > 1, na.rm = T) #5152
sum(ps$corgblup< -1 | ps$corgblup > 1, na.rm = T) #56038
ps <- ps |>
  mutate(across(coh2:vartrait, ~replace(., which(corgblup < -1 |
                                                   corgblup > 1 |
                                                   corg < -1 |
                                                   corg > 1|
                                                   (abs(corg - corgblup) > 0.2)), NA)))

sum(abs(ps$corg - ps$corgblup)>0.2, na.rm = T)
sum(is.na(ps$corg))
#240146
ntraits <- ceiling(sum(!is.na(ps$coh2))*0.01)

sel <- ps |> slice_max(coh2, n=ntraits) |>
  mutate(across(where(is.double), \(x) round(x, 4)))


for(i in 1:ntraits){
  eval(parse(text = paste(paste0("wave$", sel$wave_1[i], "_", sel$wave_2[i]),
                          "<- wave[,sel$wave_1[i]] / wave[,sel$wave_2[i]]")))
}

sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
psratios <- wave |> select(wave_1728_wave_1661:wave_381_wave_1262)
psratios <- t(psratios)
colnames(psratios) <- wave$plot_id
#distance
dist <- get_dist(psratios, method = "pearson")
# fviz_dist(dist)
hcpn <- hclust(dist, method = "average")
sub_grp <- cutree(hcpn, k = 3)
table(sub_grp)
psratios <- psratios |>
  as.data.frame() |>
  rownames_to_column() |>  left_join(sel |> select(rowname, coh2)) |>
  arrange(coh2) |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", rowname))),
         wave_2 = as.numeric(gsub(".*_","", rowname))) |>
  left_join(data.frame(rowname = names(sub_grp), group = sub_grp)) |>
  group_by(group) |> filter(!duplicated(wave_1)) |> filter(!duplicated(wave_2)) |>
  slice_max(coh2, n=3) |> select(coh2, wave_1, wave_2)
print(psratios)
psratios<- as.data.frame(psratios)
#absolute difference between wave_1 and wave_2
psratios <- psratios %>%
  mutate(diff = abs(wave_1 - wave_2))
psratio_sel <- psratios %>%
  group_by(group) %>%
  slice(which.max(diff)) %>%
  ungroup()

psratio_transform <- psratio_sel %>%
  mutate(across(c(wave_1, wave_2), ~paste("wave", ., sep = "_"))) %>%
  mutate(wave_sel= paste(wave_1, wave_2, sep = "_")) %>% select(-group,-coh2,-diff)
#write.csv(psratio_transform, "./data/psratio_transform.csv", row.names=F)
#psratio_transform<- read.csv("./data/psratio_transform.csv")
com_col <- colnames(wave)[colnames(wave) %in% psratio_transform$wave_sel]
pswave_pheno<- wave %>% select(1:9,plsr_sla_sorghum,all_of(com_col))
write.csv(pswave_pheno,"./data/pswave_pheno.csv", row.names = F)
pswave_pheno<- read.csv("./data/pswave_pheno.csv")

#plotting
l <- psratio_transform |>
  mutate(wave_1 = as.numeric(gsub("_.*", "",gsub("wave_","", wave_1))),
         wave_2 = as.numeric(gsub(".*_","", wave_2))) |> ungroup() |>
  select(wave_1, wave_2) |> mutate(label = "*") |> as.data.frame()
heatmap <-ps |> mutate(wave_1 = as.numeric(gsub("wave_","", wave_1)),
                      wave_2 = as.numeric(gsub("wave_","", wave_2))) |>
  ggplot(aes(x = wave_1, y = wave_2)) +
  geom_tile(aes(fill = coh2)) +
  scale_fill_gradient(low = "white", high = "red",na.value="black") +
  labs(x = "Wave2", y = "Wave1", title = "HeatmapPLSR Sla") +
  # scale_x_continuous(limits = c(350, 2500)) +
  # scale_y_continuous(limits = c(350, 2500)) +
  geom_label(data = l, aes(label = label))+
  theme_minimal()

ggsave(plot=heatmap,"Heatmap_plsrsla.jpeg")
