#loading packages
library(data.table)
library(tidyverse)
library(cvTools)
library(asreml)
source("./scripts/aux_function.R")
#----------------2nd step --------------------
N_blues<- read.csv("N_blues.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))
sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                    reps = 20, seed = 123)
kin<- fread('kin_additive.txt', data.table= F)
rownames(kin) <- colnames(kin)
kin<- as.matrix(kin)

#----- ST Narea --------
# ---------------------narea-running cross validation for EFMW loc---------------------
result_N <-
  crossv(
    sort = sort,
    train = N_bluesEF,
    validation = N_bluesMW,
    mytrait = "narea",
    kin=kin,
    scheme = "EFMW")

corr_N_EFMW<- data.frame(result_N$ac)
fwrite(corr_N_EFMW, "./output/corr_N_EFMW.csv")

# ---------------------narea-running cross validation for MWEF loc---------------------
result_N <-
  crossv(
    sort = sort,
    train = N_bluesMW,
    validation = N_bluesEF,
    mytrait = "narea",
    kin=kin,
    scheme = "MWEF")

corr_N_MWEF<- data.frame(result_N$ac)
fwrite(corr_N_MWEF, "./output/narea/corr_N_MWEF.csv")
#---- ST SLA ----------
S_blues<- read.csv("./data/S_blues.csv")
S_bluesEF<- S_blues %>% filter(env== "EF")
S_bluesMW<- S_blues %>% filter(env== "MW")
S_bluesEF<- S_bluesEF%>% mutate(taxa= factor(taxa))
S_bluesMW<- S_bluesMW%>% mutate(taxa= factor(taxa))
result_S <-
  crossv(
    sort = sort,
    train = S_bluesEF,
    validation = S_bluesMW,
    mytrait = "sla",
    kin=kin,
    scheme = "EFMW")

corr_S_EFMW<- data.frame(result_S$ac)
fwrite(corr_S_EFMW, "./output/corr_S_EFMW.csv")

#training mw and val EF
result_S <-
  crossv(
    sort = sort,
    train = S_bluesMW,
    validation = S_bluesEF,
    mytrait = "sla",
    kin=kin,
    scheme = "MWEF")

corr_S_MWEF<- data.frame(result_S$ac)
fwrite(corr_S_MWEF, "./output/corr_S_MWEF.csv")

# ---ST pn(plsr-narea)---------------------
pn_blues<- read.csv("./data/pn_blues.csv")
pn_bluesEF<- pn_blues %>% filter(env== "EF")
pn_bluesMW<- pn_blues %>% filter(env== "MW")
pn_bluesEF<- pn_bluesEF%>% mutate(taxa= factor(taxa))
pn_bluesMW<- pn_bluesMW%>% mutate(taxa= factor(taxa))
result_pn <-
  crossv(
    sort = sort,
    train = pn_bluesEF,
    validation = pn_bluesMW,
    mytrait = "plsr_narea",
    kin=kin,
    scheme = "EFMW")

corr_pn_EFMW<- data.frame(result_pn$ac)
fwrite(corr_pn_EFMW, "./output/corr_pn_EFMW.csv")

#training mw and val EF
result_pn <-
  crossv(
    sort = sort,
    train = pn_bluesMW,
    validation = pn_bluesEF,
    mytrait = "plsr_narea",
    kin=kin,
    scheme = "MWEF")

corr_pn_MWEF<- data.frame(result_pn$ac)
fwrite(corr_pn_MWEF, "./output/corr_pn_MWEF.csv")

# ---ST ps(plsr-sla)---------------------
ps_blues<- read.csv("./data/ps_blues.csv")
ps_bluesEF<- ps_blues %>% filter(env== "EF")
ps_bluesMW<- ps_blues %>% filter(env== "MW")
ps_bluesEF<- ps_bluesEF%>% mutate(taxa= factor(taxa))
ps_bluesMW<- ps_bluesMW%>% mutate(taxa= factor(taxa))
result_ps <-
  crossv(
    sort = sort,
    train = ps_bluesEF,
    validation = ps_bluesMW,
    mytrait = "plsr_sla",
    kin=kin,
    scheme = "EFMW")

corr_ps_EFMW<- data.frame(result_ps$ac)
fwrite(corr_ps_EFMW, "./output/corr_ps_EFMW.csv")

#training mw and val EF
result_ps <-
  crossv(
    sort = sort,
    train = ps_bluesMW,
    validation = ps_bluesEF,
    mytrait = "plsr_sla",
    kin=kin,
    scheme = "MWEF")

corr_ps_MWEF<- data.frame(result_ps$ac)
fwrite(corr_ps_MWEF, "./output/corr_ps_MWEF.csv")

#coh2 0 wave ratio was selected to use as a test synthetic trait so fitting Cv1, Cv2 for EFMW and MWEF scenario

#----------------CV2 --------------------
#------ Narea------
#------ Multi-trait T & Narea for loc EF------
acNT<-c()
fwrite(data.frame("taxa","GEBV"), "NT.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rNT<-list()
  for(i in 1:5){
    test<-N_bluesEF1
    test[test$taxa %in% sort[[j]][[i]],"narea"] <- NA
    cat(i, j,"\n")
    modelNT <- asreml(
      fixed = cbind(narea, wave_524_wave_681) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rNT[[i]] <- modelNT$predictions$pvals %>%
      filter(trait == "narea" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raN<- Reduce(rbind, rNT)
  raN <- raN %>% left_join(N_bluesMW1[,c("taxa", "narea")])
  fwrite(raN, "NT.txt", sep = "\t", append = T, col.names = F)
  acNT[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite(as.matrix(acNT), "acNT.txt", sep = "\t", col.names = F)

#------ Multi-trait S3------
acST<-c()
fwrite(data.frame("taxa","GEBV"), "ST.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rST<-list()
  for(i in 1:5){
    test<-S_bluesEF1
    test[test$taxa %in% sort[[j]][[i]],"sla"] <- NA

    modelST <- asreml(
      fixed = cbind(sla, wave_1998_wave_2482) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rST[[i]] <- modelST$predictions$pvals %>%
      filter(trait == "sla" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raS<- Reduce(rbind, rST)
  raS <- raS %>% left_join(S_bluesMW1[,c("taxa", "sla")])
  fwrite(raS, "ST.txt", sep = "\t", append = T, col.names = F)
  acST[j]<-cor(raS[,2], raS[,3], use = "complete.obs")
}
fwrite(as.matrix(acST), "acST.txt", sep = "\t", col.names = F)
#------ plsr_sla------
#------ Multi-trait plsr_sla for loc EF------

acpsT<-c()
fwrite(data.frame("taxa","GEBV"), "psT.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rpsT<-list()
  for(i in 1:5){
    test<-ps_bluesEF1
    test[test$taxa %in% sort[[j]][[i]],"plsr_sla"] <- NA
    cat(i, j,"\n")
    modelpsT <- asreml(
      fixed = cbind(plsr_sla, wave_1434_wave_2258) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rpsT[[i]] <- modelpsT$predictions$pvals %>%
      filter(trait == "plsr_sla" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raps<- Reduce(rbind, rpsT)
  raps <- raps %>% left_join(ps_bluesMW1[,c("taxa", "plsr_sla")])
  fwrite(raps, "psT.txt", sep = "\t", append = T, col.names = F)
  acpsT[j]<-cor(raps[,2], raps[,3], use = "complete.obs")
}
fwrite(as.matrix(acpsT), "acpsT.txt", sep = "\t", col.names = F)

#------ plsr_narea------
#------ Multi-trait test & plsr_narea for loc EF------

acpnT<-c()
fwrite(data.frame("taxa","GEBV"), "pnT.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rpnT<-list()
  for(i in 1:5){
    test<-pn_bluesEF1
    test[test$taxa %in% sort[[j]][[i]],"plsr_narea"] <- NA
    cat(i, j,"\n")
    modelpnT <- asreml(
      fixed = cbind(plsr_narea, wave_2201_wave_409) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rpnT[[i]] <- modelpnT$predictions$pvals %>%
      filter(trait == "plsr_narea" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  rapn<- Reduce(rbind, rpnT)
  rapn <- rapn %>% left_join(pn_bluesMW1[,c("taxa", "plsr_narea")])
  fwrite(rapn, "pnT.txt", sep = "\t", append = T, col.names = F)
  acpnT[j]<-cor(rapn[,2], rapn[,3], use = "complete.obs")
}
fwrite(as.matrix(acpnT), "acpnT.txt", sep = "\t", col.names = F)
#Cv1 for test synthetic trait with four biological traits
#------ Multi-trait N2------
acNT<-c()
fwrite(data.frame("taxa","GEBV"), "NT_CV1.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rNT<-list()
  for(i in 1:5){
    test<-N_bluesEF1
    test[test$taxa %in% sort[[j]][[i]],"narea"] <- NA
    test[test$taxa %in% sort[[j]][[i]],"wave_524_wave_681"] <- NA
    cat(i, j,"\n")
    modelNT <- asreml(
      fixed = cbind(narea, wave_524_wave_681) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rNT[[i]] <- modelNT$predictions$pvals %>%
      filter(trait == "narea" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raN<- Reduce(rbind, rNT)
  fwrite(raN, "NT_CV1.txt", sep = "\t", append = T, col.names = F)
  raN <- raN %>% left_join(N_bluesMW1[,c("taxa", "narea")])
  acNT[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite(as.matrix(acNT), "acNT_CV1.txt", sep = "\t", col.names = F)

#------ Multi-trait S2------
acST<-c()
fwrite(data.frame("taxa","GEBV"), "ST_CV1.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rST<-list()
  for(i in 1:5){
    test<-S_bluesEF1
    test[test$taxa %in% sort[[j]][[i]],"sla"] <- NA
    test[test$taxa %in% sort[[j]][[i]],"wave_1998_wave_2482"] <- NA
    modelST <- asreml(
      fixed = cbind(sla, wave_1998_wave_2482) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rST[[i]] <- modelST$predictions$pvals %>%
      filter(trait == "sla" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raS<- Reduce(rbind, rST)
  fwrite(raS, "ST_CV1.txt", sep = "\t", append = T, col.names = F)
  raS <- raS %>% left_join(S_bluesMW1[,c("taxa", "sla")])
  acST[j]<-cor(raS[,2], raS[,3], use = "complete.obs")
}
fwrite(as.matrix(acST), "acST_CV1.txt", sep = "\t", col.names = F)
#------ Multi-trait PN2------
acpnT<-c()
fwrite(data.frame("taxa","GEBV"), "pnT_CV1.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rpnT<-list()
  for(i in 1:5){
    test<-pn_bluesEF1
    test[test$taxa %in% sort[[j]][[i]],"plsr_narea"] <- NA
    test[test$taxa %in% sort[[j]][[i]],"wave_2201_wave_409"] <- NA
    modelpnT <- asreml(
      fixed = cbind(plsr_narea, wave_2201_wave_409) ~ trait,
      random =  ~  corgh(trait):vm(taxa,source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rpnT[[i]] <- modelpnT$predictions$pvals %>%
      filter(trait == "plsr_narea" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raN<- Reduce(rbind, rpnT)
  fwrite(raN, "pnT_CV1.txt", sep = "\t", append = T, col.names = F)
  raN <- raN %>% left_join(pn_bluesMW1[,c("taxa", "plsr_narea")])
  acpnT[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite(as.matrix(acpnT), "acpnT_CV1.txt", sep = "\t", col.names = F)

#------ Multi-trait ps with test synthetic trait ------
acpsT<-c()
fwrite(data.frame("taxa","GEBV"), "psT_CV1.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rpsT<-list()
  for(i in 1:5){
    test<-ps_bluesEF1
    test[test$taxa %in% sort[[j]][[i]],"plsr_sla"] <- NA
    test[test$taxa %in% sort[[j]][[i]],"wave_1434_wave_2258"] <- NA
    modelpsT <- asreml(
      fixed = cbind(plsr_sla, wave_1434_wave_2258) ~ trait,
      random =  ~  corgh(trait):vm(taxa,source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rpsT[[i]] <- modelpsT$predictions$pvals %>%
      filter(trait == "plsr_sla" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raS<- Reduce(rbind, rpsT)
  fwrite(raS, "psT_CV1.txt", sep = "\t", append = T, col.names = F)
  raS <- raS %>% left_join(ps_bluesMW1[,c("taxa", "plsr_sla")])
  acpsT[j]<-cor(raS[,2], raS[,3], use = "complete.obs")
}
fwrite(as.matrix(acpsT), "acpsT_CV1.txt", sep = "\t", col.names = F)

#multitrait model with test synthetic trait for MWEF scenario

#----------------CV2 --------------------
#------ Narea------
#------ Multi-trait model of synthetic test trait with Narea for loc mw------
acNT_mwef<-c()
fwrite(data.frame("taxa","GEBV"), "NT_mwef.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rNT_mwef<-list()
  for(i in 1:5){
    test<-N_bluesMW1
    test[test$taxa %in% sort[[j]][[i]],"narea"] <- NA
    cat(i, j,"\n")
    modelNT_mwef <- asreml(
      fixed = cbind(narea, wave_524_wave_681) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rNT_mwef[[i]] <- modelNT_mwef$predictions$pvals %>%
      filter(trait == "narea" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raN<- Reduce(rbind, rNT_mwef)
  raN <- raN %>% left_join(N_bluesEF1[,c("taxa", "narea")])
  fwrite(raN, "NT_mwef.txt", sep = "\t", append = T, col.names = F)
  acNT_mwef[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite(as.matrix(acNT_mwef), "acNT_mwef.txt", sep = "\t", col.names = F)
#------ Multi-trait model with synthetic test trait and sla------
acST<-c()
fwrite(data.frame("taxa","GEBV"), "ST_mwef.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rST<-list()
  for(i in 1:5){
    test<-S_bluesMW1
    test[test$taxa %in% sort[[j]][[i]],"sla"] <- NA

    modelST <- asreml(
      fixed = cbind(sla, wave_1998_wave_2482) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rST[[i]] <- modelST$predictions$pvals %>%
      filter(trait == "sla" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raS<- Reduce(rbind, rST)
  raS <- raS %>% left_join(S_bluesEF1[,c("taxa", "sla")])
  fwrite(raS, "ST_mwef.txt", sep = "\t", append = T, col.names = F)
  acST[j]<-cor(raS[,2], raS[,3], use = "complete.obs")
}
fwrite(as.matrix(acST), "acST_mwef.txt", sep = "\t", col.names = F)

#------ Multi-trait synthetic test trait with  plsr_narea for loc MW------

acpnT_mwef<-c()
fwrite(data.frame("taxa","GEBV"), "pnT_mwef.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rpnT_mwef<-list()
  for(i in 1:5){
    test<-pn_bluesMW1
    test[test$taxa %in% sort[[j]][[i]],"plsr_narea"] <- NA
    cat(i, j,"\n")
    modelpnT_mwef <- asreml(
      fixed = cbind(plsr_narea, wave_2201_wave_409) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rpnT_mwef[[i]] <- modelpnT_mwef$predictions$pvals %>%
      filter(trait == "plsr_narea" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  rapn<- Reduce(rbind, rpnT_mwef)
  rapn <- rapn %>% left_join(pn_bluesEF1[,c("taxa", "plsr_narea")])
  fwrite(rapn, "pnT_mwef.txt", sep = "\t", append = T, col.names = F)
  acpnT_mwef[j]<-cor(rapn[,2], rapn[,3], use = "complete.obs")
}
fwrite(as.matrix(acpnT_mwef), "acpnT_mwef.txt", sep = "\t", col.names = F)

#-----synthetic test trait with plsr_sla for loc MW------

acpsT_mwef<-c()
fwrite(data.frame("taxa","GEBV"), "psT_mwef.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rpsT_mwef<-list()
  for(i in 1:5){
    test<-ps_bluesMW1
    test[test$taxa %in% sort[[j]][[i]],"plsr_sla"] <- NA
    cat(i, j,"\n")
    modelpsT_mwef <- asreml(
      fixed = cbind(plsr_sla, wave_1434_wave_2258) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rpsT_mwef[[i]] <- modelpsT_mwef$predictions$pvals %>%
      filter(trait == "plsr_sla" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raps<- Reduce(rbind, rpsT_mwef)
  raps <- raps %>% left_join(ps_bluesEF1[,c("taxa", "plsr_sla")])
  fwrite(raps, "psT_mwef.txt", sep = "\t", append = T, col.names = F)
  acpsT_mwef[j]<-cor(raps[,2], raps[,3], use = "complete.obs")
}
fwrite(as.matrix(acpsT_mwef), "acpsT_mwef.txt", sep = "\t", col.names = F)
#CV1 model with synthetic test trait and four biological traits

#------  CV1 Multi-trait Narea and synthetic test trait------
acNT<-c()
fwrite(data.frame("taxa","GEBV"), "NT_CV1_mwef.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rNT<-list()
  for(i in 1:5){
    test<-N_bluesMW1
    test[test$taxa %in% sort[[j]][[i]],"narea"] <- NA
    test[test$taxa %in% sort[[j]][[i]],"wave_524_wave_681"] <- NA
    cat(i, j,"\n")
    modelNT <- asreml(
      fixed = cbind(narea, wave_524_wave_681) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rNT[[i]] <- modelNT$predictions$pvals %>%
      filter(trait == "narea" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raN<- Reduce(rbind, rNT)
  fwrite(raN, "NT_CV1_mwef.txt", sep = "\t", append = T, col.names = F)
  raN <- raN %>% left_join(N_bluesEF1[,c("taxa", "narea")])
  acNT[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite(as.matrix(acNT), "acNT_CV1_mwef.txt", sep = "\t", col.names = F)


#------ CV1 Multi-trait test trait and sla------
acST<-c()
fwrite(data.frame("taxa","GEBV"), "ST_CV1_mwef.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rST<-list()
  for(i in 1:5){
    test<-S_bluesMW1
    test[test$taxa %in% sort[[j]][[i]],"sla"] <- NA
    test[test$taxa %in% sort[[j]][[i]],"wave_1998_wave_2482"] <- NA
    modelST <- asreml(
      fixed = cbind(sla, wave_1998_wave_2482) ~ trait,
      random =  ~  corgh(trait):vm(taxa, source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rST[[i]] <- modelST$predictions$pvals %>%
      filter(trait == "sla" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raS<- Reduce(rbind, rST)
  fwrite(raS, "ST_CV1_mwef.txt", sep = "\t", append = T, col.names = F)
  raS <- raS %>% left_join(S_bluesEF1[,c("taxa", "sla")])
  acST[j]<-cor(raS[,2], raS[,3], use = "complete.obs")
}
fwrite(as.matrix(acST), "acST_CV1_mwef.txt", sep = "\t", col.names = F)

#------ CV1 Multi-trait test trait and plsr narea------
acpnT<-c()
fwrite(data.frame("taxa","GEBV"), "pnT_CV1_mwef.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rpnT<-list()
  for(i in 1:5){
    test<-pn_bluesMW1
    test[test$taxa %in% sort[[j]][[i]],"plsr_narea"] <- NA
    test[test$taxa %in% sort[[j]][[i]],"wave_2201_wave_409"] <- NA
    modelpnT <- asreml(
      fixed = cbind(plsr_narea, wave_2201_wave_409) ~ trait,
      random =  ~  corgh(trait):vm(taxa,source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rpnT[[i]] <- modelpnT$predictions$pvals %>%
      filter(trait == "plsr_narea" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raN<- Reduce(rbind, rpnT)
  fwrite(raN, "pnT_CV1_mwef.txt", sep = "\t", append = T, col.names = F)
  raN <- raN %>% left_join(pn_bluesEF1[,c("taxa", "plsr_narea")])
  acpnT[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite(as.matrix(acpnT), "acpnT_CV1_mwef.txt", sep = "\t", col.names = F)

#------CV1  Multi-trait test trait & plsr sla------
acpsT<-c()
fwrite(data.frame("taxa","GEBV"), "psT_CV1_mwef.txt", sep = "\t", col.names = F)
for(j in 1:length(sort)){
  rpsT<-list()
  for(i in 1:5){
    test<-ps_bluesMW1
    test[test$taxa %in% sort[[j]][[i]],"plsr_sla"] <- NA
    test[test$taxa %in% sort[[j]][[i]],"wave_1434_wave_2258"] <- NA
    modelpsT <- asreml(
      fixed = cbind(plsr_sla, wave_1434_wave_2258) ~ trait,
      random =  ~  corgh(trait):vm(taxa,source= kin, singG= "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa"))

    rpsT[[i]] <- modelpsT$predictions$pvals %>%
      filter(trait == "plsr_sla" & taxa %in% sort[[j]][[i]]) %>%
      select(taxa, predicted.value)
  }

  raS<- Reduce(rbind, rpsT)
  fwrite(raS, "psT_CV1_mwef.txt", sep = "\t", append = T, col.names = F)
  raS <- raS %>% left_join(ps_bluesEF1[,c("taxa", "plsr_sla")])
  acpsT[j]<-cor(raS[,2], raS[,3], use = "complete.obs")
}
fwrite(as.matrix(acpsT), "acpsT_CV1_mwef.txt", sep = "\t", col.names = F)
