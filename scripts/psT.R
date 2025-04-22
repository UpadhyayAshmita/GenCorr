#loading package
library(tidyverse)
library(data.table)
library(asreml)
library(fs)

#loading kinship
kin<- fread('kin_additive.txt', data.table= F)
rownames(kin) <- colnames(kin)
kin<- as.matrix(kin)
#----------------2nd step --------------------
N_blues<- read.csv("N_blues.csv")
N_bluesEF<- N_blues %>% filter(env== "EF")
N_bluesEF<- N_bluesEF%>% mutate(taxa= factor(taxa))

#creating sort
# ---------------------creating list with 5 fold and 20 reps----------------------
create_folds<- function(individuals, nfolds, reps, seed = 123){
  library(cvTools)
  library(dplyr)
  set.seed(seed)
  sort<-list()
  individuals <- as.factor(individuals)
  nl <- length(unique(individuals))
  for(a in 1:reps){
    folds <- cvFolds(nl,type ="random", K= nfolds)
    Sample<-cbind(folds$which,folds$subsets)
    cv<-split(levels(individuals)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of taxa
    sort[[a]]<-cv
  }
  return(sort)
}

sort<- create_folds( individuals= N_bluesEF$taxa, nfolds= 5,
                     reps = 20, seed = 123)

ps_blues1<- read.csv("ps_blues1.csv")
ps_bluesEF1<- ps_blues1 %>% filter(env== "EF")
ps_bluesMW1<- ps_blues1 %>% filter(env== "MW")
ps_bluesMW1<- ps_bluesMW1%>% mutate(taxa= factor(taxa))
ps_bluesEF1<- ps_bluesEF1%>% mutate(taxa= factor(taxa))

#----------------CV2 --------------------
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


