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
N_bluesMW<- N_blues %>% filter(env== "MW")
N_bluesMW<- N_bluesMW%>% mutate(taxa= factor(taxa))

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

S_blues1<- read.csv("S_blues1.csv")
S_bluesEF1<- S_blues1 %>% filter(env== "EF")
S_bluesMW1<- S_blues1 %>% filter(env== "MW")
S_bluesEF1<- S_bluesEF1%>% mutate(taxa= factor(taxa))
S_bluesMW1<- S_bluesMW1%>% mutate(taxa= factor(taxa))

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
