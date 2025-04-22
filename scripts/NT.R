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

N_blues1<- read.csv("N_blues1.csv")
N_bluesEF1<- N_blues1 %>% filter(env== "EF")
N_bluesMW1<- N_blues1 %>% filter(env== "MW")
N_bluesEF1<- N_bluesEF1%>% mutate(taxa= factor(taxa))
N_bluesMW1<- N_bluesMW1%>% mutate(taxa= factor(taxa))

#----------------CV2 --------------------
#------ Narea------
#------ Multi-trait N1 for loc EF------
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

