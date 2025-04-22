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

pn_blues1<- read.csv("pn_blues1.csv")
pn_bluesEF1<- pn_blues1 %>% filter(env== "EF")
pn_bluesMW1<- pn_blues1 %>% filter(env== "MW")
pn_bluesEF1<- pn_bluesEF1%>% mutate(taxa= factor(taxa))
pn_bluesMW1<- pn_bluesMW1%>% mutate(taxa= factor(taxa))

#----------------CV2 --------------------
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


