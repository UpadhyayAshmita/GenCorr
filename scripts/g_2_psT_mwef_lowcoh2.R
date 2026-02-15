#loading package
library(tidyverse)
library(data.table)
library(asreml)
library(fs)
library(argparse)

# ---------------- args ----------------
parser <- ArgumentParser()
parser$add_argument("--wave_col", required = TRUE,
                    help = "Wave column name (e.g. wave_1998_wave_2482)")
args <- parser$parse_args()
wave_col <- args$wave_col

# loading kinship
kin <- fread('./data/kin_additive.txt', data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

#---------------- 1) use N_blues only to define MW taxa + folds ----------------
N_blues <- read.csv("./output/N_blues.csv")
N_bluesMW <- N_blues %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# ---------------------creating list with 5 fold and 20 reps----------------------
create_folds <- function(individuals, nfolds, reps, seed = 123){
  library(cvTools)
  library(dplyr)
  set.seed(seed)
  sort <- list()
  individuals <- as.factor(individuals)
  nl <- length(unique(individuals))
  for(a in 1:reps){
    folds <- cvFolds(nl, type ="random", K = nfolds)
    Sample <- cbind(folds$which, folds$subsets)
    cv <- split(levels(individuals)[Sample[,2]], f = Sample[,1])
    sort[[a]] <- cv
  }
  return(sort)
}

# CHANGED: folds now from MW
sort <- create_folds(individuals = N_bluesMW$taxa,
                     nfolds = 5, reps = 20, seed = 123)

#---------------- 2) load PLSR-SLA data (EF/MW) ----------------
psT_blues1 <- read.csv("./output/ps_blues_lowcoh2.csv")

psT_bluesEF1 <- psT_blues1 %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
psT_bluesMW1 <- psT_blues1 %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# -------------------- CV1 --------------------
# mask plsr_sla + wave for test fold taxa
acpsT_CV1 <- c()

fwrite(data.frame(taxa = character(), GEBV = numeric()),
       "psT_CV1_MWEF.csv", sep = ",", col.names = TRUE)

for(j in 1:length(sort)){
  rpsT <- list()
  for(i in 1:5){

    test <- psT_bluesMW1
    test[test$taxa %in% sort[[j]][[i]], "plsr_sla"] <- NA
    test[test$taxa %in% sort[[j]][[i]], wave_col] <- NA

    cat("psT_CV1", i, j, "\n")

    modelpsT <- asreml(
      fixed = as.formula(paste0("cbind(plsr_sla, ", wave_col, ") ~ trait")),
      random = ~ corgh(trait):vm(taxa, source = kin, singG = "NSD"),
      residual = ~ units:corgh(trait),
      data = test,
      na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa")
    )

    rpsT[[i]] <- modelpsT$predictions$pvals %>%
      filter(trait == "plsr_sla" & taxa %in% sort[[j]][[i]]) %>%
      transmute(taxa = as.character(taxa),
                GEBV = predicted.value)
  }

  ra <- Reduce(rbind, rpsT)
  fwrite(ra, "psT_CV1_MWEF.csv", sep = ",",
         append = TRUE, col.names = FALSE)

  ra_eval <- ra %>% left_join(psT_bluesEF1[, c("taxa", "plsr_sla")])
  acpsT_CV1[j] <- cor(ra_eval$GEBV, ra_eval$plsr_sla,
                      use = "complete.obs")
}

fwrite(data.frame(rep = 1:length(acpsT_CV1), acc = acpsT_CV1),
       "acpsT_CV1_MWEF.csv", sep = ",", col.names = TRUE)

# -------------------- CV2 --------------------
# mask plsr_sla only for test fold taxa
acpsT_CV2 <- c()

fwrite(data.frame(taxa = character(), GEBV = numeric()),
       "psT_CV2_MWEF.csv", sep = ",", col.names = TRUE)

for(j in 1:length(sort)){
  rpsT <- list()
  for(i in 1:5){

    test <- psT_bluesMW1
    test[test$taxa %in% sort[[j]][[i]], "plsr_sla"] <- NA
    # wave NOT masked in CV2

    cat("psT_CV2", i, j, "\n")

    modelpsT <- asreml(
      fixed = as.formula(paste0("cbind(plsr_sla, ", wave_col, ") ~ trait")),
      random = ~ corgh(trait):vm(taxa, source = kin, singG = "NSD"),
      residual = ~ units:corgh(trait),
      data = test,
      na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa")
    )

    rpsT[[i]] <- modelpsT$predictions$pvals %>%
      filter(trait == "plsr_sla" & taxa %in% sort[[j]][[i]]) %>%
      transmute(taxa = as.character(taxa),
                GEBV = predicted.value)
  }

  ra <- Reduce(rbind, rpsT)
  fwrite(ra, "psT_CV2_MWEF.csv", sep = ",",
         append = TRUE, col.names = FALSE)

  ra_eval <- ra %>% left_join(psT_bluesEF1[, c("taxa", "plsr_sla")])
  acpsT_CV2[j] <- cor(ra_eval$GEBV, ra_eval$plsr_sla,
                      use = "complete.obs")
}

fwrite(data.frame(rep = 1:length(acpsT_CV2), acc = acpsT_CV2),
       "acpsT_CV2_MWEF.csv", sep = ",", col.names = TRUE)
