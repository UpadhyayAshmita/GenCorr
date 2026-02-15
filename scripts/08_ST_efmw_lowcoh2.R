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

#loading kinship
kin <- fread('./data/kin_additive.txt', data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

#---------------- 1) use N_blues only to define EF taxa + folds ----------------
N_blues <- read.csv("./output/N_blues.csv")
N_bluesEF <- N_blues %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))

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

sort <- create_folds(individuals = N_bluesEF$taxa, nfolds = 5,
                     reps = 20, seed = 123)

#---------------- 2) load SLA data (EF/MW) ----------------
# change this path to your lowcoh2 SLA file if needed
S_blues1 <- read.csv("./output/S_blues_lowcoh2.csv")
S_bluesEF1 <- S_blues1 %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
S_bluesMW1 <- S_blues1 %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# -------------------- CV1 --------------------
# mask sla + wave for test fold taxa
acST_CV1 <- c()

# init CSV with header
fwrite(data.frame(taxa = character(), GEBV = numeric()),
       "ST_CV1.csv", sep = ",", col.names = TRUE)

for(j in 1:length(sort)){
  rST <- list()
  for(i in 1:5){
    test <- S_bluesEF1

    test[test$taxa %in% sort[[j]][[i]], "sla"] <- NA
    test[test$taxa %in% sort[[j]][[i]], wave_col] <- NA   # CV1 masks wave too

    cat("CV1", i, j, "\n")

    modelST <- asreml(
      fixed = as.formula(paste0("cbind(sla, ", wave_col, ") ~ trait")),
      random =  ~  corgh(trait):vm(taxa, source = kin, singG = "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa")
    )

    rST[[i]] <- modelST$predictions$pvals %>%
      filter(trait == "sla" & taxa %in% sort[[j]][[i]]) %>%
      transmute(taxa = as.character(taxa), GEBV = predicted.value)
  }

  raS <- Reduce(rbind, rST)

  # append to CSV (no header)
  fwrite(raS, "ST_CV1.csv", sep = ",", append = TRUE, col.names = FALSE)

  raS_eval <- raS %>% left_join(S_bluesMW1[, c("taxa", "sla")])
  acST_CV1[j] <- cor(raS_eval$GEBV, raS_eval$sla, use = "complete.obs")
}

# write accuracy as CSV with header
fwrite(data.frame(rep = 1:length(acST_CV1), acc = acST_CV1),
       "acST_CV1.csv", sep = ",", col.names = TRUE)

# -------------------- CV2 --------------------
# mask sla only for test fold taxa (wave stays observed)
acST_CV2 <- c()

# init CSV with header
fwrite(data.frame(taxa = character(), GEBV = numeric()),
       "ST_CV2.csv", sep = ",", col.names = TRUE)

for(j in 1:length(sort)){
  rST <- list()
  for(i in 1:5){
    test <- S_bluesEF1

    test[test$taxa %in% sort[[j]][[i]], "sla"] <- NA
    # wave NOT masked in CV2

    cat("CV2", i, j, "\n")

    modelST <- asreml(
      fixed = as.formula(paste0("cbind(sla, ", wave_col, ") ~ trait")),
      random =  ~  corgh(trait):vm(taxa, source = kin, singG = "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa")
    )

    rST[[i]] <- modelST$predictions$pvals %>%
      filter(trait == "sla" & taxa %in% sort[[j]][[i]]) %>%
      transmute(taxa = as.character(taxa), GEBV = predicted.value)
  }

  raS <- Reduce(rbind, rST)

  # append to CSV (no header)
  fwrite(raS, "ST_CV2.csv", sep = ",", append = TRUE, col.names = FALSE)

  raS_eval <- raS %>% left_join(S_bluesMW1[, c("taxa", "sla")])
  acST_CV2[j] <- cor(raS_eval$GEBV, raS_eval$sla, use = "complete.obs")
}

# write accuracy as CSV with header
fwrite(data.frame(rep = 1:length(acST_CV2), acc = acST_CV2),
       "acST_CV2.csv", sep = ",", col.names = TRUE)
