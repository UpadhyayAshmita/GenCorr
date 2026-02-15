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

#---------------- 2) load PLSR-NAREA data (EF/MW) ----------------
pnT_blues1 <- read.csv("./output/pn_blues_lowcoh2.csv")

pnT_bluesEF1 <- pnT_blues1 %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
pnT_bluesMW1 <- pnT_blues1 %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# response trait column name:
response_trait <- "plsr_narea"

# -------------------- CV1 --------------------
# mask plsr_narea + wave for test fold taxa
acpnT_CV1 <- c()

# init CSV with header
fwrite(data.frame(taxa = character(), GEBV = numeric()),
       "pnT_CV1.csv", sep = ",", col.names = TRUE)

for(j in 1:length(sort)){
  rpnT <- list()
  for(i in 1:5){
    test <- pnT_bluesEF1

    test[test$taxa %in% sort[[j]][[i]], response_trait] <- NA
    test[test$taxa %in% sort[[j]][[i]], wave_col] <- NA   # CV1 masks wave too

    cat("pnT_CV1", i, j, "\n")

    modelpnT <- asreml(
      fixed = as.formula(paste0("cbind(", response_trait, ", ", wave_col, ") ~ trait")),
      random =  ~  corgh(trait):vm(taxa, source = kin, singG = "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa")
    )

    rpnT[[i]] <- modelpnT$predictions$pvals %>%
      filter(trait == response_trait & taxa %in% sort[[j]][[i]]) %>%
      transmute(taxa = as.character(taxa), GEBV = predicted.value)
  }

  ra <- Reduce(rbind, rpnT)

  # append to CSV (no header)
  fwrite(ra, "pnT_CV1.csv", sep = ",", append = TRUE, col.names = FALSE)

  ra_eval <- ra %>% left_join(pnT_bluesMW1[, c("taxa", response_trait)])
  acpnT_CV1[j] <- cor(ra_eval$GEBV, ra_eval[[response_trait]], use = "complete.obs")
}

# write accuracy as CSV with header
fwrite(data.frame(rep = 1:length(acpnT_CV1), acc = acpnT_CV1),
       "acpnT_CV1.csv", sep = ",", col.names = TRUE)

# -------------------- CV2 --------------------
# mask plsr_narea only for test fold taxa (wave stays observed)
acpnT_CV2 <- c()

# init CSV with header
fwrite(data.frame(taxa = character(), GEBV = numeric()),
       "pnT_CV2.csv", sep = ",", col.names = TRUE)

for(j in 1:length(sort)){
  rpnT <- list()
  for(i in 1:5){
    test <- pnT_bluesEF1

    test[test$taxa %in% sort[[j]][[i]], response_trait] <- NA
    # wave NOT masked in CV2

    cat("pnT_CV2", i, j, "\n")

    modelpnT <- asreml(
      fixed = as.formula(paste0("cbind(", response_trait, ", ", wave_col, ") ~ trait")),
      random =  ~  corgh(trait):vm(taxa, source = kin, singG = "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa")
    )

    rpnT[[i]] <- modelpnT$predictions$pvals %>%
      filter(trait == response_trait & taxa %in% sort[[j]][[i]]) %>%
      transmute(taxa = as.character(taxa), GEBV = predicted.value)
  }

  ra <- Reduce(rbind, rpnT)

  # append to CSV (no header)
  fwrite(ra, "pnT_CV2.csv", sep = ",", append = TRUE, col.names = FALSE)

  ra_eval <- ra %>% left_join(pnT_bluesMW1[, c("taxa", response_trait)])
  acpnT_CV2[j] <- cor(ra_eval$GEBV, ra_eval[[response_trait]], use = "complete.obs")
}

# write accuracy as CSV with header
fwrite(data.frame(rep = 1:length(acpnT_CV2), acc = acpnT_CV2),
       "acpnT_CV2.csv", sep = ",", col.names = TRUE)
