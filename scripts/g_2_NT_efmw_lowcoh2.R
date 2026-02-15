#loading package
library(tidyverse)
library(data.table)
library(asreml)
library(fs)
library(argparse)

# ---------------- args ----------------
parser <- ArgumentParser()
parser$add_argument("--wave_col", required = TRUE,
                    help = "Wave column name (e.g. wave_524_wave_681)")
args <- parser$parse_args()

wave_col <- args$wave_col

#loading kinship
kin <- fread('./data/kin_additive.txt', data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

#----------------2nd step --------------------
N_blues <- read.csv("./output/N_blues.csv")
N_bluesEF <- N_blues %>% filter(env == "EF")
N_bluesEF <- N_bluesEF %>% mutate(taxa = factor(taxa))

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

N_blues1 <- read.csv("./output/N_blues_lowcoh2.csv")
N_bluesEF1 <- N_blues1 %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
N_bluesMW1 <- N_blues1 %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# -------------------- CV1 --------------------
# mask narea + wave for test fold taxa
acNT_CV1 <- c()

# init CSV with header
fwrite(data.frame(taxa = character(), GEBV = numeric()),
       "NT_CV1.csv", sep = ",", col.names = TRUE)

for(j in 1:length(sort)){
  rNT <- list()
  for(i in 1:5){
    test <- N_bluesEF1

    test[test$taxa %in% sort[[j]][[i]], "narea"] <- NA
    test[test$taxa %in% sort[[j]][[i]], wave_col] <- NA

    cat("CV1", i, j, "\n")

    modelNT <- asreml(
      fixed = as.formula(paste0("cbind(narea, ", wave_col, ") ~ trait")),
      random =  ~  corgh(trait):vm(taxa, source = kin, singG = "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa")
    )

    rNT[[i]] <- modelNT$predictions$pvals %>%
      filter(trait == "narea" & taxa %in% sort[[j]][[i]]) %>%
      transmute(taxa = as.character(taxa), GEBV = predicted.value)
  }

  raN <- Reduce(rbind, rNT)

  # append to CSV (no header)
  fwrite(raN, "NT_CV1.csv", sep = ",", append = TRUE, col.names = FALSE)

  raN_eval <- raN %>% left_join(N_bluesMW1[, c("taxa", "narea")])
  acNT_CV1[j] <- cor(raN_eval$GEBV, raN_eval$narea, use = "complete.obs")
}

# write accuracy as CSV with header
fwrite(data.frame(rep = 1:length(acNT_CV1), acc = acNT_CV1),
       "acNT_CV1.csv", sep = ",", col.names = TRUE)

# -------------------- CV2 --------------------
# mask narea only for test fold taxa (wave stays observed)
acNT_CV2 <- c()

# init CSV with header
fwrite(data.frame(taxa = character(), GEBV = numeric()),
       "NT_CV2.csv", sep = ",", col.names = TRUE)

for(j in 1:length(sort)){
  rNT <- list()
  for(i in 1:5){
    test <- N_bluesEF1

    test[test$taxa %in% sort[[j]][[i]], "narea"] <- NA
    # wave NOT masked in CV2

    cat("CV2", i, j, "\n")

    modelNT <- asreml(
      fixed = as.formula(paste0("cbind(narea, ", wave_col, ") ~ trait")),
      random =  ~  corgh(trait):vm(taxa, source = kin, singG = "NSD"),
      residual =  ~ units:corgh(trait),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "trait:taxa")
    )

    rNT[[i]] <- modelNT$predictions$pvals %>%
      filter(trait == "narea" & taxa %in% sort[[j]][[i]]) %>%
      transmute(taxa = as.character(taxa), GEBV = predicted.value)
  }

  raN <- Reduce(rbind, rNT)

  # append to CSV (no header)
  fwrite(raN, "NT_CV2.csv", sep = ",", append = TRUE, col.names = FALSE)

  raN_eval <- raN %>% left_join(N_bluesMW1[, c("taxa", "narea")])
  acNT_CV2[j] <- cor(raN_eval$GEBV, raN_eval$narea, use = "complete.obs")
}

# write accuracy as CSV with header
fwrite(data.frame(rep = 1:length(acNT_CV2), acc = acNT_CV2),
       "acNT_CV2.csv", sep = ",", col.names = TRUE)
