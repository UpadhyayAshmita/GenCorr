#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(asreml)
library(argparse)
library(cvTools)

# Argument parsing
parser <- ArgumentParser(description = 'Run synthetic trait models in parallel')
parser$add_argument('--trait', type = 'character', help = 'Synthetic trait')
parser$add_argument('--trait_index', type = 'integer', help = 'Trait index (1, 2, or 3)')
parser$add_argument('--cv_scheme', type = 'character', help = 'CV scheme (CV1 or CV2)')
args <- parser$parse_args()

trait <- args$trait
trait_index <- args$trait_index
cv_scheme <- args$cv_scheme

wave_col <- paste0("wave_", trait)

# Load kinship matrix
kin <- fread('./data/kin_additive.txt', data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

# Load and prepare N_blues data
N_blues <- read.csv("./output/N_blues.csv")

N_bluesEF <- N_blues %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
N_bluesMW <- N_blues %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# Function to create folds (same as yours, just no library() inside)
create_folds <- function(individuals, nfolds, reps, seed = 123) {
  set.seed(seed)
  sort <- list()
  individuals <- as.factor(individuals)
  nl <- length(unique(individuals))
  for (a in 1:reps) {
    folds <- cvFolds(nl, type = "random", K = nfolds)
    Sample <- cbind(folds$which, folds$subsets)
    cv <- split(levels(individuals)[Sample[, 2]], f = Sample[, 1])
    sort[[a]] <- cv
  }
  return(sort)
}

sort <- create_folds(individuals = N_bluesEF$taxa, nfolds = 5, reps = 20, seed = 123)

# Load and prepare S_blues data
pn_blues <- read.csv("./output/pn_blues.csv")
pn_bluesEF <- pn_blues %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
pn_bluesMW <- pn_blues %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# Function to run model for a given trait and CV scheme
run_trait_model <- function(trait, trait_index, cv_scheme) {

  ac_rep <- numeric(length(sort))

  fname    <- paste0("pnW", trait_index, "_", cv_scheme, ".csv")
  ac_fname <- paste0("acpnW", trait_index, "_", cv_scheme, ".csv")

  # header that matches what we write (taxa, pred, obs)
  fwrite(data.frame(taxa=character(), pred=numeric(), obs=numeric()),
         fname)

  for (j in 1:length(sort)) {

    r_rep <- vector("list", 5)

    for (i in 1:5) {

      test_ids <- sort[[j]][[i]]

      test <- pn_bluesEF
      test[test$taxa %in% test_ids, "plsr_narea"] <- NA

      if (cv_scheme == "CV1") {
        test[test$taxa %in% test_ids, wave_col] <- NA
      }

      # FIT
      fit <- asreml(
        fixed    = as.formula(paste0("cbind(plsr_narea, ", wave_col, ") ~ trait")),
        random   = ~ corgh(trait):vm(taxa, source = kin, singG = "NSD"),
        residual = ~ units:corgh(trait),
        data     = test,
        na.action = na.method(x = "include", y = "include")
      )

      # PREDICT (this is the key fix)
      pr <- predict(fit, classify = "trait:taxa")$pvals

      r_rep[[i]] <- pr %>%
        filter(trait == "plsr_narea" & taxa %in% test_ids) %>%
        select(taxa, predicted.value) %>%
        rename(pred = predicted.value)
    }

    raS <- Reduce(rbind, r_rep)
    raS <- raS %>% left_join(pn_bluesMW[, c("taxa", "plsr_narea")], by = "taxa") %>%
      rename(obs = plsr_narea)

    fwrite(raS, fname, sep = ",", append = TRUE, col.names = FALSE)

    ac_rep[j] <- cor(raS$pred, raS$obs, use = "complete.obs")
  }

  fwrite(as.matrix(ac_rep), ac_fname, sep = ",", col.names = FALSE)
}

# Run the model
run_trait_model(trait, trait_index, cv_scheme)
