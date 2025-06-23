#!/usr/bin/env Rscript

# Load necessary packages
library(tidyverse)
library(data.table)
library(asreml)
library(fs)
library(parallel)
library(argparse)

# Argument parsing
parser <- ArgumentParser(description = 'Run synthetic trait models in parallel')
parser$add_argument('--trait', type = 'character', help = 'Synthetic trait')
parser$add_argument('--trait_index', type = 'integer', help = 'Trait index (1, 2, or 3)')
parser$add_argument('--cv_scheme', type = 'character', help = 'CV scheme')
args <- parser$parse_args()

# Capture arguments
trait <- args$trait
trait_index <- args$trait_index
cv_scheme <- args$cv_scheme

# Load kinship matrix
kin <- fread('kin_additive.txt', data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

# Load and prepare N_blues data
N_blues <- read.csv("N_blues.csv")
N_bluesEF <- N_blues %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
N_bluesMW <- N_blues %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# Function to create folds
create_folds <- function(individuals, nfolds, reps, seed = 123) {
  library(cvTools)
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

# Function to run model for a given trait and CV scheme
run_trait_model <- function(trait, trait_index, cv_scheme) {
  ac_rep <- c()
  fname <- paste0("NW", trait_index, "_", cv_scheme, ".csv")
  ac_fname <- paste0("acNW", trait_index, "_", cv_scheme, ".csv")
  fwrite(data.frame("taxa", "GEBV"), fname, sep = ",", col.names = FALSE)

  for (j in 1:length(sort)) {
    r_rep <- list()
    for (i in 1:5) {
      test <- N_bluesEF
      test[test$taxa %in% sort[[j]][[i]], "narea"] <- NA
      if (cv_scheme == "CV1") {
        test[test$taxa %in% sort[[j]][[i]], paste0("wave_", trait)] <- NA
      }

      model_rep <- asreml(
        fixed = as.formula(paste0("cbind(narea, wave_", trait, ") ~ trait")),
        random = ~ corgh(trait):vm(taxa, source = kin, singG = "NSD"),
        residual = ~ units:corgh(trait),
        data = test, na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "trait:taxa")
      )

      r_rep[[i]] <- model_rep$predictions$pvals %>%
        filter(trait == "narea" & taxa %in% sort[[j]][[i]]) %>%
        select(taxa, predicted.value)
    }

    raS <- Reduce(rbind, r_rep)
    raS <- raS %>% left_join(N_bluesMW[, c("taxa", "narea")])

    fwrite(raS, fname, sep = ",", append = TRUE, col.names = FALSE)
    ac_rep[j] <- cor(raS[, 2], raS[, 3], use = "complete.obs")
  }

  fwrite(as.matrix(ac_rep), ac_fname, sep = ",", col.names = FALSE)
}

# Run the model
run_trait_model(trait, trait_index, cv_scheme)
