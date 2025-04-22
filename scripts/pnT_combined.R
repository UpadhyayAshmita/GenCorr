#!/usr/bin/env Rscript

# Load necessary packages
library(tidyverse)
library(data.table)
library(asreml)
library(fs)
library(argparse)

# Argument parsing
parser <- ArgumentParser(description = 'Run synthetic trait models')
parser$add_argument('--trait', type = 'character', help = 'Synthetic trait')
parser$add_argument('--repetition', type = 'integer', help = 'Repetition number')
parser$add_argument('--cv_scheme', type = 'character', help = 'Cross-validation scheme')
args <- parser$parse_args()

# Ensure arguments are correctly captured
trait <- args$trait
repetition <- args$repetition
cv_scheme <- args$cv_scheme

# Load kinship matrix
kin <- fread('kin_additive.txt', data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)
# Load and prepare N_blues data
N_blues <- read.csv(paste0("N_blues_rep", repetition, ".csv"))
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

# Load and prepare pn_blues data
pn_blues <- read.csv(paste0("pn_blues1_rep", repetition, ".csv"))
pn_bluesEF <- pn_blues %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
pn_bluesMW <- pn_blues %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# Function to run model for the given trait and CV scheme
run_trait_model <- function() {
  ac_rep <- c()
  fname <- paste0("pnNT", repetition, "_", cv_scheme, "_rep", repetition, ".txt")
  ac_fname <- paste0("acpnNT", repetition, "_", cv_scheme, "_rep", repetition, ".txt")

  fwrite(data.frame("taxa", "GEBV"), fname, sep = "\t", col.names = FALSE)

  for (j in 1:length(sort)) {
    r_rep <- list()
    for (i in 1:5) {
      test <- pn_bluesEF
      test[test$taxa %in% sort[[j]][[i]], "pn"] <- NA
      if (cv_scheme == "CV1") {
        test[test$taxa %in% sort[[j]][[i]], paste0("wave_", trait)] <- NA
      }

      model_rep <- asreml(
        fixed = as.formula(paste0("cbind(pn, wave_", trait, ") ~ trait")),
        random = ~ corgh(trait):vm(taxa, source = kin, singG = "NSD"),
        residual = ~ units:corgh(trait),
        data = test, na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "trait:taxa")
      )

      r_rep[[i]] <- model_rep$predictions$pvals %>%
        filter(trait == "pn" & taxa %in% sort[[j]][[i]]) %>%
        select(taxa, predicted.value)
    }

    raS <- Reduce(rbind, r_rep)
    raS <- raS %>% left_join(pn_bluesMW[, c("taxa", "pn")])

    fwrite(raS, fname, sep = "\t", append = TRUE, col.names = FALSE)
    ac_rep[j] <- cor(raS[, 2], raS[, 3], use = "complete.obs")
  }

  fwrite(as.matrix(ac_rep), ac_fname, sep = "\t", col.names = FALSE)
}

# Run the model for the specified trait and CV scheme
run_trait_model()
