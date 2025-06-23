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
parser$add_argument('--trait1', type = 'character', required = TRUE, help = 'First synthetic trait')
parser$add_argument('--trait2', type = 'character', required = TRUE, help = 'Second synthetic trait')
parser$add_argument('--trait3', type = 'character', required = TRUE, help = 'Third synthetic trait')
parser$add_argument('--repetition', type = 'integer', required = TRUE, help = 'Repetition number')
args <- parser$parse_args()

# Capture arguments
trait1 <- args$trait1
trait2 <- args$trait2
trait3 <- args$trait3
repetition <- args$repetition

# Define traits as a list
traits <- list(trait1, trait2, trait3)

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

# Load and prepare S_blues data
S_blues <- read.csv(paste0("S_blues_rep", repetition, ".csv"))
S_bluesEF <- S_blues %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
S_bluesMW <- S_blues %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))


# Function to run model for a given trait and CV scheme
run_trait_model <- function(trait, trait_index, cv_scheme) {
  ac_rep <- c()
  fname <- paste0("SW", trait_index, "_", cv_scheme, "_rep", repetition, "_mwef.txt")
  ac_fname <- paste0("acSW", trait_index, "_", cv_scheme, "_rep", repetition, "_mwef.txt")
  fwrite(data.frame("taxa", "GEBV"), fname, sep = "\t", col.names = FALSE)
  for (j in 1:length(sort)) {
    r_rep <- list()
    for (i in 1:5) {
      test <- S_bluesMW
      test[test$taxa %in% sort[[j]][[i]], "sla"] <- NA
      if (cv_scheme == "CV1") {
        test[test$taxa %in% sort[[j]][[i]], paste0("wave_", trait)] <- NA
      }

      model_rep <- asreml(
        fixed = as.formula(paste0("cbind(sla, wave_", trait, ") ~ trait")),
        random = ~ corgh(trait):vm(taxa, source = kin, singG = "NSD"),
        residual = ~ units:corgh(trait),
        data = test, na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "trait:taxa")
      )

      r_rep[[i]] <- model_rep$predictions$pvals %>%
        filter(trait == "sla" & taxa %in% sort[[j]][[i]]) %>%
        select(taxa, predicted.value)
    }

    raS <- Reduce(rbind, r_rep)
    raS <- raS %>% left_join(S_bluesEF[, c("taxa", "sla")])

    fwrite(raS, fname, sep = "\t", append = TRUE, col.names = FALSE)
    ac_rep[j] <- cor(raS[, 2], raS[, 3], use = "complete.obs")
  }

  fwrite(as.matrix(ac_rep), ac_fname, sep = "\t", col.names = FALSE)
}

# Define CV schemes
cv_schemes <- c("CV2", "CV1")

# Function to execute the model for a given trait and CV scheme
run_parallel_cv <- function(trait) {
  mclapply(cv_schemes, function(cv_scheme) {
    run_trait_model(trait, which(traits == trait), cv_scheme)
  }, mc.cores = 2)  # Use 2 cores for the CV schemes
}

# Run each trait in parallel across different cores
mclapply(traits, function(trait) {
  run_parallel_cv(trait)
}, mc.cores = 3)  # Use 3 cores for the traits
