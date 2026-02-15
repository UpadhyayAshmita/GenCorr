#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(asreml)
library(parallel)
library(argparse)

# ---------------- args ----------------
parser <- ArgumentParser(description = "Run synthetic trait models in parallel (EF->MW and MW->EF)")
parser$add_argument("--trait1", type = "character", required = TRUE, help = "Trait 1 (WITHOUT leading 'wave_')")
parser$add_argument("--trait2", type = "character", required = TRUE, help = "Trait 2 (WITHOUT leading 'wave_')")
parser$add_argument("--trait3", type = "character", required = TRUE, help = "Trait 3 (WITHOUT leading 'wave_')")
parser$add_argument("--repetition", type = "integer", required = TRUE, help = "Repetition number")
args <- parser$parse_args()

traits <- list(args$trait1, args$trait2, args$trait3)
repetition <- args$repetition

# ---------------- paths ----------------
out_dir <- "output"
data_dir <- "data"

# ---------------- inputs ----------------
kin <- fread(file.path(data_dir, "kin_additive.txt"), data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

N_blues <- read.csv(file.path(out_dir, paste0("N_blues_rep", repetition, ".csv")))
N_bluesEF <- N_blues %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
N_bluesMW <- N_blues %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# ---------------- folds ----------------
create_folds <- function(individuals, nfolds, reps, seed = 123) {
  library(cvTools)
  set.seed(seed)
  individuals <- as.factor(individuals)
  nl <- length(unique(individuals))

  sort <- vector("list", reps)
  for (a in seq_len(reps)) {
    folds <- cvFolds(nl, type = "random", K = nfolds)
    Sample <- cbind(folds$which, folds$subsets)
    sort[[a]] <- split(levels(individuals)[Sample[, 2]], f = Sample[, 1])
  }
  sort
}

# based on EF taxa levels)
sort <- create_folds(individuals = N_bluesEF$taxa, nfolds = 5, reps = 20, seed = 123)

# ---------------- core model runner ----------------
run_trait_model <- function(trait, trait_index, cv_scheme, direction) {
  # direction:
  #   "efmw" => train on EF, evaluate vs MW narea
  #   "mwef" => train on MW, evaluate vs EF narea

  if (direction == "efmw") {
    train_df <- N_bluesEF
    eval_df  <- N_bluesMW
    suffix   <- ""
  } else if (direction == "mwef") {
    train_df <- N_bluesMW
    eval_df  <- N_bluesEF
    suffix   <- "_mwef"
  } else {
    stop("Unknown direction: ", direction)
  }

  fname    <- file.path(out_dir, paste0("NW", trait_index, "_", cv_scheme, "_rep", repetition, suffix, ".csv"))
  ac_fname <- file.path(out_dir, paste0("acNW", trait_index, "_", cv_scheme, "_rep", repetition, suffix, ".csv"))

  # Write header once (empty data frame)
  fwrite(data.frame(taxa = character(), GEBV = numeric()), fname)

  ac_rep <- numeric(length(sort))

  for (j in seq_along(sort)) {
    r_rep <- vector("list", 5)

    for (i in 1:5) {
      test <- train_df

      # Mask narea for test fold
      test[test$taxa %in% sort[[j]][[i]], "narea"] <- NA

      # CV1 also masks the wave trait for the test fold
      if (cv_scheme == "CV1") {
        test[test$taxa %in% sort[[j]][[i]], paste0("wave_", trait)] <- NA
      }

      model_rep <- asreml(
        fixed = as.formula(paste0("cbind(narea, wave_", trait, ") ~ trait")),
        random = ~ corgh(trait):vm(taxa, source = kin, singG = "NSD"),
        residual = ~ units:corgh(trait),
        data = test,
        na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "trait:taxa")
      )

      r_rep[[i]] <- model_rep$predictions$pvals %>%
        filter(trait == "narea" & taxa %in% sort[[j]][[i]]) %>%
        select(taxa, predicted.value)
    }

    raS <- Reduce(rbind, r_rep)

    # Join observed narea from the opposite environment (eval_df)
    raS <- raS %>% left_join(eval_df[, c("taxa", "narea")], by = "taxa")

    # Append GEBVs
    fwrite(
      raS %>% transmute(taxa = taxa, GEBV = predicted.value),
      fname,
      append = TRUE
    )

    # Accuracy = corr(GEBV, observed narea in eval env)
    ac_rep[j] <- cor(raS$predicted.value, raS$narea, use = "complete.obs")
  }

  fwrite(data.frame(accuracy = ac_rep), ac_fname)
}

cv_schemes <- c("CV2", "CV1")
directions <- c("efmw", "mwef")

# run both directions for each trait (still parallel)
run_parallel_for_trait <- function(trait) {
  trait_index <- which(unlist(traits) == trait)

  # 4 tasks per trait: 2 CV schemes x 2 directions
  tasks <- expand.grid(cv = cv_schemes, dir = directions, stringsAsFactors = FALSE)

  mclapply(
    seq_len(nrow(tasks)),
    function(k) {
      run_trait_model(
        trait = trait,
        trait_index = trait_index,
        cv_scheme = tasks$cv[k],
        direction = tasks$dir[k]
      )
    },
    mc.cores = 4
  )
}

mclapply(traits, run_parallel_for_trait, mc.cores = 3)
