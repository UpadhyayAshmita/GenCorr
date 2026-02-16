#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(asreml)
  library(parallel)
  library(argparse)
})

# ---------------- args ----------------
parser <- ArgumentParser(description = "Run PS synthetic trait models (EFMW + MWEF) in parallel")
parser$add_argument("--trait1", type="character", required=TRUE)
parser$add_argument("--trait2", type="character", required=TRUE)
parser$add_argument("--trait3", type="character", required=TRUE)
parser$add_argument("--repetition", type="integer", required=TRUE)
args <- parser$parse_args()

traits <- list(args$trait1, args$trait2, args$trait3)
repetition <- args$repetition

# ------------- paths (as requested) -------------
kin_path   <- file.path("./data", "kin_additive.txt")
blues_path <- file.path("./output", paste0("ps_blues_rep", repetition, ".csv"))

# ---------------- load kin ----------------
kin <- fread(kin_path, data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

# ---------------- load blues ----------------
ps_blues <- read.csv(blues_path)
ps_bluesEF <- ps_blues %>% filter(env == "EF") %>% mutate(taxa = factor(taxa))
ps_bluesMW <- ps_blues %>% filter(env == "MW") %>% mutate(taxa = factor(taxa))

# ---------------- folds----------------
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

sort <- create_folds(individuals = ps_bluesEF$taxa, nfolds = 5, reps = 20, seed = 123)

# -------- helper: robust wave column name --------
wave_colname <- function(trait_arg) {
  if (startsWith(trait_arg, "wave_")) trait_arg else paste0("wave_", trait_arg)
}

# ---------------- core runner ----------------
run_trait_model <- function(trait, trait_index, cv_scheme, direction_tag) {
  train_df <- if (direction_tag == "EFMW") ps_bluesEF else ps_bluesMW
  valid_df <- if (direction_tag == "EFMW") ps_bluesMW else ps_bluesEF

  wcol <- wave_colname(trait)

  suffix <- if (direction_tag == "MWEF") "_mwef" else ""
  fname    <- file.path("./output", paste0("psW", trait_index, "_", cv_scheme, "_rep", repetition, suffix, ".csv"))
  ac_fname <- file.path("./output", paste0("acpsW", trait_index, "_", cv_scheme, "_rep", repetition, suffix, ".csv"))

  fwrite(data.frame(taxa="taxa", GEBV="GEBV", observed="observed"),
         fname, sep = ",", col.names = FALSE)

  ac_rep <- c()

  for (j in seq_along(sort)) {
    r_rep <- vector("list", 5)

    for (i in 1:5) {
      test <- train_df
      test[test$taxa %in% sort[[j]][[i]], "ps"] <- NA

      if (cv_scheme == "CV1") {
        if (!wcol %in% names(test)) stop("Missing wave column: ", wcol)
        test[test$taxa %in% sort[[j]][[i]], wcol] <- NA
      }

      model_rep <- asreml(
        fixed = as.formula(paste0("cbind(ps, ", wcol, ") ~ trait")),
        random = ~ corgh(trait):vm(taxa, source = kin, singG = "NSD"),
        residual = ~ units:corgh(trait),
        data = test,
        na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "trait:taxa")
      )

      r_rep[[i]] <- model_rep$predictions$pvals %>%
        filter(trait == "ps" & taxa %in% sort[[j]][[i]]) %>%
        select(taxa, predicted.value)
    }

    raS <- Reduce(rbind, r_rep)
    raS <- raS %>% left_join(valid_df[, c("taxa", "ps")], by = "taxa")

    fwrite(raS, fname, sep = ",", append = TRUE, col.names = FALSE)
    ac_rep[j] <- cor(raS[, 2], raS[, 3], use = "complete.obs")
  }

  fwrite(data.frame(accuracy = ac_rep), ac_fname, sep = ",", col.names = TRUE)
}

cv_schemes <- c("CV2", "CV1")
directions <- c("EFMW", "MWEF")

run_parallel_for_trait <- function(trait) {
  trait_index <- which(unlist(traits) == trait)
  jobs <- expand.grid(cv_scheme = cv_schemes, direction = directions, stringsAsFactors = FALSE)

  mclapply(seq_len(nrow(jobs)), function(k) {
    run_trait_model(trait, trait_index, jobs$cv_scheme[k], jobs$direction[k])
  }, mc.cores = 2)
}

mclapply(unlist(traits), run_parallel_for_trait, mc.cores = 3)
