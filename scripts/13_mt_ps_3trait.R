#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(asreml)
library(argparse)
library(cvTools)

# --- 1. Handle Arguments ---
parser <- ArgumentParser()
parser$add_argument('--cv_scheme', type = 'character', help='CV1 or CV2')
parser$add_argument('--w1', type = 'character', help='Column name for wave 1')
parser$add_argument('--w2', type = 'character', help='Column name for wave 2')
args <- parser$parse_args()

cv_scheme <- args$cv_scheme
w1_col    <- args$w1
w2_col    <- args$w2

# --- 2. Settings ---
asreml.options(workspace = "4gb", pworkspace = "2gb", maxit = 50)

# --- 3. Load and Prepare GINV (Your Proven Stable Logic) ---
# Note: Using your established kinship_additive file path logic
ginv_obj <- readRDS("../synthetic_traits/data/GINV.rds")
sparse_mat <- ginv_obj$Ginv.sparse

taxa_names <- attr(sparse_mat, "rowNames")
if(is.null(taxa_names)) {
  taxa_names <- if(is.list(ginv_obj$rcn)) ginv_obj$rcn[[1]] else ginv_obj$rcn
}

attr(sparse_mat, "rowNames") <- as.character(taxa_names)
attr(sparse_mat, "INVERSE") <- TRUE

# --- 4. Load Phenotypes ---
# File: ps_blues.csv | Trait: plsr_sla
ps_data <- read.csv("./output/ps_blues.csv")

ps_EF <- ps_data %>%
  filter(env == "EF") %>%
  mutate(taxa = factor(taxa, levels = taxa_names))

ps_MW <- ps_data %>%
  filter(env == "MW") %>%
  mutate(taxa = factor(taxa, levels = taxa_names))

# --- 5. Get Starting Values ---
cat(sprintf("\nStep 1: Fitting full model for %s (plsr_sla + Waves)...\n", cv_scheme))
# Formula uses the actual column name: plsr_sla
fml_str <- paste0("cbind(plsr_sla, ", w1_col, ", ", w2_col, ") ~ trait")
fixed_fml <- as.formula(fml_str)

pre_fit <- asreml(fixed = fixed_fml,
                  random = ~ us(trait):vm(taxa, source = sparse_mat),
                  residual = ~ units:diag(trait),
                  data = ps_EF,
                  maxit = 30)

start_values <- pre_fit$vparameters[grep("trait!us", names(pre_fit$vparameters))]

# --- 6. CV Setup ---
set.seed(123)
u_indiv <- levels(ps_EF$taxa)
sort_folds <- list()
for (a in 1:20) {
  folds <- cvFolds(length(u_indiv), type = "random", K = 5)
  Sample <- cbind(folds$which, folds$subsets)
  sort_folds[[a]] <- split(u_indiv[Sample[, 2]], f = Sample[, 1])
}

# --- 7. Unified Loop ---
# Output files use 'ps' as requested
fname <- paste0("MT_ps_3Trait_", cv_scheme, "_preds.csv")
ac_fname <- paste0("MT_ps_3Trait_", cv_scheme, "_accuracy.csv")
fwrite(data.frame(rep=integer(), fold=integer(), taxa=character(), pred=numeric(), obs=numeric()), fname)

ac_rep <- numeric(20)
for (j in 1:20) {
  r_rep <- vector("list", 5)
  for (i in 1:5) {
    test_ids <- sort_folds[[j]][[i]]
    test_data <- ps_EF

    # MASKING LOGIC (Using plsr_sla from the file)
    if (cv_scheme == "CV1") {
      test_data[test_data$taxa %in% test_ids, c("plsr_sla", w1_col, w2_col)] <- NA
    } else {
      test_data[test_data$taxa %in% test_ids, "plsr_sla"] <- NA
    }

    fit <- asreml(fixed = fixed_fml,
                  random = ~ us(trait, init = start_values):vm(taxa, source = sparse_mat),
                  residual = ~ units:diag(trait),
                  data = test_data,
                  maxit = 50,
                  na.action = na.method(x = "include", y = "include"))

    # Predict trait:taxa and filter for plsr_sla
    pr <- predict(fit, classify = "trait:taxa", levels = list("taxa" = test_ids))$pvals
    r_rep[[i]] <- pr %>% filter(trait == "plsr_sla") %>%
      mutate(rep = j, fold = i) %>%
      select(rep, fold, taxa, predicted.value) %>% rename(pred = predicted.value)
    rm(fit, pr); gc()
  }

  # Join with MW observed values for plsr_sla
  raS <- bind_rows(r_rep) %>%
    left_join(ps_MW[, c("taxa", "plsr_sla")], by = "taxa") %>%
    rename(obs = plsr_sla)

  fwrite(raS, fname, append = TRUE, col.names = FALSE)
  ac_rep[j] <- cor(raS$pred, raS$obs, use = "complete.obs")
  cat(sprintf("Rep %d Complete: Acc = %.4f\n", j, ac_rep[j]))
}
fwrite(data.frame(rep=1:20, accuracy=ac_rep), ac_fname)
