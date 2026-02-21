#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(asreml)
library(argparse)
library(cvTools)

# --- 1. Handle Arguments ---
parser <- ArgumentParser()
parser$add_argument('--trait', type = 'character', help='narea, ps, pn, sla, plsr_sla, or plsr_narea')
parser$add_argument('--cv_scheme', type = 'character', help='CV1 or CV2')
parser$add_argument('--w1', type = 'character', help='Column name for wave 1')
parser$add_argument('--w2', type = 'character', help='Column name for wave 2')
parser$add_argument('--rep', type = 'integer', help='Specific repetition to run (1-20)')
args <- parser$parse_args()

trait_name <- args$trait
cv_scheme  <- args$cv_scheme
w1_col     <- args$w1
w2_col     <- args$w2
curr_rep   <- args$rep

# --- 2. Generalized Directory Setup ---
asreml.options(workspace = "4gb", pworkspace = "2gb", maxit = 50)

# EFMW Scenario: Removing MWEF tag from directories
tmp_dir <- file.path(".", paste0("TMP_", trait_name, "_EFMW_", cv_scheme))
out_dir <- file.path("./output", trait_name)

if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 3. Load GINV (Using kinship_additive logic) ---
ginv_obj <- readRDS("../synthetic_traits/data/GINV.rds")
sparse_mat <- ginv_obj$Ginv.sparse
taxa_names <- attr(sparse_mat, "rowNames")
if(is.null(taxa_names)) taxa_names <- rownames(sparse_mat)

attr(sparse_mat, "rowNames") <- as.character(taxa_names)
attr(sparse_mat, "INVERSE") <- TRUE

# --- 4. Generalized Data Loading Logic ---
if (trait_name == "plsr_sla") {
  file_name <- "ps_blues.csv"
} else if (trait_name == "plsr_narea") {
  file_name <- "pn_blues.csv"
} else if (trait_name == "sla") {
  file_name <- "S_blues.csv"
} else {
  file_name <- "N_blues.csv"
}

pheno_path <- file.path("./output", file_name)
pheno_data <- read.csv(pheno_path)

# The Handshake
pheno_data <- pheno_data %>%
  filter(taxa %in% taxa_names) %>%
  mutate(taxa = factor(as.character(taxa), levels = taxa_names))

# --- EFMW Scenario Logic ---
# Train/Folds on EF | Observe/Test on MW
data_EF <- pheno_data %>% filter(env == "EF") %>% droplevels()
data_MW <- pheno_data %>% filter(env == "MW") %>% droplevels()

# --- 5. Starting Values (from EF) ---
fml_str <- paste0("cbind(", trait_name, ", ", w1_col, ", ", w2_col, ") ~ trait")
fixed_fml <- as.formula(fml_str)

# Pre-fitting using the training environment (EF)
pre_fit <- asreml(fixed = fixed_fml,
                  random = ~ us(trait):vm(taxa, source = sparse_mat),
                  residual = ~ units:diag(trait),
                  data = data_EF, maxit = 30,
                  na.action = na.method(x = "include", y = "include"))

start_values <- pre_fit$vparameters[grep("trait!us", names(pre_fit$vparameters))]

# --- 6. CV Setup (Based on EF Taxa) ---
set.seed(123)
u_indiv <- levels(data_EF$taxa)
sort_folds <- list()
for (a in 1:20) {
  folds <- cvFolds(length(u_indiv), type = "random", K = 5)
  Sample <- cbind(folds$which, folds$subsets)
  sort_folds[[a]] <- split(u_indiv[Sample[, 2]], f = Sample[, 1])
}

# --- 7. Run Specific Repetition ---
temp_fname  <- file.path(tmp_dir, paste0("rep", curr_rep, "_preds.csv"))
temp_acname <- file.path(tmp_dir, paste0("rep", curr_rep, "_accuracy.csv"))

r_rep <- vector("list", 5)
for (i in 1:5) {
  test_ids <- sort_folds[[curr_rep]][[i]]
  test_data <- data_EF # Masking within EF

  if (cv_scheme == "CV1") {
    test_data[test_data$taxa %in% test_ids, c(trait_name, w1_col, w2_col)] <- NA
  } else {
    test_data[test_data$taxa %in% test_ids, trait_name] <- NA
  }

  fit <- asreml(fixed = fixed_fml,
                random = ~ us(trait, init = start_values):vm(taxa, source = sparse_mat),
                residual = ~ units:diag(trait),
                data = test_data, maxit = 50,
                na.action = na.method(x = "include", y = "include"))

  pr <- predict(fit, classify = "trait:taxa", levels = list("taxa" = test_ids))$pvals
  r_rep[[i]] <- pr %>% filter(trait == trait_name) %>%
    mutate(rep = curr_rep, fold = i) %>%
    select(rep, fold, taxa, predicted.value) %>% rename(pred = predicted.value)
  rm(fit, pr); gc()
}

# Validation: Join with MW observations
raS <- bind_rows(r_rep) %>%
  left_join(data_MW[, c("taxa", trait_name)], by = "taxa") %>%
  rename(obs = !!sym(trait_name))

fwrite(raS, temp_fname)
accuracy_val <- cor(raS$pred, raS$obs, use = "complete.obs")
fwrite(data.frame(rep = curr_rep, accuracy = accuracy_val), temp_acname)

# --- 8. Generalized Merge Logic ---
all_acc_files <- list.files(tmp_dir, pattern = "_accuracy.csv", full.names = TRUE)

if (length(all_acc_files) == 20) {
  cat("\nMerging EFMW results for:", trait_name, "\n")

  # Output files do not include "MWEF" tag for this scenario
  final_acc <- lapply(all_acc_files, fread) %>% bind_rows() %>% arrange(rep)
  fwrite(final_acc, file.path(out_dir, paste0("acc_mt_", trait_name, "_", cv_scheme, ".csv")))

  all_pred_files <- list.files(tmp_dir, pattern = "_preds.csv", full.names = TRUE)
  final_preds <- lapply(all_pred_files, fread) %>% bind_rows() %>% arrange(rep, fold)
  fwrite(final_preds, file.path(out_dir, paste0("MT_", trait_name, "_", cv_scheme, "_FINAL_preds.csv")))

  unlink(tmp_dir, recursive = TRUE)
}
