#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript st_narea_envcv_one_rep.R <rep_id (1-5)>", call. = FALSE)
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 5) stop("rep_id must be an integer in 1..5", call. = FALSE)

# ---- user functions ----
# Make sure this path matches where create_folds() and crossv() live
source("function/aux_function.R")

# ---- inputs ----
infile <- sprintf("./output/N_blues_rep%d.csv", rep_id)

# Load kinship ONCE per run
kin <- fread("./data/kin_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

# ---- load data ----
N_blues <- read.csv(infile)

N_bluesEF <- N_blues %>%
  filter(env == "EF") %>%
  mutate(taxa = factor(taxa))

N_bluesMW <- N_blues %>%
  filter(env == "MW") %>%
  mutate(taxa = factor(taxa))

# folds based on EF taxa (same as your code)
sort <- create_folds(
  individuals = N_bluesEF$taxa,
  nfolds = 5,
  reps = 20,
  seed = 123
)

# ---- run EF -> MW ----
res_EFMW <- crossv(
  sort = sort,
  train = N_bluesEF,
  validation = N_bluesMW,
  mytrait = "narea",
  kin = kin,
  scheme = sprintf("EFMW_rep%d", rep_id)
)

corr_N_EFMW <- data.frame(res_EFMW$ac)
fwrite(corr_N_EFMW, sprintf("./output/corr_N_EFMWrep%d.csv", rep_id))

# ---- run MW -> EF ----
res_MWEF <- crossv(
  sort = sort,
  train = N_bluesMW,
  validation = N_bluesEF,
  mytrait = "narea",
  kin = kin,
  scheme = sprintf("MWEF_rep%d", rep_id)
)

corr_N_MWEF <- data.frame(res_MWEF$ac)

# Keep your current behavior, but make sure directory exists if you want /output/narea/
# If you truly want the mixed paths exactly as-is, remove the dir creation and change path.
out_mwef_dir <- "./output"
if (!dir.exists(out_mwef_dir)) dir.create(out_mwef_dir, recursive = TRUE)

fwrite(corr_N_MWEF, sprintf("%s/corr_N_MWEFrep%d.csv", out_mwef_dir, rep_id))

message(sprintf("DONE rep %d: wrote corr_N_EFMWrep%d.csv and corr_N_MWEFrep%d.csv", rep_id, rep_id, rep_id))
