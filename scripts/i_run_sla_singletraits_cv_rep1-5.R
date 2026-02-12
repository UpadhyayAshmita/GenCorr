#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript i_run_sla_singletraits_cv_re1-5.R <rep_id (1-5)>", call. = FALSE)
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 5) stop("rep_id must be an integer in 1..5", call. = FALSE)

# ---- user functions ----
source("function/aux_function.R")  # contains create_folds() + crossv()

# ---- inputs ----
sfile <- sprintf("./output/S_blues_rep%d.csv", rep_id)
nfile <- sprintf("./output/N_blues_rep%d.csv", rep_id)

# ---- kinship ----
kin <- fread("./data/kin_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

# ---- load data ----
S_blues <- read.csv(sfile)
N_blues <- read.csv(nfile)

S_bluesEF <- S_blues %>%
  filter(env == "EF") %>%
  mutate(taxa = factor(taxa))

S_bluesMW <- S_blues %>%
  filter(env == "MW") %>%
  mutate(taxa = factor(taxa))

N_bluesEF <- N_blues %>%
  filter(env == "EF") %>%
  mutate(taxa = factor(taxa))

# ---- folds based on N_bluesEF taxa (matches your old SLA pipeline) ----
sort <- create_folds(
  individuals = N_bluesEF$taxa,
  nfolds = 5,
  reps = 20,
  seed = 123
)

# ---- accuracy output dir only (GEBV stays as-is inside crossv(): ./output/gebv_*) ----
acc_dir <- "./output/cv_singletrait"
if (!dir.exists(acc_dir)) dir.create(acc_dir, recursive = TRUE)

# ---- run EF -> MW ----
res_EFMW <- crossv(
  sort = sort,
  train = S_bluesEF,
  validation = S_bluesMW,
  mytrait = "sla",
  kin = kin,
  scheme = sprintf("EFMW_rep%d", rep_id)
)

corr_S_EFMW <- data.frame(res_EFMW$ac)
fwrite(corr_S_EFMW, file.path(acc_dir, sprintf("corr_S_EFMWrep%d.csv", rep_id)))

# ---- run MW -> EF ----
res_MWEF <- crossv(
  sort = sort,
  train = S_bluesMW,
  validation = S_bluesEF,
  mytrait = "sla",
  kin = kin,
  scheme = sprintf("MWEF_rep%d", rep_id)
)

corr_S_MWEF <- data.frame(res_MWEF$ac)
fwrite(corr_S_MWEF, file.path(acc_dir, sprintf("corr_S_MWEFrep%d.csv", rep_id)))

message(sprintf(
  "DONE rep %d: saved accuracies to %s; GEBV files written by crossv() to ./output/",
  rep_id, acc_dir
))
