#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript i_run_pn_singletraits_cv_re1-5.R <rep_id (1-5)>", call. = FALSE)
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 5) stop("rep_id must be an integer in 1..5", call. = FALSE)

# ---- user functions ----
source("function/aux_function.R")  # contains create_folds() + crossv()

# ---- inputs (all in ./output/) ----
pnfile <- sprintf("./output/pn_blues_rep%d.csv", rep_id)
nfile  <- sprintf("./output/N_blues_rep%d.csv",  rep_id)

# ---- kinship ----
kin <- fread("./data/kin_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

# ---- load data ----
pn_blues <- read.csv(pnfile)
N_blues  <- read.csv(nfile)

pn_bluesEF <- pn_blues %>%
  filter(env == "EF") %>%
  mutate(taxa = factor(taxa))

pn_bluesMW <- pn_blues %>%
  filter(env == "MW") %>%
  mutate(taxa = factor(taxa))

N_bluesEF <- N_blues %>%
  filter(env == "EF") %>%
  mutate(taxa = factor(taxa))

# ---- folds based on N_bluesEF taxa (matches your old pn pipeline) ----
sort <- create_folds(
  individuals = N_bluesEF$taxa,
  nfolds = 5,
  reps = 20,
  seed = 123
)

# accuracy output dir only
acc_dir <- "./output/cv_singletrait"
if (!dir.exists(acc_dir)) dir.create(acc_dir, recursive = TRUE)

# ---- run EF -> MW ----
res_EFMW <- crossv(
  sort = sort,
  train = pn_bluesEF,
  validation = pn_bluesMW,
  mytrait = "pn",
  kin = kin,
  scheme = sprintf("EFMW_rep%d", rep_id)
)

corr_pn_EFMW <- data.frame(res_EFMW$ac)
fwrite(corr_pn_EFMW, file.path(acc_dir, sprintf("corr_pn_EFMWrep%d.csv", rep_id)))

# ---- run MW -> EF ----
res_MWEF <- crossv(
  sort = sort,
  train = pn_bluesMW,
  validation = pn_bluesEF,
  mytrait = "pn",
  kin = kin,
  scheme = sprintf("MWEF_rep%d", rep_id)
)

corr_pn_MWEF <- data.frame(res_MWEF$ac)
fwrite(corr_pn_MWEF, file.path(acc_dir, sprintf("corr_pn_MWEFrep%d.csv", rep_id)))

message(sprintf(
  "DONE rep %d: saved accuracies to %s; GEBV files written by crossv() to ./output/",
  rep_id, acc_dir
))
