#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(cvTools)
  library(asreml)
})

source("./function/aux_function.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript scripts/run_st_cv_one.R <trait> <scheme>\n  trait in: narea,sla,plsr_narea,plsr_sla\n  scheme in: EFMW,MWEF")

trait  <- args[1]
scheme <- args[2]

# ---------- helpers ----------
read_trait_file <- function(trait){
  if (trait == "narea")      return("./output/N_blues.csv")
  if (trait == "sla")        return("./output/S_blues.csv")
  if (trait == "plsr_narea") return("./output/pn_blues.csv")
  if (trait == "plsr_sla")   return("./output/ps_blues.csv")
  stop("Unknown trait: ", trait)
}

out_file <- function(trait, scheme){
  dir.create("./output/cv_singletrait", showWarnings = FALSE, recursive = TRUE)
  file.path("./output/cv_singletrait", paste0("corr_", trait, "_", scheme, ".csv"))
}

# ---------- data ----------
df <- read.csv(read_trait_file(trait)) %>%
  mutate(
    env  = as.character(env),
    taxa = factor(taxa)
  )

dfEF <- df %>% filter(env == "EF")
dfMW <- df %>% filter(env == "MW")

# folds MUST be created on the TRAINING set taxa (to avoid mismatch)
# so we create them conditionally by scheme
make_folds <- function(train_df){
  create_folds(individuals = train_df$taxa, nfolds = 5, reps = 20, seed = 123)
}

kin <- fread("./data/kin_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

if (scheme == "EFMW") {
  train <- dfEF
  valid <- dfMW
  sort  <- make_folds(train)
} else if (scheme == "MWEF") {
  train <- dfMW
  valid <- dfEF
  sort  <- make_folds(train)
} else {
  stop("Unknown scheme: ", scheme)
}

# keep ONLY taxa present in BOTH train/valid and in kinship
common_ids <- Reduce(intersect, list(levels(train$taxa), levels(valid$taxa), rownames(kin)))
train <- train %>% filter(taxa %in% common_ids) %>% droplevels()
valid <- valid %>% filter(taxa %in% common_ids) %>% droplevels()
kin2  <- kin[common_ids, common_ids, drop = FALSE]

message("trait=", trait, " scheme=", scheme, " n_train=", nrow(train), " n_valid=", nrow(valid), " n_kin=", nrow(kin2))

# ---------- run CV ----------
res <- crossv(
  sort       = sort,
  train      = train,
  validation = valid,
  mytrait    = trait,
  kin        = kin2,
  scheme     = scheme
)

corr <- data.frame(res$ac)
fwrite(corr, out_file(trait, scheme))
message("Wrote: ", out_file(trait, scheme))
