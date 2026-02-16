#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ---------------- args ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript i_run_ps_singletraits_cv_re1-5_trycatch.R <rep_id (1-5)>", call. = FALSE)
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1 || rep_id > 5) stop("rep_id must be an integer in 1..5", call. = FALSE)

# ---------------- source only for create_folds (leave your old crossv untouched) ----------------
source("function/aux_function.R")  # provides create_folds()

# ---------------- robust crossv just for this script ----------------
crossv_try <- function(sort,
                       train,
                       validation,
                       kin,
                       mytrait,
                       scheme) {
  suppressPackageStartupMessages({
    library(asreml)
    library(dplyr)
    library(data.table)
  })

  ac <- c()
  gebv <- list()

  kin_taxa <- rownames(kin)

  # align taxa levels to kinship (critical for vm())
  train$taxa <- factor(train$taxa, levels = kin_taxa)
  validation$taxa <- factor(validation$taxa, levels = kin_taxa)

  if (any(is.na(train$taxa))) {
    bad <- unique(as.character(train$taxa[is.na(train$taxa)]))
    stop("Taxa in train not found in kin rownames: ", paste(bad, collapse = ", "))
  }
  if (any(is.na(validation$taxa))) {
    bad <- unique(as.character(validation$taxa[is.na(validation$taxa)]))
    stop("Taxa in validation not found in kin rownames: ", paste(bad, collapse = ", "))
  }

  for (j in seq_along(sort)) {
    r <- vector("list", length(sort[[j]]))

    for (i in seq_along(sort[[j]])) {
      test <- na.omit(train)
      test$taxa <- factor(test$taxa, levels = kin_taxa)

      # mask phenotype for this fold
      test[test$taxa %in% sort[[j]][[i]], mytrait] <- NA

      cat("scheme=", scheme, " rep=", j, " fold=", i,
          " n=", nrow(test),
          " n_obs=", sum(!is.na(test[[mytrait]])), "\n", sep = "")

      model <- tryCatch(
        asreml(
          fixed  = as.formula(paste0(mytrait, " ~ 1")),
          random = ~ vm(taxa, source = kin, singG = "NSD"),
          data   = test,
          na.action = na.method(x = "include"),
          predict = predict.asreml(classify = "taxa")
        ),
        error = function(e) e
      )

      if (inherits(model, "error")) {
        cat("ASReml FAILED (", scheme, ") rep=", j, " fold=", i,
            " : ", conditionMessage(model), "\n", sep = "")
        r[[i]] <- data.frame(taxa = sort[[j]][[i]], predicted.value = NA_real_)
        next
      }

      if (!isTRUE(model$converge)) model <- update.asreml(model)
      if (!isTRUE(model$converge)) model <- update.asreml(model)

      if (is.null(model$predictions$pvals)) {
        cat("ASReml returned no predictions (", scheme, ") rep=", j, " fold=", i, "\n", sep = "")
        r[[i]] <- data.frame(taxa = sort[[j]][[i]], predicted.value = NA_real_)
        next
      }

      p <- model$predictions$pvals
      out <- p[p$taxa %in% sort[[j]][[i]], 1:2]
      if (nrow(out) == 0) out <- data.frame(taxa = sort[[j]][[i]], predicted.value = NA_real_)

      r[[i]] <- out
    }

    gebv[[j]] <- Reduce(rbind, r) %>%
      left_join(validation[, c("taxa", mytrait)]) %>%
      mutate(rep = j)

    ac[j] <- suppressWarnings(cor(gebv[[j]][, 2], gebv[[j]][, 3], use = "complete.obs"))
  }

  gebv <- bind_rows(gebv)

  # keep same gebv naming/location as your original pipeline
  write.csv(gebv,
            paste0("./output/", "gebv_", mytrait, "_", scheme, ".csv"),
            row.names = FALSE)

  list(gebv = gebv, ac = ac)
}

# ---------------- inputs ----------------
psfile <- sprintf("./output/ps_blues_rep%d.csv", rep_id)
nfile  <- sprintf("./output/N_blues_rep%d.csv",  rep_id)

kin <- fread("./data/kin_additive.txt", data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

ps_blues <- read.csv(psfile)
N_blues  <- read.csv(nfile)

ps_bluesEF <- ps_blues %>%
  filter(env == "EF") %>%
  mutate(taxa = factor(taxa))

ps_bluesMW <- ps_blues %>%
  filter(env == "MW") %>%
  mutate(taxa = factor(taxa))

N_bluesEF <- N_blues %>%
  filter(env == "EF") %>%
  mutate(taxa = factor(taxa))

# folds based on N_bluesEF taxa (matches your old ps pipeline)
sort <- create_folds(
  individuals = N_bluesEF$taxa,
  nfolds = 5,
  reps = 20,
  seed = 123
)

# accuracies go to cv_singletrait (gebv stays in ./output via crossv_try write.csv)
acc_dir <- "./output/cv_singletrait"
if (!dir.exists(acc_dir)) dir.create(acc_dir, recursive = TRUE)

# ---------------- run EF -> MW ----------------
res_EFMW <- crossv_try(
  sort = sort,
  train = ps_bluesEF,
  validation = ps_bluesMW,
  mytrait = "ps",
  kin = kin,
  scheme = sprintf("EFMW_rep%d", rep_id)
)

corr_ps_EFMW <- data.frame(res_EFMW$ac)
fwrite(corr_ps_EFMW, file.path(acc_dir, sprintf("corr_ps_EFMWrep%d.csv", rep_id)))

# ---------------- run MW -> EF ----------------
res_MWEF <- crossv_try(
  sort = sort,
  train = ps_bluesMW,
  validation = ps_bluesEF,
  mytrait = "ps",
  kin = kin,
  scheme = sprintf("MWEF_rep%d", rep_id)
)

corr_ps_MWEF <- data.frame(res_MWEF$ac)
fwrite(corr_ps_MWEF, file.path(acc_dir, sprintf("corr_ps_MWEFrep%d.csv", rep_id)))

message(sprintf(
  "DONE ps rep %d: accuracies saved to %s; GEBV files written to ./output/",
  rep_id, acc_dir
))
