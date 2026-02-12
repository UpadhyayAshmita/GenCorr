library(tidyverse)
library(data.table)
library(asreml)
library(cvTools)

# Load kinship matrix
kin <- fread('./data/kin_additive.txt', data.table = FALSE)
rownames(kin) <- colnames(kin)
kin <- as.matrix(kin)

# Function to create folds
create_folds <- function(individuals, nfolds, reps, seed = 123) {
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

# Load  S_blues data
S_blues <- read.csv("./output/S_blues.csv")

taxa_levels <- rownames(kin)

S_bluesEF <- S_blues |>
  filter(env == "EF") |>
  mutate(taxa = factor(taxa, levels = taxa_levels))

S_bluesMW <- S_blues |>
  filter(env == "MW") |>
  mutate(taxa = factor(taxa, levels = taxa_levels))
sort <- create_folds(individuals = S_bluesMW$taxa, nfolds = 5, reps = 20, seed = 123)

w1_col <- "wave_1640_wave_1655"
w2_col <- "wave_738_wave_1111"
cv_scheme <- "CV1"

# use CV1 for now (always mask both ratios)
j <- 1
i <- 1
test_ids <- sort[[j]][[i]]

test <- S_bluesMW

# mask held-out taxa in MW
test[test$taxa %in% test_ids, "sla"] <- NA_real_
test[test$taxa %in% test_ids, w1_col] <- NA_real_
test[test$taxa %in% test_ids, w2_col] <- NA_real_

# fit 3-trait MT on MW
fixed_fml <- as.formula(paste0("cbind(sla, ", w1_col, ", ", w2_col, ") ~ trait"))

fit <- asreml(
  fixed     = fixed_fml,
  random    = ~ diag(trait):vm(taxa, source = kin, singG = "NSD"),
  residual  = ~ units:corgh(trait),
  data      = test,
  na.action = na.method(x = "include", y = "include"),
  workspace = 2e9   # 2 GB
)

# predict taxa effects
pr <- predict(fit, classify = "trait:taxa")$pvals

raS <- pr |>
  filter(trait == "sla", taxa %in% test_ids) |>
  transmute(taxa, pred = predicted.value) |>
  left_join(
    S_bluesEF |>
      select(taxa, sla) |>
      rename(obs = sla),
    by = "taxa"
  )

cor(raS$pred, raS$obs, use = "complete.obs")

