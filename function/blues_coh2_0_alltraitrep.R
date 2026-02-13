library(data.table)
library(dplyr)
library(tidyr)
library(janitor)
library(asreml)

prep_wave_data <- function(df) {
  df |>
    clean_names() |>
    mutate(
      name2 = factor(name2),
      taxa  = factor(taxa),
      loc   = factor(loc),
      set   = factor(set),
      block = factor(block),
      range = factor(range),
      row   = factor(row)
    ) |>
    arrange(loc, range, row)
}

detect_ratio_col <- function(df) {
  ratio_cols <- grep("^wave_\\d+_wave_\\d+$", names(df), value = TRUE)
  if (length(ratio_cols) != 1) {
    stop("Expected exactly 1 ratio column like wave_####_wave_####; found: ",
         paste(ratio_cols, collapse = ", "))
  }
  ratio_cols[[1]]
}

fit_ratio_blues <- function(df, ratio_col, loc_value) {
  fml <- stats::as.formula(paste(ratio_col, "~ set + name2"))

  m <- asreml(
    fixed = fml,
    random = ~ block,
    data = df,
    subset = loc == loc_value,
    na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
  m <- update.asreml(m)

  temp <- m$predictions$pvals[, 1:2]
  temp$wave <- ratio_col

  temp |>
    as.data.frame() |>
    pivot_wider(names_from = wave, values_from = predicted.value)
}

fit_trait_blues <- function(df, trait_col, loc_value) {
  if (!trait_col %in% names(df)) {
    stop("Trait column '", trait_col, "' not found in data.")
  }

  fml <- stats::as.formula(paste(trait_col, "~ name2 + set"))

  m <- asreml(
    fixed = fml,
    random = ~ block,
    data = df,
    subset = loc == loc_value,
    na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2")
  )
  m <- update.asreml(m)

  out <- data.frame(
    name2 = m$predictions$pvals$name2,
    value = round(m$predictions$pvals$predicted.value, 3)
  )
  names(out)[names(out) == "value"] <- trait_col
  out
}

run_blues_all_reps_oldnames <- function(
    trait_col,          # e.g. "narea", "sla", "fs_plsr_narea", "plsr_sla_sorghum"
    in_prefix,          # e.g. "Nwave", "Swave", "pnwave", "pswave"
    out_prefix,         # e.g. "N", "S", "pn", "ps"
    reps = 1:5,
    in_dir = "./output",
    out_dir = "./output"
) {
  meta <- vector("list", length(reps))

  for (i in seq_along(reps)) {
    rep <- reps[i]

    in_path <- file.path(in_dir, paste0(in_prefix, "_rep", rep, "_lowcoh2.csv"))
    df <- read.csv(in_path) |> prep_wave_data()

    ratio_col <- detect_ratio_col(df)
    meta[[i]] <- data.frame(rep = rep, ratio_col = ratio_col)

    # --- ratio BLUES ---
    ratio_mw <- fit_ratio_blues(df, ratio_col, "MW")
    ratio_ef <- fit_ratio_blues(df, ratio_col, "EF")

    fwrite(ratio_mw, file.path(out_dir, paste0(out_prefix, "ratio_bluesMW1_rep", rep, ".csv")))
    fwrite(ratio_ef, file.path(out_dir, paste0(out_prefix, "ratio_bluesEF1_rep", rep, ".csv")))

    # --- trait BLUES ---
    trait_mw <- fit_trait_blues(df, trait_col, "MW")
    trait_ef <- fit_trait_blues(df, trait_col, "EF")

    # --- join ratio + trait (same as your pipeline) ---
    blues_mw <- left_join(ratio_mw, trait_mw, by = "name2")
    blues_ef <- left_join(ratio_ef, trait_ef, by = "name2")

    fwrite(blues_mw, file.path(out_dir, paste0(out_prefix, "bluesMW1_rep", rep, ".csv")))
    fwrite(blues_ef, file.path(out_dir, paste0(out_prefix, "bluesEF1_rep", rep, ".csv")))
  }

  meta <- bind_rows(meta)
  cat("\n===== Ratio column used per rep for", trait_col, "=====\n")
  print(meta)

  invisible(meta)
}
