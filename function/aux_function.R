nareapreprocessing_allrep <- function(rep_id,
                                          breakdown_dir = "./output",
                                          phenotypes_path = "./data/phenotypes_whole.csv",
                                          sample_path = "./data/sample_per_rep.csv",
                                          out_dir = "./output",
                                          k_groups = 3) {  #k_group is for cutting the dendogram in three branch/groups

  dir.create(out_dir, showWarnings = FALSE)

  # Inputs
  breakdown_path <- file.path(breakdown_dir, paste0("narea_breakdown", rep_id, ".csv"))
  N <- data.table::fread(breakdown_path, data.table = FALSE)

  phenotypes_whole <- read.csv(phenotypes_path)
  indv_sample <- read.csv(sample_path)

  rep_col <- paste0("Rep_", rep_id)
  phenotypes <- phenotypes_whole %>%
    dplyr::filter(.data$Name2 %in% indv_sample[[rep_col]])

  wave <- phenotypes %>% janitor::clean_names()

  # QC mask (same)
  N <- N %>%
    dplyr::mutate(dplyr::across(
      coh2:vartrait,
      ~ replace(
        .,
        which(
          corgblup < -1 |
            corgblup > 1 |
            corg < -1 |
            corg > 1 |
            (abs(corg - corgblup) > 0.2)
        ),
        NA
      )
    ))

  ntraits <- ceiling(sum(!is.na(N$coh2)) * 0.01)

  sel <- N %>%
    dplyr::slice_max(coh2, n = ntraits) %>%
    dplyr::mutate(dplyr::across(where(is.double), \(x) round(x, 4)))

  # Make sure wave_1/wave_2 are clean strings like "wave_1093"
  sel$wave_1 <- gsub("\\s+", "_", as.character(sel$wave_1))
  sel$wave_2 <- gsub("\\s+", "_", as.character(sel$wave_2))

  # Create ratio columns in wave (same math)
  for (i in seq_len(nrow(sel))) {
    w1 <- sel$wave_1[i]
    w2 <- sel$wave_2[i]
    out <- paste0(w1, "_", w2)

    if (!w1 %in% names(wave) || !w2 %in% names(wave)) {
      stop(paste("Missing wave columns in `wave`:", w1, "or", w2, "at i =", i))
    }
    wave[[out]] <- wave[[w1]] / wave[[w2]]
  }

  sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)

  # Automatic range endpoints (replaces hard-coded per rep)
  first_ratio <- sel$rowname[1]
  last_ratio  <- sel$rowname[nrow(sel)]

  if (!first_ratio %in% names(wave)) stop(paste("first_ratio not found in wave:", first_ratio))
  if (!last_ratio  %in% names(wave)) stop(paste("last_ratio not found in wave:", last_ratio))

  Nratios <- wave %>%
    dplyr::select(dplyr::all_of(first_ratio):dplyr::all_of(last_ratio))

  Nratios_t <- t(Nratios)
  colnames(Nratios_t) <- wave$plot_id

  dist <- factoextra::get_dist(Nratios_t, method = "pearson")
  hcN <- hclust(dist, method = "average")
  sub_grp <- cutree(hcN, k = k_groups)

  Nratios_out <- Nratios_t %>%
    as.data.frame() %>%
    tibble::rownames_to_column("rowname") %>%
    dplyr::left_join(sel %>% dplyr::select(rowname, coh2), by = "rowname") %>%
    dplyr::arrange(coh2) %>%
    dplyr::mutate(
      wave_1 = as.numeric(gsub("_.*", "", gsub("wave_", "", rowname))),
      wave_2 = as.numeric(gsub(".*_", "", rowname))
    ) %>%
    dplyr::left_join(data.frame(rowname = names(sub_grp), group = sub_grp), by = "rowname") %>%
    dplyr::group_by(group) %>%
    dplyr::filter(!duplicated(wave_1)) %>%
    dplyr::filter(!duplicated(wave_2)) %>%
    dplyr::slice_max(coh2, n = 1) %>%
    dplyr::ungroup()

  # max diff per group
  Nratio_transform <- Nratios_out %>%
    dplyr::mutate(diff = abs(wave_1 - wave_2)) %>%
    dplyr::group_by(group) %>%
    dplyr::slice(which.max(diff)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(across(c(wave_1, wave_2), ~ paste("wave", ., sep = "_"))) %>%
    dplyr::mutate(wave_sel = paste(wave_1, wave_2, sep = "_")) %>%
    dplyr::select(-diff)

  # outputs
  write.csv(Nratios_out, file.path(out_dir, paste0("Nratios_rep", rep_id, ".csv")), row.names = FALSE)
  write.csv(Nratio_transform, file.path(out_dir, paste0("Nratio_transform_rep", rep_id, ".csv")), row.names = FALSE)

  com_col <- names(wave)[names(wave) %in% Nratio_transform$wave_sel]
  Nwave_pheno <- wave %>% dplyr::select(1:10, all_of(com_col))
  write.csv(Nwave_pheno, file.path(out_dir, paste0("Nwave_pheno_rep", rep_id, ".csv")), row.names = FALSE)

  invisible(list(
    sel = sel,
    Nratios_out = Nratios_out,
    Nratio_transform = Nratio_transform,
    Nwave_pheno = Nwave_pheno
  ))
}


slapreprocessing_allrep<- function(rep_id,
                        breakdown_dir = "./output",
                        phenotypes_path = "./data/phenotypes_whole.csv",
                        sample_path = "./data/sample_per_rep.csv",
                        out_dir = "./output",
                        k_groups = 3) {
  dir.create(out_dir, showWarnings = FALSE)

  # ---------------- inputs ----------------
  breakdown_path <- file.path(breakdown_dir, paste0("sla_breakdown", rep_id, ".csv"))
  S <- data.table::fread(breakdown_path, data.table = FALSE)

  phenotypes_whole <- read.table(phenotypes_path, header = TRUE, sep = ",")
  indv_sample <- read.csv(sample_path)

  rep_col <- paste0("Rep_", rep_id)
  phenotypes <- phenotypes_whole %>% dplyr::filter(Name2 %in% indv_sample[[rep_col]])
  wave <- phenotypes %>% janitor::clean_names()

  # ---------------- QC (same logic) ----------------
  S <- S %>%
    dplyr::mutate(dplyr::across(
      coh2:vartrait,
      ~ replace(
        .,
        which(
          corgblup < -1 |
            corgblup > 1 |
            corg < -1 |
            corg > 1 |
            (abs(corg - corgblup) > 0.2)
        ),
        NA
      )
    ))

  ntraits <- ceiling(sum(!is.na(S$coh2)) * 0.01)

  sel <- S %>%
    dplyr::slice_max(coh2, n = ntraits) %>%
    dplyr::mutate(dplyr::across(where(is.double), \(x) round(x, 4)))

  # ---------------- create ratio columns in wave ----------------
  sel$wave_1 <- gsub("\\s+", "_", as.character(sel$wave_1))
  sel$wave_2 <- gsub("\\s+", "_", as.character(sel$wave_2))

  for (i in seq_len(nrow(sel))) {
    w1 <- sel$wave_1[i]                 # e.g. "wave_1726"
    w2 <- sel$wave_2[i]                 # e.g. "wave_1363"
    out <- paste0(w1, "_", w2)          # e.g. "wave_1726_wave_1363"

    if (!w1 %in% names(wave) || !w2 %in% names(wave)) {
      stop(paste("Missing wave columns in `wave`:", w1, "or", w2, "at i =", i))
    }

    wave[[out]] <- wave[[w1]] / wave[[w2]]
  }

  sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
  sel$rowname <- gsub("\\s+", "_", sel$rowname)

  # ---------------- automatic contiguous range selection ----------------
  first_ratio <- sel$rowname[1]
  last_ratio  <- sel$rowname[nrow(sel)]

  if (!first_ratio %in% names(wave)) stop(paste("first_ratio not found in wave:", first_ratio))
  if (!last_ratio  %in% names(wave)) stop(paste("last_ratio not found in wave:", last_ratio))

  Sratios <- wave %>% dplyr::select(dplyr::all_of(first_ratio):dplyr::all_of(last_ratio))

  Sratios_t <- t(Sratios)
  colnames(Sratios_t) <- wave$plot_id

  # ---------------- clustering ----------------
  dist <- factoextra::get_dist(Sratios_t, method = "pearson")
  hcS <- hclust(dist, method = "average")
  sub_grp <- cutree(hcS, k = k_groups)

  # ---------------- select reps per cluster (same logic) ----------------
  Sratios_out <- Sratios_t %>%
    as.data.frame() %>%
    tibble::rownames_to_column("rowname") %>%
    dplyr::left_join(sel %>% dplyr::select(rowname, coh2), by = "rowname") %>%
    dplyr::arrange(coh2) %>%
    dplyr::mutate(
      wave_1 = as.numeric(gsub("_.*", "", gsub("wave_", "", rowname))),
      wave_2 = as.numeric(gsub(".*_", "", rowname))
    ) %>%
    dplyr::left_join(data.frame(rowname = names(sub_grp), group = sub_grp), by = "rowname") %>%
    dplyr::group_by(group) %>%
    dplyr::filter(!duplicated(wave_1)) %>%
    dplyr::filter(!duplicated(wave_2)) %>%
    dplyr::slice_max(coh2, n = 1) %>%
    dplyr::ungroup()

  # ---------------- choose max diff per group ----------------
  Sratio_transform <- Sratios_out %>%
    dplyr::mutate(diff = abs(wave_1 - wave_2)) %>%
    dplyr::group_by(group) %>%
    dplyr::slice(which.max(diff)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(across(c(wave_1, wave_2), ~ paste("wave", ., sep = "_"))) %>%
    dplyr::mutate(wave_sel = paste(wave_1, wave_2, sep = "_")) %>%
    dplyr::select(-diff)

  # ---------------- outputs ----------------
  write.csv(Sratios_out, file.path(out_dir, paste0("Sratios_rep", rep_id, ".csv")), row.names = FALSE)
  write.csv(Sratio_transform, file.path(out_dir, paste0("Sratio_transform_rep", rep_id, ".csv")), row.names = FALSE)

  com_col <- colnames(wave)[colnames(wave) %in% Sratio_transform$wave_sel]

  #select(1:9, sla, com_col)
  Swave_pheno <- wave %>% dplyr::select(1:9, sla, all_of(com_col))
  write.csv(Swave_pheno, file.path(out_dir, paste0("Swave_pheno_rep", rep_id, ".csv")), row.names = FALSE)

  invisible(list(
    sel = sel,
    Sratios_out = Sratios_out,
    Sratio_transform = Sratio_transform,
    Swave_pheno = Swave_pheno
  ))
}


pnpreprocessing_allrep <- function(rep_id,
                       breakdown_dir = "./output",
                       phenotypes_path = "./data/phenotypes_whole.csv",
                       sample_path = "./data/sample_per_rep.csv",
                       out_dir = "./output",
                       k_groups = 3,
                       pn_trait_col = "fs_plsr_narea") {

  dir.create(out_dir, showWarnings = FALSE)

  # ---------------- inputs ----------------
  breakdown_path <- file.path(breakdown_dir, paste0("pn_breakdown", rep_id, ".csv"))
  pn <- data.table::fread(breakdown_path, data.table = FALSE)

  phenotypes_whole <- read.table(phenotypes_path, header = TRUE, sep = ",")
  indv_sample <- read.csv(sample_path)

  rep_col <- paste0("Rep_", rep_id)
  phenotypes <- phenotypes_whole %>% dplyr::filter(Name2 %in% indv_sample[[rep_col]])
  wave <- phenotypes %>% janitor::clean_names()

  # ---------------- QC mask (same as your script) ----------------
  pn <- pn %>%
    dplyr::mutate(dplyr::across(
      coh2:vartrait,
      ~ replace(
        .,
        which(
          corgblup < -1 |
            corgblup > 1 |
            corg < -1 |
            corg > 1 |
            (abs(corg - corgblup) > 0.2)
        ),
        NA
      )
    ))

  ntraits <- ceiling(sum(!is.na(pn$coh2)) * 0.01)

  sel <- pn %>%
    dplyr::slice_max(coh2, n = ntraits) %>%
    dplyr::mutate(dplyr::across(where(is.double), \(x) round(x, 4)))

  # ---------------- create ratio columns in wave (same math as eval(parse())) ----------------
  sel$wave_1 <- gsub("\\s+", "_", as.character(sel$wave_1))
  sel$wave_2 <- gsub("\\s+", "_", as.character(sel$wave_2))

  for (i in seq_len(nrow(sel))) {
    w1 <- sel$wave_1[i]                # e.g. "wave_823"
    w2 <- sel$wave_2[i]                # e.g. "wave_726"
    out <- paste0(w1, "_", w2)         # e.g. "wave_823_wave_726"

    if (!w1 %in% names(wave) || !w2 %in% names(wave)) {
      stop(paste("Missing wave columns in `wave`:", w1, "or", w2, "at i =", i))
    }

    wave[[out]] <- wave[[w1]] / wave[[w2]]
  }

  # rowname must match ratio columns in wave
  sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
  sel$rowname <- gsub("\\s+", "_", sel$rowname)

  # ---------------- automatic contiguous range selection ----------------
  first_ratio <- sel$rowname[1]
  last_ratio  <- sel$rowname[nrow(sel)]

  if (!first_ratio %in% names(wave)) stop(paste("first_ratio not found in wave:", first_ratio))
  if (!last_ratio  %in% names(wave)) stop(paste("last_ratio not found in wave:", last_ratio))

  pnratios <- wave %>% dplyr::select(dplyr::all_of(first_ratio):dplyr::all_of(last_ratio))

  pnratios_t <- t(pnratios)
  colnames(pnratios_t) <- wave$plot_id

  # ---------------- clustering ----------------
  dist <- factoextra::get_dist(pnratios_t, method = "pearson")
  hcpn <- hclust(dist, method = "average")
  sub_grp <- cutree(hcpn, k = k_groups)

  # ---------------- representative per cluster ----------------
  pnratios_out <- pnratios_t %>%
    as.data.frame() %>%
    tibble::rownames_to_column("rowname") %>%
    dplyr::left_join(sel %>% dplyr::select(rowname, coh2), by = "rowname") %>%
    dplyr::arrange(coh2) %>%
    dplyr::mutate(
      wave_1 = as.numeric(gsub("_.*", "", gsub("wave_", "", rowname))),
      wave_2 = as.numeric(gsub(".*_", "", rowname))
    ) %>%
    dplyr::left_join(data.frame(rowname = names(sub_grp), group = sub_grp), by = "rowname") %>%
    dplyr::group_by(group) %>%
    dplyr::filter(!duplicated(wave_1)) %>%
    dplyr::filter(!duplicated(wave_2)) %>%
    dplyr::slice_max(coh2, n = 1) %>%
    dplyr::ungroup()

  # ---------------- select max diff per group ----------------
  pnratio_transform <- pnratios_out %>%
    dplyr::mutate(diff = abs(wave_1 - wave_2)) %>%
    dplyr::group_by(group) %>%
    dplyr::slice(which.max(diff)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(across(c(wave_1, wave_2), ~ paste("wave", ., sep = "_"))) %>%
    dplyr::mutate(wave_sel = paste(wave_1, wave_2, sep = "_")) %>%
    dplyr::select(-diff)

  # ---------------- outputs (match your pn filenames) ----------------
  write.csv(pnratios_out, file.path(out_dir, paste0("pnratios_rep", rep_id, ".csv")), row.names = FALSE)
  write.csv(pnratio_transform, file.path(out_dir, paste0("pnratio_transform_rep", rep_id, ".csv")), row.names = FALSE)

  # final phenotype table (same as your script)
  com_col <- colnames(wave)[colnames(wave) %in% pnratio_transform$wave_sel]

  # wave is clean_names()'d, so make sure pn trait col exists (or throw a clear error)
  if (!pn_trait_col %in% names(wave)) {
    stop(paste0(
      "Trait column `", pn_trait_col, "` not found in `wave`. ",
      "Available columns include: ", paste(head(names(wave), 30), collapse = ", "), " ..."
    ))
  }

  pnwave_pheno <- wave %>% dplyr::select(1:9, dplyr::all_of(pn_trait_col), dplyr::all_of(com_col))
  write.csv(pnwave_pheno, file.path(out_dir, paste0("pnwave_pheno_rep", rep_id, ".csv")), row.names = FALSE)

  invisible(list(
    sel = sel,
    pnratios_out = pnratios_out,
    pnratio_transform = pnratio_transform,
    pnwave_pheno = pnwave_pheno
  ))
}

pspreprocessing_allrep <- function(rep_id,
                       breakdown_dir = "./output",
                       phenotypes_path = "./data/phenotypes_whole.csv",
                       sample_path = "./data/sample_per_rep.csv",
                       out_dir = "./output",
                       k_groups = 3,
                       ps_trait_col = "plsr_sla_sorghum") {

  dir.create(out_dir, showWarnings = FALSE)

  # ---------------- inputs ----------------
  breakdown_path <- file.path(breakdown_dir, paste0("ps_breakdown", rep_id, ".csv"))
  ps <- data.table::fread(breakdown_path, data.table = FALSE)

  phenotypes_whole <- read.table(phenotypes_path, header = TRUE, sep = ",")
  indv_sample <- read.csv(sample_path)

  rep_col <- paste0("Rep_", rep_id)
  phenotypes <- phenotypes_whole %>%
    dplyr::filter(Name2 %in% indv_sample[[rep_col]])

  wave <- phenotypes %>% janitor::clean_names()

  # ---------------- QC masking ----------------
  ps <- ps %>%
    dplyr::mutate(dplyr::across(
      coh2:vartrait,
      ~ replace(
        .,
        which(
          corgblup < -1 |
            corgblup > 1 |
            corg < -1 |
            corg > 1 |
            (abs(corg - corgblup) > 0.2)
        ),
        NA
      )
    ))

  ntraits <- ceiling(sum(!is.na(ps$coh2)) * 0.01)

  sel <- ps %>%
    dplyr::slice_max(coh2, n = ntraits) %>%
    dplyr::mutate(dplyr::across(where(is.double), \(x) round(x, 4)))

  # ---------------- create ratio columns (same math as eval(parse())) ----------------
  sel$wave_1 <- gsub("\\s+", "_", as.character(sel$wave_1))
  sel$wave_2 <- gsub("\\s+", "_", as.character(sel$wave_2))

  for (i in seq_len(nrow(sel))) {
    w1 <- sel$wave_1[i]
    w2 <- sel$wave_2[i]
    out <- paste0(w1, "_", w2)

    if (!w1 %in% names(wave) || !w2 %in% names(wave)) {
      stop(paste("Missing wave columns in `wave`:", w1, "or", w2, "at i =", i))
    }

    wave[[out]] <- wave[[w1]] / wave[[w2]]
  }

  sel$rowname <- paste0(sel$wave_1, "_", sel$wave_2)
  sel$rowname <- gsub("\\s+", "_", sel$rowname)

  # ---------------- automatic contiguous selection ----------------
  first_ratio <- sel$rowname[1]
  last_ratio  <- sel$rowname[nrow(sel)]

  if (!first_ratio %in% names(wave))
    stop(paste("first_ratio not found in wave:", first_ratio))
  if (!last_ratio %in% names(wave))
    stop(paste("last_ratio not found in wave:", last_ratio))

  psratios <- wave %>%
    dplyr::select(dplyr::all_of(first_ratio):dplyr::all_of(last_ratio))

  psratios_t <- t(psratios)
  colnames(psratios_t) <- wave$plot_id

  # ---------------- clustering ----------------
  dist <- factoextra::get_dist(psratios_t, method = "pearson")
  hcps <- hclust(dist, method = "average")
  sub_grp <- cutree(hcps, k = k_groups)

  # ---------------- representative per cluster ----------------
  psratios_out <- psratios_t %>%
    as.data.frame() %>%
    tibble::rownames_to_column("rowname") %>%
    dplyr::left_join(sel %>% dplyr::select(rowname, coh2), by = "rowname") %>%
    dplyr::arrange(coh2) %>%
    dplyr::mutate(
      wave_1 = as.numeric(gsub("_.*", "", gsub("wave_", "", rowname))),
      wave_2 = as.numeric(gsub(".*_", "", rowname))
    ) %>%
    dplyr::left_join(data.frame(rowname = names(sub_grp), group = sub_grp), by = "rowname") %>%
    dplyr::group_by(group) %>%
    dplyr::filter(!duplicated(wave_1)) %>%
    dplyr::filter(!duplicated(wave_2)) %>%
    dplyr::slice_max(coh2, n = 1) %>%
    dplyr::ungroup()

  # ---------------- select max diff per group ----------------
  psratio_transform <- psratios_out %>%
    dplyr::mutate(diff = abs(wave_1 - wave_2)) %>%
    dplyr::group_by(group) %>%
    dplyr::slice(which.max(diff)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(across(c(wave_1, wave_2), ~ paste("wave", ., sep = "_"))) %>%
    dplyr::mutate(wave_sel = paste(wave_1, wave_2, sep = "_")) %>%
    dplyr::select(-diff)

  # ---------------- outputs ----------------
  write.csv(psratios_out,
            file.path(out_dir, paste0("psratios_rep", rep_id, ".csv")),
            row.names = FALSE)

  write.csv(psratio_transform,
            file.path(out_dir, paste0("psratio_transform_rep", rep_id, ".csv")),
            row.names = FALSE)

  com_col <- colnames(wave)[colnames(wave) %in% psratio_transform$wave_sel]

  if (!ps_trait_col %in% names(wave)) {
    stop(paste("Trait column not found in wave:", ps_trait_col))
  }

  pswave_pheno <- wave %>%
    dplyr::select(1:9, dplyr::all_of(ps_trait_col), dplyr::all_of(com_col))

  write.csv(pswave_pheno,
            file.path(out_dir, paste0("pswave_pheno_rep", rep_id, ".csv")),
            row.names = FALSE)

  invisible(list(
    sel = sel,
    psratios_out = psratios_out,
    psratio_transform = psratio_transform,
    pswave_pheno = pswave_pheno
  ))
}



#for processing the blues and filtering the file with name west data
process_narea_blues <- function(rep_id, kin, names_west) {

  # Read EF & MW
  NbluesEF <- read.csv(paste0("./output/NbluesEF_rep", rep_id, ".csv"))
  NbluesMW <- read.csv(paste0("./output/NbluesMW_rep", rep_id, ".csv"))

  # Combine
  Nblues_rep <- bind_rows(
    "EF" = NbluesEF,
    "MW" = NbluesMW,
    .id = "env"
  )

  # Join + mutate
  Nblues_rep <- Nblues_rep |>
    left_join(names_west %>% dplyr::select(name2, Corrected_names)) |>
    mutate(
      taxa = ifelse(
        (name2 != Corrected_names) & grepl("PI", Corrected_names),
        NA,
        Corrected_names
      )
    )

  # Filter
  N_blues <- subset(Nblues_rep, !is.na(taxa)) |>
    select(-name2, -Corrected_names)

  N_blues <- droplevels(
    N_blues[N_blues$taxa %in% rownames(kin), ]
  )

  # Write
  write.csv(
    N_blues,
    paste0("./output/N_blues_rep", rep_id, ".csv"),
    row.names = FALSE
  )
}

#process the Sla blues
process_sla_blues <- function(rep_id, kin, names_west) {
  # Read EF & MW
  sblues_ef <- read.csv(paste0("./output/SbluesEF_rep", rep_id, ".csv"))
  sblues_mw <- read.csv(paste0("./output/SbluesMW_rep", rep_id, ".csv"))

  # Combine
  sblues_rep <- bind_rows(
    EF = sblues_ef,
    MW = sblues_mw,
    .id = "env"
  )

  # Join + mutate (same logic)
  sblues_rep <- sblues_rep |>
    left_join(names_west |> dplyr::select(name2, Corrected_names)) |>
    mutate(
      taxa = ifelse(
        (name2 != Corrected_names) & grepl("PI", Corrected_names),
        NA,
        Corrected_names
      )
    )

  # Filter (same logic)
  s_blues <- subset(sblues_rep, !is.na(taxa)) |>
    select(-name2, -Corrected_names)

  s_blues <- droplevels(
    s_blues[s_blues$taxa %in% rownames(kin), ]
  )

  # Write
  write.csv(
    s_blues,
    paste0("./output/S_blues_rep", rep_id, ".csv"),
    row.names = FALSE
  )

  invisible(s_blues)
}

process_pn_blues <- function(rep_id, kin, names_west) {
  # Read EF & MW
  pnblues_ef <- read.csv(paste0("./output/pnbluesEF_rep", rep_id, ".csv"))
  pnblues_mw <- read.csv(paste0("./output/pnbluesMW_rep", rep_id, ".csv"))

  # Combine
  pnblues_rep <- bind_rows(
    EF = pnblues_ef,
    MW = pnblues_mw,
    .id = "env"
  )

  # Join + mutate (same logic)
  pnblues_rep <- pnblues_rep |>
    left_join(names_west |> dplyr::select(name2, Corrected_names)) |>
    mutate(
      taxa = ifelse(
        (name2 != Corrected_names) & grepl("PI", Corrected_names),
        NA,
        Corrected_names
      )
    )

  # Filter (same logic)
  pn_blues <- subset(pnblues_rep, !is.na(taxa)) |>
    select(-name2, -Corrected_names)

  pn_blues <- droplevels(
    pn_blues[pn_blues$taxa %in% rownames(kin), ]
  )

  # Write
  write.csv(
    pn_blues,
    paste0("./output/pn_blues_rep", rep_id, ".csv"),
    row.names = FALSE
  )

  invisible(pn_blues)
}

process_ps_blues <- function(rep_id, kin, names_west) {
  # Read EF & MW
  psblues_ef <- read.csv(paste0("./output/psbluesEF_rep", rep_id, ".csv"))
  psblues_mw <- read.csv(paste0("./output/psbluesMW_rep", rep_id, ".csv"))

  # Combine
  psblues_rep <- bind_rows(
    EF = psblues_ef,
    MW = psblues_mw,
    .id = "env"
  )

  # Join + mutate (same logic)
  psblues_rep <- psblues_rep |>
    left_join(names_west |> dplyr::select(name2, Corrected_names)) |>
    mutate(
      taxa = ifelse(
        (name2 != Corrected_names) & grepl("PI", Corrected_names),
        NA,
        Corrected_names
      )
    )

  # Filter (same logic)
  ps_blues <- subset(psblues_rep, !is.na(taxa)) |>
    select(-name2, -Corrected_names)

  ps_blues <- droplevels(
    ps_blues[ps_blues$taxa %in% rownames(kin), ]
  )

  # Write
  write.csv(
    ps_blues,
    paste0("./output/ps_blues_rep", rep_id, ".csv"),
    row.names = FALSE
  )

  invisible(ps_blues)
}






get_coh2_0_syntrait <- function(
    trait_col,
    breakdown_prefix,
    breakdown_dir = "./output",
    output_prefix = breakdown_prefix,
    output_dir = "./output",
    reps = 1:5,
    indv_sample_path = "./data/sample_per_rep.csv",
    phenotypes_path = "./data/phenotypes_whole.csv",
    keep_cols = 1:9
) {
  phenotypes_whole <- read.table(phenotypes_path, header = TRUE, sep = ",")
  indv_sample <- read.csv(indv_sample_path)

  selected_list <- vector("list", length(reps))

  for (i in seq_along(reps)) {
    rep <- reps[i]
    rep_col <- paste0("Rep_", rep)

    breakdown_path <- file.path(breakdown_dir, paste0(breakdown_prefix, "_breakdown", rep, ".csv"))
    out_path <- file.path(output_dir, paste0(output_prefix, "_rep", rep, "_lowcoh2.csv"))

    if (!file.exists(breakdown_path)) stop("Missing file: ", breakdown_path)
    if (!rep_col %in% names(indv_sample)) stop("Missing column in sample file: ", rep_col)

    B <- data.table::fread(breakdown_path, data.table = FALSE)
    if (!all(c("coh2", "wave_1", "wave_2") %in% names(B))) {
      stop("Breakdown file must have columns coh2, wave_1, wave_2: ", breakdown_path)
    }

    wave <- phenotypes_whole |>
      dplyr::filter(!(Name2 %in% indv_sample[[rep_col]])) |>
      janitor::clean_names()

    if (!trait_col %in% names(wave)) stop("Trait column not found in phenotypes: ", trait_col)

    selected <- B |>
      dplyr::mutate(
        coh2 = as.numeric(trimws(as.character(coh2))),
        abs_coh2 = abs(coh2)
      ) |>
      dplyr::slice_min(abs_coh2, n = 1, with_ties = FALSE)

    selected$rep <- rep
    selected_list[[i]] <- selected

    wave_a <- janitor::make_clean_names(selected$wave_1[1])
    wave_b <- janitor::make_clean_names(selected$wave_2[1])
    ratio_name <- paste(wave_a, wave_b, sep = "_")

    missing_wave_cols <- setdiff(c(wave_a, wave_b), names(wave))
    if (length(missing_wave_cols) > 0) {
      stop("Missing wave columns in phenotypes: ", paste(missing_wave_cols, collapse = ", "))
    }

    out <- wave |>
      dplyr::select(dplyr::all_of(keep_cols), dplyr::all_of(trait_col), dplyr::all_of(c(wave_a, wave_b))) |>
      dplyr::mutate(!!ratio_name := .data[[wave_a]] / .data[[wave_b]]) |>
      dplyr::select(-dplyr::all_of(c(wave_a, wave_b)))

    write.csv(out, out_path, row.names = FALSE)
  }

  selected_summary <- dplyr::bind_rows(selected_list) |>
    dplyr::select(rep, wave_1, wave_2, coh2)

  cat("\n===== Selected traits (coh2 closest to 0) for", breakdown_prefix, "=====\n")
  print(selected_summary)

  invisible(selected_summary)
}

#getting coh2=0 synthetic traits blue for all rep all traits
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



#cleaning the blues with name west and kin data
library(data.table)
library(dplyr)

clean_with_names_west_and_kin_all_reps <- function(
    out_prefix,                     # "N" or "S" or "pn" or "ps"
    reps = 1:5,
    out_dir = "./output",
    names_west_path = "./data/Names_WEST_SF.csv",
    kin                             # pass your kin object
) {
  # Load + prep Names_WEST exactly like your scripts
  names_west <- fread(names_west_path)
  names_west$Name2 <- gsub(" ", "", names_west$Name2)
  colnames(names_west)[colnames(names_west) == "Name2"] <- "name2"

  if (is.null(rownames(kin))) {
    stop("kin must have rownames (taxa IDs) for filtering.")
  }

  for (rep in reps) {
    ef_path <- file.path(out_dir, paste0(out_prefix, "bluesEF1_rep", rep, ".csv"))
    mw_path <- file.path(out_dir, paste0(out_prefix, "bluesMW1_rep", rep, ".csv"))

    if (!file.exists(ef_path)) stop("Missing file: ", ef_path)
    if (!file.exists(mw_path)) stop("Missing file: ", mw_path)

    blues_ef <- fread(ef_path)
    blues_mw <- fread(mw_path)

    # combine EF + MW like your script
    blues_all <- bind_rows(
      "EF" = blues_ef,
      "MW" = blues_mw,
      .id = "env"
    )

    # join Names_WEST + apply your PI rule
    blues_all <- blues_all |>
      left_join(names_west |> dplyr::select(name2, Corrected_names), by = "name2") |>
      mutate(
        taxa = ifelse(
          (name2 != Corrected_names) & grepl("PI", Corrected_names),
          NA,
          Corrected_names
        )
      )

    # remove NA taxa, drop helper columns
    out_final <- blues_all |>
      filter(!is.na(taxa)) |>
      select(-name2, -Corrected_names)

    # filter to taxa in kin rownames (exactly like your script)
    out_final <- droplevels(out_final[out_final$taxa %in% rownames(kin), ])

    # write output with the same old naming pattern
    out_path <- file.path(out_dir, paste0(out_prefix, "_blues1_rep", rep, ".csv"))
    write.csv(out_final, out_path, row.names = FALSE)
  }

  invisible(TRUE)
}


#creating sort
# ---------------------creating list with 5 fold and 20 reps----------------------
create_folds<- function(individuals, nfolds, reps, seed = 123){
  library(cvTools)
  library(dplyr)
  set.seed(seed)
  sort<-list()
  individuals <- as.factor(individuals)
  nl <- length(unique(individuals))
  for(a in 1:reps){
    folds <- cvFolds(nl,type ="random",
                     K= nfolds)
    Sample<-cbind(folds$which,folds$subsets)
    cv<-split(levels(individuals)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of taxa
    sort[[a]]<-cv
  }
  return(sort)
}

create_samples_df <- function(individuals, percent = 20, reps = 20, seed = 123) {
  set.seed(seed)  # Set a seed for reproducibility

  total_individuals <- length(individuals)
  sample_size <- ceiling(percent/100 * total_individuals)  # Calculate sample size as 20% of total

  samples_list <- vector("list", reps)  # Initialize list to store samples

  for (i in 1:reps) {
    samples_list[[i]] <- sample(individuals, sample_size)  # Randomly sample individuals without replacement
  }
  samples_df <- do.call(cbind.data.frame, samples_list)
  names(samples_df) <- paste("Rep", 1:reps, sep = "_")

  return(samples_df)
}


# ---------------------cross validation function---------------------
crossv <- function(sort,
                   train,
                   validation,
                   kin,
                   mytrait,
                   scheme
) {
  library(asreml)
  library(dplyr)
  library(data.table)
  ac <- c()
  gebv <- list()

  for(j in 1:length(sort)){
    r<-list()
    for(i in 1:5){
      test<-na.omit(train)
      test[test$taxa %in% sort[[j]][[i]], mytrait ] <- NA
      cat(i, j,"\n")
      # ---------------------fitting model---------------------

      model <- asreml(
        fixed = as.formula(paste0(mytrait, " ~ 1")),
        random =  ~ vm(taxa,source = kin, singG= "NSD" ),
        data = test, na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "taxa"))

      # Update if necessary
      if(!model$converge){ model <- update.asreml(model) }
      if(!model$converge){ model <- update.asreml(model) }

      r[[i]] <- model$predictions$pvals[model$predictions$pvals$taxa %in% sort[[j]][[i]],1:2]
    }

    gebv[[j]] <- Reduce(rbind, r)
    gebv[[j]] <- gebv[[j]] %>%
      left_join(validation[,c("taxa", mytrait)]) %>%
      mutate(rep = j)
    ac[j]<-cor(gebv[[j]][,2], gebv[[j]][,3], use = "complete.obs")
  }
  gebv<- bind_rows(gebv, .id = "rep")
  write.csv(gebv, paste0("./output/","gebv_",mytrait, "_", scheme, ".csv"),row.names=F)
  result<- list(gebv= gebv,ac= ac)
  return(result)
}
