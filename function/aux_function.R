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
