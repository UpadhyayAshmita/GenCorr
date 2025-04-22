library(data.table)
library(lme4)
library(tidyverse)
library(furrr)
library(janitor)
source("aux_functions_lme4.R")

phenotypes <- read.table("./data/phenotypes1.csv", header = T, sep = ",")

#get input variables
args <- commandArgs(trailingOnly = TRUE)
trait <- args[1]
start <- args[2]
end <- args[3]

phenotypes <- phenotypes |>
  clean_names()|>
  mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>
  arrange(loc, range, row) |>
  rename(pn = fs_plsr_narea, ps = plsr_sla_sorghum, pnm = fs_plsr_nmass, np = n_perc)

# number of reps
nrep <- psych::harmonic.mean(c(rep(1, 847), rep(16, 6), rep(17, 1)))

# --------------------------------------------------------
# --------------------- biological traits ----------------
# --------------------------------------------------------
if (trait == "narea"){
	#phenotypes["narea"] = scale(phenotypes["narea"])
	model <- lmer(narea ~ loc + set + loc:set + (1 | loc:block) + (1 | name2) + (1 | loc:name2), data = phenotypes)
	blups <- ranef(model)$name2 |>  rownames_to_column("name2")
} else if (trait == "sla"){
	#phenotypes["sla"] = scale(phenotypes["sla"])
	model <- lmer(sla ~ loc + set + loc:set + (1 | loc:block) + (1 | name2) + (1 | loc:name2), data = phenotypes)
	blups <- ranef(model)$name2 |>  rownames_to_column("name2")
} else if (trait == "pn") {
	model <- lmer(pn ~ loc + set + loc:set + (1 | loc:block) + (1 | name2) + (1 | loc:name2), data = phenotypes)
	blups <- ranef(model)$name2 |>  rownames_to_column("name2")
} else if (trait == "ps"){
	#phenotypes["ps"] = scale(phenotypes["ps"])
	model <- lmer(ps ~ loc + set + loc:set + (1 | loc:block) + (1 | name2) + (1 | loc:name2), data = phenotypes)
	blups <- ranef(model)$name2 |>  rownames_to_column("name2")
} else if (trait == "pnm") {
	model <- lmer(pnm ~ loc + set + loc:set + (1 | loc:block) + (1 | name2) + (1 | loc:name2), data = phenotypes)
	blups <- ranef(model)$name2 |>  rownames_to_column("name2")
} else if (trait == "np") {
	#phenotypes["ps"] = scale(phenotypes["ps"])
	model <- lmer(pn ~ loc + set + loc:set + (1 | loc:block) + (1 | name2) + (1 | loc:name2), data = phenotypes)
	blups <- ranef(model)$name2 |>  rownames_to_column("name2")
}
print(model)
h2 <- geth2(model)
vartrait <- varg(model)

#print(coh2(wave_1=1425, wave_2=1472, trait="np", flip=0))
#print(coh2(wave_1=2300, wave_2=1500, trait="narea", flip=0))
#print(coh2(wave_1=2300, wave_2=1500, trait="narea", flip=1))
#q()

result_df <- data.frame()

for (i in start:(end)){
  for (j in 350:2500) {
    if (i != j) {
      result <- coh2(wave_1 = i, wave_2 = j, trait = trait, flip=0)
      result_df <- rbind(result_df, result)
    }
    if (j%%50==0){
      cat(trait, i, j, "\n", sep=",")
      filename = paste("./output1/", trait, "_", start, "_", end, ".csv", sep="")
      write.table(result_df, file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
      result_df <-data.frame()
    }
  }
}

# Combine the results into a single data frame
#library(dplyr)
#combined_df <- rbind(result_list)

# View the combined data frame
#head(result_df)
cat(start, end, "done...\n", sep=",")
filename = paste("./output1/", trait, "_", start, "_", end, ".csv", sep="")
write.table(result_df, file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
