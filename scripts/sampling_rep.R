#filter data based on each Rep for 165 unique indv
phenotypes_whole <- read.table("./data/phenotypes_whole.csv", header = T, sep = ",")
Name2<- unique(phenotypes_whole$Name2)
# ---------------------creating list with random 20% sample of total Name2 with 3 reps
Name2<- unique(phenotypes_whole$Name2)
random_samples <- function(individuals, percentage, reps, seed = 123) {
  library(dplyr)
  set.seed(seed)
  individuals <- as.factor(individuals)
  n <- length(individuals)
  sample_size <- ceiling(percentage / 100 * n)
  samples_df <- data.frame(matrix(nrow = sample_size, ncol = reps))
  names(samples_df) <- paste0("Rep_", 1:reps)
  for (i in 1:reps) {
    samples_df[[i]] <- sample(individuals, sample_size)  # Randomly sample individuals
  }
  return(samples_df)
}
sorted <- random_samples(individuals= Name2, 20, 5)  # 20% sample, 5 repetitions
print(sorted)
write.csv(sorted, "./data/indv_sample.csv", row.names = F)
indv_sample<- read.csv("./data/indv_sample.csv")
phenotypes<- phenotypes_whole %>% filter(Name2 %in% indv_sample$Rep_4)
write.csv(phenotypes, "./data/phenotypes4.csv", row.names=F)
