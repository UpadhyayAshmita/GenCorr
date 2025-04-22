library(tidyverse)
library(Metrics)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Function to clean MT column and extract coh2 values
clean_data <- function(data) {
  data %>%
    # Extract the coh2 value from the MT column and store in a new column 'coh2'
    mutate(coh2 = as.numeric(sub(".*coh2=([0-9.]+).*", "\\1", MT))) %>%
    # Remove everything after S1, S2, S3 in MT column
    mutate(MT = gsub("\\(.*\\)", "", MT)) %>%
    # Trim any whitespace in MT
    mutate(MT = trimws(MT))
}

#code for whole rep1-5 facet wrapped with its complete model
#plsr-narea
#ST for plsr-narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMW.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn.ac) %>% select(-result_pn.ac)
#Multi trait cv1 for plsr_narea NW1, NW2, NW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for pnarea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apnd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.54)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.53)", Type = "CV1")
acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.54)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.53)", Type = "CV2")
# Clean pn_data_whole
pn_data_whole <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1, acpnW1, acpnW2, acpnW3)
pn_data_whole$Data <- "Complete"
pn_data_whole <- clean_data(pn_data_whole)


#replication-1
#ST for pnarea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep1.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep1.ac) %>% select(-result_pn_rep1.ac)
#Multi trait cv1 for plsr-narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for plsr-narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.69)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.68)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.69)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.69)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.68)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.69)", Type = "CV2")
# Combine all data
pn_data_rep1 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep1$Data <- "Subset-1"
pn_data_rep1<- clean_data(pn_data_rep1)

#Set 2
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep2.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep2.ac) %>% select(-result_pn_rep2.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.64)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.60)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.64)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.60)", Type = "CV2")
# Combine all data
pn_data_rep2 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep2$Data <- "Subset-2"
pn_data_rep2<- clean_data(pn_data_rep2)

#Set 3
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep3.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep3.ac) %>% select(-result_pn_rep3.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.51)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.57)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.48)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.51)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.57)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.48)", Type = "CV2")
# Combine all data
pn_data_rep3 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep3$Data <- "Subset-3"
pn_data_rep3<- clean_data(pn_data_rep3)

#Set 4
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep4.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep4.ac) %>% select(-result_pn_rep4.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
pn_data_rep4 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep4$Data <- "Subset-4"
pn_data_rep4<- clean_data(pn_data_rep4)
#Set 5
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep5.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep5.ac) %>% select(-result_pn_rep5.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.59)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.59)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV2")
# Combine all data
pn_data_rep5 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep5$Data <- "Subset-5"
pn_data_rep5<- clean_data(pn_data_rep5)
# Combine all datasets
all_data <- bind_rows(pn_data_whole, pn_data_rep1, pn_data_rep2, pn_data_rep3, pn_data_rep4, pn_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Common theme to apply to all plots for better visibility
common_theme <- theme_bw(base_size = 20) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        legend.text = element_text(size = 20),    # Adjust legend text size
        axis.text = element_text(size = 20),      # Adjust axis text size
        axis.title = element_text(size = 20))     # Adjust axis title size

# Create a ggplot with facet_wrap
plot1 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
 common_theme +
  labs(title = "PLSR-Narea",#train EF val mw
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position = "right") # Bold facet labels for clarity

#ggsave("./figure/pn_compsubsetfig1.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#MWef
#plsr-narea
#ST for plsr-narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEF.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn.ac) %>% select(-result_pn.ac)
#Multi trait CV1_mwef for plsr_narea NW1, NW2, NW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for pnarea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apnd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.54)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.53)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.54)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.53)", Type = "CV2")

# Combine all data
pn_data_mwef <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1, acpnW1, acpnW2, acpnW3)
pn_data_mwef$Data<- "Complete"

#replication 1
#ST for pnarea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep1.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep1.ac) %>% select(-result_pn_rep1.ac)
#Multi trait cv1 for plsr-narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for plsr-narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.69)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.68)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.69)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.69)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.68)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.69)", Type = "CV2")
# Combinne all data
pn_data_rep1 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep2.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep2.ac) %>% select(-result_pn_rep2.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.64)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.60)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.64)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.60)", Type = "CV2")
# Combine all data
pn_data_rep2 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep3.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep3.ac) %>% select(-result_pn_rep3.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.51)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.57)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.48)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.51)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.57)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.48)", Type = "CV2")
# Combine all data
pn_data_rep3 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep4.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep4.ac) %>% select(-result_pn_rep4.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
pn_data_rep4 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep5.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep5.ac) %>% select(-result_pn_rep5.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("acpnW1_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("acpnW2_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("acpnW3_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("acpnW1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("acpnW2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("acpnW3_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.59)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.59)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV2")
# Combine all data
pn_data_rep5 <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1,acpnW1, acpnW2, acpnW3)
pn_data_rep5$Data <- "Subset-5"


# Combine all datasets
all_data <- bind_rows(pn_data_mwef, pn_data_rep1, pn_data_rep2, pn_data_rep3, pn_data_rep4, pn_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Create a ggplot with facet_wrap
plot5 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-Narea",#(Training MW & validating on Ef)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position = "right") # Bold facet labels for clarity

#ggsave("./figure/pn_compsubsetfig1mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#plsr-sla(EFMW)
#plsr-sla
#ST for PLSR-SLA
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMW.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps.ac) %>% select(-result_ps.ac)
#Multi trait cv1 for plsr_sla NW1, NW2, NW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for psarea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait apsd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.46)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.46)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.46)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.46)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.46)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.46)", Type = "CV2")
# Combine all data
ps_data_whole <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1, acpsW1, acpsW2, acpsW3)
ps_data_whole$Data<- "Complete"


#replication 1
#ST for pnarea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep1.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep1.ac) %>% select(-result_ps_rep1.ac)
#Multi trait cv1 for plsr-sla psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for plsr-sla psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.61)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.65)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.58)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.61)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.65)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.58)", Type = "CV2")
# Combinne all data
ps_data_rep1 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep2.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep2.ac) %>% select(-result_ps_rep2.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.57)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.59)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.57)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.57)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.59)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.57)", Type = "CV2")
# Combine all data
ps_data_rep2 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep3.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep3.ac) %>% select(-result_ps_rep3.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.47)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.47)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
ps_data_rep3 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep4.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep4.ac) %>% select(-result_ps_rep4.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.61)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.50)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.61)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.50)", Type = "CV2")
# Combine all data
ps_data_rep4 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep5.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep5.ac) %>% select(-result_ps_rep5.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.50)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.47)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.50)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.47)", Type = "CV2")
# Combine all data
ps_data_rep5 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep5$Data <- "Subset-5"
# Combine all datasets
all_data <- bind_rows(ps_data_whole, ps_data_rep1, ps_data_rep2, ps_data_rep3, ps_data_rep4, ps_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Create a ggplot with facet_wrap
plot2 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-SLA", #(Training Ef & validating on MW)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/ps_compsubsetfig1.jpeg", plot = plot, width = 10, height = 8, dpi = 300)

#plsr-sla(MWEF)
#plsr sla
#ST for PLSR-SLA
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEF.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps.ac) %>% select(-result_ps.ac)
#Multi trait cv1 for plsr_sla NW1, NW2, NW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for psarea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apsd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.46)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.46)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.46)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.46)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.46)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.46)", Type = "CV2")

# Combine all data
ps_data_whole <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_whole$Data<- "Complete"

#replication 1
#ST for pnarea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep1.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep1.ac) %>% select(-result_ps_rep1.ac)
#Multi trait cv1 for plsr-sla psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for plsr-sla psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.61)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.65)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.58)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.61)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.65)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.58)", Type = "CV2")
# Combinne all data
ps_data_rep1 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep2.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep2.ac) %>% select(-result_ps_rep2.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.57)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.59)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.57)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.57)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.59)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.57)", Type = "CV2")
# Combine all data
ps_data_rep2 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep3.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep3.ac) %>% select(-result_ps_rep3.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.47)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.47)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
ps_data_rep3 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep4.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep4.ac) %>% select(-result_ps_rep4.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.61)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.50)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.61)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.50)", Type = "CV2")
# Combine all data
ps_data_rep4 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep5.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep5.ac) %>% select(-result_ps_rep5.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("acpsW1_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("acpsW2_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("acpsW1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("acpsW2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("acpsW3_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.50)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.47)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.50)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.47)", Type = "CV2")
# Combine all data
ps_data_rep5 <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_rep5$Data <- "Subset-5"


# Combine all datasets
all_data <- bind_rows(ps_data_whole, ps_data_rep1, ps_data_rep2, ps_data_rep3, ps_data_rep4, ps_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Create a ggplot with facet_wrap
plot6 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-SLA",#(Training MW & validating on Ef)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position = "right") # Bold facet labels for clarity
#ggsave("./figure/ps_compsubsetfig1mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#EFMW- narea
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMW.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N.ac) %>% select(-result_N.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2= 0.49)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2= 0.49)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2= 0.49)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2= 0.49)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV2")
# Combine all data
N_data_whole <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1,acNW2, acNW3)#acNT_CV1,acNT
N_data_whole$Data <- "Complete"

#replication 1
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep1.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep1.ac) %>% select(-result_N_rep1.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.63)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.63)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.63)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.63)", Type = "CV2")
# Combine all data
N_data_rep1 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep1$Data <- "Subset-1"

#replication 2
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep2.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep2.ac) %>% select(-result_N_rep2.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.66)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.58)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.55)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.66)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.58)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.55)", Type = "CV2")
# Combine all data
N_data_rep2 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep2$Data <- "Subset-2"
#replication 3
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep3.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep3.ac) %>% select(-result_N_rep3.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.52)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.56)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.51)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.52)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.56)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.51)", Type = "CV2")
# Combine all data
N_data_rep3 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep3$Data <- "Subset-3"
#replication 4
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep4.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep4.ac) %>% select(-result_N_rep4.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV2")
# Combine all data
N_data_rep4 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep4$Data <- "Subset-4"
#replication 5
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep5.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep5.ac) %>% select(-result_N_rep5.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.48)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.48)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
N_data_rep5 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep5$Data <- "Subset-5"
all_data <- bind_rows(N_data_whole, N_data_rep1, N_data_rep2, N_data_rep3, N_data_rep4, N_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))
# Create a ggplot with facet_wrap
# Create the ggplot with the adjusted legend order
plot3 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "Narea",#(Training Ef & validating on MW)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/Narea_compsubsetfig1.jpeg", plot = plot, width = 10, height = 8, dpi = 300)




#MWEF narea
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEF.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N.ac) %>% select(-result_N.ac)
#Multi trait CV1_mwef for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("acNW1_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("acNW2_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("acNW3_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("acNW1_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("acNW2_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("acNW3_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2= 0.49)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2= 0.49)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2= 0.49)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2= 0.49)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV2")

# Combine all data
N_data_mwef_whole <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1, acNW1, acNW2, acNW3)
N_data_mwef_whole$Data <- "Complete"

#Set 1 mwef
#ST for narea
corr_N_MWEF<- read.csv("./output/corr_N_MWEFrep1.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep1.ac) %>% select(-result_N_rep1.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.63)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.63)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.63)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.63)", Type = "CV2")
# Combine all data
N_data_rep1 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep1$Data <- "Subset-1"
#Set 2
#ST for narea
corr_N_MWEF<- read.csv("./output/corr_N_MWEFrep2.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep2.ac) %>% select(-result_N_rep2.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.66)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.58)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.55)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.66)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.58)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.55)", Type = "CV2")
# Combine all data
N_data_rep2 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep2$Data <- "Subset-2"
#Set 3
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep3.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep3.ac) %>% select(-result_N_rep3.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.52)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.56)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.51)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.52)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.56)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.51)", Type = "CV2")
# Combine all data
N_data_rep3 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep3$Data <- "Subset-3"

#Set 4
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep4.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep4.ac) %>% select(-result_N_rep4.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV2")
# Combine all data
N_data_rep4 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep4$Data <- "Subset-4"
#Set 5
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep5.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep5.ac) %>% select(-result_N_rep5.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.48)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.48)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
N_data_rep5 <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1, acNW2, acNW3)
N_data_rep5$Data <- "Subset-5"

# Combine all datasets
all_data <- bind_rows(N_data_mwef_whole, N_data_rep1, N_data_rep2, N_data_rep3, N_data_rep4, N_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Create a ggplot with facet_wrap
plot7 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "Narea",#(Training MW & validating on EF)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position = "right") # Bold facet labels for clarity

#ggsave("./figure/Narea_compsubsetfig1mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#sla-efmw
#ST for SLA EF
corr_S_EFMW<- read.csv("./output/corr_S_EFMW.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S.ac) %>% select(-result_S.ac)
#Multi trait cv1 for sla SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for sla SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2= 0.57)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2= 0.54)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2= 0.57)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2= 0.54)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV2")

# Combine all data
S_data_whole <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_whole$Data <- "Complete"


#replication 1
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep1.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep1.ac) %>% select(-result_S_rep1.ac)
#Multi trait cv1 for sla SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for sla SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.63)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.65)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.63)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.63)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.65)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.63)", Type = "CV2")
# CombiSe all data
S_data_rep1 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep1$Data <- "Subset-1"

#replication 2
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep2.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep2.ac) %>% select(-result_S_rep2.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.64)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.60)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.64)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.60)", Type = "CV2")
# Combine all data
S_data_rep2 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep2$Data <- "Subset-2"
#replication 3
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep3.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep3.ac) %>% select(-result_S_rep3.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.51)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.57)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.48)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.51)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.57)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.48)", Type = "CV2")
# Combine all data
S_data_rep3 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep3$Data <- "Subset-3"
#replication 4
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep4.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep4.ac) %>% select(-result_S_rep4.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
S_data_rep4 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep5.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep5.ac) %>% select(-result_S_rep5.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.59)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.59)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV2")
# Combine all data
S_data_rep5 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep5$Data <- "Subset-5"

# Combine all datasets
all_data <- bind_rows(S_data_whole, S_data_rep1, S_data_rep2, S_data_rep3, S_data_rep4, S_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Create a ggplot with facet_wrap
plot4 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme+
  labs(title =NULL,#(Training Ef & validating on MW)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/sla_compsubsetfig1.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#mwef-sla
#ST for SLA MWEF
corr_S_MWEF<- read.csv("./output/corr_S_MWEF.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S.ac) %>% select(-result_S.ac)
#Multi trait cv1 for sla SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for sla SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2= 0.57)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2= 0.54)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2= 0.57)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2= 0.54)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV2")

# Combine all data
S_data_mwef <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1, acSW1, acSW2, acSW3)
S_data_mwef$Data <- "Complete"

#ST for sla
corr_S_MWEF<- read.csv("./output/corr_S_MWEFrep1.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep1.ac) %>% select(-result_S_rep1.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.63)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.63)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.63)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.63)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.63)", Type = "CV2")
# Combine all data
S_data_rep1 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep1$Data <- "Subset-1"

#replication 2
#ST for narea
corr_S_MWEF<- read.csv("./output/corr_S_MWEFrep2.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep2.ac) %>% select(-result_S_rep2.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.66)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.58)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.55)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.66)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.58)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.55)", Type = "CV2")
# Combine all data
S_data_rep2 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep2$Data <- "Subset-2"


#replication 3
#ST for sla
corr_S_MWEF<- read.csv("./output/sla/corr_S_MWEFrep3.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep3.ac) %>% select(-result_S_rep3.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.52)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.56)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.51)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.52)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.56)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.51)", Type = "CV2")
# Combine all data
S_data_rep3 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep3$Data <- "Subset-3"
#replication 4
#ST for sla
corr_S_MWEF<- read.csv("./output/sla/corr_S_MWEFrep4.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep4.ac) %>% select(-result_S_rep4.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.50)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.52)", Type = "CV2")
# Combine all data
S_data_rep4 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep4$Data <- "Subset-4"
#replication 5
#ST for narea
corr_S_MWEF<- read.csv("./output/sla/corr_S_MWEFrep5.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep5.ac) %>% select(-result_S_rep5.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.48)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.48)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.48)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
S_data_rep5 <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_rep5$Data <- "Subset-5"

# Combine all dataset
all_data <- bind_rows(S_data_mwef, S_data_rep1, S_data_rep2, S_data_rep3, S_data_rep4, S_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Create a ggplot with facet_wrap
plot8 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "SLA",#(Training MW & validating on Ef)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position = "right") # Bold facet labels for clarity

#ggsave("./figure/sla_compsubsetfig1mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#making  four traits efmw model into one figure
#Remove the x-axis title for the first three plots
plot1_no_legend <- plot1 + theme(legend.position = "none", axis.title.y= element_blank())
plot2_no_legend <- plot2 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank())
plot3_no_legend <- plot3 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank())
plot4_no_legend <- plot4 + theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y= element_blank())  # Keep x-axis title for the last plot

# Extract the shared legend from one of the plots (e.g., plot1)
shared_legend <- get_legend(plot1)

# Combine the four plots into one figure, placing them in one column (4 rows)
combined_plot <- plot_grid(plot3_no_legend, plot4_no_legend,
                           plot2_no_legend, plot1_no_legend,
                           labels = c("a", "b", "c", "d"),  # Label each subplot
                           label_size =  20,  # Label size
                           ncol = 1,  # Arrange in one column (four rows)
                           align = "v")  # Align vertically
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5,size= 22)+
    theme(plot.margin= margin(r=10)),  # Add shared y-axis label
  combined_plot,
  ncol = 2,
  rel_widths = c(0.05, 0.5)  # Adjust the relative width of the label vs. the plots
)
# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(combined_with_y_axis, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust the width ratio

# Display the final combined plot
print(final_plot)
ggsave("./figure/firstfour.jpeg", plot = final_plot, width = 15, height = 20, dpi = 300)

#for mwef plot
#Remove the x-axis title for the first three plots
plot5_no_legend <- plot5 + theme(legend.position = "none", axis.title.y= element_blank())
plot6_no_legend <- plot6 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank())
plot7_no_legend <- plot7 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank())
plot8_no_legend <- plot8 + theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y= element_blank())  # Keep x-axis title for the last plot

# Extract the shared legend from one of the plots (e.g., plot1)
shared_legend <- get_legend(plot5)

# Combine the four plots into one figure, placing them in one column (4 rows)
combined_plot <- plot_grid(plot7_no_legend, plot8_no_legend,
                           plot6_no_legend, plot5_no_legend,
                           labels = c("a", "b", "c", "d"),  # Label each subplot
                           label_size = 20,  # Label size
                           ncol = 1,  # Arrange in one column (four rows)
                           align = "v")  # Align vertically
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5, size= 22)+
    theme(plot.margin= margin(r=10)),  # Add shared y-axis label
  combined_plot,
  ncol = 2,
  rel_widths = c(0.05, 0.5)  # Adjust the relative width of the label vs. the plots
)
# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(combined_with_y_axis, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust the width ratio

# Display the final combined plot
print(final_plot)
ggsave("./figure/lastfour_mwef.pdf", plot = final_plot, width = 15, height = 20, dpi = 300)


#figure 4 i.e single trait and coh2 0 for each trait complete and 5 subsets
#plsr-narea
#ST for plsr-narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMW.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn.ac) %>% select(-result_pn.ac)
acpnT_CV1<- read.csv("./output/acpnT_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("./output/acpnT.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apnd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1<- acpnT_CV1 %>% mutate(MT= "S0(coh2= 0)", Type= "CV1")
acpnT<- acpnT %>% mutate (MT= "S0(coh2= 0)", Type ="CV2")
# Combine all data

pn_data_whole <- bind_rows(acpn, acpnT_CV1,acpnT)
pn_data_whole$Data<- "Complete"


#replication-1
#ST for pnarea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep1.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep1.ac) %>% select(-result_pn_rep1.ac)
acpnT_CV1<- read.csv("acpnNT1_CV1_rep1.txt", sep= ",", header= FALSE) %>% mutate(acc= V1) %>% select(-V1)
#multitrait Cv2 for plsr-narea pnW1, pnW2, pnW3
acpnT<- read.csv("acpnNT1_CV2_rep1.txt", sep= ",", header= FALSE) %>% mutate(acc= V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep1 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep1$Data <- "Subset-1"

#replication 2
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep2.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep2.ac) %>% select(-result_pn_rep2.ac)
acpnT_CV1<- read.csv("acpnNT2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("acpnNT2_CV2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep2 <- bind_rows(acpn,acpnT_CV1, acpnT)
pn_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep3.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep3.ac) %>% select(-result_pn_rep3.ac)
acpnT_CV1<- read.csv("acpnNT3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("acpnNT3_CV2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep3 <- bind_rows(acpn,acpnT_CV1, acpnT)
pn_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep4.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep4.ac) %>% select(-result_pn_rep4.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("acpnNT4_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("acpnNT4_CV2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep4 <-bind_rows(acpn,acpnT_CV1, acpnT)
pn_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep5.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep5.ac) %>% select(-result_pn_rep5.ac)
acpnT_CV1<- read.csv("acpnNT5_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("acpnNT5_CV2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep5 <- bind_rows(acpn,acpnT_CV1, acpnT)
pn_data_rep5$Data <- "Subset-5"

# Combine all datasets
all_data <- bind_rows(pn_data_whole, pn_data_rep1, pn_data_rep2, pn_data_rep3, pn_data_rep4, pn_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))
common_theme <- theme_bw(base_size = 20) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        legend.text = element_text(size = 16),    # Adjust legend text size
        axis.text = element_text(size = 18),      # Adjust axis text size (increased from 16)
        axis.title = element_text(size = 20))
# Create a ggplot with facet_wrap
plot9 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
 common_theme +
  labs(title = "PLSR-Narea",#(Training Ef & validating on MW)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position = "right") # Bold facet labels for clarity

#ggsave("./figure/pn_compsubsetfig3_efmw.jpeg", plot = plot, width = 10, height = 8, dpi = 300)


#MWef
#plsr-narea
#ST for plsr-narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEF.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn.ac) %>% select(-result_pn.ac)
acpnT_CV1<- read.csv("./output/acpnT_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("./output/acpnT_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apnd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_mwef_whole <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_mwef_whole$Data<- "Complete"

#replication 1
#ST for pnarea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep1.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep1.ac) %>% select(-result_pn_rep1.ac)
#Multi trait cv1 for plsr-narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("acpnT1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for plsr-narea pnW1, pnW2, pnW3
acpnT<- read.csv("acpnT1_CV2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combinne all data
pn_data_rep1 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep2.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep2.ac) %>% select(-result_pn_rep2.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("acpnT2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnT<- read.csv("acpnT2_CV2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep2 <-bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep3.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep3.ac) %>% select(-result_pn_rep3.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("acpnT3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnT<- read.csv("acpnT3_CV2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep3 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep4.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep4.ac) %>% select(-result_pn_rep4.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("acpnT4_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnT<- read.csv("acpnT4_CV2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep4 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep5.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep5.ac) %>% select(-result_pn_rep5.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("acpnT5_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnT<- read.csv("acpnT5_CV2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
pn_data_rep5 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep5$Data <- "Subset-5"


# Combine all datasets
all_data <- bind_rows(pn_data_mwef_whole, pn_data_rep1, pn_data_rep2, pn_data_rep3, pn_data_rep4, pn_data_rep5)

all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))


# Create a ggplot with facet_wrap
plot13 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-Narea",#(Training MW & validating on Ef)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/pn_compsubsetfig3_mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#plsr-sla(EFMW)
#plsr-sla
#ST for PLSR-SLA
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMW.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps.ac) %>% select(-result_ps.ac)
#Multi trait cv1 for plsr_sla NW1, NW2, NW3
acpsT_CV1<- read.csv("./output/acpsT_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for psarea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait apsd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT<- acpsT %>%  mutate(MT= "S0(coh2=0)", Type="CV2")
# Combine all data
ps_data_whole <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_whole$Data<- "Complete"


#replication 1
#ST for pnarea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep1.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep1.ac) %>% select(-result_ps_rep1.ac)
#Multi trait cv1 for plsr-sla psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for plsr-sla psW1, psW2, psW3
acpsT<- read.csv("acpsT1_CV2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


# Add labels for multi-trait aSd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combinne all data
ps_data_rep1 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep2.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep2.ac) %>% select(-result_ps_rep2.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("acpsT2_CV2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_rep2 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep3.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep3.ac) %>% select(-result_ps_rep3.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("acpsT3_CV2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_rep3 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep4.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep4.ac) %>% select(-result_ps_rep4.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT4_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("acpsT4_CV2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_rep4 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep5.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep5.ac) %>% select(-result_ps_rep5.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT5_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("acpsT5_CV2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_rep5 <- bind_rows(acps, acpsT, acpsT_CV1)
ps_data_rep5$Data <- "Subset-5"
# Combine all datasets
all_data <- bind_rows(ps_data_whole, ps_data_rep1, ps_data_rep2, ps_data_rep3, ps_data_rep4, ps_data_rep5)
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))
# Create a ggplot with facet_wrap
plot10 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "PLSR-SLA",#(Training EF & validating on MW)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/ps_compsubsetfig3_efmw.jpeg", plot = plot, width = 10, height = 8, dpi = 300)

#plsr-sla(MWEF)
#plsr sla
#ST for PLSR-SLA
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEF.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps.ac) %>% select(-result_ps.ac)
#Multi trait cv1 for plsr_sla NW1, NW2, NW3
acpsT_CV1<- read.csv("./output/acpsT_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for psarea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apsd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_mwef_whole <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_mwef_whole$Data<- "Complete"

#replication 1
#ST for pnarea
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEFrep1.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep1.ac) %>% select(-result_ps_rep1.ac)
#Multi trait cv1 for plsr-sla psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for plsr-sla psW1, psW2, psW3
acpsT<- read.csv("acpsT1_CV2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combinne all data
ps_data_rep1 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEFrep2.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep2.ac) %>% select(-result_ps_rep2.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("acpsT2_CV2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_rep2 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEFrep3.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep3.ac) %>% select(-result_ps_rep3.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("acpsT3_CV2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_rep3 <- bind_rows(acps, acpsT_CV1,acpsT)
ps_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEFrep4.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep4.ac) %>% select(-result_ps_rep4.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT4_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("acpsT4_CV2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")

acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_rep4 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEFrep5.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep5.ac) %>% select(-result_ps_rep5.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("acpsT5_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("acpsT5_CV2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
ps_data_rep5 <- bind_rows(acps, acpsT_CV1,acpsT)
ps_data_rep5$Data <- "Subset-5"


# Combine all datasets
all_data <- bind_rows(ps_data_mwef_whole, ps_data_rep1, ps_data_rep2, ps_data_rep3, ps_data_rep4, ps_data_rep5)
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))


# Create a ggplot with facet_wrap
plot14 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "PLSR-SLA",#(Training MW & validating on EF)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/ps_compsubsetfig3_mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)


#EFMW- narea
#ST for narea
corr_N_EFMW<- read.csv("./output/narea/corr_N_EFMW.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N.ac) %>% select(-result_N.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2= 0)", Type = "CV1")
acNT<- acNT %>%  mutate(MT= "S0(coh2= 0)",Type= "CV2")
# Combine all data
N_data_whole <- bind_rows(acN,acNT_CV1,acNT)#acNT_CV1,acNT
N_data_whole$Data <- "Complete"

#replication 1
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep1.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep1.ac) %>% select(-result_N_rep1.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT1_CV2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep1 <- bind_rows(acN,acNT_CV1,acNT)
N_data_rep1$Data <- "Subset-1"

#Set 2
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep2.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep2.ac) %>% select(-result_N_rep2.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT2_CV2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep2 <- bind_rows(acN,acNT_CV1,acNT)
N_data_rep2$Data <- "Subset-2"
#Set 3
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep3.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep3.ac) %>% select(-result_N_rep3.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT3_CV2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep3 <- bind_rows(acN,acNT_CV1,acNT)
N_data_rep3$Data <- "Subset-3"

#Set 4
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep4.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep4.ac) %>% select(-result_N_rep4.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT4_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT4_CV2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep4 <- bind_rows(acN,acNT_CV1,acNT)
N_data_rep4$Data <- "Subset-4"
#Set 5
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep5.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep5.ac) %>% select(-result_N_rep5.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT5_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT5_CV2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep5 <-bind_rows(acN,acNT_CV1,acNT)
N_data_rep5$Data <- "Subset-5"

# Combine all datasets
all_data <- bind_rows(N_data_whole, N_data_rep1, N_data_rep2, N_data_rep3, N_data_rep4, N_data_rep5)
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))


# Create a ggplot with facet_wrap
plot11 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Narea",#(Training EF & validating on MW)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/Narea_compsubsetfig3_efmw.jpeg", plot = plot, width = 10, height = 8, dpi = 300)




#MWEF narea
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEF.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N.ac) %>% select(-result_N.ac)
#Multi trait CV1_mwef for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2= 0)", Type = "CV1")

acNT <- acNT %>% mutate(MT = "S0(coh2= 0)", Type = "CV2")

# Combine all data
N_data_mwef_whole <- bind_rows(acN,acNT_CV1,acNT)
N_data_mwef_whole$Data <- "Complete"
#replication 1 mwef
#ST for narea
corr_N_MWEF<- read.csv("./output/corr_N_MWEFrep1.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep1.ac) %>% select(-result_N_rep1.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT1_CV2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep1 <- bind_rows(acN, acNT, acNT_CV1)
N_data_rep1$Data <- "Subset-1"
#Subset 2
#ST for narea
corr_N_MWEF<- read.csv("./output/corr_N_MWEFrep2.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep2.ac) %>% select(-result_N_rep2.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT2_CV2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep2 <- bind_rows(acN,acNT, acNT_CV1)
N_data_rep2$Data <- "Subset-2"
#Subset 3
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep3.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep3.ac) %>% select(-result_N_rep3.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT3_CV2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep3 <- bind_rows(acN,acNT,acNT_CV1)
N_data_rep3$Data <- "Subset-3"

#Subset 4
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep4.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep4.ac) %>% select(-result_N_rep4.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT4_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT4_CV2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")

acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep4 <- bind_rows(acN,acNT, acNT_CV1)
N_data_rep4$Data <- "Subset-4"
#Subset 5
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep5.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep5.ac) %>% select(-result_N_rep5.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("acNT5_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("acNT5_CV2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
N_data_rep5 <- bind_rows(acN,acNT,acNT_CV1)
N_data_rep5$Data <- "Subset-5"

# Combine all dataSubsets
all_data <- bind_rows(N_data_mwef_whole, N_data_rep1, N_data_rep2, N_data_rep3, N_data_rep4, N_data_rep5)
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))



# Create a ggplot with facet_wrap
plot15 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Narea",#(Training MW & validating on EF)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position = "right") # Bold facet labels for clarity

#ggsave("./figure/Narea_compsubsetfig3_mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#sla-efmw
#ST for SLA EF
corr_S_EFMW<- read.csv("./output/corr_S_EFMW.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S.ac) %>% select(-result_S.ac)
#Multi trait cv1 for sla SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for sla SW1, SW2, SW3
acST<- read.csv("./output/acST.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2= 0)", Type = "CV1")
acST<- acST %>%  mutate(MT= "S0(coh2= 0)",Type= "CV2")

# Combine all data
S_data_whole <- bind_rows(acS,acST_CV1,acST)
S_data_whole$Data <- "Complete"


#Subset 1
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep1.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep1.ac) %>% select(-result_S_rep1.ac)
acST_CV1<- read.csv("acST1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for sla SW1, SW2, SW3
acST<- read.csv("acST1_CV2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# CombiSe all data
S_data_rep1 <- bind_rows(acS,acST_CV1,acST)
S_data_rep1$Data <- "Subset-1"

#Subset 2
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep2.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep2.ac) %>% select(-result_S_rep2.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("acST2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST2_CV2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_rep2 <- bind_rows(acS,acST_CV1,acST)
S_data_rep2$Data <- "Subset-2"
#Subset 3
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep3.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep3.ac) %>% select(-result_S_rep3.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("acST3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST3_CV2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_rep3 <- bind_rows(acS,acST, acST_CV1)
S_data_rep3$Data <- "Subset-3"

#Subset 4
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep4.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep4.ac) %>% select(-result_S_rep4.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("acST4_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST4_CV2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_rep4 <- bind_rows(acS,acST_CV1,acST)
S_data_rep4$Data <- "Subset-4"

#Subset 5
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep5.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep5.ac) %>% select(-result_S_rep5.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("acST5_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST5_CV2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_rep5 <- bind_rows(acS, acST, acST_CV1,)
S_data_rep5$Data <- "Subset-5"

# Combine all dataSubsets
all_data <- bind_rows(S_data_whole, S_data_rep1, S_data_rep2, S_data_rep3, S_data_rep4, S_data_rep5)
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Create a ggplot with facet_wrap
plot12 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "SLA",#(Training EF & validating on MW)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right")
#ggsave("./figure/sla_compsubsetfig3_efmw.jpeg", plot = plot, width = 10, height = 8, dpi = 300)



#mwef-sla
#ST for SLA MWEF
corr_S_MWEF<- read.csv("./output/corr_S_MWEF.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S.ac) %>% select(-result_S.ac)
#Multi trait cv1 for sla SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for sla SW1, SW2, SW3
acST<- read.csv("./output/acST_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")

acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_mwef_whole <- bind_rows(acS,acST,acST_CV1)
S_data_mwef_whole$Data <- "Complete"

#rep-1
#ST for sla
corr_S_MWEF<- read.csv("./output/corr_S_MWEFrep1.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep1.ac) %>% select(-result_S_rep1.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("acST1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST1_CV2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")

acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_rep1 <- bind_rows(acS,acST_CV1,acST)
S_data_rep1$Data <- "Subset-1"

#Subset 2
#ST for narea
corr_S_MWEF<- read.csv("./output/corr_S_MWEFrep2.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep2.ac) %>% select(-result_S_rep2.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("acST2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST2_CV2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")

acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_rep2 <- bind_rows(acS,acST_CV1,acST)
S_data_rep2$Data <- "Subset-2"


#Subset 3
#ST for sla
corr_S_MWEF<- read.csv("./output/sla/corr_S_MWEFrep3.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep3.ac) %>% select(-result_S_rep3.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("acST3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST3_CV2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")

acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")


# Combine all data
S_data_rep3 <-bind_rows(acS,acST_CV1,acST)
S_data_rep3$Data <- "Subset-3"
#Subset 4
#ST for sla
corr_S_MWEF<- read.csv("./output/sla/corr_S_MWEFrep4.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep4.ac) %>% select(-result_S_rep4.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("acST4_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST4_CV2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")
acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_rep4 <- bind_rows(acS,acST_CV1,acST)
S_data_rep4$Data <- "Subset-4"
#Subset 5
#ST for narea
corr_S_MWEF<- read.csv("./output/sla/corr_S_MWEFrep5.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep5.ac) %>% select(-result_S_rep5.ac)
acST_CV1<- read.csv("acST5_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST5_CV2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S0(coh2=0)", Type = "CV1")

acST <- acST %>% mutate(MT = "S0(coh2=0)", Type = "CV2")

# Combine all data
S_data_rep5 <- bind_rows(acS,acST,acST_CV1)
S_data_rep5$Data <- "Subset-5"

# Combine all dataSubset
all_data <- bind_rows(S_data_mwef_whole, S_data_rep1, S_data_rep2, S_data_rep3, S_data_rep4, S_data_rep5)
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

# Create a ggplot with facet_wrap
plot16 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "SLA",#(Training MW & validating on EF)",
       x = NULL,
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/sla_compsubsetfig3_mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)

#for mwef four figure in one plot
shared_legend <- get_legend(plot13)

# Remove legends from individual plots
plot13_no_legend <- plot9 + theme(legend.position = "none")
plot14_no_legend <- plot10 + theme(legend.position = "none")
plot15_no_legend <- plot11 + theme(legend.position = "none")
plot1_no_legend <- plot12 + theme(legend.position = "none")

combined_plot <- plot_grid(plot9_no_legend, plot10_no_legend,
                           plot11_no_legend, plot12_no_legend,
                           labels = c("a", "b", "c", "d"),  # Label each subplot
                           label_size = 14,  # Label size
                           ncol = 2,  # Arrange in a 2x2 grid
                           nrow = 2)
final_plot <- plot_grid(combined_plot, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust width ratio

print(final_plot)

ggsave("./figure/zerocoh2vsST_efmw_completevssubsets.jpeg", plot = final_plot, width = 20, height = 10, dpi = 300)

#Remove the x-axis title for the first three plots
plot9_no_legend <- plot9 + theme(legend.position = "none", axis.title.y= element_blank())
plot10_no_legend <- plot10 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())
plot11_no_legend <- plot11 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())
plot12_no_legend <- plot12 + theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())  # Keep x-axis title for the last plot

# Extract the shared legend from one of the plots (e.g., plot1)
shared_legend <- get_legend(plot9)

# Combine the four plots into one figure, placing them in one column (4 rows)
combined_plot <- plot_grid(plot11_no_legend, plot12_no_legend,
                           plot10_no_legend, plot9_no_legend,
                           labels = c("a", "b", "c", "d"),  # Label each subplot
                           label_size = 20,  # Label size
                           ncol = 1,  # Arrange in one column (four rows)
                           align = "v")  # Align vertically
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5),  # Add shared y-axis label
  combined_plot,
  ncol = 2,
  rel_widths = c(0.05, 1)  # Adjust the relative width of the label vs. the plots
)
# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(combined_with_y_axis, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust the width ratio

# Display the final combined plot
print(final_plot)
ggsave("./figure/zerocoh2vsST_efmw_completevssubsets.pdf", plot = final_plot, width = 10, height = 15, dpi = 300)

#Remove the x-axis title for the first three plots
plot13_no_legend <- plot13 + theme(legend.position = "none", axis.title.y= element_blank())
plot14_no_legend <- plot14 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())
plot15_no_legend <- plot15 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())
plot16_no_legend <- plot16 + theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())  # Keep x-axis title for the last plot

# Extract the shared legend from one of the plots (e.g., plot1)
shared_legend <- get_legend(plot13)

# Combine the four plots into one figure, placing them in one column (4 rows)
combined_plot <- plot_grid(plot15_no_legend, plot16_no_legend,
                           plot14_no_legend, plot13_no_legend,
                           labels = c("a", "b", "c", "d"),  # Label each subplot
                           label_size = 18,  # Label size
                           ncol = 1,  # Arrange in one column (four rows)
                           align = "v")  # Align vertically
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5),  # Add shared y-axis label
  combined_plot,
  ncol = 2,
  rel_widths = c(0.05, 1)  # Adjust the relative width of the label vs. the plots
)
# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(combined_with_y_axis, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust the width ratio

# Display the final combined plot
print(final_plot)
ggsave("./figure/zerocoh2vsST_mwef_completevssubsets.pdf", plot = final_plot, width = 10, height = 15, dpi = 300)



#for fig 2 i.e all complete model only
#plsr-narea
#ST for plsr-narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMW.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn.ac) %>% select(-result_pn.ac)
#Multi trait cv1 for plsr_narea NW1, NW2, NW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for pnarea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apnd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.54)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.53)", Type = "CV1")
acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.54)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.53)", Type = "CV2")
# Combine all data
pn_data_whole <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1, acpnW1, acpnW2, acpnW3)
pn_data_whole$Data<- "PLSR-Narea"

#MWef
#plsr-narea
#ST for plsr-narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEF.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn.ac) %>% select(-result_pn.ac)
#Multi trait CV1_mwef for plsr_narea NW1, NW2, NW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for pnarea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apnd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnW1_CV1 <- acpnW1_CV1 %>% mutate(MT = "S1(coh2=0.54)", Type = "CV1")
acpnW2_CV1 <- acpnW2_CV1 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV1")
acpnW3_CV1 <- acpnW3_CV1 %>% mutate(MT = "S3(coh2=0.53)", Type = "CV1")

acpnW1 <- acpnW1 %>% mutate(MT = "S1(coh2=0.54)", Type = "CV2")
acpnW2 <- acpnW2 %>% mutate(MT = "S2(coh2=0.53)", Type = "CV2")
acpnW3 <- acpnW3 %>% mutate(MT = "S3(coh2=0.53)", Type = "CV2")

# Combine all data
pn_data_mwef <- bind_rows(acpn, acpnW1_CV1, acpnW2_CV1, acpnW3_CV1, acpnW1, acpnW2, acpnW3)
pn_data_mwef$Data<- "PLSR-Narea"

#plsr-sla(EFMW)
#plsr-sla
#ST for PLSR-SLA
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMW.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps.ac) %>% select(-result_ps.ac)
#Multi trait cv1 for plsr_sla NW1, NW2, NW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for psarea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait apsd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.46)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.46)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.46)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.46)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.46)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.46)", Type = "CV2")
# Combine all data
ps_data_whole <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1, acpsW1, acpsW2, acpsW3)
ps_data_whole$Data<- "PLSR-SLA"

#ST for PLSR-SLA
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEF.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps.ac) %>% select(-result_ps.ac)
#Multi trait cv1 for plsr_sla NW1, NW2, NW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for psarea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apsd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsW1_CV1 <- acpsW1_CV1 %>% mutate(MT = "S1(coh2=0.46)", Type = "CV1")
acpsW2_CV1 <- acpsW2_CV1 %>% mutate(MT = "S2(coh2=0.46)", Type = "CV1")
acpsW3_CV1 <- acpsW3_CV1 %>% mutate(MT = "S3(coh2=0.46)", Type = "CV1")

acpsW1 <- acpsW1 %>% mutate(MT = "S1(coh2=0.46)", Type = "CV2")
acpsW2 <- acpsW2 %>% mutate(MT = "S2(coh2=0.46)", Type = "CV2")
acpsW3 <- acpsW3 %>% mutate(MT = "S3(coh2=0.46)", Type = "CV2")

# Combine all data
ps_data_mwef <- bind_rows(acps, acpsW1_CV1, acpsW2_CV1, acpsW3_CV1,acpsW1, acpsW2, acpsW3)
ps_data_mwef$Data<- "PLSR-SLA"

#EFMW- narea
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMW.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N.ac) %>% select(-result_N.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")
# Combine all data
N_data_whole <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1,acNW1,acNW2, acNW3)#acNT_CV1,acNT
N_data_whole$Data <- "Narea"


#MWEF narea
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEF.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N.ac) %>% select(-result_N.ac)
#Multi trait CV1_mwef for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("acNW1_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("acNW2_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("acNW3_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("acNW1_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("acNW2_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("acNW3_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2= 0.49)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2= 0.49)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2= 0.49)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2= 0.49)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV2")

# Combine all data
N_data_mwef <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1, acNW1, acNW2, acNW3)
N_data_mwef$Data <- "Narea"


#sla-efmw
#ST for SLA EF
corr_S_EFMW<- read.csv("./output/corr_S_EFMW.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S.ac) %>% select(-result_S.ac)
#Multi trait cv1 for sla SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for sla SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2= 0.57)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2= 0.54)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2= 0.57)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2= 0.54)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV2")

# Combine all data
S_data_whole <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_whole$Data <- "SLA"
#mwef-sla
#ST for SLA MWEF
corr_S_MWEF<- read.csv("./output/corr_S_MWEF.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S.ac) %>% select(-result_S.ac)
#Multi trait cv1 for sla SW1, SW2, SW3
acSW1_CV1<- read.csv("./output/acSW1_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2_CV1<- read.csv("./output/acSW2_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3_CV1<- read.csv("./output/acSW3_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for sla SW1, SW2, SW3
acSW1<- read.csv("./output/acSW1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW2<- read.csv("./output/acSW2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acSW3<- read.csv("./output/acSW3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2= 0.57)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2= 0.54)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2= 0.57)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2= 0.54)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2= 0.49)", Type = "CV2")

# Combine all data
S_data_mwef <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1, acSW1, acSW2, acSW3)
S_data_mwef$Data <- "SLA"

# Combine all dataset
all_data <- bind_rows(pn_data_whole, ps_data_whole,N_data_whole,S_data_whole)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Narea","SLA","PLSR-SLA","PLSR-Narea"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))
ps_data_whole$Type <- factor(ps_data_whole$Type, levels = c("ST", "CV1", "CV2"))

# Create individual plots for each trait (Narea, SLA, PLSR-SLA, PLSR-Narea)

# Create individual plots for each trait (Narea, SLA, PLSR-SLA, PLSR-Narea)
plot_pn <- ggplot(subset(pn_data_whole, Data == "PLSR-Narea"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-Narea", x = "Model Type", y = "Accuracy") +
  scale_y_continuous(breaks = seq(0, max(pn_data_whole$acc), by = 0.02))+# Remove individual axis titles
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")
ggsave("./figure/plsrnarea_seminar.jpeg", plot = plot_pn, width = 12, height = 7 , dpi = 300)

plot_ps <- ggplot(subset(ps_data_whole, Data == "PLSR-SLA"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-SLA", x = "Model Type", y = "Accuracy") +
  scale_y_continuous(breaks = seq(0, max(ps_data_whole$acc), by = 0.02))+# Remove individual axis titles
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")
ggsave("./figure/plsrsla_seminar.jpeg", plot = plot_ps, width = 12, height = 7 , dpi = 300)


plot1 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-Narea",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, max(all_data$acc), by = 0.03))+

  # Adjust theme for clarity
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

plot_sla <- ggplot(subset(all_data, Data == "SLA"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "SLA", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

plot_plsr_sla <- ggplot(subset(all_data, Data == "PLSR-SLA"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "PLSR-SLA", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

plot_plsr_narea <- ggplot(subset(all_data, Data == "PLSR-Narea"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "PLSR-Narea", x = "Model Type", y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Extract the shared legend from one of the plots
shared_legend <- get_legend(plot_narea + theme(legend.position = "right"))

# Combine the individual plots into a single figure without x/y axis titles or legend
combined_plots <- plot_grid(
  plot_narea, plot_sla, plot_plsr_sla, plot_plsr_narea,
  labels = c("a", "b", "c", "d"),  # Add labels
  label_size = 12,  # Label size
  ncol = 1,  # Arrange in one column (4 rows)
  align = "v"
)

# Add the shared y-axis label ("Accuracy")
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5),  # Add shared y-axis label
  combined_plots,
  ncol = 2,  # Arrange label and plots in two columns
  rel_widths = c(0.04, 1)  # Adjust the relative width of the label vs. the plots
)

# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(
  combined_with_y_axis, shared_legend,
  ncol = 2,  # Place the legend beside the plot grid
  rel_widths = c(4, 0.4)  # Adjust the width ratio
)

# Add the shared x-axis label ("Model Type") at the bottom
final_with_x_axis <- plot_grid(
  final_plot,
  ggdraw() + draw_label("Model Type", vjust = 1),  # Add shared x-axis label
  ncol = 1,  # Arrange vertically
  rel_heights = c(1, 0.05)  # Adjust relative height of the plot and the x-axis label
)

print(final_with_x_axis)

ggsave("./figure/completemodel_fig2efmw.pdf", plot = final_plot, width = 7, height = 10 , dpi = 300)


# Combine all data set
all_data <- bind_rows(pn_data_mwef, ps_data_mwef,N_data_mwef,S_data_mwef)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Narea","SLA","PLSR-SLA","PLSR-Narea"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))
# Create individual plots for each trait (Narea, SLA, PLSR-SLA, PLSR-Narea)
plot_narea <- ggplot(subset(all_data, Data == "Narea"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Narea", x = NULL, y = NULL) +  # Remove individual axis titles
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

plot_sla <- ggplot(subset(all_data, Data == "SLA"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "SLA", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

plot_plsr_sla <- ggplot(subset(all_data, Data == "PLSR-SLA"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "PLSR-SLA", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

plot_plsr_narea <- ggplot(subset(all_data, Data == "PLSR-Narea"), aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "PLSR-Narea", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Extract the shared legend from one of the plots
shared_legend <- get_legend(plot_narea + theme(legend.position = "right"))

# Combine the individual plots into a single figure without x/y axis titles or legend
combined_plots <- plot_grid(
  plot_narea, plot_sla, plot_plsr_sla, plot_plsr_narea,
  labels = c("a", "b", "c", "d"),  # Add labels
  label_size = 12,  # Label size
  ncol = 1,  # Arrange in one column (4 rows)
  align = "v"
)

# Add the shared y-axis label ("Accuracy")
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5),  # Add shared y-axis label
  combined_plots,
  ncol = 2,  # Arrange label and plots in two columns
  rel_widths = c(0.04, 1)  # Adjust the relative width of the label vs. the plots
)

# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(
  combined_with_y_axis, shared_legend,
  ncol = 2,  # Place the legend beside the plot grid
  rel_widths = c(4, 0.4)  # Adjust the width ratio
)

# Add the shared x-axis label ("Model Type") at the bottom
final_with_x_axis <- plot_grid(
  final_plot,
  ggdraw() + draw_label("Model Type", vjust = 1),  # Add shared x-axis label
  ncol = 1,  # Arrange vertically
  rel_heights = c(1, 0.05)  # Adjust relative height of the plot and the x-axis label
)

print(final_with_x_axis)

ggsave("./figure/completemodel_fig2mwef.pdf", plot = final_with_x_axis, width = 7, height = 10 , dpi = 300)


# Create a ggplot with facet_wrap
plotB <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Complete model(Training MW & validating on EF)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

ggsave("./figure/completemodel_fig2mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)

# Extract the shared legend from one of the plots (e.g., plot1)
shared_legend <- get_legend(plotA)

# Remove legends from individual plots
plotA_no_legend <- plotA+ theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())
plotB_no_legend <- plotB + theme(legend.position = "none")


combined_plot <- plot_grid(plotA_no_legend, plotB_no_legend,
                           labels = c("a", "b"),  # Label each subplot
                           label_size = 18,  # Label size
                           ncol = 1,  # Arrange in a 2x2 grid
                           align="v")
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5),  # Add shared y-axis label
  combined_plot,
  ncol = 2,
  rel_widths = c(0.05, 1)  # Adjust the relative width of the label vs. the plots
)
# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(combined_with_y_axis, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust the width ratio

# Display the final combined plot
print(final_plot)
ggsave("./figure/firstfour.jpeg", plot = final_plot, width = 15, height = 20, dpi = 300)
combined_with_y_axis<- plot
final_plot <- plot_grid(combined_plot, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust width ratio

print(final_plot)

ggsave("./figure/allcompletemodel_ef&mw.jpeg", plot = final_plot, width = 20, height = 10, dpi = 300)



#making  four traits efmw model into one figure
#Remove the x-axis title for the first three plots
plot1_no_legend <- plot1 + theme(legend.position = "none", axis.title.y= element_blank())
plot2_no_legend <- plot2 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank())
plot3_no_legend <- plot3 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank())
plot4_no_legend <- plot4 + theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y= element_blank())  # Keep x-axis title for the last plot

# Extract the shared legend from one of the plots (e.g., plot1)
shared_legend <- get_legend(plot1)

# Combine the four plots into one figure, placing them in one column (4 rows)
combined_plot <- plot_grid(plot3_no_legend, plot4_no_legend,
                           plot2_no_legend, plot1_no_legend,
                           labels = c("a", "b", "c", "d"),  # Label each subplot
                           label_size = 18,  # Label size
                           ncol = 1,  # Arrange in one column (four rows)
                           align = "v")  # Align vertically
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.35, hjust = 0.5,size= 22)+
    theme(plot.margin= margin(r=20)),  # Add shared y-axis label
  combined_plot,
  ncol = 2,
  rel_widths = c(0.05, 1.5)  # Adjust the relative width of the label vs. the plots
)
# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(combined_with_y_axis, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust the width ratio

# Display the final combined plot
print(final_plot)
ggsave("./figure/firstfour.pdf", plot = final_plot, width = 15, height = 20, dpi = 300)
#for mwef plot
#Remove the x-axis title for the first three plots
plot5_no_legend <- plot5 + theme(legend.position = "none", axis.title.y= element_blank())
plot6_no_legend <- plot6 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank())
plot7_no_legend <- plot7 + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank())
plot8_no_legend <- plot8 + theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y= element_blank())  # Keep x-axis title for the last plot

# Extract the shared legend from one of the plots (e.g., plot1)
shared_legend <- get_legend(plot5)

# Combine the four plots into one figure, placing them in one column (4 rows)
combined_plot <- plot_grid(plot7_no_legend, plot8_no_legend,
                           plot6_no_legend, plot5_no_legend,
                           labels = c("a", "b", "c", "d"),  # Label each subplot
                           label_size = 18,  # Label size
                           ncol = 1,  # Arrange in one column (four rows)
                           align = "v")  # Align vertically
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5),  # Add shared y-axis label
  combined_plot,
  ncol = 2,
  rel_widths = c(0.05, 0.5)  # Adjust the relative width of the label vs. the plots
)
# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(combined_with_y_axis, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust the width ratio

# Display the final combined plot
print(final_plot)
ggsave("./figure/lastfour_mwef.pdf", plot = final_plot, width = 15, height = 20, dpi = 300)
# Narea Plot
narea_data <- data.frame(
  Target_Trait = rep("Narea", 9),
  Synthetic_Trait = rep(c("s1", "s2", "s3"), each = 3),
  Metric = rep(c("h2", "coh2", "corg"), 3),
  Value = c(0.61, 0.49, 0.82, 0.61, 0.49, 0.76, 0.61, 0.49, 0.83)
)


# SLA Plot
SLA_data <- data.frame(
  Target_Trait = rep("SLA", 9),
  Synthetic_Trait = rep(c("s1", "s2", "s3"), each = 3),
  Metric = rep(c("h2", "coh2", "corg"), 3),
  Value = c(0.55, 0.47, 0.79, 0.58, 0.50, 0.81, 0.59, 0.48, 0.80)
)
# PLSR-Narea Plot
PLSRNarea_data <- data.frame(
  Target_Trait = rep("PLSR-Narea", 9),
  Synthetic_Trait = rep(c("s1", "s2", "s3"), each = 3),
  Metric = rep(c("h2", "coh2", "corg"), 3),
  Value = c(0.68, 0.54, 0.82, 0.68, 0.53, 0.88, 0.68, 0.53, 0.79)
)

# PLSR-SLA Plot
PLSRSLA_data <- data.frame(
  Target_Trait = rep("PLSR-SLA", 9),
  Synthetic_Trait = rep(c("s1", "s2", "s3"), each = 3),
  Metric = rep(c("h2", "coh2", "corg"), 3),
  Value = c(0.63, 0.46, 0.79, 0.63, 0.46, 0.78, 0.63, 0.46, 0.72))

# Ensure Metric is a factor with the correct order
narea_data$Metric <- factor(narea_data$Metric, levels = c("h2", "corg", "coh2"))
SLA_data$Metric <- factor(SLA_data$Metric, levels = c("h2", "corg", "coh2"))
PLSRNarea_data$Metric <- factor(PLSRNarea_data$Metric, levels = c("h2", "corg", "coh2"))
PLSRSLA_data$Metric <- factor(PLSRSLA_data$Metric, levels = c("h2", "corg", "coh2"))

  # Define common theme with consistent border size
common_theme <- theme_bw(base_size = 26) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        legend.text = element_text(size = 36),    # Adjust legend text size
        #axis.text.x = element_text(size = 26),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 30),      # Adjust axis text size
        axis.text.y = element_text(size=30))     # Adjust axis title size

# Create individual plots with titles and no legends
narea_plot <- ggplot(narea_data, aes(x = Synthetic_Trait, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.4, color = "black") +
  geom_hline(yintercept = 0.34, linetype = "dotted", color = "black", size = 1) +
  labs(title = "Narea", x= "Synthetic Traits", y= "Value") +
  scale_fill_manual(values = c("h2" = "#0072B2", "coh2" = "#E69F00", "corg" = "#009E73")) +
  ylim(0, 1) +
  common_theme + theme(legend.position = "right")

SLA_plot <- ggplot(SLA_data, aes(x = Synthetic_Trait, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.4, color = "black") +
  geom_hline(yintercept = 0.32, linetype = "dotted", color = "black", size = 1) +
  labs(title = "SLA", x= "Synthetic Traits", y= "Value") +
  ylim(0, 1) +
  scale_fill_manual(values = c("h2" = "#0072B2", "coh2" = "#E69F00", "corg" = "#009E73")) +
  common_theme + theme(legend.position = "right")

PLSRNarea_plot <- ggplot(PLSRNarea_data, aes(x = Synthetic_Trait, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.4, color = "black") +
  geom_hline(yintercept = 0.4, linetype = "dotted", color = "black", size = 1) +
  labs(title = "PLSR-Narea", x= "Synthetic Traits", y= "Value") +
  ylim(0, 1) +
  scale_fill_manual(values = c("h2" = "#0072B2", "coh2" = "#E69F00", "corg" = "#009E73")) +
  common_theme + theme(legend.position = "right")

PLSRSLA_plot <- ggplot(PLSRSLA_data, aes(x = Synthetic_Trait, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.4, color = "black") +
  geom_hline(yintercept = 0.21, linetype = "dotted", color = "black", size = 1) +
  labs(title = "PLSR-SLA", x = "Synthetic Trait", y = "Value") +
  scale_fill_manual(values = c("h2" = "#0072B2", "coh2" = "#E69F00", "corg" = "#009E73")) +
  ylim(0, 1) +
  common_theme + theme(legend.position = "right")


#Remove the x-axis title for the first three plots
plotnarea_no_legend <- narea_plot + theme(legend.position = "none", axis.title.y= element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank())
plotsla_no_legend <- SLA_plot + theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())
plotplsrsla_no_legend <- PLSRSLA_plot+ theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y= element_blank(),axis.text.x = element_blank())
plotplsrnarea_no_legend <- PLSRNarea_plot + theme(legend.position = "none",axis.title.y = element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank())  # Keep x-axis title for the last plot

# Extract the shared legend from one of the plots (e.g., plot1)
shared_legend <- get_legend(narea_plot)
# Combine the four plots into one figure, placing them in one column (4 rows)
combined_plot <- plot_grid(plotnarea_no_legend, plotsla_no_legend,
                           plotplsrsla_no_legend, plotplsrnarea_no_legend,
                           labels = c("a", "b", "c", "d"),  # Label each subplot
                           label_size =  30,  # Label size
                           ncol = 1,  # Arrange in one column (four rows)
                           align = "v")  # Align vertically
combined_with_y_axis <- plot_grid(
  ggdraw() + draw_label("Value", angle = 90, vjust = 0.5, hjust = 0.5,size= 36)+
    theme(plot.margin= margin(r=10)),  # Add shared y-axis label
  combined_plot,
  ncol = 2,
  rel_widths = c(0.05, 0.8)  # Adjust the relative width of the label vs. the plots
)
# Combine the plot grid with the shared legend on the right
final_plot <- plot_grid(combined_with_y_axis, shared_legend,
                        ncol = 2,  # Place the legend beside the plot grid
                        rel_widths = c(4, 0.4))  # Adjust the width ratio

# Display the final combined plot
print(final_plot)
ggsave("./figure/fig.jpeg", plot = final_plot, width = 18, height = 10, dpi = 300)



