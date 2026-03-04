library(tidyverse)
library(cowplot)
library(scales)


#h2 comparison
# =========================
# 1) PUT YOUR NEW SYNTHETIC TRAIT NAMES HERE
# (these appear on x-axis)
# =========================
synthetic_labels_narea <- c(
  "wave_1511_wave_2314\n(S1)",
  "wave_1039_wave_405\n(S2)",
  "wave_910_wave_731\n(S3)"
)

synthetic_labels_sla <- c(
  "wave_1683_wave_1666\n(S1)",
  "wave_1640_wave_1655\n(S2)",
  "wave_738_wave_1111\n(S3)"
)

synthetic_labels_plsr_sla <- c(
  "wave_2242_wave_1540\n(S1)",
  "wave_2156_wave_1385\n(S2)",
  "wave_393_wave_779\n(S3)"
)

synthetic_labels_plsr_narea <- c(
  "wave_1253_wave_376\n(S1)",
  "wave_867_wave_732\n(S2)",
  "wave_1339_wave_2266\n(S3)"
)

# order must be: S1(h2,corg,coh2), S2(h2,corg,coh2), S3(h2,corg,coh2)
# =========================
narea_values <- c(
  0.61, 0.89, 0.49,
  0.61, 0.85, 0.49,
  0.61, 0.76, 0.49
)

sla_values <- c(
  0.62, 0.98, 0.52,
  0.57, 0.84, 0.50,
  0.62, 0.76, 0.49
)

plsr_sla_values <- c(
  0.63, 0.76, 0.47,
  0.63, 0.80, 0.47,
  0.63, 0.75, 0.47
)

plsr_narea_values <- c(
  0.68, 0.88, 0.54,
  0.68, 0.79, 0.54,
  0.68, 0.82, 0.54
)

# dotted reference lines (match your sample)
ref_narea <- 0.34
ref_sla <- 0.33
ref_plsr_sla <- 0.28
ref_plsr_narea <- 0.38

# =========================
# 3) BUILD LONG DATA (no changes needed)
# =========================
make_trait_df <- function(trait, syn_labels, values) {
  expand_grid(
    Target_Trait = trait,
    Synthetic_Trait = syn_labels,
    Metric = c("h2", "corg", "coh2")
  ) %>%
    mutate(Value = values,
           Metric = factor(Metric, levels = c("h2", "corg", "coh2")))
}

narea_data <- make_trait_df("Narea", synthetic_labels_narea, narea_values)
sla_data <- make_trait_df("SLA", synthetic_labels_sla, sla_values)
plsr_sla_data <- make_trait_df("PLSR-SLA", synthetic_labels_plsr_sla, plsr_sla_values)
plsr_narea_data <- make_trait_df("PLSR-Narea", synthetic_labels_plsr_narea, plsr_narea_values)

# =========================
# 4) STYLE (match sample)
# =========================
fill_cols <- c("h2" = "#0072B2", "corg" = "#009E73", "coh2" = "#E69F00")

common_theme <- theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 11, hjust = 0),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# =========================
# 5) PLOTS (one legend)
# =========================
plot_bar <- function(df, title, ref_line, show_x = FALSE, show_legend = FALSE) {
  ggplot(df, aes(x = Synthetic_Trait, y = Value, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.6),
             width = 0.55, color = "black", linewidth = 0.3) +
    geom_hline(yintercept = ref_line, linetype = "dotted", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = fill_cols, name = "Metric") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = title,
      x = if (show_x) "Synthetic Traits" else NULL,
      y = NULL
    ) +
    common_theme +
    theme(
      legend.position = if (show_legend) "right" else "none",
      axis.title.x = if (show_x) element_text() else element_blank(),
      axis.text.x = if (show_x) element_text() else element_text()
    )
}

p1 <- plot_bar(narea_data, "Narea", ref_narea, show_x = FALSE, show_legend = FALSE)
p2 <- plot_bar(sla_data, "SLA", ref_sla, show_x = FALSE, show_legend = FALSE)
p3 <- plot_bar(plsr_sla_data, "PLSR-SLA", ref_plsr_sla, show_x = FALSE, show_legend = FALSE)
p4 <- plot_bar(plsr_narea_data, "PLSR-Narea", ref_plsr_narea, show_x = TRUE, show_legend = FALSE)

# shared legend
leg <- get_legend(plot_bar(narea_data, "Narea", ref_narea, show_x = FALSE, show_legend = TRUE))

# stack + labels
panels <- plot_grid(p1, p2, p3, p4,
                    ncol = 1, labels = c("a","b","c","d"),
                    label_size = 12, align = "v")

# shared y label
with_y <- plot_grid(
  ggdraw() + draw_label("Value", angle = 90, vjust = 0.5),
  panels,
  ncol = 2,
  rel_widths = c(0.06, 1)
)

final_plot <- plot_grid(with_y, leg, ncol = 2, rel_widths = c(1, 0.18))

print(final_plot)

ggsave("./figure/rerun_h2comparisoncompletemodel.pdf", plot = final_plot,
       width = 7, height = 10, dpi = 300)


# =========================
# COHERITABILITY VALUES
# =========================
coh2_narea <- c(S0 = 0, S1 = 0.49, S2 = 0.49, S3 = 0.49)
coh2_sla   <- c(S0 = 0, S1 = 0.52, S2 = 0.50, S3 = 0.49)
coh2_pn    <- c(S0 = 0, S1 = 0.54, S2 = 0.54, S3 = 0.54)
coh2_ps    <- c(S0 = 0, S1 = 0.47, S2 = 0.47, S3 = 0.47)

# =========================
# X labels: force 2 decimals, except S0 = (0)
# =========================
lab_from_coh2 <- function(x) {
  fmt <- function(v) sprintf("%.2f", v)
  c(
    ST = "ST",
    S0 = "S0\n(0)",
    S1 = paste0("S1\n(", fmt(x["S1"]), ")"),
    S2 = paste0("S2\n(", fmt(x["S2"]), ")"),
    S3 = paste0("S3\n(", fmt(x["S3"]), ")")
  )
}
lab_N  <- lab_from_coh2(coh2_narea)
lab_S  <- lab_from_coh2(coh2_sla)
lab_PN <- lab_from_coh2(coh2_pn)
lab_PS <- lab_from_coh2(coh2_ps)

# =========================
# READERS (mixed formats)
# =========================
read_acc_two_col <- function(path) {
  read.csv(path) %>% select(acc)
}
read_acc_one_col <- function(path) {
  read.csv(path, header = FALSE) %>%
    transmute(acc = suppressWarnings(as.numeric(V1))) %>%
    filter(!is.na(acc))
}

# =========================
# COLORS (match sample)
# =========================
fill_cols <- c(ST = "#F8766D", CV1 = "#00BA38", CV2 = "#619CFF")

# =========================
# THEME to match the "paper" plot:
# - rectangle border around each panel
# - no extra axis titles inside panels
# - minimal grid
# =========================
paper_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0, size = 10),
    # rectangle border around each panel
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    # keep some space
    plot.margin = margin(5, 5, 5, 5)
  )

# =========================
# ---- BUILD DATA (same as before) ----
# If your existing N_data/S_data/PN_data/PS_data are already created,
# you can SKIP this whole section and start from "PLOTS" below.
# =========================

N_data <- bind_rows(
  read.csv("./output/cv_singletrait/corr_narea_EFMW.csv") %>%
    transmute(acc = res.ac, MT="ST", Type="ST"),
  read_acc_two_col("./output/acNT_CV1.csv") %>% mutate(MT="S0", Type="CV1"),
  read_acc_two_col("./output/acNT_CV2.csv") %>% mutate(MT="S0", Type="CV2"),
  read_acc_one_col("./output/acNW1_CV1.csv") %>% mutate(MT="S1", Type="CV1"),
  read_acc_one_col("./output/acNW2_CV1.csv") %>% mutate(MT="S2", Type="CV1"),
  read_acc_one_col("./output/acNW3_CV1.csv") %>% mutate(MT="S3", Type="CV1"),
  read_acc_one_col("./output/acNW1_CV2.csv") %>% mutate(MT="S1", Type="CV2"),
  read_acc_one_col("./output/acNW2_CV2.csv") %>% mutate(MT="S2", Type="CV2"),
  read_acc_one_col("./output/acNW3_CV2.csv") %>% mutate(MT="S3", Type="CV2")
) %>% mutate(
  MT = factor(MT, levels = c("ST","S0","S1","S2","S3")),
  Type = factor(Type, levels = c("ST","CV1","CV2"))
)

S_data <- bind_rows(
  read.csv("./output/cv_singletrait/corr_sla_EFMW.csv") %>%
    transmute(acc = res.ac, MT="ST", Type="ST"),
  read_acc_two_col("./output/acST_CV1.csv") %>% mutate(MT="S0", Type="CV1"),
  read_acc_two_col("./output/acST_CV2.csv") %>% mutate(MT="S0", Type="CV2"),
  read_acc_one_col("./output/acSW1_CV1.csv") %>% mutate(MT="S1", Type="CV1"),
  read_acc_one_col("./output/acSW2_CV1.csv") %>% mutate(MT="S2", Type="CV1"),
  read_acc_one_col("./output/acSW3_CV1.csv") %>% mutate(MT="S3", Type="CV1"),
  read_acc_one_col("./output/acSW1_CV2.csv") %>% mutate(MT="S1", Type="CV2"),
  read_acc_one_col("./output/acSW2_CV2.csv") %>% mutate(MT="S2", Type="CV2"),
  read_acc_one_col("./output/acSW3_CV2.csv") %>% mutate(MT="S3", Type="CV2")
) %>% mutate(
  MT = factor(MT, levels = c("ST","S0","S1","S2","S3")),
  Type = factor(Type, levels = c("ST","CV1","CV2"))
)

PN_data <- bind_rows(
  read.csv("./output/cv_singletrait/corr_plsr_narea_EFMW.csv") %>%
    transmute(acc = res.ac, MT="ST", Type="ST"),
  read_acc_two_col("./output/acpnT_CV1.csv") %>% mutate(MT="S0", Type="CV1"),
  read_acc_two_col("./output/acpnT_CV2.csv") %>% mutate(MT="S0", Type="CV2"),
  read_acc_one_col("./output/acpnW1_CV1.csv") %>% mutate(MT="S1", Type="CV1"),
  read_acc_one_col("./output/acpnW2_CV1.csv") %>% mutate(MT="S2", Type="CV1"),
  read_acc_one_col("./output/acpnW3_CV1.csv") %>% mutate(MT="S3", Type="CV1"),
  read_acc_one_col("./output/acpnW1_CV2.csv") %>% mutate(MT="S1", Type="CV2"),
  read_acc_one_col("./output/acpnW2_CV2.csv") %>% mutate(MT="S2", Type="CV2"),
  read_acc_one_col("./output/acpnW3_CV2.csv") %>% mutate(MT="S3", Type="CV2")
) %>% mutate(
  MT = factor(MT, levels = c("ST","S0","S1","S2","S3")),
  Type = factor(Type, levels = c("ST","CV1","CV2"))
)

PS_data <- bind_rows(
  read.csv("./output/cv_singletrait/corr_plsr_sla_EFMW.csv") %>%
    transmute(acc = res.ac, MT="ST", Type="ST"),
  read_acc_two_col("./output/acpsT_CV1.csv") %>% mutate(MT="S0", Type="CV1"),
  read_acc_two_col("./output/acpsT_CV2.csv") %>% mutate(MT="S0", Type="CV2"),
  read_acc_one_col("./output/acpsW1_CV1.csv") %>% mutate(MT="S1", Type="CV1"),
  read_acc_one_col("./output/acpsW2_CV1.csv") %>% mutate(MT="S2", Type="CV1"),
  read_acc_one_col("./output/acpsW3_CV1.csv") %>% mutate(MT="S3", Type="CV1"),
  read_acc_one_col("./output/acpsW1_CV2.csv") %>% mutate(MT="S1", Type="CV2"),
  read_acc_one_col("./output/acpsW2_CV2.csv") %>% mutate(MT="S2", Type="CV2"),
  read_acc_one_col("./output/acpsW3_CV2.csv") %>% mutate(MT="S3", Type="CV2")
) %>% mutate(
  MT = factor(MT, levels = c("ST","S0","S1","S2","S3")),
  Type = factor(Type, levels = c("ST","CV1","CV2"))
)

# =========================
# PLOTS (NO "acc" labels, NO "MT", panel borders)
# =========================
p1 <- ggplot(N_data, aes(MT, acc, fill = Type)) +
  geom_boxplot(width = 0.6, outlier.size = 0.7) +
  scale_fill_manual(values = fill_cols, name = "Type") +
  scale_x_discrete(labels = lab_N) +
  labs(title = "Narea") +
  paper_theme + theme(legend.position = "none")

p2 <- ggplot(S_data, aes(MT, acc, fill = Type)) +
  geom_boxplot(width = 0.6, outlier.size = 0.7) +
  scale_fill_manual(values = fill_cols, name = "Type") +
  scale_x_discrete(labels = lab_S) +
  labs(title = "SLA") +
  paper_theme + theme(legend.position = "none")

p3 <- ggplot(PN_data, aes(MT, acc, fill = Type)) +
  geom_boxplot(width = 0.6, outlier.size = 0.7) +
  scale_fill_manual(values = fill_cols, name = "Type") +
  scale_x_discrete(labels = lab_PN) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + # 2 decimals for panel c
  labs(title = "PLSR-Narea") +
  paper_theme + theme(legend.position = "none")

p4 <- ggplot(PS_data, aes(MT, acc, fill = Type)) +
  geom_boxplot(width = 0.6, outlier.size = 0.7) +
  scale_fill_manual(values = fill_cols, name = "Type") +
  scale_x_discrete(labels = lab_PS) +
  labs(title = "PLSR-SLA") +
  paper_theme + theme(legend.position = "none")

# Shared legend
leg <- get_legend(
  ggplot(N_data, aes(MT, acc, fill = Type)) +
    geom_boxplot() +
    scale_fill_manual(values = fill_cols, name = "Type") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right")
)

# Stack + tags ad
panels <- plot_grid(
  p1, p2, p3, p4,
  ncol = 1,
  labels = c("a", "b", "c", "d"),
  align = "v",
  label_size = 12
)

# Add ONE shared y-axis label (Accuracy) on the left
with_y <- plot_grid(
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5),
  panels,
  ncol = 2,
  rel_widths = c(0.06, 1)
)

# Add legend on right
final_plot <- plot_grid(
  with_y, leg,
  ncol = 2,
  rel_widths = c(1, 0.18)
)

# Add ONE shared x-axis label (Model Type) at bottom
final_plot2 <- plot_grid(
  final_plot,
  ggdraw() + draw_label("Model Type", vjust = 1),
  ncol = 1,
  rel_heights = c(1, 0.06)
)

print(final_plot2)

ggsave("./figure/rerun_allcompletemodel_comparison.pdf",
       plot = final_plot2, width = 7, height = 10, dpi = 300) #this gives us complete mdoel all trait comparison




#now the script is for comparing subset and complete model ST and MT models with two synthetic traits
#ST vs MT 2 trait vs MT one trait
library(dplyr)
library(ggplot2)
library(scales)

# =========================
# Narea (ST vs MT)
# =========================

acN <- read.csv("./output/corr_N_EFMW.csv")
st_narea <- acN %>%
  mutate(accuracy = result_N.ac, model = "ST", cv = "ST") %>%
  select(model, cv, accuracy)

mt_narea_cv1 <- read.csv("./output/acc_mt_narea_CV1.csv") %>%
  mutate(model = "MT(S1+S2)", cv = "CV1")
mt_narea_cv2 <- read.csv("./output/MT_Narea_3Trait_CV2_FINAL_accuracy.csv") %>%
  mutate(model = "MT(S1+S2)", cv = "CV2")

if ("acc" %in% names(mt_narea_cv1)) mt_narea_cv1 <- mt_narea_cv1 %>% mutate(accuracy = acc)
if ("acc" %in% names(mt_narea_cv2)) mt_narea_cv2 <- mt_narea_cv2 %>% mutate(accuracy = acc)

mt_narea_cv1 <- mt_narea_cv1 %>% select(model, cv, accuracy)
mt_narea_cv2 <- mt_narea_cv2 %>% select(model, cv, accuracy)

n_s1_cv1 <- read.csv("./output/acNW1_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S1)", cv = "CV1") %>% select(model, cv, accuracy)
n_s2_cv1 <- read.csv("./output/acNW2_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S2)", cv = "CV1") %>% select(model, cv, accuracy)
n_s3_cv1 <- read.csv("./output/acNW3_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S3)", cv = "CV1") %>% select(model, cv, accuracy)

n_s1_cv2 <- read.csv("./output/acNW1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S1)", cv = "CV2") %>% select(model, cv, accuracy)
n_s2_cv2 <- read.csv("./output/acNW2.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S2)", cv = "CV2") %>% select(model, cv, accuracy)
n_s3_cv2 <- read.csv("./output/acNW3.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S3)", cv = "CV2") %>% select(model, cv, accuracy)

narea_plot_df <- bind_rows(
  st_narea,
  mt_narea_cv1, mt_narea_cv2,
  n_s1_cv1, n_s2_cv1, n_s3_cv1,
  n_s1_cv2, n_s2_cv2, n_s3_cv2
) %>%
  mutate(
    cv = factor(cv, levels = c("ST", "CV1", "CV2")),
    model = factor(model, levels = c("ST", "MT(S1+S2)", "MT(S1)", "MT(S2)", "MT(S3)"))
  )

# =========================
# SLA (ST vs MT)
# =========================

acS <- read.csv("./output/corr_S_EFMW.csv")
st_sla <- acS %>%
  mutate(accuracy = result_S.ac, model = "ST", cv = "ST") %>%
  select(model, cv, accuracy)

mt_sla_cv1 <- read.csv("./output/acc_mt_sla_CV1.csv") %>%
  mutate(model = "MT(S1+S2)", cv = "CV1")
mt_sla_cv2 <- read.csv("./output/MT_SLA_3Trait_CV2_accuracy.csv") %>%
  mutate(model = "MT(S1+S2)", cv = "CV2")

if ("acc" %in% names(mt_sla_cv1)) mt_sla_cv1 <- mt_sla_cv1 %>% mutate(accuracy = acc)
if ("acc" %in% names(mt_sla_cv2)) mt_sla_cv2 <- mt_sla_cv2 %>% mutate(accuracy = acc)

mt_sla_cv1 <- mt_sla_cv1 %>% select(model, cv, accuracy)
mt_sla_cv2 <- mt_sla_cv2 %>% select(model, cv, accuracy)

s_s1_cv1 <- read.csv("./output/acSW1_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S1)", cv = "CV1") %>% select(model, cv, accuracy)
s_s2_cv1 <- read.csv("./output/acSW2_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S2)", cv = "CV1") %>% select(model, cv, accuracy)
s_s3_cv1 <- read.csv("./output/acSW3_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S3)", cv = "CV1") %>% select(model, cv, accuracy)

s_s1_cv2 <- read.csv("./output/acSW1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S1)", cv = "CV2") %>% select(model, cv, accuracy)
s_s2_cv2 <- read.csv("./output/acSW2.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S2)", cv = "CV2") %>% select(model, cv, accuracy)
s_s3_cv2 <- read.csv("./output/acSW3.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S3)", cv = "CV2") %>% select(model, cv, accuracy)

sla_plot_df <- bind_rows(
  st_sla,
  mt_sla_cv1, mt_sla_cv2,
  s_s1_cv1, s_s2_cv1, s_s3_cv1,
  s_s1_cv2, s_s2_cv2, s_s3_cv2
) %>%
  mutate(
    cv = factor(cv, levels = c("ST", "CV1", "CV2")),
    model = factor(model, levels = c("ST", "MT(S1+S2)", "MT(S1)", "MT(S2)", "MT(S3)"))
  )

# =========================
# PS (PLSR-SLA) (ST vs MT)
# =========================

acps <- read.csv("./output/corr_ps_EFMW.csv") %>%
  mutate(accuracy = result_ps.ac) %>%
  select(-result_ps.ac)

st_ps <- acps %>%
  mutate(model = "ST", cv = "ST") %>%
  select(model, cv, accuracy)

mt_ps_cv1 <- read.csv("./output/acc_mt_ps_CV1.csv") %>%
  mutate(model = "MT(S1+S2)", cv = "CV1")
mt_ps_cv2 <- read.csv("./output/MT_ps_3Trait_CV2_accuracy.csv") %>%
  mutate(model = "MT(S1+S2)", cv = "CV2")

if ("acc" %in% names(mt_ps_cv1)) mt_ps_cv1 <- mt_ps_cv1 %>% mutate(accuracy = acc)
if ("acc" %in% names(mt_ps_cv2)) mt_ps_cv2 <- mt_ps_cv2 %>% mutate(accuracy = acc)

mt_ps_cv1 <- mt_ps_cv1 %>% select(model, cv, accuracy)
mt_ps_cv2 <- mt_ps_cv2 %>% select(model, cv, accuracy)

ps_s1_cv1 <- read.csv("./output/acpsW1_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S1)", cv = "CV1") %>% select(model, cv, accuracy)
ps_s2_cv1 <- read.csv("./output/acpsW2_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S2)", cv = "CV1") %>% select(model, cv, accuracy)
ps_s3_cv1 <- read.csv("./output/acpsW3_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S3)", cv = "CV1") %>% select(model, cv, accuracy)

ps_s1_cv2 <- read.csv("./output/acpsW1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S1)", cv = "CV2") %>% select(model, cv, accuracy)
ps_s2_cv2 <- read.csv("./output/acpsW2.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S2)", cv = "CV2") %>% select(model, cv, accuracy)
ps_s3_cv2 <- read.csv("./output/acpsW3.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S3)", cv = "CV2") %>% select(model, cv, accuracy)

ps_plot_df <- bind_rows(
  st_ps,
  mt_ps_cv1, mt_ps_cv2,
  ps_s1_cv1, ps_s2_cv1, ps_s3_cv1,
  ps_s1_cv2, ps_s2_cv2, ps_s3_cv2
) %>%
  mutate(
    cv = factor(cv, levels = c("ST", "CV1", "CV2")),
    model = factor(model, levels = c("ST", "MT(S1+S2)", "MT(S1)", "MT(S2)", "MT(S3)"))
  )

# =========================
# PN (PLSR-Narea) (ST vs MT)
# =========================

acpn <- read.csv("./output/corr_pn_EFMW.csv") %>%
  mutate(accuracy = result_pn.ac) %>%
  select(-result_pn.ac)

st_pn <- acpn %>%
  mutate(model = "ST", cv = "ST") %>%
  select(model, cv, accuracy)

mt_pn_cv1 <- read.csv("./output/acc_mt_pn_CV1.csv") %>%
  mutate(model = "MT(S1+S2)", cv = "CV1")
mt_pn_cv2 <- read.csv("./output/acc_mt_pn_CV2.csv") %>%
  mutate(model = "MT(S1+S2)", cv = "CV2")

if ("acc" %in% names(mt_pn_cv1)) mt_pn_cv1 <- mt_pn_cv1 %>% mutate(accuracy = acc)
if ("acc" %in% names(mt_pn_cv2)) mt_pn_cv2 <- mt_pn_cv2 %>% mutate(accuracy = acc)

mt_pn_cv1 <- mt_pn_cv1 %>% select(model, cv, accuracy)
mt_pn_cv2 <- mt_pn_cv2 %>% select(model, cv, accuracy)

pn_s1_cv1 <- read.csv("./output/acpnW1_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S1)", cv = "CV1") %>% select(model, cv, accuracy)
pn_s2_cv1 <- read.csv("./output/acpnW2_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S2)", cv = "CV1") %>% select(model, cv, accuracy)
pn_s3_cv1 <- read.csv("./output/acpnW3_CV1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S3)", cv = "CV1") %>% select(model, cv, accuracy)

pn_s1_cv2 <- read.csv("./output/acpnW1.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S1)", cv = "CV2") %>% select(model, cv, accuracy)
pn_s2_cv2 <- read.csv("./output/acpnW2.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S2)", cv = "CV2") %>% select(model, cv, accuracy)
pn_s3_cv2 <- read.csv("./output/acpnW3.txt", header = FALSE) %>%
  mutate(accuracy = V1, model = "MT(S3)", cv = "CV2") %>% select(model, cv, accuracy)

pn_plot_df <- bind_rows(
  st_pn,
  mt_pn_cv1, mt_pn_cv2,
  pn_s1_cv1, pn_s2_cv1, pn_s3_cv1,
  pn_s1_cv2, pn_s2_cv2, pn_s3_cv2
) %>%
  mutate(
    cv = factor(cv, levels = c("ST", "CV1", "CV2")),
    model = factor(model, levels = c("ST", "MT(S1+S2)", "MT(S1)", "MT(S2)", "MT(S3)"))
  )

# =========================
# Final 4-panel stacked plot (facet, legend right, strip top-left)
# =========================

all4_df <- bind_rows(
  narea_plot_df %>% mutate(panel = "Narea"),
  sla_plot_df   %>% mutate(panel = "SLA"),
  pn_plot_df    %>% mutate(panel = "PLSR-Narea"),
  ps_plot_df    %>% mutate(panel = "PLSR-SLA")
) %>%
  mutate(
    panel = factor(
      panel,
      levels = c("Narea", "SLA", "PLSR-Narea", "PLSR-SLA")  # <- NEW ORDER
    )
  )
p_all <- ggplot(all4_df, aes(x = model, y = accuracy, fill = cv)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 4),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  facet_wrap(~panel, ncol = 1, scales = "free_y", strip.position = "top") +
  labs(x = "Model Type", y = "Accuracy", fill = "Type") +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 11, hjust = 0),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.spacing = unit(0.7, "lines")
  )

print(p_all)

ggsave("./figure/all4_traits_ST_vs_MT.png", p_all, width = 5.8, height = 10, dpi = 300)


#for subset model i.e 5 reps for each traits

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
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for plsr-narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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

plot1 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-Narea",
       x = "Model Type",
       y = "Predictive ability") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, max(all_data$acc), by = 0.03))+

  # Adjust theme for clarity
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")


#plsr-sla(EFMW)
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
ps_data_whole<- clean_data(ps_data_whole)

#replication 1
#ST for pnarea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep1.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep1.ac) %>% select(-result_ps_rep1.ac)
#Multi trait cv1 for plsr-sla psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for plsr-sla psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep1<- clean_data(ps_data_rep1)
#replication 2
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep2.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep2.ac) %>% select(-result_ps_rep2.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("acpsW3_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
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
ps_data_rep2<- clean_data(ps_data_rep2)
#replication 3
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep3.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep3.ac) %>% select(-result_ps_rep3.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep3<- clean_data(ps_data_rep3)

#replication 4
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep4.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep4.ac) %>% select(-result_ps_rep4.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep4<- clean_data(ps_data_rep4)
#replication 5
#ST for narea
corr_ps_EFMW<- read.csv("./output/./output/corr_ps_EFMWrep5.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep5.ac) %>% select(-result_ps_rep5.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep5<- clean_data(ps_data_rep5)
# Combine all datasets
all_data <- bind_rows(ps_data_whole, ps_data_rep1, ps_data_rep2, ps_data_rep3, ps_data_rep4, ps_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

plot2 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-SLA",
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
N_data_whole$Data <- "Complete"
N_data_whole<- clean_data(N_data_whole)
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
N_data_rep1<- clean_data(N_data_rep1)
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
N_data_rep2<- clean_data(N_data_rep2)

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
N_data_rep3<- clean_data(N_data_rep3)

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
N_data_rep4<- clean_data(N_data_rep4)

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
N_data_rep5<- clean_data(N_data_rep5)

all_data <- bind_rows(N_data_whole, N_data_rep1, N_data_rep2, N_data_rep3, N_data_rep4, N_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))
# Create the ggplot with the adjusted legend order
plot3 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "Narea",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, max(all_data$acc), by = 0.035),
                     labels = scales::number_format(accuracy = 0.01))+

  # Adjust theme for clarity
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")
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
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.57)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.54)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.57)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.54)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")

# Combine all data
S_data_whole <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1,acSW1, acSW2, acSW3)
S_data_whole$Data <- "Complete"
S_data_whole<- clean_data(S_data_whole)
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
S_data_rep1<- clean_data(S_data_rep1)
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
S_data_rep2<- clean_data(S_data_rep2)
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
S_data_rep3<- clean_data(S_data_rep3)

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
S_data_rep4<- clean_data(S_data_rep4)

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
S_data_rep5<- clean_data(S_data_rep5)

# Combine all datasets
all_data <- bind_rows(S_data_whole, S_data_rep1, S_data_rep2, S_data_rep3, S_data_rep4, S_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

plot4 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "SLA",
       x = "Model Type",
       y = "Accuracy") +
  facet_grid(~ Data, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, max(all_data$acc), by = 0.035),
                     labels = scales::number_format(accuracy = 0.01))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")


#making  four traits efmw model into one figure
#Remove the x-axis title for the first three plots
plot1_no_legend <- plot1 + theme(legend.position = "none", axis.title.y= element_blank(),axis.text.x=element_blank(), axis.title.x=element_blank())
plot2_no_legend <- plot2 + theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank())
plot3_no_legend <- plot3 + theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x=element_blank(), axis.title.y=element_blank())
plot4_no_legend <- plot4 + theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y= element_blank(), axis.text.x=element_blank())  # Keep x-axis title for the last plot

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
  ggdraw() + draw_label("Accuracy", angle = 90, vjust = 0.5, hjust = 0.5,size= 24)+
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
#print(final_plot)
ggsave("./figure/firstfour.pdf", plot = final_plot, width = 20, height = 20, dpi = 300)






# for MWEF of last four fig
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
pn_data_mwef<- clean_data(pn_data_mwef)
#replication 1
#ST for pnarea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep1.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep1.ac) %>% select(-result_pn_rep1.ac)
#Multi trait cv1 for plsr-narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for plsr-narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
pn_data_rep1<- clean_data(pn_data_rep1)

#replication 2
#ST for narea
corr_pn_MWEF<- read.csv("./output/./output/corr_pn_MWEFrep2.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep2.ac) %>% select(-result_pn_rep2.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
#replication 3
#ST for narea
corr_pn_MWEF<- read.csv("./output/./output/corr_pn_MWEFrep3.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep3.ac) %>% select(-result_pn_rep3.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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


#replication 4
#ST for narea
corr_pn_MWEF<- read.csv("./output/./output/corr_pn_MWEFrep4.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep4.ac) %>% select(-result_pn_rep4.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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

#replication 5
#ST for narea
corr_pn_MWEF<- read.csv("./output/./output/corr_pn_MWEFrep5.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep5.ac) %>% select(-result_pn_rep5.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnW1_CV1<- read.csv("./output/acpnW1_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2_CV1<- read.csv("./output/acpnW2_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3_CV1<- read.csv("./output/acpnW3_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnW1<- read.csv("./output/acpnW1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW2<- read.csv("./output/acpnW2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnW3<- read.csv("./output/acpnW3_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
all_data <- bind_rows(pn_data_mwef, pn_data_rep1, pn_data_rep2, pn_data_rep3, pn_data_rep4, pn_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

plot5 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-Narea",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") +

  # Add annotations for coh2 values for each facet with dynamic y position and label inside a box
  geom_label(data = all_data %>% filter(MT %in% c("S1", "S2", "S3"), !is.na(coh2)) %>% distinct(MT, Data, .keep_all = TRUE),
             aes(label = round(coh2, 2), x = MT, group = Data,
                 y = case_when(
                   MT == "S1" ~ max(acc) + 0.05,  # Position for S1
                   MT == "S2" ~ max(acc) + 0.07,  # Higher position for S2
                   MT == "S3" ~ max(acc) + 0.09   # Even higher for S3
                 )),
             color = "black",        # Text color
             fill = "white",         # Box fill color
             size = 6,               # Text size
             fontface = "bold",      # Bold text
             label.size = 0.1) +     # Border size around the box

  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")



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
ps_data_whole<- clean_data(ps_data_whole)
#replication 1
#ST for pnarea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep1.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep1.ac) %>% select(-result_ps_rep1.ac)
#Multi trait cv1 for plsr-sla psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for plsr-sla psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep1<- clean_data(ps_data_rep1)
#replication 2
#ST for narea
corr_ps_EFMW<- read.csv("./output/./output/corr_ps_EFMWrep2.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep2.ac) %>% select(-result_ps_rep2.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep2<- clean_data(ps_data_rep2)
#replication 3
#ST for narea
corr_ps_EFMW<- read.csv("./output/./output/corr_ps_EFMWrep3.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep3.ac) %>% select(-result_ps_rep3.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep3<- clean_data(ps_data_rep3)
#replication 4
#ST for narea
corr_ps_EFMW<- read.csv("./output/./output/corr_ps_EFMWrep4.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep4.ac) %>% select(-result_ps_rep4.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep4<- clean_data(ps_data_rep4)
#replication 5
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep5.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep5.ac) %>% select(-result_ps_rep5.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsW1_CV1<- read.csv("./output/acpsW1_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2_CV1<- read.csv("./output/acpsW2_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3_CV1<- read.csv("./output/acpsW3_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsW1<- read.csv("./output/acpsW1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW2<- read.csv("./output/acpsW2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpsW3<- read.csv("./output/acpsW3_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

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
ps_data_rep5<- clean_data(ps_data_rep5)

# Combine all datasets
all_data <- bind_rows(ps_data_whole, ps_data_rep1, ps_data_rep2, ps_data_rep3, ps_data_rep4, ps_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

plot6 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "PLSR-SLA",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") +

  # Add annotations for coh2 values for each facet with dynamic y position and label inside a box
  geom_label(data = all_data %>% filter(MT %in% c("S1", "S2", "S3"), !is.na(coh2)) %>% distinct(MT, Data, .keep_all = TRUE),
             aes(label = round(coh2, 2), x = MT, group = Data,
                 y = case_when(
                   MT == "S1" ~ max(acc) + 0.05,  # Position for S1
                   MT == "S2" ~ max(acc) + 0.07,  # Higher position for S2
                   MT == "S3" ~ max(acc) + 0.09   # Even higher for S3
                 )),
             color = "black",        # Text color
             fill = "white",         # Box fill color
             size = 6,               # Text size
             fontface = "bold",      # Bold text
             label.size = 0.1) +     # Border size around the box

  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

#MWEF narea
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEF.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N.ac) %>% select(-result_N.ac)
#Multi trait CV1_mwef for narea NW1, NW2, NW3
acNW1_CV1<- read.csv("./output/acNW1_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2_CV1<- read.csv("./output/acNW2_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3_CV1<- read.csv("./output/acNW3_CV1_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNW1<- read.csv("./output/acNW1_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW2<- read.csv("./output/acNW2_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acNW3<- read.csv("./output/acNW3_CV2_mwef.csv",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNW1_CV1 <- acNW1_CV1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV1")
acNW2_CV1 <- acNW2_CV1 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV1")
acNW3_CV1 <- acNW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acNW1 <- acNW1 %>% mutate(MT = "S1(coh2=0.49)", Type = "CV2")
acNW2 <- acNW2 %>% mutate(MT = "S2(coh2=0.49)", Type = "CV2")
acNW3 <- acNW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")

# Combine all data
N_data_mwef_whole <- bind_rows(acN, acNW1_CV1, acNW2_CV1, acNW3_CV1, acNW1, acNW2, acNW3)
N_data_mwef_whole$Data <- "Complete"
N_data_mwef_whole<- clean_data(N_data_mwef_whole)
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
N_data_rep1<- clean_data(N_data_rep1)
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
N_data_rep2<- clean_data(N_data_rep2)

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
N_data_rep3<- clean_data(N_data_rep3)

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
N_data_rep4<- clean_data(N_data_rep4)

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
N_data_rep5<- clean_data(N_data_rep5)

# Combine all datasets
all_data <- bind_rows(N_data_mwef_whole, N_data_rep1, N_data_rep2, N_data_rep3, N_data_rep4, N_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))
plot7 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "Narea",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") +

  # Add annotations for coh2 values for each facet with dynamic y position and label inside a box
  geom_label(data = all_data %>% filter(MT %in% c("S1", "S2", "S3"), !is.na(coh2)) %>% distinct(MT, Data, .keep_all = TRUE),
             aes(label = round(coh2, 2), x = MT, group = Data,
                 y = case_when(
                   MT == "S1" ~ max(acc) + 0.05,  # Position for S1
                   MT == "S2" ~ max(acc) + 0.07,  # Higher position for S2
                   MT == "S3" ~ max(acc) + 0.09   # Even higher for S3
                 )),
             color = "black",        # Text color
             fill = "white",         # Box fill color
             size = 6,               # Text size
             fontface = "bold",      # Bold text
             label.size = 0.1) +     # Border size around the box

  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")


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
acSW1_CV1 <- acSW1_CV1 %>% mutate(MT = "S1(coh2=0.57)", Type = "CV1")
acSW2_CV1 <- acSW2_CV1 %>% mutate(MT = "S2(coh2=0.54)", Type = "CV1")
acSW3_CV1 <- acSW3_CV1 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV1")

acSW1 <- acSW1 %>% mutate(MT = "S1(coh2=0.57)", Type = "CV2")
acSW2 <- acSW2 %>% mutate(MT = "S2(coh2=0.54)", Type = "CV2")
acSW3 <- acSW3 %>% mutate(MT = "S3(coh2=0.49)", Type = "CV2")

# Combine all data
S_data_mwef <- bind_rows(acS, acSW1_CV1, acSW2_CV1, acSW3_CV1, acSW1, acSW2, acSW3)
S_data_mwef$Data <- "Complete"
S_data_mwef<- clean_data(S_data_mwef)
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
S_data_rep1<- clean_data(S_data_rep1)
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
S_data_rep2<- clean_data(S_data_rep2)


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
S_data_rep3<- clean_data(S_data_rep3)

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
S_data_rep4<- clean_data(S_data_rep4)

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
S_data_rep5<- clean_data(S_data_rep5)

# Combine all dataset
all_data <- bind_rows(S_data_mwef, S_data_rep1, S_data_rep2, S_data_rep3, S_data_rep4, S_data_rep5)

# Adjust facet labels for better clarity
all_data$Data <- factor(all_data$Data, levels = c("Complete", "Subset-1", "Subset-2", "Subset-3", "Subset-4", "Subset-5"))
all_data$Type <- factor(all_data$Type, levels = c("ST", "CV1", "CV2"))

plot8 <- ggplot(all_data, aes(x = MT, y = acc, fill = Type)) +
  geom_boxplot() +
  common_theme +
  labs(title = "SLA",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") +

  # Add annotations for coh2 values for each facet with dynamic y position and label inside a box
  geom_label(data = all_data %>% filter(MT %in% c("S1", "S2", "S3"), !is.na(coh2)) %>% distinct(MT, Data, .keep_all = TRUE),
             aes(label = round(coh2, 2), x = MT, group = Data,
                 y = max(acc) + 0.05),  # Dynamic y position above the boxplot
             color = "red",          # Text color
             fill = "white",           # Box fill color
             size = 6,                 # Text size
             fontface = "bold",        # Bold text
             label.size = 0.5) +       # Border size around the box

  # Adjust theme for clarity
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

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



#result 3rd part
#figure 4 i.e single trait and coh2 0 for each trait complete and 5 subsets
#plsr-narea
#ST for plsr-narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMW.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn.ac) %>% select(-result_pn.ac)
acpnT_CV1<- read.csv("./output/acpnT_CV1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("./output/acpnT.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apnd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1<- acpnT_CV1 %>% mutate(MT= "S-0", Type= "CV1")
acpnT<- acpnT %>% mutate (MT= "S-0", Type ="CV2")
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
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
pn_data_rep1 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep2.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep2.ac) %>% select(-result_pn_rep2.ac)
acpnT_CV1<- read.csv("./output/acpnNT2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("./output/acpnNT2_CV2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
pn_data_rep2 <- bind_rows(acpn,acpnT_CV1, acpnT)
pn_data_rep2$Data <- "Subset-2"
#replication 3
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep3.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep3.ac) %>% select(-result_pn_rep3.ac)
acpnT_CV1<- read.csv("./output/acpnNT3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("./output/acpnNT3_CV2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
pn_data_rep3 <- bind_rows(acpn,acpnT_CV1, acpnT)
pn_data_rep3$Data <- "Subset-3"

#replication 4
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep4.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep4.ac) %>% select(-result_pn_rep4.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("./output/acpnNT4_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("./output/acpnNT4_CV2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
pn_data_rep4 <-bind_rows(acpn,acpnT_CV1, acpnT)
pn_data_rep4$Data <- "Subset-4"
#replication 5
#ST for narea
corr_pn_EFMW<- read.csv("./output/corr_pn_EFMWrep5.csv")
acpn <- corr_pn_EFMW %>% mutate(acc = result_pn_rep5.ac) %>% select(-result_pn_rep5.ac)
acpnT_CV1<- read.csv("./output/acpnNT5_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("./output/acpnNT5_CV2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

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
        legend.text = element_text(size = 20),    # Adjust legend text size
        axis.text = element_text(size = 20),      # Adjust axis text size (increased from 16)
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
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT<- acpsT %>%  mutate(MT= "S-0", Type="CV2")
# Combine all data
ps_data_whole <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_whole$Data<- "Complete"


#replication 1
#ST for pnarea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep1.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep1.ac) %>% select(-result_ps_rep1.ac)
#Multi trait cv1 for plsr-sla psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for plsr-sla psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT1_CV2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


# Add labels for multi-trait aSd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combinne all data
ps_data_rep1 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep2.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep2.ac) %>% select(-result_ps_rep2.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT2_CV2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
ps_data_rep2 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep3.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep3.ac) %>% select(-result_ps_rep3.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT3_CV2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
ps_data_rep3 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep4.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep4.ac) %>% select(-result_ps_rep4.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT4_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT4_CV2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
ps_data_rep4 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_ps_EFMW<- read.csv("./output/corr_ps_EFMWrep5.csv")
acps <- corr_ps_EFMW %>% mutate(acc = result_ps_rep5.ac) %>% select(-result_ps_rep5.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT5_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT5_CV2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

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
  common_theme +
  labs(title = "PLSR-SLA",#(Training EF & validating on MW)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

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
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT<- acNT %>%  mutate(MT= "S-0",Type= "CV2")
# Combine all data
N_data_whole <- bind_rows(acN,acNT_CV1,acNT)#acNT_CV1,acNT
N_data_whole$Data <- "Complete"

#replication 1
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep1.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep1.ac) %>% select(-result_N_rep1.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT1_CV2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_rep1 <- bind_rows(acN,acNT_CV1,acNT)
N_data_rep1$Data <- "Subset-1"

#Set 2
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep2.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep2.ac) %>% select(-result_N_rep2.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT2_CV2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_rep2 <- bind_rows(acN,acNT_CV1,acNT)
N_data_rep2$Data <- "Subset-2"
#Set 3
#ST for narea
corr_N_EFMW<- read.csv("./output/./output/corr_N_EFMWrep3.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep3.ac) %>% select(-result_N_rep3.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT3_CV2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_rep3 <- bind_rows(acN,acNT_CV1,acNT)
N_data_rep3$Data <- "Subset-3"

#Set 4
#ST for narea
corr_N_EFMW<- read.csv("./output/./output/corr_N_EFMWrep4.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep4.ac) %>% select(-result_N_rep4.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT4_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT4_CV2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_rep4 <- bind_rows(acN,acNT_CV1,acNT)
N_data_rep4$Data <- "Subset-4"
#Set 5
#ST for narea
corr_N_EFMW<- read.csv("./output/corr_N_EFMWrep5.csv")
acN <- corr_N_EFMW %>% mutate(acc = result_N_rep5.ac) %>% select(-result_N_rep5.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT5_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT5_CV2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

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
  common_theme +
  labs(title = "Narea",#(Training EF & validating on MW)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/Narea_compsubsetfig3_efmw.jpeg", plot = plot, width = 10, height = 8, dpi = 300)
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
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acST<- acST %>%  mutate(MT= "S-0",Type= "CV2")

# Combine all data
S_data_whole <- bind_rows(acS,acST_CV1,acST)
S_data_whole$Data <- "Complete"


#Subset 1
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep1.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep1.ac) %>% select(-result_S_rep1.ac)
acST_CV1<- read.csv("./output/acST1_CV1_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for sla SW1, SW2, SW3
acST<- read.csv("./output/acST1_CV2_rep1.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

# CombiSe all data
S_data_rep1 <- bind_rows(acS,acST_CV1,acST)
S_data_rep1$Data <- "Subset-1"

#Subset 2
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep2.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep2.ac) %>% select(-result_S_rep2.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST2_CV1_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("./output/acST2_CV2_rep2.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
S_data_rep2 <- bind_rows(acS,acST_CV1,acST)
S_data_rep2$Data <- "Subset-2"
#Subset 3
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep3.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep3.ac) %>% select(-result_S_rep3.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST3_CV1_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("./output/acST3_CV2_rep3.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
S_data_rep3 <- bind_rows(acS,acST, acST_CV1)
S_data_rep3$Data <- "Subset-3"

#Subset 4
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep4.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep4.ac) %>% select(-result_S_rep4.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST4_CV1_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("./output/acST4_CV2_rep4.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
S_data_rep4 <- bind_rows(acS,acST_CV1,acST)
S_data_rep4$Data <- "Subset-4"

#Subset 5
#ST for narea
corr_S_EFMW<- read.csv("./output/corr_S_EFMWrep5.csv")
acS <- corr_S_EFMW %>% mutate(acc = result_S_rep5.ac) %>% select(-result_S_rep5.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST5_CV1_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("./output/acST5_CV2_rep5.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

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
  common_theme +
  labs(title = "SLA",#(Training EF & validating on MW)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right")
#ggsave("./figure/sla_compsubsetfig3_efmw.jpeg", plot = plot, width = 10, height = 8, dpi = 300)


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


#MWEF plot
#plsr-narea
#ST for plsr-narea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEF.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn.ac) %>% select(-result_pn.ac)
acpnT_CV1<- read.csv("./output/acpnT_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
acpnT<- read.csv("./output/acpnT_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apnd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
pn_data_mwef_whole <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_mwef_whole$Data<- "Complete"

#replication 1
#ST for pnarea
corr_pn_MWEF<- read.csv("./output/corr_pn_MWEFrep1.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep1.ac) %>% select(-result_pn_rep1.ac)
#Multi trait cv1 for plsr-narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("./output/acpnT1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for plsr-narea pnW1, pnW2, pnW3
acpnT<- read.csv("./output/acpnT1_CV2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combinne all data
pn_data_rep1 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_pn_MWEF<- read.csv("./output/./output/corr_pn_MWEFrep2.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep2.ac) %>% select(-result_pn_rep2.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("./output/acpnT2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnT<- read.csv("./output/acpnT2_CV2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
pn_data_rep2 <-bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_pn_MWEF<- read.csv("./output/./output/corr_pn_MWEFrep3.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep3.ac) %>% select(-result_pn_rep3.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("./output/acpnT3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnT<- read.csv("./output/acpnT3_CV2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
pn_data_rep3 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_pn_MWEF<- read.csv("./output/./output/corr_pn_MWEFrep4.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep4.ac) %>% select(-result_pn_rep4.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("./output/acpnT4_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnT<- read.csv("./output/acpnT4_CV2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
pn_data_rep4 <- bind_rows(acpn,acpnT_CV1,acpnT)
pn_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_pn_MWEF<- read.csv("./output/./output/corr_pn_MWEFrep5.csv")
acpn <- corr_pn_MWEF %>% mutate(acc = result_pn_rep5.ac) %>% select(-result_pn_rep5.ac)
#Multi trait cv1 for narea pnW1, pnW2, pnW3
acpnT_CV1<- read.csv("./output/acpnT5_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
#multitrait Cv2 for narea pnW1, pnW2, pnW3
acpnT<- read.csv("./output/acpnT5_CV2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)
# Add labels for multi-trait and type
acpn <- acpn %>% mutate(MT = "ST", Type = "ST")
acpnT_CV1 <- acpnT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpnT <- acpnT %>% mutate(MT = "S-0", Type = "CV2")

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

#PLSR_SLA mwef
corr_ps_MWEF<- read.csv("./output/./output/corr_ps_MWEF.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps.ac) %>% select(-result_ps.ac)
#Multi trait cv1 for plsr_sla NW1, NW2, NW3
acpsT_CV1<- read.csv("./output/./output/acpsT_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for psarea psW1, psW2, psW3
acpsT<- read.csv("./output/./output/acpsT_CV1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait apsd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
ps_data_mwef_whole <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_mwef_whole$Data<- "Complete"

#replication 1
#ST for pnarea
corr_ps_MWEF<- read.csv("./output/./output/corr_ps_MWEFrep1.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep1.ac) %>% select(-result_ps_rep1.ac)
#Multi trait cv1 for plsr-sla psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for plsr-sla psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT1_CV2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait aSd type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combinne all data
ps_data_rep1 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep1$Data <- "Subset-1"
#replication 2
#ST for narea
corr_ps_MWEF<- read.csv("./output/./output/corr_ps_MWEFrep2.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep2.ac) %>% select(-result_ps_rep2.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT2_CV2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
ps_data_rep2 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep2$Data <- "Subset-2"

#replication 3
#ST for narea
corr_ps_MWEF<- read.csv("./output/./output/corr_ps_MWEFrep3.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep3.ac) %>% select(-result_ps_rep3.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT3_CV2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
ps_data_rep3 <- bind_rows(acps, acpsT_CV1,acpsT)
ps_data_rep3$Data <- "Subset-3"


#replication 4
#ST for narea
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEFrep4.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep4.ac) %>% select(-result_ps_rep4.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT4_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT4_CV2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")

acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
ps_data_rep4 <- bind_rows(acps,acpsT_CV1,acpsT)
ps_data_rep4$Data <- "Subset-4"

#replication 5
#ST for narea
corr_ps_MWEF<- read.csv("./output/corr_ps_MWEFrep5.csv")
acps <- corr_ps_MWEF %>% mutate(acc = result_ps_rep5.ac) %>% select(-result_ps_rep5.ac)
#Multi trait cv1 for narea psW1, psW2, psW3
acpsT_CV1<- read.csv("./output/acpsT5_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea psW1, psW2, psW3
acpsT<- read.csv("./output/acpsT5_CV2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acps <- acps %>% mutate(MT = "ST", Type = "ST")
acpsT_CV1 <- acpsT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acpsT <- acpsT %>% mutate(MT = "S-0", Type = "CV2")

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
  common_theme +
  labs(title = "PLSR-SLA",#(Training MW & validating on EF)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/ps_compsubsetfig3_mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)


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
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")

acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_mwef_whole <- bind_rows(acN,acNT_CV1,acNT)
N_data_mwef_whole$Data <- "Complete"
#replication 1 mwef
#ST for narea
corr_N_MWEF<- read.csv("./output/corr_N_MWEFrep1.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep1.ac) %>% select(-result_N_rep1.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT1_CV2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_rep1 <- bind_rows(acN, acNT, acNT_CV1)
N_data_rep1$Data <- "Subset-1"
#Subset 2
#ST for narea
corr_N_MWEF<- read.csv("./output/corr_N_MWEFrep2.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep2.ac) %>% select(-result_N_rep2.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT2_CV2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_rep2 <- bind_rows(acN,acNT, acNT_CV1)
N_data_rep2$Data <- "Subset-2"
#Subset 3
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep3.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep3.ac) %>% select(-result_N_rep3.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT3_CV2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_rep3 <- bind_rows(acN,acNT,acNT_CV1)
N_data_rep3$Data <- "Subset-3"

#Subset 4
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep4.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep4.ac) %>% select(-result_N_rep4.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT4_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT4_CV2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")

acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
N_data_rep4 <- bind_rows(acN,acNT, acNT_CV1)
N_data_rep4$Data <- "Subset-4"
#Subset 5
#ST for narea
corr_N_MWEF<- read.csv("./output/narea/corr_N_MWEFrep5.csv")
acN <- corr_N_MWEF %>% mutate(acc = result_N_rep5.ac) %>% select(-result_N_rep5.ac)
#Multi trait cv1 for narea NW1, NW2, NW3
acNT_CV1<- read.csv("./output/acNT5_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea NW1, NW2, NW3
acNT<- read.csv("./output/acNT5_CV2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acN <- acN %>% mutate(MT = "ST", Type = "ST")
acNT_CV1 <- acNT_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acNT <- acNT %>% mutate(MT = "S-0", Type = "CV2")

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
  common_theme +
  labs(title = "Narea",#(Training MW & validating on EF)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position = "right") # Bold facet labels for clarity

#ggsave("./figure/Narea_compsubsetfig3_mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)

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
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")

acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
S_data_mwef_whole <- bind_rows(acS,acST,acST_CV1)
S_data_mwef_whole$Data <- "Complete"

#rep-1
#ST for sla
corr_S_MWEF<- read.csv("./output/corr_S_MWEFrep1.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep1.ac) %>% select(-result_S_rep1.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST1_CV1_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("./output/acST1_CV2_rep1_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")

acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
S_data_rep1 <- bind_rows(acS,acST_CV1,acST)
S_data_rep1$Data <- "Subset-1"

#Subset 2
#ST for narea
corr_S_MWEF<- read.csv("./output/corr_S_MWEFrep2.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep2.ac) %>% select(-result_S_rep2.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST2_CV1_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("./output/acST2_CV2_rep2_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")

acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
S_data_rep2 <- bind_rows(acS,acST_CV1,acST)
S_data_rep2$Data <- "Subset-2"


#Subset 3
#ST for sla
corr_S_MWEF<- read.csv("./output/sla/corr_S_MWEFrep3.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep3.ac) %>% select(-result_S_rep3.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST3_CV1_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("acST3_CV2_rep3_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")

acST <- acST %>% mutate(MT = "S-0", Type = "CV2")


# Combine all data
S_data_rep3 <-bind_rows(acS,acST_CV1,acST)
S_data_rep3$Data <- "Subset-3"
#Subset 4
#ST for sla
corr_S_MWEF<- read.csv("./output/./output/sla/corr_S_MWEFrep4.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep4.ac) %>% select(-result_S_rep4.ac)
#Multi trait cv1 for narea SW1, SW2, SW3
acST_CV1<- read.csv("./output/acST4_CV1_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)


#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("./output/acST4_CV2_rep4_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")
acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

# Combine all data
S_data_rep4 <- bind_rows(acS,acST_CV1,acST)
S_data_rep4$Data <- "Subset-4"
#Subset 5
#ST for narea
corr_S_MWEF<- read.csv("./output/./output/sla/corr_S_MWEFrep5.csv")
acS <- corr_S_MWEF %>% mutate(acc = result_S_rep5.ac) %>% select(-result_S_rep5.ac)
acST_CV1<- read.csv("./output/acST5_CV1_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

#multitrait Cv2 for narea SW1, SW2, SW3
acST<- read.csv("./output/acST5_CV2_rep5_mwef.txt",sep=",",header= FALSE)%>% mutate(acc = V1) %>% select(-V1)

# Add labels for multi-trait and type
acS <- acS %>% mutate(MT = "ST", Type = "ST")
acST_CV1 <- acST_CV1 %>% mutate(MT = "S-0", Type = "CV1")

acST <- acST %>% mutate(MT = "S-0", Type = "CV2")

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
  common_theme +
  labs(title = "SLA",#(Training MW & validating on EF)",
       x = "Model Type",
       y = "Accuracy") +
  facet_wrap(~ Data, ncol = 6, scales = "free_x") + # Facet and allow free x-axis scaling
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add panel borders
        strip.background = element_rect(fill = "lightgray"), # Differentiate the facets
        strip.text = element_text(face = "bold"),
        legend.position="right") # Bold facet labels for clarity

#ggsave("./figure/sla_compsubsetfig3_mwef.jpeg", plot = plot, width = 10, height = 8, dpi = 300)
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


