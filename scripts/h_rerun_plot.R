library(tidyverse)
library(cowplot)
library(scales)

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
       plot = final_plot2, width = 7, height = 10, dpi = 300)

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

# =========================
# 2) PUT YOUR NEW VALUES HERE
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


