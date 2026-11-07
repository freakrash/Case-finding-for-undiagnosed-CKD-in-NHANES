# Figure 1: decision curve (net benefit vs. threshold), C1-unflagged headline stratum, validation.
# Figure 2: subgroup calibration (O/E ratio by race/ethnicity), MICE-based, with reference line at 1.

source("R/00_config.R")
library(ggplot2)

dca <- read.csv(file.path(TABLES_DIR, "decision_curve.csv"))
dca_long <- reshape(dca, direction = "long", varying = c("model","test_all","test_none","c1","c2"),
                     v.names = "net_benefit", times = c("Model","Test-all","Test-none","Rule C1","Rule C2"),
                     timevar = "strategy")
dca_long <- dca_long[dca_long$strategy != "Rule C1", ]  # C1 == test-none in this stratum, redundant on the plot

p1 <- ggplot(dca_long, aes(x = threshold, y = net_benefit, color = strategy, linetype = strategy)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = DCA_PRIMARY, linetype = "dotted", color = "grey40") +
  annotate("text", x = DCA_PRIMARY + 0.004, y = max(dca_long$net_benefit) * 0.95,
           label = "primary\nthreshold (5%)", hjust = 0, size = 3, color = "grey40") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Model" = "#1b5e8a", "Test-all" = "#c0392b",
                                 "Test-none" = "grey50", "Rule C2" = "#e29b1a")) +
  scale_linetype_manual(values = c("Model" = "solid", "Test-all" = "solid",
                                    "Test-none" = "dashed", "Rule C2" = "solid")) +
  labs(x = "Decision threshold (probability of undiagnosed CKD)", y = "Net benefit",
       title = "Decision-curve analysis: model-directed testing vs. alternatives",
       subtitle = "C1-unflagged headline stratum, temporal external validation (complete-case)",
       color = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(FIGURES_DIR, "figure1_decision_curve.png"), p1, width = 7, height = 5, dpi = 300)
cat("[14] Figure 1 saved.\n")

race_mi <- read.csv(file.path(TABLES_DIR, "subgroup_oe_race_mi.csv"))
race_mi$level <- factor(race_mi$level, levels = race_mi$level[order(race_mi$oe)])

p2 <- ggplot(race_mi, aes(x = oe, y = level)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = events_unw), color = "#1b5e8a") +
  geom_text(aes(label = sprintf("n=%d, %d events", n_unw, events_unw)), hjust = -0.15, size = 3, color = "grey30") +
  scale_size_continuous(range = c(2, 7), guide = "none") +
  scale_x_continuous(limits = c(0, max(race_mi$oe) * 1.6)) +
  labs(x = "Observed / Expected ratio (>1 = model under-predicts risk)", y = NULL,
       title = "Subgroup calibration by race/ethnicity",
       subtitle = "C1-unflagged headline stratum, validation, MICE-based predictions (primary)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(FIGURES_DIR, "figure2_subgroup_calibration.png"), p2, width = 7.5, height = 4.5, dpi = 300)
cat("[14] Figure 2 saved.\n")
