#!/usr/bin/env Rscript
# Test how close assigned properties to true values
library(tidyverse)
library(devtools)
load_all()

# Import dms data
s_dms <- read_csv("analysis/starr_ace2_spike.csv") %>%
  select(position = site_SARS2, wt=wildtype, mut=mutant, score=expr_avg) %>%
  deep_mutational_scan("Starr S", annotate = TRUE, study = "Starr et al. 2020", gene = "S")

# Import and summarise measured results
columns <- cols(uniprot = col_character(), name = col_character(), position = col_double(), wt = col_character(),
  mut = col_character(), sift_score = col_double(), template = col_character(),
  relative_surface_accessibility = col_double(), foldx_ddg = col_double(), ptm = col_character(),
  int_uniprot = col_character(), int_name = col_character(), int_template = col_character(),
  interaction_energy = col_double(), diff_interaction_energy = col_double(), diff_interface_residues = col_integer(),
  freq = col_double(), mut_escape_mean = col_double(), mut_escape_max = col_double(), annotation = col_character()
)
s_measure <- read_tsv("analysis/sars_mutfunc.tsv", col_types = columns) %>%
  filter(name == "s") %>%
  group_by(position, wt) %>%
  summarise(mean_sift = mean(log10(sift_score + 0.00001), na.rm = TRUE),
            total_energy = clamp(mean(foldx_ddg, na.rm = TRUE), -10, 10),
            all_atom_rel = mean(relative_surface_accessibility, na.rm = TRUE),
            .groups = "drop")

# Calculate hexagonal binned landscape summary
hex_pred <- select(landscape_properties(s_dms, bins = 20), position, wt, mean_sift, total_energy, all_atom_rel)

# Combine predictions and plot
s_pred <- describe_clusters(s_dms, full = TRUE) %>%
  select(position, wt, mean_sift, total_energy, all_atom_rel) %>%
  bind_rows(pred_subtype = ., pred_hex = hex_pred, measure = s_measure, .id = "type") %>%
  pivot_longer(c(mean_sift, total_energy, all_atom_rel)) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  pivot_longer(starts_with("pred"), names_to = "method", values_to = "pred", names_prefix = "pred_")

lms <- filter(s_pred, !is.na(measure), !is.infinite(measure), !is.na(pred), !is.infinite(pred)) %>%
  group_by(name, method) %>%
  group_modify(~broom::glance(lm(measure~pred, data = .x)))

p <- ggplot(s_pred, aes(x = measure, y = pred, colour = method)) +
  facet_wrap(~name, nrow = 1, scales = "free") +
  geom_point(shape = 20) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggpubr::theme_pubclean()
