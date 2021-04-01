#!/usr/bin/env Rscript
# Test how close assigned properties to true values
library(tidyverse)
library(devtools)
library(hexbin)
load_all()

### Spike predictions ###
# Import data
s_dms <- read_csv("analysis/starr_ace2_spike.csv") %>%
  select(position = site_SARS2, wt=wildtype, mut=mutant, score=expr_avg) %>%
  deep_mutational_scan("Starr S", annotate = TRUE, study = "Starr et al. 2020", gene = "S")

s_mutfunc <- read_tsv("analysis/sars_mutfunc.tsv") %>%
  mutate(sift_score = log10(sift_score + 0.0001)) %>%
  mutate(across(total_energy:entropy_complex, clamp, lower = -10, upper = 10)) %>%
  group_by(position, wt) %>%
  summarise(across(-mut, mean, na.rm = TRUE), .groups = "drop") %>%
  rename(mean_sift = sift_score, all_atom_rel = relative_surface_accessibility)

pred_aa <- group_by(deep_landscape, wt) %>%
  summarise(across(c(mean_sift:energy_ionisation, all_atom_rel), mean, na.rm = TRUE)) %>%
  pivot_longer(-wt, values_to = "pred_wt")

# Combine predictions and plot
s_pred <- bind_rows(pred_subtype = describe_clusters(s_dms, full = TRUE),
                    pred_hex = landscape_properties(s_dms, bins = 20),
                    measure = s_mutfunc, .id = "type") %>%
  select(-name, -cluster, -global_cluster_freq, - group, -description, -notes, -one_of(amino_acids), -mean_score) %>%
  pivot_longer(c(-type, -position, -wt)) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  left_join(pred_aa, by = c("wt", "name")) %>%
  pivot_longer(starts_with("pred"), names_to = "method", values_to = "pred", names_prefix = "pred_")

lms_spike <- filter(s_pred, !is.na(measure), !is.infinite(measure), !is.na(pred), !is.infinite(pred)) %>%
  group_by(name, method) %>%
  group_modify(~broom::glance(lm(measure~pred, data = .x)))

p_spike <- ggplot(s_pred, aes(x = measure, y = pred, colour = method)) +
  facet_wrap(~name, nrow = 5, scales = "free") +
  geom_point(shape = 20, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggpubr::theme_pubclean()

### Whole dataset predictions
dl <- mutate(deep_landscape,
             hex = hexbin(deepscanscape::deep_landscape$umap1, deepscanscape::deep_landscape$umap2, xbins = 20, IDs = TRUE)@cID)

pred_subtype <- group_by(deep_landscape, cluster) %>%
  summarise(across(c(mean_sift:energy_ionisation, all_atom_rel), mean, na.rm = TRUE)) %>%
  pivot_longer(-cluster, values_to = "pred_subtype")

pred_hex <- group_by(dl, hex) %>%
  summarise(across(c(mean_sift:energy_ionisation, all_atom_rel), mean, na.rm = TRUE)) %>%
  pivot_longer(-hex, values_to = "pred_hex")

dl_preds <- select(dl, cluster, wt, hex, mean_sift:energy_ionisation, all_atom_rel) %>%
  pivot_longer(c(-cluster, -wt, -hex), values_to = "measure") %>%
  left_join(pred_aa, by = c("wt", "name")) %>%
  left_join(pred_hex, by = c("hex", "name")) %>%
  left_join(pred_subtype, by = c("cluster", "name")) %>%
  pivot_longer(starts_with("pred"), names_to = "method", values_to = "pred", names_prefix = "pred_") %>%
  drop_na()

lms_dl <- filter(dl_preds, !is.na(measure), !is.infinite(measure), !is.na(pred), !is.infinite(pred)) %>%
  group_by(name, method) %>%
  group_modify(~broom::glance(lm(measure~pred, data = .x)))

p_dl <- ggplot(dl_preds, aes(x = measure, y = pred, colour = method)) +
  facet_wrap(~name, nrow = 5, scales = "free") +
  geom_point(shape = 20, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggpubr::theme_pubclean()
