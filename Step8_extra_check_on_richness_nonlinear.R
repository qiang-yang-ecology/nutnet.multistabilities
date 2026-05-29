# ---- Nonlinear directional richness response vs enrichment intensity (quadratic) ----

library(lme4)
library(broom.mixed)
library(ggplot2)
library(dplyr)
library(ggeffects)
library(gridExtra)
library(cowplot)
library(grid)
library(gtable)
library(patchwork)

## ----------------------------------------------------------
## 3. Quadratic mixed models for N and P (site_code as random effect)
## ----------------------------------------------------------

# N model: lrr ~ (Ri_sc + Ri_sc^2) * year_trt + (1 | site_code)
m_LRR_N_quad <- lmer(
  lrr ~ (Ri_sc + I(Ri_sc^2)) * year_trt + (1 | site_code),
  data = df.lrr.dir.N
)

# P model: lrr ~ (Ri_sc + Ri_sc^2) * year_trt + (1 | site_code)
m_LRR_P_quad <- lmer(
  lrr ~ (Ri_sc + I(Ri_sc^2)) * year_trt + (1 | site_code),
  data = df.lrr.dir.P
)

# Extract fixed-effect tables
result.N.quad <- broom.mixed::tidy(m_LRR_N_quad, effects = "fixed")
result.P.quad <- broom.mixed::tidy(m_LRR_P_quad, effects = "fixed")

## ----------------------------------------------------------
## 4. Save coefficient tables of nonlinear models into one PDF page
## ----------------------------------------------------------

# Clean up table appearance
tbl_N_quad <- result.N.quad %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

tbl_P_quad <- result.P.quad %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Convert to table grobs
g_N_quad <- tableGrob(tbl_N_quad, rows = NULL)
g_P_quad <- tableGrob(tbl_P_quad, rows = NULL)

# Titles
title_N_quad <- textGrob("Nonlinear model results for nitrogen (N)",
                         gp = gpar(fontsize = 14, fontface = "bold"))
title_P_quad <- textGrob("Nonlinear model results for phosphorus (P)",
                         gp = gpar(fontsize = 14, fontface = "bold"))

# Attach titles
panel_N_quad <- gtable::gtable_add_rows(
  g_N_quad,
  heights = grobHeight(title_N_quad) + unit(4, "pt"),
  pos = 0
)
panel_N_quad <- gtable::gtable_add_grob(
  panel_N_quad, title_N_quad, 1, 1, 1, ncol(panel_N_quad)
)

panel_P_quad <- gtable::gtable_add_rows(
  g_P_quad,
  heights = grobHeight(title_P_quad) + unit(4, "pt"),
  pos = 0
)
panel_P_quad <- gtable::gtable_add_grob(
  panel_P_quad, title_P_quad, 1, 1, 1, ncol(panel_P_quad)
)

# Combine vertically
final_tables_quad <- grid.arrange(panel_N_quad, panel_P_quad,
                                  ncol = 1, heights = c(1, 1))

# Save PDF
pdf(file = paste0(getwd(), "/myfigure/LRR_nonlin_model_tables_N_and_P_richness.pdf"),
    width = 7, height = 7)
grid.draw(final_tables_quad)
dev.off()

## ----------------------------------------------------------
## 5. Tipping-curve prediction figure:
##    x-axis = scaled log intensity (Ri_sc), colours = Year 1/5/10
## ----------------------------------------------------------

# Choose representative years
years_sel <- c(1, 5, 10)

years_sel_N <- years_sel[years_sel %in% unique(df.lrr.dir.N$year_trt)]
years_sel_P <- years_sel[years_sel %in% unique(df.lrr.dir.P$year_trt)]

year_term_N <- paste0("year_trt [", paste(years_sel_N, collapse = ","), "]")
year_term_P <- paste0("year_trt [", paste(years_sel_P, collapse = ","), "]")

# Predictions along Ri_sc for N
pred_N_tip <- ggpredict(
  m_LRR_N_quad,
  terms = c("Ri_sc [all]", year_term_N)
) %>%
  as_tibble() %>%
  mutate(
    Time = factor(
      group,
      levels = as.character(years_sel_N),
      labels = paste("Year", years_sel_N)
    )
  )

# Predictions along Ri_sc for P
pred_P_tip <- ggpredict(
  m_LRR_P_quad,
  terms = c("Ri_sc [all]", year_term_P)
) %>%
  as_tibble() %>%
  mutate(
    Time = factor(
      group,
      levels = as.character(years_sel_P),
      labels = paste("Year", years_sel_P)
    )
  )

# Colours for years
col_vec_time <- c(
  "Year 1"  = "#18A1D1",  # red
  "Year 5"  = "#96C594",  # green
  "Year 10" = "#E361BE"   # blue
)
fill_vec_time <- col_vec_time

# Plot function: LRR in richness vs Ri_sc (scaled log intensity)
plot_tipping_richness <- function(pred_data, nutrient_label) {
  ggplot(pred_data,
         aes(x = x,                  # Ri_sc
             y = predicted,
             colour = Time,
             group = Time)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_line(size = 1.3) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Time),
                alpha = 0.18, colour = NA) +
    scale_x_continuous(
      name   = paste0("Scaled log ", nutrient_label, " enrichment intensity"),
      breaks = pretty(range(pred_data$x), n = 5)
    ) +
    scale_colour_manual(values = col_vec_time) +
    scale_fill_manual(values  = fill_vec_time) +
    labs(
      y      = "Predicted LRR in richness",
      colour = "Time since treatment",
      fill   = "Time since treatment"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text       = element_text(colour = "black"),
      legend.position = "bottom",
      plot.title      = element_blank()
    )
}

p_N_rich_tipping <- plot_tipping_richness(pred_N_tip, "N")
p_P_rich_tipping <- plot_tipping_richness(pred_P_tip, "P")

# Combine N and P panels with shared legend
p_rich_tipping_both <-
  (p_N_rich_tipping | p_P_rich_tipping) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.tag        = element_text(face = "bold")
  )

p_rich_tipping_both

# Save to PDF
pdf(
  file = paste0(getwd(),
                "/myfigure/Nonlinear_tipping_curves_richness_LRR_vs_scaled_log_intensity_N_and_P.pdf"),
  width = 8.5,
  height = 5
)
print(p_rich_tipping_both)
dev.off()