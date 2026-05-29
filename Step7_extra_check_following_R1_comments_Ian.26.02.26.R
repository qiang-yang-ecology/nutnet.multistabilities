# ---- Directional biomass response vs enrichment intensity (Ian's check) ----
rm(list = ls())
gc()
# necessary R packages
Packages <- c(
  "tidyverse", "magrittr", "grid", "gridExtra", "rnaturalearth", "rnaturalearthdata",
  "ggcorrplot", "scales", "sp", "sf", "broom", "nlme", "scico", "vegan", "dlookr", "GGally",
  "PerformanceAnalytics"
)
pacman::p_load(char = Packages)
filter <- dplyr::filter
# functions with conflicting status (presence in different packages used here)
select <- dplyr::select
summarise <- dplyr::summarise
# working path
main.path <- getwd()
data.path <- paste0(main.path, "/mydata")
figure.path <- paste0(main.path, "/myfigure")
set.seed(1234)

load(file = paste0(data.path, "/df.tidy.stabilities.Rdata"))
load(file = paste0(data.path, "/df.tidy.Rdata"))

## ----------------------------------------------------------
## 1. Recreate df.tidy.mass and control selection (as in Step2 Rmd)
## ----------------------------------------------------------

# df.tidy must already be loaded from df.tidy.Rdata
# Only plots with biomass data
df.tidy.mass <- df.tidy %>%
  filter(data.availability != "only coverage data") %>%
  select(site_code:trt, biomass.data) %>%
  rename(data = biomass.data)

# All potential control plots (unnested biomass)
df.controls <- df.tidy.mass %>%
  unnest(cols = data) %>%
  ungroup() %>%
  filter(trt == "Control") %>%
  select(site_code, block, plot, year_trt) %>%
  distinct()

set.seed(1234)
df.controls.selected <- df.controls %>%
  arrange(site_code, block, desc(year_trt)) %>%
  group_by(site_code, block) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-year_trt)

# Attach control biomass to each block
df.tidy.mass.control <- df.tidy.mass %>%
  left_join(df.controls.selected, .)   # same as in Rmd

df.tidy.mass.control <- df.tidy.mass.control %>%
  select(site_code, block, plot, data) %>%
  rename(control.plot = plot, control.data = data)

df.tidy.mass.2 <- df.tidy.mass %>%
  left_join(df.tidy.mass.control) %>%
  select(site_code, block, plot, control.plot, trt, data, control.data)

## ----------------------------------------------------------
## 2. Compute mean directional LRR per plot for +N and +P
## ----------------------------------------------------------

# the total biomass of the control and the treatment N, P and their log response ratio (LRR) for each year after treatment (year_trt >= 1).
df.lrr.dir <- df.tidy.mass.2 %>%
  filter(trt %in% c("N", "P")) %>%
  mutate(
    results = map2(
      .x = data,
      .y = control.data,
      .f = function(dat1, dat2) {
        # Combine treatment and control; keep overlapping years and year_trt >= 1
        dat <- dat1 %>%
          mutate(g = "experimental") %>%
          bind_rows(dat2 %>% mutate(g = "control")) %>%
          select(year_trt, g, mass) %>%
          filter(
            year_trt %in% intersect(dat1$year_trt, dat2$year_trt),
            year_trt >= 1
          )
        
        # Total biomass per year per group
        dat.nested <- dat %>%
          group_by(year_trt, g) %>%
          summarise(total_mass = sum(mass), .groups = "drop")
        
        # Log response ratio per year
        lRR.total.mass <- dat.nested %>%
          pivot_wider(names_from = g, values_from = total_mass) %>%
          mutate(lrr = log(experimental / control))
        
        return(lRR.total.mass)
      }
    )
  ) %>%
  select(site_code, block, plot, trt, results) %>%
  unnest(cols = results)

df.lrr.dir <- na.omit(df.lrr.dir) # removed one record
df.lrr.dir <- df.lrr.dir %>%
  filter(is.finite(lrr)) # removed records with infinite LRR (zero control biomass)

# now add nutrient addition rates (N_added_proportion, P_added_proportion) from df.tidy for each plot, to calculate enrichment intensity (Ri) later.
df.lrr.dir.N <- df.lrr.dir %>%
  filter(trt == "N") %>%
  left_join(
    df.tidy %>%
      select(site_code, block, plot, trt, N_added_proportion) %>%
      distinct(),
    by = c("site_code", "block", "plot", "trt")
  ) %>% 
  mutate(
    Ri_raw = N_added_proportion,
    Ri     = log(Ri_raw),
    Ri_sc  = scale(Ri)[, 1]
  ) 

m_LRR_N <- lmer(
  lrr ~ Ri_sc * year_trt + (1 | site_code),
  data = df.lrr.dir.N
)

df.lrr.dir.P <- df.lrr.dir %>%
  filter(trt == "P") %>%
  left_join(
    df.tidy %>%
      select(site_code, block, plot, trt, P_added_proportion) %>%
      distinct(),
    by = c("site_code", "block", "plot", "trt")
  ) %>% 
  mutate(
    Ri_raw = P_added_proportion,
    Ri     = log(Ri_raw),
    Ri_sc  = scale(Ri)[, 1]
  ) 

m_LRR_P <- lmer(
  lrr ~ Ri_sc * year_trt + (1 | site_code),
  data = df.lrr.dir.P
)

result.N <- broom.mixed::tidy(m_LRR_N, effects = "fixed")
result.P <- broom.mixed::tidy(m_LRR_P, effects = "fixed")

# save as pdf
library(gridExtra)
library(cowplot)

# Clean up table appearance a bit (optional)
tbl_N <- result.N %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

tbl_P <- result.P %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

# Convert to table grobs
g_N <- tableGrob(tbl_N, rows = NULL)
g_P <- tableGrob(tbl_P, rows = NULL)

# Add titles above tables
title_N <- textGrob("Model results for nitrogen (N)", gp = gpar(fontsize = 14, fontface = "bold"))
title_P <- textGrob("Model results for phosphorus (P)", gp = gpar(fontsize = 14, fontface = "bold"))

# Combine table + title for each
panel_N <- gtable::gtable_add_rows(g_N, heights = grobHeight(title_N) + unit(4, "pt"), pos = 0)
panel_N <- gtable::gtable_add_grob(panel_N, title_N, 1, 1, 1, ncol(panel_N))

panel_P <- gtable::gtable_add_rows(g_P, heights = grobHeight(title_P) + unit(4, "pt"), pos = 0)
panel_P <- gtable::gtable_add_grob(panel_P, title_P, 1, 1, 1, ncol(panel_P))

# Combine the two panels vertically
final_tables <- grid.arrange(panel_N, panel_P, ncol = 1, heights = c(1, 1))

pdf(paste0(getwd(), "/myfigure/LRR_model_tables_N_and_P.pdf.pdf"), width = 7, height = 7)
grid.draw(final_tables)
dev.off()


# ------------------------------------------------------------------------
# make prediction figure
# ------------------------------------------------------------------------
library(ggeffects)
library(patchwork)
pred_N <- ggpredict(
  m_LRR_N,
  terms = c("year_trt [all]", "Ri_sc [-1, 0, 1]")
) %>%
  as_tibble() %>%
  mutate(
    Intensity = case_when(
      group == "-1" ~ "Low intensity",
      group == "0"  ~ "Medium intensity",
      group == "1"  ~ "High intensity"
    ),
    Intensity = factor(
      Intensity,
      levels = c("Low intensity", "Medium intensity", "High intensity")
    )
  )

pred_P <- ggpredict(
  m_LRR_P,
  terms = c("year_trt [all]", "Ri_sc [-1, 0, 1]")
) %>%
  as_tibble() %>%
  mutate(
    Intensity = case_when(
      group == "-1" ~ "Low intensity",
      group == "0"  ~ "Medium intensity",
      group == "1"  ~ "High intensity"
    ),
    Intensity = factor(
      Intensity,
      levels = c("Low intensity", "Medium intensity", "High intensity")
    )
  )

col_vec  <- c(
  "Low intensity"    = "#FDE725FF",
  "Medium intensity" = "#21908CFF",
  "High intensity"   = "#440154FF"
)
fill_vec <- col_vec

plot_pred <- function(pred_data) {
  ggplot(pred_data,
         aes(x = x, y = predicted,
             colour = Intensity, group = Intensity)) +
    geom_line(size = 1.25) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Intensity),
                alpha = 0.18, colour = NA) +
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15)) +
    scale_colour_manual(values = col_vec) +
    scale_fill_manual(values = fill_vec) +
    labs(
      x = "Year",
      y = "Predicted LRR in total biomass",
      title = NULL,
      colour = "Intensity",
      fill   = "Intensity"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text       = element_text(color = "black"),
      legend.position = "bottom",
      plot.title      = element_blank()
    )
}

# Plots WITH legends (do NOT remove legend here)
p_N_with_legend <- plot_pred(pred_N)
p_P_with_legend <- plot_pred(pred_P)
p_lrr_pred <- p_N_with_legend +p_P_with_legend+
  plot_annotation(
    tag_levels = 'a' # This adds A, B, C...
  )+
  plot_layout(guides = 'collect') & theme(plot.tag = element_text(face = 'bold'),legend.position = 'bottom')

pdf(paste0(getwd(), "/myfigure/Predicted effects of nutrient-addition intensity and time on directional biomass responses (LRR) under N and P enrichment.pdf"), width = 7, height = 5)
print(p_lrr_pred)
dev.off()

