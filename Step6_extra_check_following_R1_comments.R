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

# ---- first check the three outliers on Fig.S4 ----

pdf(paste0(getwd(), "/myfigure/correlation.two.temporal.functional.variability.mass.pdf"), width = 5.5, height = 5)
# Fit a linear model
model <- lm(func_var_det ~ func_var, data = df.tidy.stabilities)

# Extract R-squared and p-value
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]

# Create the plot
plot(
  df.tidy.stabilities$func_var, 
  df.tidy.stabilities$func_var_det, 
  pch = 1, col = "gray",
  xlab = "Functional Variability", 
  ylab = "Detrended Functional Variability", 
)

# Add regression line
abline(model, col = "black", lwd = 2)

# Add R-squared and p-value in the upper left corner
text(
  x = -1.5, 
  y = 3,
  labels = paste0("R² = ", round(r_squared, 3), "\n", "P = ", format.pval(p_value, digits = 3)),
  adj = c(0, 1),  # Align text to top-left corner
  cex = 1  # Slightly increase text size
)

dev.off()

plot(
  df.tidy.stabilities$func_var, 
  df.tidy.stabilities$func_var_det, 
  pch = 1, col = "blue",
  xlab = "Functional Variability", 
  ylab = "Detrended Functional Variability", 
)

# Add regression line
abline(model, col = "red", lwd = 2)

# Add R-squared and p-value in the upper left corner
text(
  x = -1, 
  y = 2.5,
  labels = paste0("R² = ", round(r_squared, 3), "\n", "P = ", format.pval(p_value, digits = 3)),
  adj = c(0, 1),  # Align text to top-left corner
  cex = 1  # Slightly increase text size
)
dev.off()

# ---- comments on legumes ----
#### check Sally's comments on the potential increases in the relative abundance of legumes under Phosphoros addition
# first select plots that have compostional stability values
df_plot_selected <- df.tidy.stabilities %>%
  filter(!is.na(comp_change)) %>% 
  select(site_code, block, trt, plot) %>% 
  distinct()
# also need to add the coverage data for the control
df_plot_control <- df_plot_selected %>% 
  select(site_code, block) %>% 
  distinct() %>% 
  left_join(
    df.tidy %>% 
      filter(trt == "Control") %>% 
      select(site_code, block, plot, trt) %>% 
      distinct()
  )
df_plot_selected <- df_plot_selected %>% 
  bind_rows(df_plot_control)

df_plot_selected <- df.tidy %>% 
  select(site_code, block, plot, trt, coverage.data) %>% 
  semi_join(df_plot_selected)
# the species list
species_list <- df_plot_selected %>% 
  select(coverage.data) %>%
  unnest(cols = coverage.data) %>%
  select(Taxon) %>% 
  unlist() %>% 
  unique()
# using Taxonstand to find legumes
library(rgbif)
family_backbone <- sapply(species_list, function(x){
  out <- name_backbone(name = x)$family
  ifelse(is.null(out), "NA", out)
})
legume_list1 <- species_list[which(family_backbone == "Fabaceae")]

# there are a few species that could not be assigned family with the function name_backbone

# I manually checked them and assigned it to legume or not
df_manual_check <- read_csv(paste0(data.path, "/legumes_status_a_few_unclear_species.csv"))
legume_list2 <- df_manual_check %>% 
  filter(legume == "legume") %>% 
  pull(name)
legume_list <- union(legume_list1, legume_list2)
# then calculate the relative abundance of legumes in each plot at each year

cal_legume_rel <- function(dt) {
  dt %>% 
    mutate(
      is_legume = if_else(Taxon %in% legume_list, "legume", "non_legume")
    ) %>% 
    group_by(year_trt, is_legume) %>% 
    summarise(cover = sum(max_cover, na.rm = TRUE), .groups = "drop") %>% 
    
    # Make sure both legume and non_legume exist for every year_trt
    complete(
      year_trt,
      is_legume = c("legume", "non_legume"),
      fill = list(cover = 0)
    ) %>% 
    
    pivot_wider(
      names_from = is_legume,
      values_from = cover,
      values_fill = 0
    ) %>% 
    
    mutate(
      total_cover = legume + non_legume,
      rel_abundance_legume = if_else(
        total_cover > 0,
        legume / total_cover,
        NA_real_
      )
    ) %>% 
    select(year_trt, rel_abundance_legume)
}
df_legumes <- df_plot_selected %>% 
  select(site_code, block, plot, trt, coverage.data) %>% 
  semi_join(df_plot_selected) %>% 
  mutate(coverage.data = map(coverage.data, cal_legume_rel))
df_legumes <- df_legumes %>% 
  unnest(cols = coverage.data)

library(glmmTMB)
m <- glmmTMB(
  rel_abundance_legume ~ trt * year_trt +
    (1 | site_code),
  data = df_legumes,
  ziformula = ~ 1, # considering zero-inflation
  family = beta_family(link = "logit")
)
# now make the figure
# a coefficient figure
library(broom.mixed)
# Extract clean coefficient table
coefs <- broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE)

# Remove intercept (optional)
coefs2 <- coefs %>% 
  filter(term != "(Intercept)")

# Pretty names
coefs2$term <- recode(coefs2$term,
                      "trtK"           = "K vs Control",
                      "trtN"           = "N vs Control",
                      "trtP"           = "P vs Control",
                      "year_trt"       = "Year (Control)",
                      "trtK:year_trt"  = "K × Year",
                      "trtN:year_trt"  = "N × Year",
                      "trtP:year_trt"  = "P × Year"
)
# assign significance
coefs2 %<>%
  mutate(sig.sign = map2_chr(estimate, p.value, function(x, p) {
    if (x < 0 & p < 0.05) {
      out <- "negative"
    } else if (x > 0 & p < 0.05) {
      out <- "positive"
    } else {
      out <- "non-significant"
    }
    return(out)
  }))

p_legume_a <- ggplot(coefs2, aes(estimate, term)) +
  geom_point(size = 3, aes(colour = sig.sign)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig.sign), height = 0.15) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw(base_size = 14) +
  labs(x = "Coefficient (logit scale)", y = "",
       colour = "Significance: ")+
  theme(axis.text = element_text(size=12, color = "black"),
        legend.position = "bottom")
# a prediction figure
library(ggeffects)
# Generate marginal effects
pred <- ggpredict(m, terms = c("year_trt", "trt"))
pred <- as_tibble(pred)

pred <- pred %>%
  mutate(
    trt_lab = factor(
      group,
      levels = c("Control", "N", "P", "K"),        # original names
      labels = c("Control", "+N", "+P", "+K")      # what you want to show
    )
  )
library(RColorBrewer)
dark2_cols <- brewer.pal(3, "Dark2")
col_vec <- c(
  "Control" = "grey40",
  "+N"      = dark2_cols[1],
  "+P"      = dark2_cols[2],
  "+K"      = dark2_cols[3]
)

p_legume_b <- ggplot(pred) +
  # main fitted lines (4)
  geom_line(aes(x = x, y = predicted,
                colour = trt_lab, group = trt_lab),
            size = 1.5) +
  # lower CI (4 dashed)
  geom_line(aes(x = x, y = conf.low,
                colour = trt_lab, group = trt_lab),
            linetype = "dashed", size = 0.9) +
  # upper CI (4 dashed)
  geom_line(aes(x = x, y = conf.high,
                colour = trt_lab, group = trt_lab),
            linetype = "dashed", size = 0.9) +
  theme_bw(base_size = 14) +
  labs(
    x = "Year",
    y = "Predicted legume share",
    colour = "Treatment: "
  ) +
  scale_colour_manual(values = col_vec) +
  theme(
    plot.title = element_blank(),
    legend.position = "bottom",
    axis.text = element_text(size=12, color = "black")
  )

library(patchwork)
p_legume <- p_legume_a + p_legume_b + 
  plot_layout(ncol = 2, widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 16))

pdf(paste0(getwd(), "/myfigure/legume.rel.abundance.pdf"), width = 11.5, height = 6)
p_legume
dev.off()

# later we found that the legume relative abundance between treatments differ at year0
# here use a boxplot to check this reflect some baseline difference
df_legumes_year0 <- df_legumes %>% 
  filter(year_trt == 0)

m0 <- lmer(rel_abundance_legume ~ trt + (1 | site_code), data = df_legumes_year0)
pairs_df <- pairs(emmeans(m0, ~ trt)) |> 
  as.data.frame() |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))%>% 
  rename(
    Comparison = contrast,
    Estimate   = estimate,
    SE         = SE,
    df         = df,
    t.ratio     = t.ratio,
    p.value     = p.value
  )

grid.table(
  pairs_df,
  rows = NULL
)
dev.off()

# ---- how the change in n, p, k affect species richness and the ratio of dominant species ----



# identify the dominant species, with control for low richness plot
df <- df.tidy %>%
  semi_join(df_plot_selected)
df <- df %>%
  select(site_code, block, plot, trt, coverage.data) %>%
  unnest(cols = coverage.data) %>%
  select(site_code, block, plot, trt, year_trt, species = Taxon, cover = max_cover)
dom_ids <- df %>%
  group_by(site_code, block, plot, trt, species) %>%
  summarise(mean_cover = mean(cover, na.rm = TRUE), .groups = "drop") %>%
  group_by(site_code, block, plot, trt) %>%
  arrange(desc(mean_cover), .by_group = TRUE) %>%
  mutate(
    richness_plot = n(),  # number of species in this plot × trt
    
    # your rule:
    # if richness <= 2 -> 1 dominant
    # if 3–4          -> 2 dominants
    # if >= 5         -> 3 dominants
    n_dom = case_when(
      richness_plot <= 2 ~ 1L,
      richness_plot <= 4 ~ 2L,
      TRUE               ~ 3L
    ),
    
    is_dominant = row_number() <= n_dom
  )
# calculate the relative share of the dominant species in each plot at each year_trt
# also calculate the species richness of each plot at each year_trt
df_dom <- df %>%
  left_join(dom_ids %>% select(site_code, block, plot, trt, species, is_dominant),
            by = c("site_code", "block", "plot", "trt", "species")) %>%
  mutate(is_dominant = ifelse(is.na(is_dominant), FALSE, is_dominant)) %>%
  group_by(site_code, block, plot, trt, year_trt) %>%
  mutate(
    total = sum(cover, na.rm = TRUE),
    rel   = cover / total
  ) %>%
  summarise(
    dom_share = sum(rel[is_dominant], na.rm = TRUE),
    richness  = n_distinct(species[cover > 0]),
    .groups = "drop"
  )


## ---- RICHNESS MODEL ----

m_rich <- glmmTMB(
  richness ~ trt * year_trt +
    (1 | site_code),
  data   = df_dom,
  family = nbinom2(link = "log")
)

summary(m_rich)

## ---- COEFFICIENT FIGURE (RICHNESS) ----

coefs_rich <- broom.mixed::tidy(m_rich, effects = "fixed", conf.int = TRUE)

coefs_rich2 <- coefs_rich %>% 
  filter(term != "(Intercept)")

coefs_rich2$term <- dplyr::recode(
  coefs_rich2$term,
  "trtK"          = "K vs Control",
  "trtN"          = "N vs Control",
  "trtP"          = "P vs Control",
  "year_trt"      = "Year (Control)",
  "trtK:year_trt" = "K × Year",
  "trtN:year_trt" = "N × Year",
  "trtP:year_trt" = "P × Year"
)

coefs_rich2 %<>%
  mutate(sig.sign = purrr::map2_chr(estimate, p.value, function(x, p) {
    if (x < 0 & p < 0.05) {
      "negative"
    } else if (x > 0 & p < 0.05) {
      "positive"
    } else {
      "non-significant"
    }
  }))

p_rich_a <- ggplot(coefs_rich2, aes(estimate, term)) +
  geom_point(size = 3, aes(colour = sig.sign)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig.sign),
                 height = 0.15) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_bw(base_size = 14) +
  labs(
    x = "Coefficient (log scale)",
    y = "",
    colour = "Significance:"
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom"
  )

## ---- PREDICTION FIGURE (RICHNESS) ----

pred_rich <- ggpredict(m_rich, terms = c("year_trt", "trt")) %>%
  as_tibble() %>%
  mutate(
    trt_lab = factor(
      group,
      levels = c("Control", "N", "P", "K"),
      labels = c("Control", "+N", "+P", "+K")
    )
  )

dark2_cols <- brewer.pal(3, "Dark2")
col_vec <- c(
  "Control" = "grey40",
  "+N"      = dark2_cols[1],
  "+P"      = dark2_cols[2],
  "+K"      = dark2_cols[3]
)

p_rich_b <- ggplot(pred_rich) +
  geom_line(aes(x = x, y = predicted,
                colour = trt_lab, group = trt_lab),
            size = 1.5) +
  geom_line(aes(x = x, y = conf.low,
                colour = trt_lab, group = trt_lab),
            linetype = "dashed", size = 0.9) +
  geom_line(aes(x = x, y = conf.high,
                colour = trt_lab, group = trt_lab),
            linetype = "dashed", size = 0.9) +
  theme_bw(base_size = 14) +
  labs(
    x = "Year",
    y = "Predicted species richness",
    colour = "Treatment:"
  ) +
  scale_colour_manual(values = col_vec) +
  theme(
    plot.title = element_blank(),
    legend.position = "bottom",
    axis.text = element_text(size = 12, color = "black")
  )

p_rich <- p_rich_a + p_rich_b +
  plot_layout(ncol = 2, widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 16))

pdf(paste0(getwd(), "/myfigure/richness.change.pdf"),
    width = 11.5, height = 6)
p_rich
dev.off()


## ---- PREP DOMINANT SHARE FOR BETA MODEL ----

# Smithson & Verkuilen transform or simple clipping to (0,1):
n_obs <- nrow(df_dom)
df_dom <- df_dom %>%
  mutate(
    dom_share_beta = (dom_share * (n_obs - 1) + 0.5) / n_obs
    # alternatively:
    # dom_share_beta = pmin(pmax(dom_share, 1e-4), 1 - 1e-4)
  )

## ---- DOMINANT SHARE MODEL ----

m_dom <- glmmTMB(
  dom_share_beta ~ trt * year_trt +
    (1 | site_code),
  data      = df_dom,
  ziformula = ~ 1,  # allow extra zeros if present
  family    = beta_family(link = "logit")
)

summary(m_dom)

## ---- COEFFICIENT FIGURE (DOM SHARE) ----

coefs_dom <- broom.mixed::tidy(m_dom, effects = "fixed", conf.int = TRUE)

coefs_dom2 <- coefs_dom %>%
  filter(term != "(Intercept)")

coefs_dom2$term <- dplyr::recode(
  coefs_dom2$term,
  "trtK"                  = "K vs Control",
  "trtN"                  = "N vs Control",
  "trtP"                  = "P vs Control",
  "year_trt"       = "Year (Control)",
  "trtK:year_trt"  = "K × Year",
  "trtN:year_trt"  = "N × Year",
  "trtP:year_trt"  = "P × Year"
)

coefs_dom2 %<>%
  mutate(sig.sign = purrr::map2_chr(estimate, p.value, function(x, p) {
    if (x < 0 & p < 0.05) {
      "negative"
    } else if (x > 0 & p < 0.05) {
      "positive"
    } else {
      "non-significant"
    }
  }))

p_dom_a <- ggplot(coefs_dom2, aes(estimate, term)) +
  geom_point(size = 3, aes(colour = sig.sign)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, colour = sig.sign),
                 height = 0.15) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 14) +
  labs(
    x = "Coefficient (logit scale)",
    y = "",
    colour = "Significance:"
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "bottom"
  )

## ---- PREDICTION FIGURE (DOM SHARE) ----

pred_dom <- ggpredict(m_dom, terms = c("year_trt", "trt")) %>%
  as_tibble() %>%
  mutate(
    # back to original year for x-axis:
    x_year = x + mean(df_dom$year_trt, na.rm = TRUE),
    trt_lab = factor(
      group,
      levels = c("Control", "N", "P", "K"),
      labels = c("Control", "+N", "+P", "+K")
    )
  )

p_dom_b <- ggplot(pred_dom) +
  geom_line(aes(x = x_year, y = predicted,
                colour = trt_lab, group = trt_lab),
            size = 1.5) +
  geom_line(aes(x = x_year, y = conf.low,
                colour = trt_lab, group = trt_lab),
            linetype = "dashed", size = 0.9) +
  geom_line(aes(x = x_year, y = conf.high,
                colour = trt_lab, group = trt_lab),
            linetype = "dashed", size = 0.9) +
  theme_bw(base_size = 14) +
  labs(
    x = "Year",
    y = "Predicted dominant share",
    colour = "Treatment:"
  ) +
  scale_colour_manual(values = col_vec) +
  theme(
    plot.title = element_blank(),
    legend.position = "bottom",
    axis.text = element_text(size = 12, color = "black")
  )

p_dom <- p_dom_a + p_dom_b +
  plot_layout(ncol = 2, widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 16))

pdf(paste0(getwd(), "/myfigure/dominant.share.change.pdf"),
    width = 11.5, height = 6)
p_dom
dev.off()


# ---- address the question about the N:P ratio ----
# I initially considered inlcuding both perturbation intensity and the NP ratio in the model,
# but the strong collinearity between them diluted the effect of perturbaiton intensity and make the results not reliable
# I therefore only conducted a supplementary analysis here focusing on NP ratio only
# stabilities of P treatment
library(lme4)
library(lmerTest)
df_P_stabilities <- df.tidy.stabilities %>%
  filter(trt == "P") 
# add base line N and P concentration
df_P_stabilities <- df_P_stabilities %>%
  left_join(
    df.tidy %>%
      select(site_code, block, trt, plot, P_added_proportion, pct_N, ppm_P) %>%
      distinct(),
    by = c("site_code", "block", "plot", "trt")
  )
# ratio of N:P
df_P_stabilities <- df_P_stabilities %>%
  mutate(
    N_P_ratio = log(pct_N / ppm_P),
    logNP_sc = scale(N_P_ratio)[, 1],
    timescale = scale(end_year)[, 1],
  )
# new model including N:P ratio as moderator
# model for fucn_var_det
m_P_func_var_det <- lmer(
  func_var_det ~ 
    logNP_sc       * timescale +
    (1 | site_code),
  data = df_P_stabilities
)
summary(m_P_func_var_det)

# model for comp_var
# one plot with comp_var as Inf, remove it
df_P_stabilities_2 <- df_P_stabilities %>% 
  filter(!is.infinite(comp_var))
m_P_comp_var <- lmer(
  comp_var ~ 
    logNP_sc       * timescale +
    (1 | site_code),
  data = df_P_stabilities_2 %>% na.omit()
)
summary(m_P_comp_var)

# model for func_change
df_P_stabilities_3 <- df_P_stabilities %>% 
  filter(!is.infinite(func_change))
m_P_func_change <- lmer(
  func_change ~ 
    logNP_sc       * timescale +
    (1 | site_code),
  data = df_P_stabilities_3)
summary(m_P_func_change)

# model for comp_change
df_stabilities_4 <- df_P_stabilities %>%
  filter(!is.infinite(comp_change))
m_P_comp_change <- lmer(
  comp_change ~ 
    logNP_sc       * timescale +
    (1 | site_code),
  data = df_stabilities_4)
summary(m_P_comp_change)

# make a figure showing the results as point-errorbar


## 1. 写一个小工具函数：抽取 fixed effects + 置信区间，并标记 component ----
get_coefs <- function(model, comp_label) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(component = comp_label)
}

## 2. 从四个模型分别抽系数并绑定在一起 ----
coefs_all <- bind_rows(
  get_coefs(m_P_func_var_det, "func_var_det"),
  get_coefs(m_P_comp_var,      "comp_var"),
  get_coefs(m_P_func_change,   "func_change"),
  get_coefs(m_P_comp_change,   "comp_change")
)

## 3. 美化 term 名称、component 名称，并标记显著性方向 ----
coefs_all <- coefs_all %>%
  mutate(
    term_pretty = recode(
      term,
      "logNP_sc"             = "N:P ratio",
      "timescale"            = "Years",
      "logNP_sc:timescale"   = "N:P ratio × Years"
    ),
    # 统一 term 顺序（行的顺序）
    term_pretty = factor(
      term_pretty,
      levels = c(
        "N:P ratio",
        "Years",
        "N:P ratio × Years"
      )
    ),
    # component 名称与顺序（列的顺序）
    component = factor(
      component,
      levels = c("func_var_det", "comp_var", "func_change", "comp_change"),
      labels = c("func_var_det", "comp_var", "func_change", "comp_change")
    ),
    # 显著性和方向：红（正）、蓝（负）、灰（不显著）
    sig.sign = case_when(
      estimate > 0 & p.value < 0.05 ~ "positive",
      estimate < 0 & p.value < 0.05 ~ "negative",
      TRUE                          ~ "non-significant"
    ),
    sig.sign = factor(
      sig.sign,
      levels = c("negative", "non-significant", "positive")
    )
  )

## 4. 颜色与之前 legume / dominance / richness 的风格保持一致 ----
col_vec_sig <- c(
  "negative"        = "blue",
  "non-significant" = "grey60",
  "positive"        = "red"
)

## 5. 画图：facet_wrap 列为 component，行是 predictors（term_pretty） ----
p_P_NP_coefs <- ggplot(coefs_all,
                       aes(x = estimate, y = term_pretty, colour = sig.sign)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.15) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  facet_wrap(~ component, nrow = 1, scales = "free_x") +
  scale_colour_manual(
    values = col_vec_sig,
    name   = "Significance:"
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "The Coefficient Estimate",
    y = ""
  ) +
  theme(
    axis.text.y      = element_text(size = 13, colour = "black"),
    axis.text.x      = element_text(size = 11, colour = "black"),
    legend.text    = element_text(size = 13),
    strip.text.x     = element_text(face = "bold", size = 13),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

pdf(paste0(getwd(), "/myfigure/P_NP_stabilities_coefs.pdf"),
    width = 9.9, height = 4.5)
p_P_NP_coefs
dev.off()


