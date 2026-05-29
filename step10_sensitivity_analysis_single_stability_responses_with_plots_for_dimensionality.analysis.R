# set up working environment ##############################
rm(list = ls())
gc()
# necessary R packages
Packages <- c(
  "tidyverse", "magrittr", "grid", "gridExtra",
  "scales", "broom", "nlme", "sjPlot", "patchwork", "lme4", "lmerTest"
)
pacman::p_load(char = Packages)
# functions with conflicting status (presence in different packages used here)
select <- dplyr::select
summarise <- dplyr::summarise

# working path
main.path <- getwd()
data.path <- paste0(main.path, "/mydata")
figure.path <- paste0(main.path, "/myfigure")
set.seed(1234)

# prepare data #################
# here i will select the plots that is only used in dimensionality analysis

load(paste0(data.path, "/df.tidy.stabilities.Rdata"))
load(paste0(data.path, "/df.tidy.Rdata"))
df.full <- df.tidy.stabilities %>%
  left_join(df.tidy %>% select(-contains("data")))
rm(df.tidy, df.tidy.stabilities)

df.full %<>%
  select(site_code:end_year, matches("_added_proportion|change|var"))
sta.names <- colnames(df.full) %>% str_subset("change|var")

#plots in dimensionality analysis
df_plots_dimensionality <- read_csv(file = "plots_grouping_according_to_perturbation.csv")

df.full <- df.full %>%
  semi_join(df_plots_dimensionality, by = c("site_code", "block", "plot"))



# scale the stability components here
# some of the columns include Inf and NA, so make a new scale function
new.scale <- function(x) {
  x[is.infinite(x)] <- NA
  y <- scale(x)
}
df.full <- df.full %>%
  mutate(
    func_var = new.scale(func_var),
    func_var_det = new.scale(func_var_det),
    comp_var = new.scale(comp_var),
    func_change = new.scale(func_change),
    comp_change = new.scale(comp_change)
  )

# data for the treatments N, P, and K
df.meta.N.P.K <- data.frame(trt = c("N", "P", "K")) %>%
  mutate(pert.name = paste0(trt, "_added_proportion")) %>%
  expand_grid(sta.names)

df.N.P.K.nested <- df.meta.N.P.K %>%
  mutate(
    data = pmap(list(trt, pert.name, sta.names), function(x1, x2, x3) {
      output <- df.full %>%
        filter(trt == x1) %>%
        select(one_of(c("site_code", "block", "plot", "end_year", x2, x3))) %>%
        set_colnames(c("site_code", "block", "plot", "end_year", "perturbation_intensity", "stability"))
    })
  )

# note that there will be missing data in the stability column
# this is due to the availability of cover data or biomass data for some plots
# remove the NA values
df.N.P.K.nested <- df.N.P.K.nested %>%
  mutate(data = map(data, na.omit))
# I also noted that there are infinite values in the stability column
#  remove them
df.N.P.K.nested <- df.N.P.K.nested %>%
  mutate(data = map(data, function(df1) df1 %>% filter(!is.infinite(stability))))

# perform the regression##############################

df.models.N.P.K <- df.N.P.K.nested %>%
  mutate(mod.out = map(data, function(dat) {
    # log transform perturbation intensity
    dat %<>% mutate(log.pert = log(perturbation_intensity))
    dat %<>% na.omit()
    mod.dat <- dat
    mod.dat %<>% mutate(across(c(end_year, log.pert), scale))
    formula_ <- stability ~ log.pert + end_year + log.pert:end_year + (1 | site_code)
    mod <- lmer(formula = formula_, data = mod.dat)
    return(mod)
  }))



# organize model coefficient estimates and p values ###############
  
df.models.N.P.K.coeffficients <- df.models.N.P.K %>%
  mutate(mod.coefficient = map(mod.out, function(mod) {
    # coefficient estimate from the summary function
    output <- summary(mod)$coefficients
    coeff <- output %>% as.data.frame()
    coeff %<>% mutate(predictor = rownames(output))
    coeff %<>% select(predictor, Estimate)
    # confidence interval from lme4::confint.merMod()
    # profile confidence intervals
    cfs <- confint.merMod(mod)
    cfs %<>% as.data.frame()
    cfs %<>% mutate(predictor = rownames(cfs))
    cfs %<>% rename(lcf = `2.5 %`, ucf = `97.5 %`)
    cfs %<>% filter(!str_detect(predictor, ".sig"))
    # P value from anova (lmerTest borrow from pkgtest)
    df.p.value <- anova(mod) %>% as.data.frame()
    df.p.value %<>% select(NumDF:`Pr(>F)`)
    df.p.value %<>% mutate(predictor = rownames(df.p.value))
    
    # add significance sign
    coeff <- coeff %>%
      left_join(cfs) %>%
      left_join(df.p.value)
    coeff %<>% filter(predictor != "(Intercept)")
    coeff %<>%
      mutate(sig.sign = map2_chr(Estimate, `Pr(>F)`, function(x, p) {
        if (x < 0 & p < 0.05) {
          out <- "negative"
        } else if (x > 0 & p < 0.05) {
          out <- "positive"
        } else {
          out <- "non-significant"
        }
        return(out)
      }))
    return(coeff)
  }))


df.models.N.P.K.coeffficients %<>%
  mutate(model.ID = 1:nrow(df.models.N.P.K.coeffficients)) %>%
  select(-pert.name, -data, -mod.out) %>%
  select(model.ID, everything()) %>%
  unnest(cols = mod.coefficient)


##### visualize model coefficient estimates ############################
  
df.coefficients <- df.models.N.P.K.coeffficients %>%
  mutate(predictor = str_replace_all(predictor, "end_year", "years")) %>%
  mutate(predictor = str_replace_all(predictor, "log.pert", "perturbation"))
stabilities <- df.coefficients %>%
  select(sta.names) %>%
  distinct() %>%
  pull()
df.coefficients %<>%
  mutate(sta.names = as.factor(sta.names)) %>%
  mutate(sta.names = fct_relevel(sta.names, stabilities)) %>%
  mutate(trt = as.factor(trt)) %>%
  mutate(trt = fct_relevel(trt, c("N", "P", "K"))) %>%
  mutate(predictor = as.factor(predictor)) %>%
  mutate(predictor = fct_relevel(predictor, c("years", "perturbation", "perturbation:years")))

p1 <- df.coefficients %>%
  # filter(str_detect(sta.names, "functional")) %>%
  ggplot(aes(x = predictor, y = Estimate, color = sig.sign)) +
  geom_point() +
  geom_errorbar(aes(ymin = lcf, ymax = ucf), width = .2) +
  facet_grid(trt ~ sta.names) +
  coord_flip() +
  geom_hline(yintercept = 0, col = "grey", linetype = "dashed", size = 0.3) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_linedraw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(color = "black", size = 8),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, color = "black", angle = 90),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold", color = "blue"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  xlab("Predictors") +
  ylab("The Coefficient Estimate")

trt.labs <- c("+N", "+P", "+K")
names(trt.labs) <- c("N", "P", "K")

# Create the plot
(p1 <- p1 + facet_grid(
  trt ~ sta.names,
  labeller = labeller(trt = trt.labs)
))



pdf(paste0(figure.path, "/lmm.coefficients_site.as.random.intercepts_sensitity_analysis.pdf"), width = 6, height = 4)
print(p1)
dev.off()


