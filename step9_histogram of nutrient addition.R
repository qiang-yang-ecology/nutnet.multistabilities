df.tidy.N <- filter(df.tidy, trt=="N")
df.tidy.P <- filter(df.tidy, trt=="P")




# make the six plots in one page
# move the infinite values if there are any
par(mfrow=c(3,2))
hist(df.tidy.N$N_added_proportion, xlab = "Value", main = "N added proportion")
hist(df.tidy.P$P_added_proportion, xlab = "Value", main = "P added proportion")
hist(log(df.tidy.N$N_added_proportion), xlab = "Value", main = "log N added proportion")
hist(log(df.tidy.P$P_added_proportion), xlab = "Value", main = "log P added proportion")
hist(scale(log(df.tidy.N$N_added_proportion)), xlab = "Value", main = "scaled log N added proportion")
hist(scale(log(df.tidy.P$P_added_proportion)), xlab = "Value", main = "scaled log P added proportion")


# later Ian asked me to send the table of the original nutrient levels as well as the perturbation intensity and a map

# the original table
load(paste0(data.path, "/df.tidy.stabilities.Rdata"))
load(paste0(data.path, "/df.tidy.Rdata"))
df_nutrients <- df.tidy %>%
  select(site_code, block, plot, trt, latitude, longitude,
         pct_N, ppm_P, ppm_K, contains("proportion"),  ) %>%
  distinct()

df.tidy.N <- filter(df_nutrients, trt=="N") %>% select(-ppm_P, -ppm_K, -P_added_proportion, -K_added_proportion)
df.tidy.P <- filter(df_nutrients, trt=="P") %>% select(-pct_N, -ppm_K, -N_added_proportion, -K_added_proportion)
df.tidy.K <- filter(df_nutrients, trt=="K") %>% select(-pct_N, -ppm_P, -N_added_proportion, -P_added_proportion)
# remove infinite values (if any)
df.tidy.N <- df.tidy.N %>%
  filter(is.finite(pct_N)|is.finite(N_added_proportion))
df.tidy.P <- df.tidy.P %>%
  filter(is.finite(ppm_P)|is.finite(P_added_proportion))
df.tidy.K <- df.tidy.K %>%
  filter(is.finite(ppm_K)|is.finite(K_added_proportion))

df.tidy.N$pct_N %>% range()
df.tidy.P$ppm_P %>% range()
df.tidy.K$ppm_K %>% range()

# the table of the original nutrient levels as well as the perturbation intensity
df.tidy.N <- df.tidy.N %>%
  mutate(
    Ri_raw = N_added_proportion,
    Ri     = log(Ri_raw),
    Ri_sc  = scale(Ri)[, 1]
  )
df.tidy.P <- df.tidy.P %>%
  mutate(
    Ri_raw = P_added_proportion,
    Ri     = log(Ri_raw),
    Ri_sc  = scale(Ri)[, 1]
  )
df.tidy.K <- df.tidy.K %>%
  mutate(
    Ri_raw = K_added_proportion,
    Ri     = log(Ri_raw),
    Ri_sc  = scale(Ri)[, 1]
  )
# save(the three tables in one excel file
# each table will be a sheet
# short file name
write.xlsx(list(N = df.tidy.N, P = df.tidy.P, K = df.tidy.K), file = paste0(data.path, "/original_nutrient_levels_and_perturbation_intensity.xlsx"))

# make maps
# calculate site level mean of original nutrient level and perturbation intensity
df.tidy.N.site <- df.tidy.N %>%
  group_by(site_code) %>%
  summarise(
    latitude = mean(latitude),
    longitude = mean(longitude),
    pct_N = mean(pct_N, na.rm = TRUE),
    Ri_sc = mean(Ri_sc, na.rm = TRUE)
  ) %>% 
  ungroup()

df.tidy.P.site <- df.tidy.P %>%
  group_by(site_code) %>%
  summarise(
    latitude = mean(latitude),
    longitude = mean(longitude),
    ppm_P = mean(ppm_P, na.rm = TRUE),
    Ri_sc = mean(Ri_sc, na.rm = TRUE)
  ) %>% 
  ungroup()

df.tidy.K.site <- df.tidy.K %>%
  group_by(site_code) %>%
  summarise(
    latitude = mean(latitude),
    longitude = mean(longitude),
    ppm_K = mean(ppm_K, na.rm = TRUE),
    Ri_sc = mean(Ri_sc, na.rm = TRUE)
  ) %>% 
  ungroup()

# make global map with color gradients as original nutrient levels and perturbation intensity
# for N, P, K using colors from Dark2 palette in RColorBrewer

library(sf)
library(rnaturalearth)
library(RColorBrewer)

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Convert site data to an sf object
sites_sf_N <- st_as_sf(df.tidy.N.site, 
                     coords = c("longitude", "latitude"), crs = 4326)
sites_sf_P <- st_as_sf(df.tidy.P.site, 
                       coords = c("longitude", "latitude"), crs = 4326)
sites_sf_K <- st_as_sf(df.tidy.K.site, 
                       coords = c("longitude", "latitude"), crs = 4326)

# Create the map for N perturbation intensity (Ri_sc)
p.map.N.Ri_sc <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray40", linewidth = 0.1) +
  geom_sf(
    data = sites_sf_N,
    aes(fill = Ri_sc),
    color = "white",
    size = 2,
    shape = 21,
    stroke = 0.2
  ) +
  scale_fill_gradient(
    low = "#E6F5EE",
    high = "#1B9E77",
    name = "N enrichment"
  ) +
  theme_minimal() +
  labs(title = "N perturbation intensity (Ri_sc)") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )

# Create the map for P perturbation intensity (Ri_sc)
p.map.P.Ri_sc <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray40", linewidth = 0.1) +
  geom_sf(
    data = sites_sf_P,
    aes(fill = Ri_sc),
    color = "white",
    size = 2,
    shape = 21,
    stroke = 0.2
  ) +
  scale_fill_gradient(
    low = "#FEE8D8",
    high = "#D95F02",
    name = "P enrichment"
  ) +
  theme_minimal() +
  labs(title = "P perturbation intensity (Ri_sc)") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )

# Create the map for K perturbation intensity (Ri_sc)
p.map.K.Ri_sc <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray40", linewidth = 0.1) +
  geom_sf(
    data = sites_sf_K,
    aes(fill = Ri_sc),
    color = "white",
    size = 2,
    shape = 21,
    stroke = 0.2
  ) +
  scale_fill_gradient(
    low = "#ECEAF6",
    high = "#7570B3",
    name = "K enrichment"
  ) +
  theme_minimal() +
  labs(title = "K perturbation intensity (Ri_sc)") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )

# basal N map
p.map.N.basal <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray40", linewidth = 0.1) +
  geom_sf(
    data = sites_sf_N,
    aes(fill = pct_N),
    color = "white",
    size = 2,
    shape = 21,
    stroke = 0.2
  ) +
  scale_fill_gradient(
    low = "#E6F5EE",
    high = "#1B9E77",
    name = "Basal N (ppt)"
  ) +
  theme_minimal() +
  labs(title = "Basal N level") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )
# basal P
p.map.P.basal <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray40", linewidth = 0.1) +
  geom_sf(
    data = sites_sf_P,
    aes(fill = ppm_P),
    color = "white",
    size = 2,
    shape = 21,
    stroke = 0.2
  ) +
  scale_fill_gradient(
    low = "#FEE8D8",
    high = "#D95F02",
    name = "Basal P (ppm)"
  ) +
  theme_minimal() +
  labs(title = "Basal P level") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )
# basal K
p.map.K.basal <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray40", linewidth = 0.1) +
  geom_sf(
    data = sites_sf_K,
    aes(fill = ppm_K),
    color = "white",
    size = 2,
    shape = 21,
    stroke = 0.2
  ) +
  scale_fill_gradient(
    low = "#ECEAF6",
    high = "#7570B3",
    name = "Basal K (ppm)"
  ) +
  theme_minimal() +
  labs(title = "Basal K level") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )

pdf(paste0(getwd(), "/myfigure/Maps of original nutrient levels and perturbation intensity.pdf"), width = 8, height = 6)
print(p.map.N.Ri_sc)
print(p.map.P.Ri_sc)
print(p.map.K.Ri_sc)
print(p.map.N.basal)
print(p.map.P.basal)
print(p.map.K.basal)
dev.off()