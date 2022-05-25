---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/25 at 16:08"
output: 
  html_document:
    keep_md: yes
  github_document:
    always_allow_html: true
---

# Goals

The goal of the statistics are to investigate the potential effects of primary production and advection on mesopelagic backscatter for each Arctic region. Therefore, I gathered remote sensing data of sea ice concentration and ocean colour, and modelled advection from Copernicus.

Since remote sensing of primary production (chlorophyll *a*) can be biased by high concentration of CDOM (chromophoric dissolved organic matter) and the photoadptation of Arctic phytoplankton (Lewis et al. 2016), I use to variables as proxy for primary production. First, I use the satellite-derived monthly averaged chlrophyll *a* concentration from CMEMS (doi: [10.48670/moi-00066](https://doi.org/10.48670/moi-00066)), and I calculate the open water days from sea ice concentration from the EUMETSAT OSI SAF remote sensing dataset (doi: [10.24381/cds.3cd8b812](https://doi.org/10.24381/cds.3cd8b812)) which I use as a proxy for primary production. For advection I use data from the Copernicus Marine Environmental Monitoring Service Arctic Ocean Physics Reanalysis product which gives monthly average of current velocity in the Arctic Ocean.

# Package loading


```r
# Load packages
library(tidyverse)    # Tidy code
library(cowplot)      # Plots on a grid
library(raster)       # Data gridding
library(sf)           # Spatial data
library(rgdal)        # Read shapefiles
library(ggpubr)       # Deal with stats
library(ggfortify)    # Plotting glm
library(RColorBrewer) # Diverging colour palettes
library(cmocean)      # Oceanographic colour palettes
library(moments)      # Overlay distributions
library(ggcorrplot)   # Correlation plots
library(kableExtra)   # Pretty tables
library(mgcv)         # Fit GAM
library(gratia)       # Visualise GAM
library(visreg)       # Visualise GAM
library(itsadug)      # Predict GAM
library(MetBrewer)    # Pretty colour palettes
library(MuMIn)        # AIC weights
library(DT)           # Interactive table
library(rstatix)      # Pipe-friendly stats
library(broom)        # LM
# Custom figure theme
theme_set(theme_bw())
theme_update(axis.text = element_text(size = 9),
             axis.title = element_text(size = 9),
             strip.text.x = element_text(size = 9, face = "plain", hjust = 0.5),
             strip.background = element_rect(colour = "transparent",
                                             fill = "transparent"),
             legend.title = element_text(size = 9),
             legend.margin = margin(0, 0, 0, 0),
             legend.box.margin = margin(0, 0, -8, 0),
             panel.grid = element_blank(), 
             plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
# Suppress summarise() warning
options(dplyr.summarise.inform = F)
```


```r
# Laea projection
cell_res <- 25
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931")
projection(arctic_laea) <- gsub("units=m", "units=km", projection(arctic_laea)) 
res(arctic_laea) <- c(cell_res, cell_res) 

# Coastline shapefiles reprojected to laea
coast_10m_laea <- readOGR("data/bathy/ne_10m_land.shp", verbose = F) %>% 
  spTransform(CRSobj = crs("EPSG:4326")) %>% 
  crop(extent(-180, 180, 0, 90)) %>%
  spTransform(CRSobj = crs(arctic_laea)) %>% 
  fortify() %>% 
  rename(xc = long, yc = lat)

# Gridded acoustic, CTD, and sea ice data
load(paste0("data/acoustics/SA_grids_", cell_res, "km.RData")) 
load(paste0("data/remote_sensing/physics_grids_", cell_res, "km.RData"))
load(paste0("data/remote_sensing/chl_grids_", cell_res, "km.RData"))
```

# Data preparation

I combine backscatter and environmental data which were gridded on the Lambert Azimuthal Equal Area grid (EPSG:6931) with a cell resolution of 25 x 25 km. Since I do not have a lot of samples from the the Beaufort Sea and West Arctic Ocean, I combine those two regions due to their geographical proximity.


```r
# Gridded SA
SA_laea <- SA_grid_laea %>%  
  dplyr::select(-lat, -lon) 
# Remote sensing: physics reanalysis
phy_laea <- phy_grid_laea %>% 
  dplyr::select(year, area, xc, yc, cell_res, depth, velocity) %>%
  # Find missing depth values
  complete(depth, nesting(year, area, xc, yc, cell_res), fill = list(velocity = NA)) %>%
  arrange(year, area, xc, yc, cell_res, depth) %>%
  mutate(depth = factor(depth))

stat_laea <- SA_laea %>%
  left_join(., phy_laea, by = c("year", "xc", "yc", "area", "cell_res")) %>%
  left_join(., chl_grid_laea, by = c("year", "xc", "yc", "area", "cell_res")) 
```

There are NAs in the remote sensing dataset.


```r
plot_grid(
  stat_laea %>% 
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = chl)) +
    scale_colour_cmocean(name = "algae", na.value = "red") + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank()),
  stat_laea %>% 
    filter(depth == 380) %>%
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = velocity)) +
    scale_colour_cmocean(name = "speed", na.value = "red") + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank()),
  ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_statistics_files/figure-html/map-NA-1.png)<!-- -->

This is because the coverage of acoustic data and remote sensing products slightly differ close to the coasts. For sea ice concentration (and open water duration) the NAs are close to land while for ocean colour they are in ice covered areas. I use the average value of neighbouring cells (24 cells around the pixel of interest) to fill missing values.


```r
stat_laea <- stat_laea %>%
  rowwise() %>%
  mutate( 
    # Ocean colour
    xc_v = if_else(is.na(velocity) == T, xc, NaN),
    yc_v = if_else(is.na(velocity) == T, yc, NaN),
    year_v = if_else(is.na(velocity) == T, year, NaN), 
    velocity = if_else(
      is.na(velocity) == T,
      mean(pull(subset(phy_grid_laea, 
                       xc >= xc_v - 2 * cell_res & xc <= xc_v + 2 * cell_res & 
                         yc >= yc_v - 2 * cell_res & yc <= yc_v + 2 * cell_res & 
                         year == year_v,
                       select = velocity),
                velocity),
           na.rm = T),
      velocity),
    # Ocean colour
    xc_chl = if_else(is.na(chl) == T, xc, NaN),
    yc_chl = if_else(is.na(chl) == T, yc, NaN),
    year_chl = if_else(is.na(chl) == T, year, NaN), 
    chl = if_else(
      is.na(chl) == T,
      mean(pull(subset(chl_grid_laea, 
                       xc >= xc_chl - 2 * cell_res & xc <= xc_chl + 2 * cell_res & 
                         yc >= yc_chl - 2 * cell_res & yc <= yc_chl + 2 * cell_res & 
                         year == year_chl,
                       select = chl),
                chl),
           na.rm = T),
      chl),
    # Convert veloctity to cm/s
    velocity = velocity * 100, 
    IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                IHO_area == "West Arctic Ocean" ~ "WAO_BF",
                                IHO_area == "Beaufort Sea" ~ "WAO_BF",
                                IHO_area == "The Northwestern Passages" ~ "CAA",
                                IHO_area == "Baffin Bay" ~ "BB",
                                IHO_area == "Davis Strait" ~ "DS"),
                      levels = c("WAO_BF", "CAA", "BB", "DS", "EAO"))) %>%
  ungroup() %>%
  dplyr::select(-xc_chl, -yc_chl, -year_chl)
```

Check if I replace NAs correctly.


```r
plot_grid(
  stat_laea %>% 
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = chl)) +
    scale_colour_cmocean(name = "algae", na.value = "red") + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank()),
  stat_laea %>%
    filter(depth == 380) %>%
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = velocity)) +
    scale_colour_cmocean(name = "speed", na.value = "red", limits = c(0, 6)) + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()),
  ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_statistics_files/figure-html/map-NA-fix-1.png)<!-- -->

There are still one missing chl values, so this cell is excluded from further analyses.


```r
stat_laea <- stat_laea %>%
  filter(is.na(chl) == F & is.na(velocity) == F)
```


```r
plot_grid(
  stat_laea %>% 
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = chl)) +
    scale_colour_cmocean(name = "algae", na.value = "red") + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank()),
  stat_laea %>%
    filter(depth == 380) %>%
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = velocity)) +
    scale_colour_cmocean(name = "speed", na.value = "red", limits = c(0, 6)) + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()),
  ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_statistics_files/figure-html/map-chl-NA-1.png)<!-- -->

Check the total number of observations per group.


```r
stat_laea %>%
  filter(depth == 318) %>%
  group_by(IHO_area) %>% 
  summarise(n = n())
```

```
## # A tibble: 5 × 2
##   IHO_area     n
##   <fct>    <int>
## 1 WAO_BF      23
## 2 CAA         33
## 3 BB          64
## 4 DS          28
## 5 EAO         45
```

# Data exploration

Maps of all variables.


```r
# IHO Regions
stat_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  geom_point(aes(col = IHO_area)) +
  ggtitle("regions") +
  facet_wrap(~ year) +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(legend.position = "top",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-SA-int-1.png)<!-- -->


```r
# Temperature
chl_grid_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = chl)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("chl (mg m-3)", name = "algae", limits = c(0,5)) +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("chlorophyll") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(legend.position = "top",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-chl-1.png)<!-- -->


```r
# Velocity
phy_grid_laea %>% 
  filter(depth == 318) %>%
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = velocity * 100)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("velo (cm/s)", name = "speed", limits = c(0, 10)) +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Velocity at 318 m depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(legend.position = "top",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-velocity318-1.png)<!-- -->


```r
# Velocity
phy_grid_laea %>% 
  filter(depth == 380) %>%
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = velocity * 100)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("velo (cm/s)", name = "speed", limits = c(0, 10)) +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Velocity at 380 m depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(legend.position = "top",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-velocity380-1.png)<!-- -->

# Spatial and interannual variability


```r
SA_diff <- SA_grid_laea %>%
  dplyr::select(year, IHO_area, NASC_int, SA_int, NASC_int, CM) %>%
  mutate(IHO_area = 
           factor(
             case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                       IHO_area == "West Arctic Ocean" ~ "WAO_BF",
                       IHO_area == "Beaufort Sea" ~ "WAO_BF",
                       IHO_area == "The Northwestern Passages" ~ "CAA",
                       IHO_area == "Baffin Bay" ~ "BB",
                       IHO_area == "Davis Strait" ~ "DS"),
             levels = c("WAO_BF", "CAA", "BB", "DS", "EAO")),
         year = factor(year))
```

## All areas

I check whether there are inter-annual and spatial variability in S\~A\~ and the centre of mass across years and areas.


```r
# Visualize data
plot_grid(SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = NASC_int)),
          SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = CM)))
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-plot-1.png)<!-- -->


```r
# Normality: QQ plot
plot_grid(ggqqplot(SA_diff, "NASC_int", facet.by = "year"),
          ggqqplot(SA_diff, "CM", facet.by = "year"),
          ggqqplot(SA_diff, "NASC_int", facet.by = "IHO_area"),
          ggqqplot(SA_diff, "CM", facet.by = "IHO_area"), 
          ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/normality-check-1.png)<!-- -->

I used a non parametric Kruskal-Wallis test to test for interannual and spatial variability because the variance is not homogenous and not very normal.


```r
bind_rows(kruskal_test(SA_diff, NASC_int ~ IHO_area),
          kruskal_test(SA_diff, NASC_int ~ year),
          kruskal_test(SA_diff, CM ~ IHO_area),
          kruskal_test(SA_diff, CM ~ year)) %>%
  mutate(var = c("IHO_area", "year", "IHO_area", "year"))
```

```
## # A tibble: 4 × 7
##   .y.          n statistic    df        p method         var     
##   <chr>    <int>     <dbl> <int>    <dbl> <chr>          <chr>   
## 1 NASC_int   199     43.2      4 9.48e- 9 Kruskal-Wallis IHO_area
## 2 NASC_int   199     17.1      2 1.97e- 4 Kruskal-Wallis year    
## 3 CM         199     69.7      4 2.58e-14 Kruskal-Wallis IHO_area
## 4 CM         199      3.50     2 1.74e- 1 Kruskal-Wallis year
```

Mesopelagic backscatter S\~A\~ varied significantly between areas (*H* = 29.792952, *p* \< 0.001) and years (*H* = 29.792952, *p* = 0.011). The center of mass varied significantly between areas (*H* = 61.185994, *p* \< 0.001) but not between years (*H* = 4.164743, *p* = 0.125).

## Within areas

I check whether there are inter-annual variability in S\~A\~ across years within areas. This will be used to decide whether a random effect is needed in the GAM. I used a non parametric Kruskal-Wallis test to test for interannual variability within areas.


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(IHO_area) %>%
  kruskal_test(NASC_int ~ year)
```

```
## # A tibble: 5 × 7
##   IHO_area .y.          n statistic    df        p method        
## * <fct>    <chr>    <int>     <dbl> <int>    <dbl> <chr>         
## 1 WAO_BF   NASC_int    23    11.6       2 0.00307  Kruskal-Wallis
## 2 CAA      NASC_int    33    15.6       2 0.000403 Kruskal-Wallis
## 3 BB       NASC_int    65     5.74      2 0.0567   Kruskal-Wallis
## 4 DS       NASC_int    28     0.430     2 0.807    Kruskal-Wallis
## 5 EAO      NASC_int    50     8.05      2 0.0178   Kruskal-Wallis
```

There were some interannual differences in SA within some region (WAO_BF - H = 6.4670868, *p* = 0.0394; CAA - H = 11.8260406, *p* = 0.0027; BB - H = 5.0981791, *p* = 0.0782; DS - H = 0.3282051, *p* = 0.8490; EAO - H = 4.9779498, *p* = 0.0830).

# LM - Latitude - S\~A\~ int

A number of studies show decreasing mesopelagic backscatter with increasing latitude. To test this hypothesis, I use a linear regression. Because there are some high mesopelagic backscatter values in the East Arctic Ocean I transformed the nautical area scattering coefficient into nautical area scattering strength in dB.

## Data preparation

First I prepare the data.


```r
SA_lm <- SA_grid_laea %>%
  mutate(IHO = 
           factor(
             case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                       IHO_area == "West Arctic Ocean" ~ "WAO",
                       IHO_area == "Beaufort Sea" ~ "BF",
                       IHO_area == "The Northwestern Passages" ~ "CAA",
                       IHO_area == "Baffin Bay" ~ "BB",
                       IHO_area == "Davis Strait" ~ "DS"),
             levels = c("WAO", "BF", "CAA", "BB", "DS", "EAO")))
```

Plot the relationship I want to test.


```r
SA_lm %>%
  ggplot(aes(x = lat, y = SA_int)) +
  geom_point() + 
  geom_smooth(aes(x = lat, y = SA_int, group = IHO, col = IHO), method = "lm") +
  geom_smooth(aes(x = lat, y = SA_int), method = "lm", col = "red") 
```

```
## `geom_smooth()` using formula 'y ~ x'
## `geom_smooth()` using formula 'y ~ x'
```

![](PanArctic_DSL_statistics_files/figure-html/lm-plot-1.png)<!-- -->

## All areas

Fit linear regression.


```r
LM_all <- lm(SA_int ~ lat, data = SA_grid_laea)
summary(LM_all)
```

```
## 
## Call:
## lm(formula = SA_int ~ lat, data = SA_grid_laea)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -20.963  -4.861  -0.018   4.878  35.060 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  53.4543     9.3752   5.702 4.29e-08 ***
## lat          -0.4865     0.1249  -3.894 0.000135 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 8.627 on 197 degrees of freedom
## Multiple R-squared:  0.07149,	Adjusted R-squared:  0.06677 
## F-statistic: 15.17 on 1 and 197 DF,  p-value: 0.0001347
```

Check residuals.


```r
par(mfrow = c(2, 2))
plot(LM_all)
```

![](PanArctic_DSL_statistics_files/figure-html/residuals-lm-1.png)<!-- -->

## LM per area

Fit linear regression for each Arctic region.


```r
# Fit models
LM_area <- SA_lm %>%
  nest_by(IHO) %>%
  # Fit linear model for each arctic region
  mutate(mod = list(lm(SA_int ~ lat, data = data))) %>%
  # Find coefficients, summary (AIC, r2, etc) and fit
  summarise(coefs = list(tidy(mod)),
            summary_mod = list(glance(mod)),
            fit = list(augment(mod, interval = "confidence"))) %>%
  ungroup()

# Extract coefficients
LM_area_coefs <- LM_area %>%
  dplyr::select(IHO, coefs) %>% 
  unnest(cols = c(coefs)) %>%
  mutate(term = if_else(term == "(Intercept)", "itcpt", "lat")) %>%
  rename(est = estimate, 
         sd = std.error,
         p_val = p.value) %>%
  pivot_wider(names_from = term, values_from = c(est, sd, statistic, p_val))

# Extract fit
LM_area_fit <- LM_area %>%
  dplyr::select(IHO, fit) %>% 
  unnest(cols = c(fit))

# Combine data
LM_area_res <- left_join(LM_area_coefs, LM_area_fit)
```

```
## Joining, by = "IHO"
```

Plot regressions. We can see that mesopelagic backscatter decreases with increasing latitude in all regions.


```r
LM_area_res %>%
  ggplot() + 
  geom_point(aes(x = lat, y = SA_int, col = IHO)) +
  # geom_abline(aes(slope = est_lat, intercept = est_itcpt, col = IHO)) + 
  geom_line(aes(x = lat, y = .fitted, col = IHO), size = 0.8) + 
  geom_ribbon(aes(x = lat, ymin = .lower, ymax = .upper, group = IHO), alpha = 0.1)
```

![](PanArctic_DSL_statistics_files/figure-html/lm-plot-area-1.png)<!-- -->

Save data


```r
# Save data
save(LM_area_coefs, LM_area_fit, LM_area_res, 
     file = "data/statistics/LM_results.RData")
```

# HGAM - Gamma NASC

## Data preparation

I select data at 318 m depth because it fits well with the DSL diurnal centre of mass. It also prevent from removing too much data like the deeper depth levels would do.


```r
SA_df <- stat_laea %>%
  filter(depth == 318) %>% 
  dplyr::select(year, xc, yc, IHO_area, NASC_int, SA_int, velocity, chl) %>%
  group_by(IHO_area) %>%
  mutate(year = factor(year)) %>%
  ungroup() %>%
  rename(IHO = IHO_area, 
         v = velocity) 
```

Plot data.


```r
plot_grid(SA_df %>%
            ggplot(aes(x = v, y = NASC_int, col = IHO)) +
            geom_point() +
            geom_smooth(method = "gam", se = F, col = "grey20") +
            facet_wrap(~ IHO, scales = "free") +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = chl, y = NASC_int, col = IHO)) +
            geom_point() +
            geom_smooth(method = "gam", se = F, col = "grey20") +
            facet_wrap(~ IHO, scales = "free") +
            theme(legend.position = "none"),
          ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/plot-data-area-1.png)<!-- -->

Check spearman correlations.


```r
# Compute Spearman correlation matrix
corr <- SA_df %>% 
  dplyr::select(-year, -xc, -yc, -IHO) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr, type = "lower", lab = T, title = "corr")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-1.png)<!-- -->

I do not expect regional functional responses to derive from a global functional responses because I expect differences in species assemblages and environmental conditions from one region of the Arctic to the other. So I fit models S and I following the nomenclature proposed by Pedersen et al. 2019.

## Model fitting

In model S the random intercept is already included in the `s(bs = "fs")` term.


```r
# Model G
GAM_G <- gam(NASC_int ~ 
               # First order effects
               s(v, k = 5, bs = "tp") + # Global smoothers
               s(chl, k = 5, bs = "tp") + # Global smoothers
               s(IHO, bs = "re") + # Random effect
               s(year, bs = "re") + # Random effect
               te(xc, yc, by = year, k = 5) +
               # Second order effects
               s(IHO, year, bs = "re"), # Random slope
            data = SA_df, family = Gamma(link = "log"), method = "ML")
# Model S
GAM_S <- gam(NASC_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                te(xc, yc, by = year, k = 5) +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "ML")
# Model I
GAM_I <- gam(NASC_int ~
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                te(xc, yc, by = year, k = 5) +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, by = IHO, k = 5, bs = "tp") +
                s(v, by = IHO, k = 5, bs = "tp"),
             data = SA_df, family = Gamma(link = "log"), method = "ML")
```

Check model metrics


```r
# Summary metrics
GAM_AIC <- AIC(GAM_G, GAM_S, GAM_I) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAM_G", "GAM_S", "GAM_I"),
           reml = round(c(GAM_G$gcv.ubre,
                          GAM_S$gcv.ubre, 
                          GAM_I$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAM_G$deviance / GAM_G$null.deviance)) * 100,
             (1 - (GAM_S$deviance / GAM_S$null.deviance)) * 100,
             (1 - (GAM_I$deviance / GAM_I$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAM_G)$r.sq,
                        summary(GAM_S)$r.sq, 
                        summary(GAM_I)$r.sq),
                      2)) %>%
  full_join(., GAM_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(AIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-d24d7d755b66a89c325b" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d24d7d755b66a89c325b">{"x":{"filter":"none","vertical":false,"data":[["GAM_G","GAM_I","GAM_S"],[36.062,32.211,31.597],[77.45,74.41,73.56],[0.22,0.19,0.17],[1147.88,1127.53,1142.16],[2249.56,2269.546,2276.319],[0,19.99,26.76],[0.99995,5e-05,0]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Model G summary.


```r
summary(GAM_G)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(v, k = 5, bs = "tp") + s(chl, k = 5, bs = "tp") + 
##     s(IHO, bs = "re") + s(year, bs = "re") + te(xc, yc, by = year, 
##     k = 5) + s(IHO, year, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.2485     0.1827   23.25   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df      F p-value    
## s(v)               1.000e+00  1.000  0.743   0.390    
## s(chl)             1.000e+00  1.000  0.032   0.857    
## s(IHO)             3.272e-05  4.000  0.000   0.630    
## s(year)            5.474e-05  2.000  0.000   0.530    
## te(xc,yc):year2015 5.774e+00  6.679 10.385  <2e-16 ***
## te(xc,yc):year2016 8.378e+00  9.590  8.246  <2e-16 ***
## te(xc,yc):year2017 1.444e+01 15.792 13.451  <2e-16 ***
## s(IHO,year)        2.460e-04 14.000  0.000   0.883    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.218   Deviance explained = 77.4%
## -ML = 1147.9  Scale est. = 1.3846    n = 193
```

Model S summary.


```r
summary(GAM_S)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + s(IHO, bs = "re") + te(xc, yc, 
##     by = year, k = 5) + s(year, IHO, bs = "re") + s(chl, IHO, 
##     bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.4459     0.2268    19.6   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df     F  p-value    
## s(year)            1.215e+00  2.000 1.679  0.05510 .  
## s(IHO)             2.320e-05  4.000 0.000  0.46136    
## te(xc,yc):year2015 3.000e+00  3.000 5.042  0.00228 ** 
## te(xc,yc):year2016 3.000e+00  3.000 7.900 5.89e-05 ***
## te(xc,yc):year2017 4.122e+00  4.585 8.380 2.64e-06 ***
## s(year,IHO)        4.900e-05 14.000 0.000  0.81495    
## s(chl,IHO)         1.179e+01 24.000 5.278  < 2e-16 ***
## s(v,IHO)           3.014e+00 24.000 0.292  0.05308 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.169   Deviance explained = 73.6%
## -ML = 1142.2  Scale est. = 1.7166    n = 193
```

Model I summary.


```r
summary(GAM_I)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + s(IHO, bs = "re") + te(xc, yc, 
##     by = year, k = 5) + s(year, IHO, bs = "re") + s(chl, by = IHO, 
##     k = 5, bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.5155     0.1839   24.56   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df      F  p-value    
## s(year)            2.842e-04  2.000  0.000 0.165517    
## s(IHO)             2.667e-06  4.000  0.000 0.373439    
## te(xc,yc):year2015 5.408e+00  6.232  4.482 0.000302 ***
## te(xc,yc):year2016 4.979e+00  5.612  5.319 9.86e-05 ***
## te(xc,yc):year2017 3.000e+00  3.000 11.007 1.73e-06 ***
## s(year,IHO)        6.271e-05 14.000  0.000 0.438548    
## s(chl):IHOWAO_BF   1.000e+00  1.000  0.028 0.867572    
## s(chl):IHOCAA      1.000e+00  1.000  0.085 0.771152    
## s(chl):IHOBB       1.000e+00  1.000  0.525 0.469540    
## s(chl):IHODS       1.000e+00  1.000  8.498 0.004052 ** 
## s(chl):IHOEAO      3.942e+00  3.998 29.560  < 2e-16 ***
## s(v):IHOWAO_BF     1.000e+00  1.000  2.470 0.117936    
## s(v):IHOCAA        1.000e+00  1.000  2.204 0.139598    
## s(v):IHOBB         2.883e+00  3.369  4.137 0.004550 ** 
## s(v):IHODS         1.000e+00  1.000  8.970 0.003170 ** 
## s(v):IHOEAO        1.000e+00  1.000  0.954 0.330198    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.192   Deviance explained = 74.4%
## -ML = 1127.5  Scale est. = 1.6457    n = 193
```

## Model checking

### Basis size and residual distribution

First I check the basis size k. `k-indexes` are \> 1 or close to 1 so the basis size is large enough. The residual plot look good too.


```r
appraise(GAM_S, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_S)
```

```
##                    k'          edf   k-index p-value
## s(year)             3 1.215169e+00        NA      NA
## s(IHO)              5 2.320265e-05        NA      NA
## te(xc,yc):year2015 24 3.000214e+00 0.7310080  0.0025
## te(xc,yc):year2016 24 3.000256e+00 0.7417886  0.0000
## te(xc,yc):year2017 24 4.122050e+00 0.7442888  0.0050
## s(year,IHO)        15 4.900103e-05        NA      NA
## s(chl,IHO)         25 1.178885e+01 0.9442870  0.8500
## s(v,IHO)           25 3.013697e+00 1.0095995  0.9850
```


```r
appraise(GAM_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_I)
```

```
##                    k'          edf   k-index p-value
## s(year)             3 2.842306e-04        NA      NA
## s(IHO)              5 2.667030e-06        NA      NA
## te(xc,yc):year2015 24 5.407951e+00 0.7881492  0.0550
## te(xc,yc):year2016 24 4.979189e+00 0.7798776  0.0350
## te(xc,yc):year2017 24 3.000055e+00 0.8020384  0.0525
## s(year,IHO)        15 6.270713e-05        NA      NA
## s(chl):IHOWAO_BF    4 1.000003e+00 0.9370173  0.7850
## s(chl):IHOCAA       4 1.000001e+00 0.9370173  0.8200
## s(chl):IHOBB        4 1.000041e+00 0.9370173  0.8225
## s(chl):IHODS        4 1.000003e+00 0.9370173  0.8375
## s(chl):IHOEAO       4 3.941687e+00 0.9370173  0.8075
## s(v):IHOWAO_BF      4 1.000012e+00 1.0024945  0.9825
## s(v):IHOCAA         4 1.000010e+00 1.0024945  0.9725
## s(v):IHOBB          4 2.882656e+00 1.0024945  0.9750
## s(v):IHODS          4 1.000005e+00 1.0024945  0.9825
## s(v):IHOEAO         4 1.000025e+00 1.0024945  0.9600
```

### Resiudals against covariates


```r
GAM_S_resid <- bind_cols(SA_df, residuals.gam(GAM_S)) %>%
  rename(resid = "...9")
plot_grid(GAM_S_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_S_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_S_resid %>%
            ggplot(aes(x = xc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "xc") +
            geom_point(), 
          GAM_S_resid %>%
            ggplot(aes(x = yc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "yc") +
            geom_point(), 
          GAM_S_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_S_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS-residuals-covariates-1.png)<!-- -->


```r
GAM_I_resid <- bind_cols(SA_df, residuals.gam(GAM_I)) %>%
  rename(resid = "...9")
plot_grid(GAM_I_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_I_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_I_resid %>%
            ggplot(aes(x = xc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "xc") +
            geom_point(), 
          GAM_I_resid %>%
            ggplot(aes(x = yc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "yc") +
            geom_point(), 
          GAM_I_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_I_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI-residuals-covariates-1.png)<!-- -->

### Autotcorrelation


```r
par(mfrow = c(1, 2))
acf(resid(GAM_S), lag.max = 20, main = "ACF")
pacf(resid(GAM_S), lag.max = 20, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS-acf-pacf-1.png)<!-- -->


```r
par(mfrow = c(1, 2))
acf(resid(GAM_I), lag.max = 20, main = "ACF")
pacf(resid(GAM_I), lag.max = 20, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI-acf-pacf-1.png)<!-- -->

## Term selection

I turn on the double penalty (`select = TRUE`) and check the covariates that has been shrunk.

### Model S


```r
# Model S
GAM_S_p <- gam(NASC_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                te(xc, yc, by = year, k = 5) +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "REML", 
             select = T)
summary(GAM_S_p)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + s(IHO, bs = "re") + te(xc, yc, 
##     by = year, k = 5) + s(year, IHO, bs = "re") + s(chl, IHO, 
##     bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    4.517      0.381   11.86   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df     F  p-value    
## s(year)            1.713e+00      2 4.106  0.00182 ** 
## s(IHO)             2.714e-04      4 0.000  0.69611    
## te(xc,yc):year2015 5.416e+00     24 1.217 2.81e-05 ***
## te(xc,yc):year2016 5.260e+00     24 1.231 2.08e-06 ***
## te(xc,yc):year2017 7.317e+00     24 1.602 9.21e-07 ***
## s(year,IHO)        3.732e-04     14 0.000  0.50413    
## s(chl,IHO)         1.052e+01     24 4.136  < 2e-16 ***
## s(v,IHO)           1.647e+00     24 0.103  0.13498    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.17   Deviance explained = 75.8%
## -REML =   1155  Scale est. = 1.5823    n = 193
```

`s(IHO)` and `s(year,IHO)` can be dropped. Refit model without those terms.


```r
# Model S
GAM_S2 <- gam(NASC_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                # s(IHO, bs = "re") +
                te(xc, yc, by = year, k = 5) +
                # Second order random effects
                # s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "ML")
summary(GAM_S2)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + te(xc, yc, by = year, k = 5) + 
##     s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.4458     0.2268    19.6   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                       edf Ref.df     F  p-value    
## s(year)             1.215  2.000 1.679  0.05510 .  
## te(xc,yc):year2015  3.000  3.001 5.041  0.00228 ** 
## te(xc,yc):year2016  3.000  3.000 7.900 5.89e-05 ***
## te(xc,yc):year2017  4.122  4.585 8.380 2.64e-06 ***
## s(chl,IHO)         11.789 24.000 5.278  < 2e-16 ***
## s(v,IHO)            3.015 24.000 0.292  0.05309 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.169   Deviance explained = 73.6%
## -ML = 1142.2  Scale est. = 1.7166    n = 193
```

Compare models.


```r
# Summary metrics
GAMS_AIC <- AIC(GAM_S, GAM_S2) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAM_S", "GAM_S2"),
           reml = round(c(GAM_S$gcv.ubre,
                          GAM_S2$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAM_S$deviance / GAM_S$null.deviance)) * 100,
             (1 - (GAM_S2$deviance / GAM_S2$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAM_S)$r.sq,
                        summary(GAM_S2)$r.sq),
                      2)) %>%
  full_join(., GAMS_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(AIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-7e92b08fca7c3901fd4b" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-7e92b08fca7c3901fd4b">{"x":{"filter":"none","vertical":false,"data":[["GAM_S","GAM_S2"],[31.597,31.599],[73.56,73.56],[0.17,0.17],[1142.16,1142.16],[2276.319,2276.321],[0,0],[0.50025,0.49975]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAM_S3 <- gam(NASC_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                # s(IHO, bs = "re") +
                te(xc, yc, by = year, k = 5) +
                # Second order random effects
                # s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
              data = SA_df, family = Gamma(link = "log"), method = "REML")
summary(GAM_S3)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + te(xc, yc, by = year, k = 5) + 
##     s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.4504     0.2915   15.27   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df     F  p-value    
## s(year)            8.276e-05  2.000 0.000 0.447663    
## te(xc,yc):year2015 8.741e+00 10.008 3.238 0.000796 ***
## te(xc,yc):year2016 6.395e+00  7.076 5.419 1.76e-05 ***
## te(xc,yc):year2017 5.071e+00  5.816 8.331  < 2e-16 ***
## s(chl,IHO)         1.073e+01 24.000 4.965  < 2e-16 ***
## s(v,IHO)           3.634e+00 24.000 0.394 0.006962 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.16   Deviance explained = 77.5%
## -REML = 1125.1  Scale est. = 1.4403    n = 193
```

Check residuals.


```r
appraise(GAM_S3)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS3-residuals-covariates-1.png)<!-- -->

```r
GAM_S3_resid <- bind_cols(SA_df, residuals.gam(GAM_S3)) %>%
  rename(resid = "...9")
plot_grid(GAM_S3_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_S3_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_S3_resid %>%
            ggplot(aes(x = xc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "xc") +
            geom_point(), 
          GAM_S3_resid %>%
            ggplot(aes(x = yc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "yc") +
            geom_point(), 
          GAM_S3_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_S3_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS3-residuals-covariates-2.png)<!-- -->

### Model I


```r
# Model I
GAM_I_p <- gam(NASC_int ~
                 # First order effects
                 # s(chl, bs = "tp", k = 5) +
                 # s(v, bs = "tp", k = 5) +
                 s(year, bs = "re") +
                 s(IHO, bs = "re") +
                 te(xc, yc, by = year, k = 5) +
                 # Second order random effects
                 s(year, IHO, bs = "re") +
                 # Second order functional effects
                 s(chl, by = IHO, k = 5, bs = "tp") +
                 s(v, by = IHO, k = 5, bs = "tp"),
               data = SA_df, family = Gamma(link = "log"), method = "REML",
               select = T)
summary(GAM_I_p)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + s(IHO, bs = "re") + te(xc, yc, 
##     by = year, k = 5) + s(year, IHO, bs = "re") + s(chl, by = IHO, 
##     k = 5, bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.9043     0.1336   36.71   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df      F  p-value    
## s(year)            5.382e-05      2  0.000 0.490526    
## s(IHO)             1.417e-01      4  0.028 0.214304    
## te(xc,yc):year2015 8.453e+00     24  4.578  < 2e-16 ***
## te(xc,yc):year2016 6.438e+00     24  2.319  < 2e-16 ***
## te(xc,yc):year2017 6.476e+00     23  1.181 0.000415 ***
## s(year,IHO)        2.103e-04     14  0.000 0.510483    
## s(chl):IHOWAO_BF   1.156e-03      4  0.000 0.386585    
## s(chl):IHOCAA      3.285e-04      4  0.000 0.372492    
## s(chl):IHOBB       5.820e-04      4  0.000 0.452178    
## s(chl):IHODS       1.410e-04      4  0.000 0.910677    
## s(chl):IHOEAO      3.912e+00      4 28.274  < 2e-16 ***
## s(v):IHOWAO_BF     1.692e-04      4  0.000 0.937982    
## s(v):IHOCAA        4.365e-04      4  0.000 0.490336    
## s(v):IHOBB         8.814e-01      4  1.518 0.007757 ** 
## s(v):IHODS         7.511e-01      4  0.621 0.063632 .  
## s(v):IHOEAO        2.705e-04      4  0.000 0.681137    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.19   Deviance explained = 76.4%
## -REML = 1148.7  Scale est. = 1.6139    n = 193
```

`s(year)`, `s(IHO)`, and `s(year,IHO)` can be dropped. Refit model without those terms.


```r
# Model I
GAM_I2 <- gam(NASC_int ~
                 # First order effects
                 # s(chl, bs = "tp", k = 5) +
                 # s(v, bs = "tp", k = 5) +
                 # s(year, bs = "re") +
                 # s(IHO, bs = "re") +
                 te(xc, yc, by = year, k = 5) +
                 # Second order random effects
                 # s(year, IHO, bs = "re") +
                 # Second order functional effects
                 s(chl, by = IHO, k = 5, bs = "tp") +
                 s(v, by = IHO, k = 5, bs = "tp"),
               data = SA_df, family = Gamma(link = "log"), method = "ML")
summary(GAM_I2)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ te(xc, yc, by = year, k = 5) + s(chl, by = IHO, k = 5, 
##     bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.5155     0.1838   24.56   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                      edf Ref.df      F  p-value    
## te(xc,yc):year2015 5.408  6.232  4.482 0.000302 ***
## te(xc,yc):year2016 4.979  5.612  5.320 9.86e-05 ***
## te(xc,yc):year2017 3.000  3.000 11.007 1.73e-06 ***
## s(chl):IHOWAO_BF   1.000  1.000  0.028 0.867615    
## s(chl):IHOCAA      1.000  1.000  0.085 0.771175    
## s(chl):IHOBB       1.000  1.000  0.525 0.469530    
## s(chl):IHODS       1.000  1.000  8.498 0.004052 ** 
## s(chl):IHOEAO      3.942  3.998 29.559  < 2e-16 ***
## s(v):IHOWAO_BF     1.000  1.000  2.471 0.117928    
## s(v):IHOCAA        1.000  1.000  2.204 0.139597    
## s(v):IHOBB         2.883  3.369  4.137 0.004550 ** 
## s(v):IHODS         1.000  1.000  8.970 0.003170 ** 
## s(v):IHOEAO        1.000  1.000  0.954 0.330194    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.192   Deviance explained = 74.4%
## -ML = 1127.5  Scale est. = 1.6457    n = 193
```
Compare models


```r
# Summary metrics
GAMI_AIC <- AIC(GAM_I, GAM_I2) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAM_I", "GAM_I2"),
           reml = round(c(GAM_I$gcv.ubre,
                          GAM_I2$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAM_I$deviance / GAM_I$null.deviance)) * 100,
             (1 - (GAM_I2$deviance / GAM_I2$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAM_I)$r.sq,
                        summary(GAM_I2)$r.sq),
                      2)) %>%
  full_join(., GAMI_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(AIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-067e3094c50318a08661" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-067e3094c50318a08661">{"x":{"filter":"none","vertical":false,"data":[["GAM_I","GAM_I2"],[32.211,32.21],[74.41,74.41],[0.19,0.19],[1127.53,1127.53],[2269.546,2269.546],[0,0],[0.5,0.5]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAM_I3 <- gam(NASC_int ~ 
                 # First order effects
                 # s(chl, bs = "tp", k = 5) +
                 # s(v, bs = "tp", k = 5) +
                 # s(year, bs = "re") +
                 # s(IHO, bs = "re") +
                 te(xc, yc, by = year, k = 5) +
                 # Second order random effects
                 # s(year, IHO, bs = "re") +
                 # Second order functional effects
                 s(chl, by = IHO, k = 5, bs = "tp") +
                 s(v, by = IHO, k = 5, bs = "tp"),
              data = SA_df, family = Gamma(link = "log"), method = "REML")
summary(GAM_I3)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ te(xc, yc, by = year, k = 5) + s(chl, by = IHO, k = 5, 
##     bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.7818     0.1987   24.07   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                      edf Ref.df      F  p-value    
## te(xc,yc):year2015 8.740 10.150  4.906 3.66e-06 ***
## te(xc,yc):year2016 7.532  8.701  6.121 9.61e-07 ***
## te(xc,yc):year2017 4.175  4.707 11.002  < 2e-16 ***
## s(chl):IHOWAO_BF   1.000  1.001  0.611  0.43521    
## s(chl):IHOCAA      2.037  2.312  2.968  0.04873 *  
## s(chl):IHOBB       1.078  1.145  0.101  0.78141    
## s(chl):IHODS       1.000  1.000  0.956  0.32985    
## s(chl):IHOEAO      3.928  3.995 35.066  < 2e-16 ***
## s(v):IHOWAO_BF     1.000  1.000  2.051  0.15410    
## s(v):IHOCAA        1.748  2.111  1.497  0.27587    
## s(v):IHOBB         3.115  3.567  4.936  0.00109 ** 
## s(v):IHODS         1.000  1.001  4.022  0.04660 *  
## s(v):IHOEAO        1.000  1.000  4.281  0.04020 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.177   Deviance explained = 78.9%
## -REML = 1108.4  Scale est. = 1.2519    n = 193
```

Check residuals.


```r
appraise(GAM_I3)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI3-residuals-covariates-1.png)<!-- -->

```r
GAM_I3_resid <- bind_cols(SA_df, residuals.gam(GAM_I3)) %>%
  rename(resid = "...9")
plot_grid(GAM_I3_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_I3_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_I3_resid %>%
            ggplot(aes(x = xc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "xc") +
            geom_point(), 
          GAM_I3_resid %>%
            ggplot(aes(x = yc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "yc") +
            geom_point(), 
          GAM_I3_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_I3_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI3-residuals-covariates-2.png)<!-- -->

## Model predictions

I predict the model to get the smooths in the response scale (NASC). I need to create a specific data frame for this with the factors of interest `IHO` and `year`, and covariates `v` and `chl`.


```r
# Find median, min, and max of SA_df
val <- SA_df %>%
  dplyr::select(IHO, year, xc, yc, v, chl) %>%
  group_by(IHO) %>%
  summarise(median_xc = median(xc),
            median_yc = median(yc),
            min_v = min(v),
            max_v = max(v),
            median_v = median(v),
            min_chl = min(chl),
            max_chl = max(chl),
            median_chl = median(chl))
# Resolution for predictions
res = 50
# IHO areas
IHO_area <- levels(SA_df$IHO)
# Empty data frame that we populate with new values
SA_new <- data.frame()

for (i in IHO_area) {
  # Select data from which we build the 
  val_tmp <- val %>% 
    filter(IHO == i)
  # Temporary data frame for each region
  SA_tmp <- data.frame(
    IHO = rep(rep(rep(i, res), 3), 2),
    year = rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2),
    var = c(rep("v", res * 3), rep("chl", res * 3)),
    # Velocity
    xc = rep(rep(rep(subset(val, IHO == i)$median_xc, res), 3), 2),
    yc = rep(rep(rep(subset(val, IHO == i)$median_yc, res), 3), 2),
    v = c(rep(seq(subset(val, IHO == i)$min_v,
                  subset(val, IHO == i)$max_v,
                  length.out = res), 3),
          rep(subset(val, IHO == "WAO")$median_v, res * 3)),
    # Chlorophyll
    chl = c(rep(subset(val, IHO == i)$median_chl, res * 3),
            rep(seq(subset(val, IHO == i)$min_chl, 
                    subset(val, IHO == i)$max_chl,
                    length.out = res), 3)))
  # Append data
  SA_new <- bind_rows(SA_new, SA_tmp)
}
# Reorder factors
SA_new <- SA_new %>%
  mutate(IHO = factor(IHO, levels = c("WAO_BF", "CAA", "BB", "DS", "EAO")))
# Remove unused variables
rm(val_tmp, SA_tmp, IHO_area, res)
```

Get GAM predictions (`GAM_S` and `GAM_I`) for the new data. I then calculate the average functional response for all years.


```r
ilink_S <- family(GAM_S3)$linkinv # Get link function
pred_S <- predict(GAM_S3, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilink_S(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilink_S(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilink_S(fit)) %>% # Calculate fit
  group_by(IHO, var, v, chl) %>%
  summarise(NASC_fit = mean(fitted), # Calculate average functional response
            NASC_lwr_ci = mean(lwr_ci),
            NASC_upr_ci = mean(upr_ci)) %>%
  mutate(SA_fit = 10 * log10(NASC_fit), # Convert to SA
         SA_lwr_ci = 10 * log10(NASC_lwr_ci),
         SA_upr_ci = 10 * log10(NASC_upr_ci))
```


```r
ilink_I <- family(GAM_I3)$linkinv # Get link function
pred_I <- predict(GAM_I3, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilink_I(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilink_I(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilink_I(fit)) %>% # Calculate fit
  group_by(IHO, var, v, chl) %>%
  summarise(NASC_fit = mean(fitted), # Calculate average functional response
            NASC_lwr_ci = mean(lwr_ci),
            NASC_upr_ci = mean(upr_ci)) %>%
  mutate(SA_fit = 10 * log10(NASC_fit), # Convert to SA
         SA_lwr_ci = 10 * log10(NASC_lwr_ci),
         SA_upr_ci = 10 * log10(NASC_upr_ci))
```


```r
# Save data
save(pred_S, pred_I, GAM_S3, GAM_I3, file = "data/statistics/GAM_results.RData")
```

## Model visualization

Plot to see if predictions worked well.


```r
plot_grid(
  pred_S %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = NASC_fit)) +
    geom_ribbon(aes(x = v, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model S"),
  pred_I %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = NASC_fit)) +
    geom_ribbon(aes(x = v, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model I"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/velo-ggplot-1.png)<!-- -->


```r
plot_grid(
  pred_S %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = NASC_fit)) +
    geom_ribbon(aes(x = chl, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model S"),
  pred_I %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = NASC_fit)) +
    geom_ribbon(aes(x = chl, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model I"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/chl-ggplot-1.png)<!-- -->


```r
GAMS_velo <- pred_S %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO),
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "Model S",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAMS_chl <- pred_S %>%
  filter(var == "chl") %>%
  ggplot() +
  geom_line(aes(x = chl, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), 
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  labs(x = expression("Chlorophyll (mg m"^-3*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "none")
# Combine plots
plot_grid(GAMS_velo, GAMS_chl, ncol = 1, rel_heights = c(1,0.8))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS-pretty-ggplot-1.png)<!-- -->

```r
# Remove unused variables
rm(GAMS_velo, GAMS_chl)
```


```r
GAMI_velo <- pred_I %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), 
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "Model I",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAMI_chl <- pred_I %>%
  filter(var == "chl") %>%
  ggplot() +
  geom_line(aes(x = chl, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO),
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  labs(x = expression("Chlorophyll (mg m"^-3*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "none")
# Combine plots
plot_grid(GAMI_velo, GAMI_chl, ncol = 1, rel_heights = c(1,0.8))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI-pretty-ggplot-1.png)<!-- -->

```r
# Remove unused variables
rm(GAMI_velo, GAMI_chl)
```


# HGAM - Gaussian SA

## Model fitting

In model S the random intercept is already included in the `s(bs = "fs")` term.


```r
# Model G
GAMSA_G <- gam(SA_int ~ 
                 # First order effects
                 s(v, k = 5, bs = "tp") + # Global smoothers
                 s(chl, k = 5, bs = "tp") + # Global smoothers
                 s(IHO, bs = "re") + # Random effect
                 s(year, bs = "re") + # Random effect
                 te(xc, yc, by = year, k = 5) +
                 # Second order effects
                 s(IHO, year, bs = "re"), # Random slope
               data = SA_df, family = "gaussian", method = "ML")
# Model S
GAMSA_S <- gam(SA_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                te(xc, yc, by = year, k = 5) +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
               data = SA_df, family = "gaussian", method = "ML")
# Model I
GAMSA_I <- gam(SA_int ~
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                te(xc, yc, by = year, k = 5) +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, by = IHO, k = 5, bs = "tp") +
                s(v, by = IHO, k = 5, bs = "tp"),
             data = SA_df, family = "gaussian", method = "ML")
```

Check model metrics


```r
# Summary metrics
GAMSA_AIC <- AIC(GAMSA_G, GAMSA_S, GAMSA_I) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAMSA_G", "GAMSA_S", "GAMSA_I"),
           reml = round(c(GAMSA_G$gcv.ubre,
                          GAMSA_S$gcv.ubre, 
                          GAMSA_I$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAMSA_G$deviance / GAMSA_G$null.deviance)) * 100,
             (1 - (GAMSA_S$deviance / GAMSA_S$null.deviance)) * 100,
             (1 - (GAMSA_I$deviance / GAMSA_I$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAMSA_G)$r.sq,
                        summary(GAMSA_S)$r.sq, 
                        summary(GAMSA_I)$r.sq),
                      2)) %>%
  full_join(., GAMSA_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(AIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-aedd47d937920ed646c4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-aedd47d937920ed646c4">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_S","GAMSA_G"],[26.205,27.596,16.896],[51.06,50.73,40.83],[0.44,0.45,0.36],[636.67,648.48,648.48],[1302.264,1306.327,1320.271],[0,4.06,18.01],[0.88397,0.11592,0.00011]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Model G summary.


```r
summary(GAMSA_G)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(v, k = 5, bs = "tp") + s(chl, k = 5, bs = "tp") + 
##     s(IHO, bs = "re") + s(year, bs = "re") + te(xc, yc, by = year, 
##     k = 5) + s(IHO, year, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  17.5167     0.8028   21.82   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df      F  p-value    
## s(v)               1.000e+00  1.000  1.662  0.19900    
## s(chl)             1.000e+00  1.000  9.204  0.00278 ** 
## s(IHO)             2.359e-05  4.000  0.000  0.67473    
## s(year)            8.623e-01  2.000  1.183  0.06448 .  
## te(xc,yc):year2015 3.000e+00  3.000 10.708 2.30e-06 ***
## te(xc,yc):year2016 3.000e+00  3.000  7.993 5.02e-05 ***
## te(xc,yc):year2017 4.678e+00  5.556  6.047 1.85e-05 ***
## s(IHO,year)        1.270e-04 14.000  0.000  0.70697    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.363   Deviance explained = 40.8%
## -ML = 648.48  Scale est. = 49.707    n = 193
```

Model S summary.


```r
summary(GAMSA_S)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(IHO, bs = "re") + te(xc, yc, 
##     by = year, k = 5) + s(year, IHO, bs = "re") + s(chl, IHO, 
##     bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.672      1.088   15.32   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df     F  p-value    
## s(year)            1.323e+00  2.000 2.542 0.018192 *  
## s(IHO)             5.076e-06  4.000 0.000 0.676843    
## te(xc,yc):year2015 3.000e+00  3.000 8.087 4.59e-05 ***
## te(xc,yc):year2016 3.641e+00  3.966 4.955 0.000606 ***
## te(xc,yc):year2017 3.958e+00  4.438 5.912 9.85e-05 ***
## s(year,IHO)        2.052e-05 14.000 0.000 0.699658    
## s(chl,IHO)         8.264e+00 23.000 1.350 0.000133 ***
## s(v,IHO)           1.315e+00 24.000 0.091 0.149941    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.445   Deviance explained = 50.7%
## -ML = 648.48  Scale est. = 43.322    n = 193
```

Model I summary.


```r
summary(GAMSA_I)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(IHO, bs = "re") + te(xc, yc, 
##     by = year, k = 5) + s(year, IHO, bs = "re") + s(chl, by = IHO, 
##     k = 5, bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  16.4840     0.8814    18.7   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df      F  p-value    
## s(year)            8.229e-05  2.000  0.000 0.216653    
## s(IHO)             1.593e-06  4.000  0.000 0.368506    
## te(xc,yc):year2015 3.000e+00  3.000 10.637 2.52e-06 ***
## te(xc,yc):year2016 4.646e+00  5.240  4.549 0.000642 ***
## te(xc,yc):year2017 3.000e+00  3.000  7.426 0.000106 ***
## s(year,IHO)        2.284e-05 14.000  0.000 0.435740    
## s(chl):IHOWAO_BF   1.000e+00  1.000  0.052 0.820121    
## s(chl):IHOCAA      1.000e+00  1.000  0.132 0.717123    
## s(chl):IHOBB       1.000e+00  1.000  1.043 0.308609    
## s(chl):IHODS       1.000e+00  1.000  5.814 0.016973 *  
## s(chl):IHOEAO      3.756e+00  3.965  8.829 5.15e-06 ***
## s(v):IHOWAO_BF     1.000e+00  1.000  2.231 0.137142    
## s(v):IHOCAA        1.000e+00  1.000  1.278 0.259896    
## s(v):IHOBB         1.000e+00  1.000  6.318 0.012888 *  
## s(v):IHODS         1.000e+00  1.000  3.835 0.051839 .  
## s(v):IHOEAO        1.000e+00  1.000  0.021 0.885045    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.443   Deviance explained = 51.1%
## -ML = 636.67  Scale est. = 43.52     n = 193
```

## Model checking

### Basis size and residual distribution

First I check the basis size k. `k-indexes` are \> 1 or close to 1 so the basis size is large enough. The residual plot look good too.


```r
appraise(GAMSA_S, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAMSA_S)
```

```
##                    k'          edf   k-index p-value
## s(year)             3 1.322630e+00        NA      NA
## s(IHO)              5 5.076413e-06        NA      NA
## te(xc,yc):year2015 24 3.000053e+00 0.8880814  0.0450
## te(xc,yc):year2016 24 3.640706e+00 0.8317125  0.0025
## te(xc,yc):year2017 24 3.958002e+00 0.8282778  0.0025
## s(year,IHO)        15 2.052434e-05        NA      NA
## s(chl,IHO)         25 8.264432e+00 1.1088913  0.9400
## s(v,IHO)           25 1.314614e+00 1.0645371  0.8075
```


```r
appraise(GAMSA_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAMSA_I)
```

```
##                    k'          edf   k-index p-value
## s(year)             3 8.228487e-05        NA      NA
## s(IHO)              5 1.593109e-06        NA      NA
## te(xc,yc):year2015 24 3.000010e+00 0.9682036  0.2575
## te(xc,yc):year2016 24 4.646085e+00 0.8601250  0.0100
## te(xc,yc):year2017 24 3.000030e+00 0.8857877  0.0275
## s(year,IHO)        15 2.283484e-05        NA      NA
## s(chl):IHOWAO_BF    4 1.000001e+00 1.1096623  0.9375
## s(chl):IHOCAA       4 1.000001e+00 1.1096623  0.9400
## s(chl):IHOBB        4 1.000005e+00 1.1096623  0.9275
## s(chl):IHODS        4 1.000002e+00 1.1096623  0.9075
## s(chl):IHOEAO       4 3.755849e+00 1.1096623  0.9325
## s(v):IHOWAO_BF      4 1.000007e+00 1.0917812  0.8800
## s(v):IHOCAA         4 1.000003e+00 1.0917812  0.8925
## s(v):IHOBB          4 1.000018e+00 1.0917812  0.8625
## s(v):IHODS          4 1.000004e+00 1.0917812  0.8750
## s(v):IHOEAO         4 1.000008e+00 1.0917812  0.8850
```

### Resiudals against covariates


```r
GAMSA_S_resid <- bind_cols(SA_df, residuals.gam(GAMSA_S)) %>%
  rename(resid = "...9")
plot_grid(GAMSA_S_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAMSA_S_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAMSA_S_resid %>%
            ggplot(aes(x = xc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "xc") +
            geom_point(), 
          GAMSA_S_resid %>%
            ggplot(aes(x = yc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "yc") +
            geom_point(), 
          GAMSA_S_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAMSA_S_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS-residuals-covariates-1.png)<!-- -->


```r
GAMSA_I_resid <- bind_cols(SA_df, residuals.gam(GAMSA_I)) %>%
  rename(resid = "...9")
plot_grid(GAMSA_I_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAMSA_I_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAMSA_I_resid %>%
            ggplot(aes(x = xc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "xc") +
            geom_point(), 
          GAMSA_I_resid %>%
            ggplot(aes(x = yc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "yc") +
            geom_point(), 
          GAMSA_I_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAMSA_I_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-residuals-covariates-1.png)<!-- -->

### Autotcorrelation


```r
par(mfrow = c(1, 2))
acf(resid(GAMSA_S), lag.max = 20, main = "ACF")
pacf(resid(GAMSA_S), lag.max = 20, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS-acf-pacf-1.png)<!-- -->


```r
par(mfrow = c(1, 2))
acf(resid(GAMSA_I), lag.max = 20, main = "ACF")
pacf(resid(GAMSA_I), lag.max = 20, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-acf-pacf-1.png)<!-- -->

## Term selection

I turn on the double penalty (`select = TRUE`) and check the covariates that has been shrunk.

### Model S


```r
# Model S
GAMSA_S_p <- gam(SA_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                te(xc, yc, by = year, k = 5) +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = "gaussian", method = "REML", 
             select = T)
summary(GAMSA_S_p)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(IHO, bs = "re") + te(xc, yc, 
##     by = year, k = 5) + s(year, IHO, bs = "re") + s(chl, IHO, 
##     bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.226      1.896    8.56 6.57e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df     F  p-value    
## s(year)            1.7487378      2 5.645 0.000854 ***
## s(IHO)             0.0002176      4 0.000 0.195671    
## te(xc,yc):year2015 2.5346439     24 0.565 0.000365 ***
## te(xc,yc):year2016 3.6425880     23 0.518 0.001275 ** 
## te(xc,yc):year2017 4.3428539     24 0.839 2.02e-05 ***
## s(year,IHO)        0.0001539     14 0.000 0.530778    
## s(chl,IHO)         8.4335642     24 1.420 4.08e-05 ***
## s(v,IHO)           1.7800930     24 0.139 0.022921 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.449   Deviance explained = 51.4%
## -REML = 661.81  Scale est. = 43.005    n = 193
```

`s(IHO)` and `s(year,IHO)` can be dropped. Refit the model without those terms.


```r
# Model S
GAMSA_S2 <- gam(SA_int ~ 
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  s(year, bs = "re") +
                  # s(IHO, bs = "re") +
                  te(xc, yc, by = year, k = 5) +
                  # Second order random effects
                  # s(year, IHO, bs = "re") +
                  # Second order functional effects
                  s(chl, IHO, bs = "fs", k = 5) +
                  s(v, IHO, bs = "fs", k = 5),
                data = SA_df, family = "gaussian", method = "ML")
summary(GAMSA_S2)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + te(xc, yc, by = year, k = 5) + 
##     s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.672      1.088   15.32   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                      edf Ref.df     F  p-value    
## s(year)            1.323  2.000 2.542 0.018191 *  
## te(xc,yc):year2015 3.000  3.000 8.086 4.59e-05 ***
## te(xc,yc):year2016 3.641  3.966 4.955 0.000606 ***
## te(xc,yc):year2017 3.958  4.438 5.912 9.85e-05 ***
## s(chl,IHO)         8.264 24.000 1.294 0.000133 ***
## s(v,IHO)           1.315 24.000 0.091 0.149949    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.445   Deviance explained = 50.7%
## -ML = 648.48  Scale est. = 43.322    n = 193
```

Compare models.


```r
# Summary metrics
GAMSAS_AIC <- AIC(GAMSA_S, GAMSA_S2) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAMSA_S", "GAMSA_S2"),
           reml = round(c(GAMSA_S$gcv.ubre,
                          GAMSA_S2$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAMSA_S$deviance / GAMSA_S$null.deviance)) * 100,
             (1 - (GAMSA_S2$deviance / GAMSA_S2$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAMSA_S)$r.sq,
                        summary(GAMSA_S2)$r.sq),
                      2)) %>%
  full_join(., GAMSAS_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(AIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-9a0fac7b1fbabae53fd4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-9a0fac7b1fbabae53fd4">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_S","GAMSA_S2"],[27.596,27.596],[50.73,50.73],[0.45,0.45],[648.48,648.48],[1306.327,1306.328],[0,0],[0.50012,0.49988]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAMSA_S3 <- gam(SA_int ~ 
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  s(year, bs = "re") +
                  # s(IHO, bs = "re") +
                  te(xc, yc, by = year, k = 5) +
                  # Second order random effects
                  # s(year, IHO, bs = "re") +
                  # Second order functional effects
                  s(chl, IHO, bs = "fs", k = 5) +
                  s(v, IHO, bs = "fs", k = 5),
                data = SA_df, family = "gaussian", method = "REML")
summary(GAMSA_S3)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + te(xc, yc, by = year, k = 5) + 
##     s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.424      1.252   13.12   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                      edf Ref.df     F  p-value    
## s(year)            1.377  2.000 2.112 0.024524 *  
## te(xc,yc):year2015 3.000  3.001 7.844  6.3e-05 ***
## te(xc,yc):year2016 6.254  7.317 3.647 0.001134 ** 
## te(xc,yc):year2017 4.325  4.941 5.445 0.000119 ***
## s(chl,IHO)         8.260 24.000 1.304 0.000118 ***
## s(v,IHO)           1.454 24.000 0.100 0.141636    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.463   Deviance explained = 53.2%
## -REML = 618.65  Scale est. = 41.907    n = 193
```

Check residuals.


```r
appraise(GAMSA_S3)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS3-residuals-covariates-1.png)<!-- -->

```r
GAMSA_S3_resid <- bind_cols(SA_df, residuals.gam(GAMSA_S3)) %>%
  rename(resid = "...9")
plot_grid(GAMSA_S3_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAMSA_S3_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAMSA_S3_resid %>%
            ggplot(aes(x = xc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "xc") +
            geom_point(), 
          GAMSA_S3_resid %>%
            ggplot(aes(x = yc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "yc") +
            geom_point(), 
          GAMSA_S3_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAMSA_S3_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS3-residuals-covariates-2.png)<!-- -->

### Model I


```r
# Model I
GAMSA_I_p <- gam(SA_int ~
                   # First order effects
                   # s(chl, bs = "tp", k = 5) +
                   # s(v, bs = "tp", k = 5) +
                   s(year, bs = "re") +
                   s(IHO, bs = "re") +
                   te(xc, yc, by = year, k = 5) +
                   # Second order random effects
                   s(year, IHO, bs = "re") +
                   # Second order functional effects
                   s(chl, by = IHO, k = 5, bs = "tp") +
                   s(v, by = IHO, k = 5, bs = "tp"),
                data = SA_df, family = "gaussian", method = "REML",
                select = T)
summary(GAMSA_I_p)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(IHO, bs = "re") + te(xc, yc, 
##     by = year, k = 5) + s(year, IHO, bs = "re") + s(chl, by = IHO, 
##     k = 5, bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   17.565      1.661   10.57   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df     F  p-value    
## s(year)            1.748e+00      2 5.380  0.00108 ** 
## s(IHO)             2.312e-05      4 0.000  0.32062    
## te(xc,yc):year2015 4.095e+00     24 1.314 1.79e-05 ***
## te(xc,yc):year2016 4.549e+00     24 0.687  0.00210 ** 
## te(xc,yc):year2017 3.983e+00     24 1.100 1.20e-05 ***
## s(year,IHO)        3.134e-05     14 0.000  0.63537    
## s(chl):IHOWAO_BF   8.132e-01      4 1.246  0.00755 ** 
## s(chl):IHOCAA      3.124e-05      4 0.000  0.49773    
## s(chl):IHOBB       6.387e-05      4 0.000  0.46524    
## s(chl):IHODS       2.135e-05      4 0.000  0.64873    
## s(chl):IHOEAO      3.699e+00      4 8.244 2.18e-06 ***
## s(v):IHOWAO_BF     2.471e-05      4 0.000  0.93766    
## s(v):IHOCAA        2.623e-05      4 0.000  0.52777    
## s(v):IHOBB         8.206e-01      4 1.146  0.01694 *  
## s(v):IHODS         5.479e-01      4 0.248  0.15781    
## s(v):IHOEAO        2.922e-05      4 0.000  0.93533    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.47   Deviance explained = 52.6%
## -REML = 658.45  Scale est. = 41.354    n = 193
```

`s(IHO)`, and `s(year,IHO)` can be dropped. Refit model without those terms.


```r
# Model I
GAMSA_I2 <- gam(SA_int ~
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  s(year, bs = "re") +
                  # s(IHO, bs = "re") +
                  te(xc, yc, by = year, k = 5) +
                  # Second order random effects
                  # s(year, IHO, bs = "re") +
                  # Second order functional effects
                  s(chl, by = IHO, k = 5, bs = "tp") +
                  s(v, by = IHO, k = 5, bs = "tp"),
                data = SA_df, family = "gaussian", method = "ML")
summary(GAMSA_I2)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + te(xc, yc, by = year, k = 5) + 
##     s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, k = 5, 
##     bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  16.4841     0.8814    18.7   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                          edf Ref.df      F  p-value    
## s(year)            0.0001465  2.000  0.000 0.216649    
## te(xc,yc):year2015 3.0000133  3.000 10.637 2.52e-06 ***
## te(xc,yc):year2016 4.6460880  5.240  4.549 0.000642 ***
## te(xc,yc):year2017 3.0000927  3.000  7.425 0.000106 ***
## s(chl):IHOWAO_BF   1.0000011  1.000  0.052 0.820115    
## s(chl):IHOCAA      1.0000010  1.000  0.132 0.717117    
## s(chl):IHOBB       1.0000058  1.000  1.043 0.308601    
## s(chl):IHODS       1.0000019  1.000  5.813 0.016976 *  
## s(chl):IHOEAO      3.7558492  3.965  8.829 5.15e-06 ***
## s(v):IHOWAO_BF     1.0000075  1.000  2.231 0.137145    
## s(v):IHOCAA        1.0000024  1.000  1.278 0.259897    
## s(v):IHOBB         1.0000102  1.000  6.318 0.012888 *  
## s(v):IHODS         1.0000042  1.000  3.835 0.051840 .  
## s(v):IHOEAO        1.0000077  1.000  0.021 0.885051    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.443   Deviance explained = 51.1%
## -ML = 636.67  Scale est. = 43.52     n = 193
```

Compare models


```r
# Summary metrics
GAMSAI_AIC <- AIC(GAMSA_I, GAMSA_I2) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAMSA_I", "GAMSA_I2"),
           reml = round(c(GAMSA_I$gcv.ubre,
                          GAMSA_I2$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAMSA_I$deviance / GAMSA_I$null.deviance)) * 100,
             (1 - (GAMSA_I2$deviance / GAMSA_I2$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAMSA_I)$r.sq,
                        summary(GAMSA_I2)$r.sq),
                      2)) %>%
  full_join(., GAMSAI_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(AIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-47be6de6478a296ea568" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-47be6de6478a296ea568">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_I2"],[26.205,26.205],[51.06,51.06],[0.44,0.44],[636.67,636.67],[1302.264,1302.264],[0,0],[0.5,0.5]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAMSA_I3 <- gam(SA_int ~
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  s(year, bs = "re") +
                  # s(IHO, bs = "re") +
                  te(xc, yc, by = year, k = 5) +
                  # Second order random effects
                  # s(year, IHO, bs = "re") +
                  # Second order functional effects
                  s(chl, by = IHO, k = 5, bs = "tp") +
                  s(v, by = IHO, k = 5, bs = "tp"),
                data = SA_df, family = "gaussian", method = "REML")
summary(GAMSA_I3)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + te(xc, yc, by = year, k = 5) + 
##     s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, k = 5, 
##     bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   17.616      1.139   15.46   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                      edf Ref.df      F  p-value    
## s(year)            1.032  2.000  1.040 0.125611    
## te(xc,yc):year2015 3.000  3.000 10.006 4.40e-06 ***
## te(xc,yc):year2016 3.759  4.158  4.890 0.000861 ***
## te(xc,yc):year2017 4.311  5.068  5.237 0.000150 ***
## s(chl):IHOWAO_BF   1.000  1.000  0.716 0.398640    
## s(chl):IHOCAA      1.920  2.205  2.022 0.127781    
## s(chl):IHOBB       1.001  1.001  1.661 0.199101    
## s(chl):IHODS       1.000  1.000  1.954 0.164030    
## s(chl):IHOEAO      3.748  3.961  9.058 4.83e-06 ***
## s(v):IHOWAO_BF     1.000  1.000  3.198 0.075551 .  
## s(v):IHOCAA        1.154  1.286  0.222 0.620204    
## s(v):IHOBB         1.449  1.761  4.160 0.041544 *  
## s(v):IHODS         1.000  1.000  3.710 0.055792 .  
## s(v):IHOEAO        1.000  1.000  0.040 0.842651    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.463   Deviance explained = 53.7%
## -REML = 590.84  Scale est. = 41.952    n = 193
```

Check residuals.


```r
appraise(GAMSA_I3)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI3-residuals-covariates-1.png)<!-- -->

```r
GAMSA_I3_resid <- bind_cols(SA_df, residuals.gam(GAMSA_I3)) %>%
  rename(resid = "...9")
plot_grid(GAMSA_I3_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAMSA_I3_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAMSA_I3_resid %>%
            ggplot(aes(x = xc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "xc") +
            geom_point(), 
          GAMSA_I3_resid %>%
            ggplot(aes(x = yc, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "yc") +
            geom_point(), 
          GAMSA_I3_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAMSA_I3_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI3-residuals-covariates-2.png)<!-- -->

## Model predictions

I predict the model to get the smooths in the response scale (NASC). I need to create a specific data frame for this with the factors of interest `IHO` and `year`, and covariates `v` and `chl`.


```r
# Find median, min and max values of SA_df
val <- SA_df %>%
  dplyr::select(IHO, year, xc, yc, v, chl) %>%
  group_by(IHO) %>%
  summarise(median_xc = median(xc),
            median_yc = median(yc),
            min_v = min(v),
            max_v = max(v),
            median_v = median(v),
            min_chl = min(chl),
            max_chl = max(chl),
            median_chl = median(chl))
# Resolution for predictions
res = 50
# IHO areas
IHO_area <- levels(SA_df$IHO)
# Empty data frame that we populate with new values
SA_new <- data.frame()

for (i in IHO_area) {
  # Select data from which we build the 
  val_tmp <- val %>% 
    filter(IHO == i)
  # Temporary data frame for each region
  SA_tmp <- data.frame(
    IHO = rep(rep(rep(i, res), 3), 2),
    year = rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2),
    var = c(rep("v", res * 3), rep("chl", res * 3)),
    # Velocity
    xc = rep(rep(rep(subset(val, IHO == i)$median_xc, res), 3), 2),
    yc = rep(rep(rep(subset(val, IHO == i)$median_yc, res), 3), 2),
    v = c(rep(seq(subset(val, IHO == i)$min_v,
                  subset(val, IHO == i)$max_v,
                  length.out = res), 3),
          rep(subset(val, IHO == "WAO")$median_v, res * 3)),
    # Chlorophyll
    chl = c(rep(subset(val, IHO == i)$median_chl, res * 3),
            rep(seq(subset(val, IHO == i)$min_chl, 
                    subset(val, IHO == i)$max_chl,
                    length.out = res), 3)))
  # Append data
  SA_new <- bind_rows(SA_new, SA_tmp)
}
# Reorder factors
SA_new <- SA_new %>%
  mutate(IHO = factor(IHO, levels = c("WAO_BF", "CAA", "BB", "DS", "EAO")))
```

Get GAM predictions (`GAM_S` and `GAM_I`) for the new data. I then calculate the average functional response for all years.


```r
ilinkSA_S <- family(GAMSA_S3)$linkinv # Get link function
predSA_S <- predict(GAMSA_S3, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilinkSA_S(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilinkSA_S(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilinkSA_S(fit)) %>% # Calculate fit
  group_by(IHO, var, v, chl) %>%
  summarise(SA_fit = mean(fitted), # Calculate average functional response
            SA_lwr_ci = mean(lwr_ci),
            SA_upr_ci = mean(upr_ci)) 
```


```r
ilinkSA_I <- family(GAMSA_I3)$linkinv # Get link function
predSA_I <- predict(GAMSA_I3, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilinkSA_I(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilinkSA_I(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilinkSA_I(fit)) %>% # Calculate fit
  group_by(IHO, var, v, chl) %>%
  summarise(SA_fit = mean(fitted), # Calculate average functional response
            SA_lwr_ci = mean(lwr_ci),
            SA_upr_ci = mean(upr_ci)) 
```


```r
# Save data
save(predSA_S, predSA_I, GAMSA_S2, GAMSA_I2, file = "data/statistics/GAMSA_results.RData")
```

## Model visualization

Plot to see if predictions worked well.


```r
plot_grid(
  predSA_S %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = SA_fit)) +
    geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = SA_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model SA S"),
  predSA_I %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = SA_fit)) +
    geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = SA_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model SA I"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/veloSA-ggplot-1.png)<!-- -->


```r
plot_grid(
  predSA_S %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = SA_fit)) +
    geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = SA_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model SA S"),
  predSA_I %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = SA_fit)) +
    geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = SA_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model SA I"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/chlSA-ggplot-1.png)<!-- -->


```r
GAMSAS_velo <- predSA_S %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO),
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "Model S",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAMSAS_chl <- predSA_S %>%
  filter(var == "chl") %>%
  ggplot() +
  geom_line(aes(x = chl, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), 
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  labs(x = expression("Chlorophyll (mg m"^-3*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "none")
# Combine plots
plot_grid(GAMSAS_velo, GAMSAS_chl, ncol = 1, rel_heights = c(1,0.8))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS-pretty-ggplot-1.png)<!-- -->


```r
GAMSAI_velo <- predSA_I %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), 
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "Model I",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAMSAI_chl <- predSA_I %>%
  filter(var == "chl") %>%
  ggplot() +
  geom_line(aes(x = chl, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO),
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  labs(x = expression("Chlorophyll (mg m"^-3*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "none")
# Combine plots
plot_grid(GAMSAI_velo, GAMSAI_chl, ncol = 1, rel_heights = c(1,0.8))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-pretty-ggplot-1.png)<!-- -->
