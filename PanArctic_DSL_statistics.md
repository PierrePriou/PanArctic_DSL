---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/24 at 14:31"
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
cell_res <- 50 
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
load("data/acoustics/SA_grids.RData") # Acoustic data
load("data/remote_sensing/physics_grids.RData") # Modelled physics data 
load("data/remote_sensing/seaice_grids.RData") # Remote sensing sea ice data
load("data/remote_sensing/chl_grid.RData") # Remote sensing sea ice data
```

# Data preparation

I combine backscatter and environmental data which were gridded on the Lambert Azimuthal Equal Area grid (EPSG:6931). Since I do not have a lot of samples from the the Beaufort Sea and West Arctic Ocean, I combine those two regions due to their geographical proximity.


```r
# Gridded SA
SA_laea <- SA_grid_laea %>%  
  dplyr::select(-lat, -lon) 
# Remote sensing: physics reanalysis
phy_laea <- phy_grid_laea %>% 
  dplyr::select(year, area, xc, yc, cell_res, depth, velocity)

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
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = velocity)) +
    scale_colour_viridis_c("v", na.value = "red") + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank()),
  ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_statistics_files/figure-html/map-NA-1.png)<!-- -->

This is because the coverage of acoustic data and remote sensing products slightly differ close to the coasts. For sea ice concentration (and open water duration) the NAs are close to land while for ocean colour they are in ice covered areas. I use the average value of neighbouring cells (8 cells around the pixel of interest) to fill missing values.


```r
stat_laea <- stat_laea %>%
  rowwise() %>%
  mutate( 
    # Ocean colour
    xc_chl = if_else(is.na(chl) == T, xc, NaN),
    yc_chl = if_else(is.na(chl) == T, yc, NaN),
    year_chl = if_else(is.na(chl) == T, year, NaN), 
    chl = if_else(
      is.na(chl) == T,
      mean(pull(subset(chl_grid_laea, 
                       xc >= xc_chl - cell_res & xc <= xc_chl + cell_res & 
                         yc >= yc_chl - cell_res & yc <= yc_chl + cell_res & 
                         year == year_chl,
                       select = chl),
                chl),
           na.rm = T),
      chl),
    # Convert veloctity to cm/s
    velocity = velocity * 100, 
    IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                IHO_area == "West Arctic Ocean" ~ "WAO",
                                IHO_area == "Beaufort Sea" ~ "BF",
                                IHO_area == "The Northwestern Passages" ~ "CAA",
                                IHO_area == "Baffin Bay" ~ "BB",
                                IHO_area == "Davis Strait" ~ "DS"),
                      levels = c("WAO", "BF", "CAA", "BB", "DS", "EAO"))) %>%
  ungroup() %>%
  dplyr::select(-xc_chl, -yc_chl, -year_chl)
```

Check if I replace NAs correctly.


```r
stat_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  geom_point(aes(col = chl)) +
  scale_colour_cmocean(name = "algae", na.value = "red") + 
  facet_wrap(~ year) +
  coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
  theme(legend.position = "top",
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-NA-fix-1.png)<!-- -->

There are still missing chl values, so these cells are excluded from further analyses.


```r
stat_laea <- stat_laea %>%
  filter(is.na(chl) == F & is.na(velocity) == F)
```


```r
stat_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  geom_point(aes(col = chl)) +
  scale_colour_cmocean(name = "algae", na.value = "red") + 
  facet_wrap(~ year) +
  coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
  theme(legend.position = "top",
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-chl-NA-1.png)<!-- -->

Check the total number of observations per group.


```r
stat_laea %>%
  group_by(IHO_area) %>% 
  summarise(n = n())
```

```
## # A tibble: 6 × 2
##   IHO_area     n
##   <fct>    <int>
## 1 WAO         35
## 2 BF          71
## 3 CAA        101
## 4 BB         258
## 5 DS         167
## 6 EAO        157
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

![](PanArctic_DSL_statistics_files/figure-html/map-velocity-1.png)<!-- -->

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
## 1 NASC_int   139     29.8      4 5.39e- 6 Kruskal-Wallis IHO_area
## 2 NASC_int   139      9.06     2 1.08e- 2 Kruskal-Wallis year    
## 3 CM         139     61.2      4 1.63e-12 Kruskal-Wallis IHO_area
## 4 CM         139      4.16     2 1.25e- 1 Kruskal-Wallis year
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
##   IHO_area .y.          n statistic    df      p method        
## * <fct>    <chr>    <int>     <dbl> <int>  <dbl> <chr>         
## 1 WAO_BF   NASC_int    17     6.47      2 0.0394 Kruskal-Wallis
## 2 CAA      NASC_int    21    11.8       2 0.0027 Kruskal-Wallis
## 3 BB       NASC_int    48     5.10      2 0.0782 Kruskal-Wallis
## 4 DS       NASC_int    25     0.328     2 0.849  Kruskal-Wallis
## 5 EAO      NASC_int    28     4.98      2 0.083  Kruskal-Wallis
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
## -20.827  -4.714  -0.059   4.627  34.498 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  50.9658    10.8055   4.717 5.84e-06 ***
## lat          -0.4550     0.1449  -3.140  0.00207 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 8.236 on 137 degrees of freedom
## Multiple R-squared:  0.06715,	Adjusted R-squared:  0.06034 
## F-statistic: 9.861 on 1 and 137 DF,  p-value: 0.002068
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
  geom_point(aes(x = lat, y = SA_int)) +
  # geom_abline(aes(slope = est_lat, intercept = est_itcpt, col = IHO)) + 
  geom_line(aes(x = lat, y = .fitted, col = IHO)) + 
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
  dplyr::select(year, xc, yc, IHO_area, NASC_int, velocity, chl) %>%
  group_by(IHO_area) %>%
  mutate(
    year = factor(year),
    NASC_int_n = (NASC_int - min(NASC_int)) / (max(NASC_int) - min(NASC_int))) %>%
  ungroup() %>%
  rename(IHO = IHO_area, 
         v = velocity) 
```

Plot data.


```r
SA_df %>%
  ggplot(aes(x = v, y = NASC_int, col = IHO)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, col = "grey20") +
  facet_wrap(~ IHO, scales = "free") +
  theme(legend.position = "none")
```

```
## `geom_smooth()` using formula 'y ~ s(x, bs = "cs")'
```

```
## Warning: Computation failed in `stat_smooth()`:
## x has insufficient unique values to support 10 knots: reduce k.
```

![](PanArctic_DSL_statistics_files/figure-html/plot-data-area-1.png)<!-- -->

```r
SA_df %>%
  ggplot(aes(x = chl, y = NASC_int, col = IHO)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, col = "grey20") +
  facet_wrap(~ IHO, scales = "free") +
  theme(legend.position = "none")
```

```
## `geom_smooth()` using formula 'y ~ s(x, bs = "cs")'
```

```
## Warning: Computation failed in `stat_smooth()`:
## x has insufficient unique values to support 10 knots: reduce k.
```

![](PanArctic_DSL_statistics_files/figure-html/plot-data-area-2.png)<!-- -->

Check spearman correlations.


```r
# Compute Spearman correlation matrix
corr <- SA_df %>% 
  dplyr::select(-year, -xc, -yc, -IHO, -NASC_int_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr, type = "lower", lab = T, title = "corr")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-1.png)<!-- -->

# Velocity and chlorophyll

I expect regional functional responses to derive from a global functional responses because I assume that mesopelagic species would react in a common way to changes in primary production and advection. I also expect differences in regional functional responses as species assemblages and environmental conditions differ from one region of the Arctic to the other. So I fit models GS and GI following the nomenclature of Pedersen et al. 2019.

## Model fitting

In model GS the random intercept is already included in the `s(bs = "fs")` term.


```r
# Model G
GAM_G <- gam(NASC_int ~ 
               # First order effects
               s(v, k = 5, bs = "tp") + # Global smoothers
               s(chl, k = 5, bs = "tp") + # Global smoothers
               s(IHO, bs = "re") + # Random effect
               s(year, bs = "re") + # Random effect
               # Second order effects
               s(IHO, year, bs = "re"), # Random slope
            data = SA_df, family = Gamma(link = "log"), method = "REML")
# Model GS
GAM_GS <- gam(NASC_int ~ 
                # First order effects
                s(chl, bs = "tp", k = 5) +
                s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "REML")
# Model GI
GAM_GI <- gam(NASC_int ~
                # First order effects
                s(chl, bs = "tp", k = 5) +
                s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, by = IHO, k = 5, bs = "tp") + 
                s(v, by = IHO, k = 5, bs = "tp"),
             data = SA_df, family = Gamma(link = "log"), method = "REML")
```

Check model metrics


```r
# Summary metrics
GAM_AIC <- AIC(GAM_G, GAM_GS, GAM_GI) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAM_G", "GAM_GS", "GAM_GI"),
           reml = round(c(GAM_G$gcv.ubre,
                          GAM_GS$gcv.ubre, 
                          GAM_GI$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAM_G$deviance / GAM_G$null.deviance)) * 100,
             (1 - (GAM_GS$deviance / GAM_GS$null.deviance)) * 100,
             (1 - (GAM_GI$deviance / GAM_GI$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAM_G)$r.sq,
                        summary(GAM_GS)$r.sq, 
                        summary(GAM_GI)$r.sq),
                      2)) %>%
  full_join(., GAM_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-9cd53a50288d835157f4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-9cd53a50288d835157f4">{"x":{"filter":"none","vertical":false,"data":[["GAM_GI","GAM_GS","GAM_G"],[33.588,33.592,21.715],[76.11,74.13,64.3],[0.13,0,0.02],[770.99,791.17,804.33],[1558.636,1570.085,1595.621],[0,11.45,36.99],[0.99675,0.00325,0]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
##     s(IHO, bs = "re") + s(year, bs = "re") + s(IHO, year, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.6225     0.6799   6.799 5.24e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df      F p-value  
## s(v)        2.2938  2.777  2.614  0.0649 .
## s(chl)      2.1867  2.669  1.873  0.1533  
## s(IHO)      2.5742  5.000 10.871  0.4260  
## s(year)     0.7108  2.000  6.937  0.4660  
## s(IHO,year) 9.1286 15.000  7.754  0.1088  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0181   Deviance explained = 64.3%
## -REML = 804.33  Scale est. = 2.5421    n = 131
```

Model GS summary.


```r
summary(GAM_GS)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(chl, bs = "tp", k = 5) + s(v, bs = "tp", k = 5) + 
##     s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, bs = "re") + 
##     s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    4.758      0.382   12.46   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                   edf Ref.df     F  p-value    
## s(chl)      1.0008306  1.001 3.664  0.05824 .  
## s(v)        1.0000741  1.000 1.595  0.20950    
## s(year)     0.1458971  2.000 0.122  0.46403    
## s(IHO)      0.0002262  5.000 0.000  0.66812    
## s(year,IHO) 8.7967462 15.000 2.437 6.32e-05 ***
## s(chl,IHO)  8.2834373 26.000 4.790 9.82e-06 ***
## s(v,IHO)    6.2736572 28.000 0.912  0.00282 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.00361   Deviance explained = 74.1%
## -REML = 791.17  Scale est. = 1.5328    n = 131
```

Model GI summary.


```r
summary(GAM_GI)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(chl, bs = "tp", k = 5) + s(v, bs = "tp", k = 5) + 
##     s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, bs = "re") + 
##     s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, k = 5, 
##     bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.7988     0.3027   15.86   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                     edf    Ref.df     F  p-value    
## s(chl)        1.001e+00 1.001e+00 0.125 0.724419    
## s(v)          1.000e+00 1.000e+00 3.023 0.085074 .  
## s(year)       2.498e-04 2.000e+00 0.000 0.532313    
## s(IHO)        3.176e-05 5.000e+00 0.000 0.818937    
## s(year,IHO)   7.376e+00 1.500e+01 1.883 6.54e-05 ***
## s(chl):IHOWAO 1.000e+00 1.000e+00 8.142 0.005228 ** 
## s(chl):IHOBF  1.437e+00 1.623e+00 1.699 0.267389    
## s(chl):IHOCAA 2.825e-05 4.849e-05 0.004 0.999653    
## s(chl):IHOBB  2.822e+00 3.255e+00 1.983 0.105609    
## s(chl):IHODS  1.000e+00 1.000e+00 0.003 0.960333    
## s(chl):IHOEAO 2.848e+00 3.252e+00 6.297 0.000747 ***
## s(v):IHOWAO   1.000e+00 1.000e+00 0.207 0.650387    
## s(v):IHOBF    5.931e-05 1.156e-04 0.015 0.998969    
## s(v):IHOCAA   1.000e+00 1.000e+00 2.123 0.148141    
## s(v):IHOBB    2.272e+00 2.715e+00 2.218 0.116333    
## s(v):IHODS    1.000e+00 1.000e+00 7.973 0.005698 ** 
## s(v):IHOEAO   3.597e+00 3.882e+00 3.243 0.009333 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Rank: 82/84
## R-sq.(adj) =  0.134   Deviance explained = 76.1%
## -REML = 770.99  Scale est. = 1.2855    n = 131
```

## Model checking

### Basis size and residual distribution

First I check the basis size k. `k-indexes` are \> 1 or close to 1 so the basis size is large enough. The residual plot look good too.


```r
appraise(GAM_GS, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGS-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_GS)
```

```
##             k'         edf   k-index p-value
## s(chl)       4 1.000830579 0.9420324  0.7250
## s(v)         4 1.000074135 0.8460029  0.2200
## s(year)      3 0.145897063        NA      NA
## s(IHO)       6 0.000226213        NA      NA
## s(year,IHO) 18 8.796746199        NA      NA
## s(chl,IHO)  30 8.283437284 0.9420324  0.7425
## s(v,IHO)    30 6.273657151 0.8460029  0.2950
```


```r
appraise(GAM_GI, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_GI)
```

```
##               k'          edf   k-index p-value
## s(chl)         4 1.000930e+00 0.9368622  0.6900
## s(v)           4 1.000091e+00 0.8582598  0.2825
## s(year)        3 2.497932e-04        NA      NA
## s(IHO)         6 3.175713e-05        NA      NA
## s(year,IHO)   18 7.376439e+00        NA      NA
## s(chl):IHOWAO  4 1.000009e+00 0.9368622  0.6500
## s(chl):IHOBF   4 1.436755e+00 0.9368622  0.6600
## s(chl):IHOCAA  4 2.824578e-05 0.9368622  0.6925
## s(chl):IHOBB   4 2.822030e+00 0.9368622  0.7050
## s(chl):IHODS   4 1.000023e+00 0.9368622  0.6500
## s(chl):IHOEAO  4 2.847750e+00 0.9368622  0.6800
## s(v):IHOWAO    4 1.000041e+00 0.8582598  0.2625
## s(v):IHOBF     4 5.930684e-05 0.8582598  0.2650
## s(v):IHOCAA    4 1.000050e+00 0.8582598  0.2550
## s(v):IHOBB     4 2.272183e+00 0.8582598  0.2550
## s(v):IHODS     4 1.000218e+00 0.8582598  0.2875
## s(v):IHOEAO    4 3.597326e+00 0.8582598  0.3375
```

### Resiudals against covariates


```r
GAM_GS_resid <- bind_cols(SA_df, residuals.gam(GAM_GS)) %>%
  rename(resid = "...9")
plot_grid(GAM_GS_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_GS_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_GS_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_GS_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGS-residuals-covariates-1.png)<!-- -->


```r
GAM_GI_resid <- bind_cols(SA_df, residuals.gam(GAM_GI)) %>%
  rename(resid = "...9")
plot_grid(GAM_GI_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_GI_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_GI_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_GI_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/veloGI-residuals-covariates-1.png)<!-- -->

## Term selection

I turn on the double penalty (`select = TRUE`) and check the covariates that has been shrunk.

### Model GS


```r
# Model GS
GAM_GS_p <- gam(NASC_int ~ 
                # First order effects
                s(chl, bs = "tp", k = 5) +
                s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "REML", 
             select = T)
summary(GAM_GS_p)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(chl, bs = "tp", k = 5) + s(v, bs = "tp", k = 5) + 
##     s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, bs = "re") + 
##     s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.6726     0.3765   12.41   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                   edf Ref.df     F  p-value    
## s(chl)      0.7896126      4 0.874 7.93e-07 ***
## s(v)        0.6567804      4 0.203  0.10692    
## s(year)     0.4011720      2 0.681  0.46356    
## s(IHO)      0.0002687      5 0.000  0.58679    
## s(year,IHO) 8.7827645     15 2.597  0.00127 ** 
## s(chl,IHO)  8.5683335     27 4.544 2.34e-05 ***
## s(v,IHO)    5.3192768     29 0.737  0.00330 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0332   Deviance explained = 73.9%
## -REML = 793.03  Scale est. = 1.5258    n = 131
```

The global smoother `s(v)` is close to 0 suggesting that there is no global pattern with regard to advection.The random effects `s(IHO)` and `s(year)` also have low edf but I decide to keep them in the model as their interaction `s(year,IHO)` is significant. I refit the model without the `s(v)` term.


```r
# Model GS
GAM_GS2 <- gam(NASC_int ~ 
                # First order effects
                s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "REML")
summary(GAM_GS2)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(chl, bs = "tp", k = 5) + s(year, bs = "re") + s(IHO, 
##     bs = "re") + s(year, IHO, bs = "re") + s(chl, IHO, bs = "fs", 
##     k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.7906     0.3876   12.36   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                   edf Ref.df     F  p-value    
## s(chl)      1.0003586      1 4.015 0.047599 *  
## s(year)     0.3979259      2 0.672 0.457078    
## s(IHO)      0.0001493      5 0.000 0.597781    
## s(year,IHO) 8.7196921     15 2.593 0.001123 ** 
## s(chl,IHO)  8.4360022     26 4.884 1.66e-05 ***
## s(v,IHO)    5.8823787     29 1.245 0.000167 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0279   Deviance explained = 73.8%
## -REML =  791.5  Scale est. = 1.5042    n = 131
```

Compare models.


```r
GAM_GS_AIC <- AIC(GAM_GS, GAM_GS2) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAM_GS", "GAM_GS2"),
           reml = round(c(GAM_GS$gcv.ubre,
                          GAM_GS2$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAM_GS$deviance / GAM_GS$null.deviance)) * 100,
             (1 - (GAM_GS2$deviance / GAM_GS2$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAM_GS)$r.sq,
                        summary(GAM_GS2)$r.sq),
                      2)) %>%
  full_join(., GAM_GS_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 10)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-845cc66a51e691162fd6" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-845cc66a51e691162fd6">{"x":{"filter":"none","vertical":false,"data":[["GAM_GS2","GAM_GS"],[32.126,33.592],[73.82,74.13],[0.03,0],[791.5,791.17],[1568.88,1570.085],[0,1.2],[0.6462280583,0.3537719417]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

`GAM_GS` has better metics. So I keep that model.

Check residuals.


```r
appraise(GAM_GS2)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGS2-residuals-covariates-1.png)<!-- -->

```r
GAM_GS2_resid <- bind_cols(SA_df, residuals.gam(GAM_GS2)) %>%
  rename(resid = "...9")
plot_grid(GAM_GS2_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_GS2_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_GS2_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_GS2_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGS2-residuals-covariates-2.png)<!-- -->

### Model GI


```r
# Model GI
GAM_GI_p <- gam(NASC_int ~
                  # First order effects
                  s(chl, bs = "tp", k = 5) +
                  s(v, bs = "tp", k = 5) +
                  s(year, bs = "re") +
                  s(IHO, bs = "re") +
                  # Second order random effects
                  s(year, IHO, bs = "re") +
                  # Second order functional effects
                  s(chl, by = IHO, k = 5, bs = "tp") + 
                  s(v, by = IHO, k = 5, bs = "tp"),
                data = SA_df, family = Gamma(link = "log"), method = "REML",
                select = T)
summary(GAM_GI_p)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(chl, bs = "tp", k = 5) + s(v, bs = "tp", k = 5) + 
##     s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, bs = "re") + 
##     s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, k = 5, 
##     bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.6699     0.2485   18.79   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                     edf Ref.df      F  p-value    
## s(chl)        0.0002481      4  0.000  0.65570    
## s(v)          1.7734232      4  2.521  0.01460 *  
## s(year)       0.0002346      2  0.000  0.61464    
## s(IHO)        0.0002470      5  0.000  0.35976    
## s(year,IHO)   9.1698294     15  2.356 2.39e-05 ***
## s(chl):IHOWAO 0.9127738      3 26.563  0.00102 ** 
## s(chl):IHOBF  0.7858752      4  4.239  0.03096 *  
## s(chl):IHOCAA 0.0002570      4  0.000  0.75611    
## s(chl):IHOBB  0.5906950      4  0.441  0.19831    
## s(chl):IHODS  0.0001718      4  0.000  0.92001    
## s(chl):IHOEAO 2.5420467      4 21.952  < 2e-16 ***
## s(v):IHOWAO   0.0001933      4  0.000  0.92336    
## s(v):IHOBF    0.0004501      4  0.000  0.70011    
## s(v):IHOCAA   0.0006383      4  0.000  0.47251    
## s(v):IHOBB    0.0007327      4  0.000  0.32383    
## s(v):IHODS    1.3401252      4  3.574  0.00831 ** 
## s(v):IHOEAO   2.6366373      4  5.270  0.00016 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.151   Deviance explained =   74%
## -REML = 787.11  Scale est. = 1.3632    n = 131
```

The global smoother `s(chl)` is close to 0 so I refit the model without it. The random effects `s(IHO)` and `s(year)` also have low edf but I decide to keep them in the model as their interaction `s(year,IHO)` is significant. I refit the model without the `s(v,year)` and `s(chl,year)` terms.


```r
# Model GI
GAM_GI2 <- gam(NASC_int ~
                 # First order effects
                 # s(chl, bs = "tp", k = 5) +
                 s(v, bs = "tp", k = 5) +
                 s(year, bs = "re") +
                 s(IHO, bs = "re") +
                 # Second order random effects
                 s(year, IHO, bs = "re") +
                 # Second order functional effects
                 s(chl, by = IHO, k = 5, bs = "tp") + 
                 s(v, by = IHO, k = 5, bs = "tp"),
               data = SA_df, family = Gamma(link = "log"), method = "REML")
summary(GAM_GI2)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(v, bs = "tp", k = 5) + s(year, bs = "re") + s(IHO, 
##     bs = "re") + s(year, IHO, bs = "re") + s(chl, by = IHO, k = 5, 
##     bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.9628     0.3045    16.3   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                     edf    Ref.df      F  p-value    
## s(v)          1.000e+00 1.000e+00  0.221 0.639154    
## s(year)       2.777e-04 2.000e+00  0.000 0.578322    
## s(IHO)        3.145e-05 5.000e+00  0.000 0.861597    
## s(year,IHO)   7.286e+00 1.500e+01  1.836 0.000105 ***
## s(chl):IHOWAO 1.000e+00 1.000e+00 15.437 0.000154 ***
## s(chl):IHOBF  1.542e+00 1.746e+00  2.895 0.059977 .  
## s(chl):IHOCAA 1.000e+00 1.000e+00  0.325 0.570053    
## s(chl):IHOBB  1.000e+00 1.000e+00  1.626 0.205040    
## s(chl):IHODS  1.000e+00 1.000e+00  0.223 0.637410    
## s(chl):IHOEAO 2.836e+00 3.240e+00 15.062  < 2e-16 ***
## s(v):IHOWAO   1.000e+00 1.000e+00  1.126 0.291058    
## s(v):IHOBF    1.000e+00 1.000e+00  2.156 0.144995    
## s(v):IHOCAA   3.719e-05 7.217e-05  0.000 0.500000    
## s(v):IHOBB    2.433e+00 2.879e+00  3.403 0.024289 *  
## s(v):IHODS    1.000e+00 1.000e+00  1.386 0.241740    
## s(v):IHOEAO   3.584e+00 3.875e+00  4.455 0.001623 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Rank: 79/80
## R-sq.(adj) =  0.151   Deviance explained = 74.9%
## -REML = 772.48  Scale est. = 1.3598    n = 131
```

Compare models.


```r
GAM_GI_AIC <- AIC(GAM_GI, GAM_GI2) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAM_GI", "GAM_GI2"),
           reml = round(c(GAM_GI$gcv.ubre,
                          GAM_GI2$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAM_GI$deviance / GAM_GI$null.deviance)) * 100,
             (1 - (GAM_GI2$deviance / GAM_GI2$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAM_GI)$r.sq,
                        summary(GAM_GI2)$r.sq),
                      2)) %>%
  full_join(., GAM_GI_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 10)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-cee52cbc79d98a2840d9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-cee52cbc79d98a2840d9">{"x":{"filter":"none","vertical":false,"data":[["GAM_GI","GAM_GI2"],[33.588,31.494],[76.11,74.92],[0.13,0.15],[770.99,772.48],[1558.636,1561.392],[0,2.76],[0.7986695994,0.2013304006]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

`GAM_GI` has better metics. So I keep that model.

Check residuals.


```r
appraise(GAM_GI2)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGI2-residuals-covariates-1.png)<!-- -->

```r
GAM_GI2_resid <- bind_cols(SA_df, residuals.gam(GAM_GI2)) %>%
  rename(resid = "...9")
plot_grid(GAM_GS2_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_GS2_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_GS2_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_GS2_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGI2-residuals-covariates-2.png)<!-- -->

## Model predictions

I predict the model to get the smooths in the response scale (NASC). I need to create a specific data frame for this with the factors of interest `IHO` and `year`, and covariates `v` and `chl`.


```r
val <- SA_df %>%
  dplyr::select(IHO, year, v, chl) %>%
  group_by(IHO) %>%
  summarise(min_v = min(v),
            max_v = max(v),
            median_v = median(v),
            min_chl = min(chl),
            max_chl = max(chl),
            median_chl = median(chl))
val
```

```
## # A tibble: 6 × 7
##   IHO   min_v max_v median_v min_chl max_chl median_chl
##   <fct> <dbl> <dbl>    <dbl>   <dbl>   <dbl>      <dbl>
## 1 WAO   0.575  4.47    3.09   0.111    0.235      0.186
## 2 BF    0.779  5.39    2.15   0.181    0.857      0.483
## 3 CAA   0.473  3.39    0.874  0.234    0.806      0.405
## 4 BB    0.119  4.07    1.81   0.283    2.54       1.01 
## 5 DS    0.370  2.74    1.49   0.319    0.685      0.413
## 6 EAO   0.450  5.20    2.07   0.0892   4.06       1.35
```


```r
# Resolution for predictions
res = 20
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
  if (i == last(IHO_area)){
    print("New data frame built.")
  }
}
```

```
## [1] "New data frame built."
```

```r
# Reorder factors
SA_new <- SA_new %>%
  mutate(IHO = factor(IHO, levels = c("WAO", "BF", "CAA", "BB", "DS", "EAO")))
```

Get GAM predictions (`GAM_GS` and `GAM_GI`) for the new data. I then calculate the average functional response for all years.


```r
ilink_GS <- family(GAM_GS)$linkinv # Get link function
pred_GS <- predict(GAM_GS, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilink_GS(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilink_GS(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilink_GS(fit)) %>% # Calculate fit
  group_by(IHO, var, v, chl) %>%
  summarise(NASC_fit = mean(fitted), # Calculate average functional response
            NASC_lwr_ci = mean(lwr_ci),
            NASC_upr_ci = mean(upr_ci)) %>%
  mutate(SA_fit = 10 * log10(NASC_fit), # Convert to SA
         SA_lwr_ci = 10 * log10(NASC_lwr_ci),
         SA_upr_ci = 10 * log10(NASC_upr_ci))
```


```r
ilink_GI <- family(GAM_GI)$linkinv # Get link function
pred_GI <- predict(GAM_GI, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilink_GI(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilink_GI(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilink_GI(fit)) %>% # Calculate fit
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
save(pred_GS, pred_GI, GAM_GS, GAM_GI, file = "data/statistics/GAM_results.RData")
```

## Model visualization

Plot to see if predictions worked well.


```r
plot_grid(
  pred_GS %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = NASC_fit)) +
    geom_ribbon(aes(x = v, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model GS"),
  pred_GI %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = NASC_fit)) +
    geom_ribbon(aes(x = v, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model GI"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/velo-ggplot-1.png)<!-- -->


```r
plot_grid(
  pred_GS %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = NASC_fit)) +
    geom_ribbon(aes(x = chl, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model GS"),
  pred_GI %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = NASC_fit)) +
    geom_ribbon(aes(x = chl, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("Model GI"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/chl-ggplot-1.png)<!-- -->


```r
GAMGS_velo <- pred_GS %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "Model GS",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAMGS_chl <- pred_GS %>%
  filter(var == "chl") %>%
  ggplot() +
  geom_line(aes(x = chl, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), alpha = 0.1) +
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
plot_grid(GAMGS_velo, GAMGS_chl, ncol = 1, rel_heights = c(1,0.8))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGS-pretty-ggplot-1.png)<!-- -->


```r
GAMGI_velo <- pred_GI %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "Model GI",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAMGI_chl <- pred_GI %>%
  filter(var == "chl") %>%
  ggplot() +
  geom_line(aes(x = chl, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), alpha = 0.1) +
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
plot_grid(GAMGI_velo, GAMGI_chl, ncol = 1, rel_heights = c(1,0.8))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGI-pretty-ggplot-1.png)<!-- -->
