---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/06/10 at 16:53"
output: 
  html_document:
    code_folding: show
    keep_md: yes
    toc: true
    toc_float: true
    toc_collapsed: true
  github_document:
    always_allow_html: true
---

# Goals

The goal of the statistics are to investigate the potential effects of advection on mesopelagic backscatter for each Arctic region. I use data from the Copernicus Marine Environmental Monitoring Service Arctic Ocean Physics Reanalysis product which gives monthly average of current velocity in the Arctic Ocean.

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
             strip.text.x = element_text(size = 9,
                                         face = "plain",
                                         hjust = 0.5),
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
cell_res <- 150
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931")
projection(arctic_laea) <- gsub("units=m",
                                "units=km",
                                projection(arctic_laea))
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
```

# Data preparation

I combine backscatter and environmental data which were gridded on the Lambert Azimuthal Equal Area grid (EPSG:6931) with a cell resolution of 150 x 150 km. Since I do not have a lot of samples from the the Beaufort Sea and West Arctic Ocean, I combine those two regions due to their geographical proximity.


```r
# Gridded SA
SA_laea <- SA_grid_laea %>%  
  dplyr::select(-lat, -lon) 
# Remote sensing: physics reanalysis
phy_laea <- phy_grid_laea %>% 
  dplyr::select(year, area, xc, yc, cell_res, depth, v_mean) %>%
  # Find missing depth values
  complete(depth, nesting(year, area, xc, yc, cell_res), 
           fill = list(v_mean = NA)) %>%
  arrange(year, area, xc, yc, cell_res, depth) %>%
  mutate(depth = factor(depth))

stat_laea <- SA_laea %>%
  left_join(., phy_laea, by = c("year", "xc", "yc", "area", "cell_res"))
```

There is 1 NA in the physics dataset.


```r
plot_grid(
  stat_laea %>%
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = LME)) +
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()),
  stat_laea %>% 
    filter(depth == 318) %>%
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = v_mean)) +
    scale_colour_cmocean(name = "speed", na.value = "red") + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank()),
  ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_statistics_files/figure-html/map-NA-1.png)<!-- -->


```r
stat_laea <- stat_laea %>%
  mutate(v = v_mean * 100,
         LME = 
           factor(case_when(LME == "Beaufort Sea LME" ~ "BF",
                            LME == "Central Arctic LME" ~ "CAO",
                            LME == "Baffin Bay LME" ~ "BB",
                            LME == "Barents Sea LME" ~ "SV",
                            LME == "Northern Canadian Archipelago LME" ~ "NAR"),
                  levels = c("CAO", "BF", "NAR", "BB", "SV")))

stat_laea %>%
  filter(depth == 318) %>%
  group_by(LME) %>% 
  summarise(n = n())
```

```
## # A tibble: 5 × 2
##   LME       n
##   <fct> <int>
## 1 CAO       4
## 2 BF       20
## 3 NAR       3
## 4 BB       33
## 5 SV       13
```
Because we sampled at the border of the CAO and there is only a limited amount of data in this region we combine CAO data based on their proximity to other LME (Barents Sea and Beaufort Sea). I exclude data from Nares Strait because we do not have enough points for a regression


```r
stat_laea <- stat_laea %>%
  mutate(LME = as.character(LME),
         LME = factor(if_else(LME == "CAO" & yc <= 0, "SV",
                      if_else(LME == "CAO" & yc > 0, "BF", LME)),
                      # if_else(LME == "NAR", "BB", LME))),
                      levels = c("BF", "BB", "SV", "NAR"))) %>%
  filter(is.na(v) == F & LME != "NAR") %>%
  droplevels()

stat_laea %>%
  filter(depth == 318) %>%
  group_by(LME) %>% 
  summarise(n = n())
```

```
## # A tibble: 3 × 2
##   LME       n
##   <fct> <int>
## 1 BF       22
## 2 BB       32
## 3 SV       15
```


```r
plot_grid(
  stat_laea %>%
    filter(depth == 318) %>%
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = LME)) +
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()),
  stat_laea %>% 
    filter(depth == 318) %>%
    ggplot(aes(x = xc,  y = yc)) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    geom_point(aes(col = v_mean)) +
    scale_colour_cmocean(name = "speed", na.value = "red") + 
    facet_wrap(~ year) +
    coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank()),
  ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_statistics_files/figure-html/map-NA2-1.png)<!-- -->

# Data exploration

Maps of all variables.


```r
# LME Regions
stat_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  geom_point(aes(col = LME)) +
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
# Velocity
phy_grid_laea %>% 
  filter(depth == 318) %>%
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = v_mean * 100)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("velo (cm/s)", name = "speed", limits = c(0, 8)) +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Velocity at 318 m depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(legend.position = "top",
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-velocity318-1.png)<!-- -->

# Spatial and interannual variability


```r
SA_diff <- SA_grid_laea %>%
  mutate(LME = case_when(LME == "Beaufort Sea LME" ~ "BF",
                         LME == "Central Arctic LME" ~ "CAO",
                         LME == "Baffin Bay LME" ~ "BB",
                         LME == "Barents Sea LME" ~ "SV",
                         LME == "Northern Canadian Archipelago LME" ~
                           "NAR"),
         LME = factor(if_else(LME == "CAO" & yc <= 0, "SV",
                      if_else(LME == "CAO" & yc > 0, "BF", LME)),
                      levels = c("BF", "NAR", "BB", "SV")),
         year = factor(year)) %>%
  dplyr::select(year, LME, NASC_int, SA_int, CM)
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
          ggqqplot(SA_diff, "NASC_int", facet.by = "LME"),
          ggqqplot(SA_diff, "CM", facet.by = "LME"), 
          ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/normality-check-1.png)<!-- -->

I used a non parametric Kruskal-Wallis test to test for interannual and spatial variability because the variance is not homogenous and not very normal.


```r
bind_rows(kruskal_test(SA_diff, NASC_int ~ LME),
          kruskal_test(SA_diff, NASC_int ~ year),
          kruskal_test(SA_diff, CM ~ LME),
          kruskal_test(SA_diff, CM ~ year)) %>%
  mutate(var = c("LME", "year", "LME", "year"))
```

```
## # A tibble: 4 × 7
##   .y.          n statistic    df       p method         var  
##   <chr>    <int>     <dbl> <int>   <dbl> <chr>          <chr>
## 1 NASC_int    72     13.0      3 0.00454 Kruskal-Wallis LME  
## 2 NASC_int    72     12.3      2 0.00214 Kruskal-Wallis year 
## 3 CM          72     16.0      3 0.00113 Kruskal-Wallis LME  
## 4 CM          72      2.64     2 0.267   Kruskal-Wallis year
```

Mesopelagic backscatter S\~A\~ varied significantly between areas (*H* = 29.792952, *p* \< 0.001) and years (*H* = 29.792952, *p* = 0.011). The center of mass varied significantly between areas (*H* = 61.185994, *p* \< 0.001) but not between years (*H* = 4.164743, *p* = 0.125).

## Within areas

I check whether there are inter-annual variability in S\~A\~ across years within areas. This will be used to decide whether a random effect is needed in the GAM. I used a non parametric Kruskal-Wallis test to test for interannual variability within areas.

### NASC


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(LME) %>%
  kruskal_test(NASC_int ~ year)
```

```
## # A tibble: 4 × 7
##   LME   .y.          n statistic    df        p method        
## * <fct> <chr>    <int>     <dbl> <int>    <dbl> <chr>         
## 1 BF    NASC_int    21     14.1      2 0.000884 Kruskal-Wallis
## 2 NAR   NASC_int     3      0        1 1        Kruskal-Wallis
## 3 BB    NASC_int    33      5.12     2 0.0775   Kruskal-Wallis
## 4 SV    NASC_int    15      4.10     2 0.129    Kruskal-Wallis
```

### CM


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(LME) %>%
  kruskal_test(NASC_int ~ year)
```

```
## # A tibble: 4 × 7
##   LME   .y.          n statistic    df        p method        
## * <fct> <chr>    <int>     <dbl> <int>    <dbl> <chr>         
## 1 BF    NASC_int    21     14.1      2 0.000884 Kruskal-Wallis
## 2 NAR   NASC_int     3      0        1 1        Kruskal-Wallis
## 3 BB    NASC_int    33      5.12     2 0.0775   Kruskal-Wallis
## 4 SV    NASC_int    15      4.10     2 0.129    Kruskal-Wallis
```

# LM - Latitude - S\~A\~ int

A number of studies show decreasing mesopelagic backscatter with increasing latitude. To test this hypothesis, I use a linear regression. Because there are some high mesopelagic backscatter values in the East Arctic Ocean I transformed the nautical area scattering coefficient into nautical area scattering strength in dB.

## Data preparation

First I prepare the data.


```r
SA_lm <- SA_grid_laea %>%
  mutate(LME = case_when(LME == "Beaufort Sea LME" ~ "BF",
                         LME == "Central Arctic LME" ~ "CAO",
                         LME == "Baffin Bay LME" ~ "BB",
                         LME == "Barents Sea LME" ~ "SV",
                         LME == "Northern Canadian Archipelago LME" ~
                           "NAR"),
         LME = factor(if_else(LME == "CAO" & yc <= 0, "SV",
                      if_else(LME == "CAO" & yc > 0, "BF", LME)),
                      levels = c("BF", "NAR", "BB", "SV")))
```

Plot the relationship I want to test.


```r
SA_lm %>%
  ggplot(aes(x = lat, y = SA_int)) +
  geom_point() + 
  geom_smooth(aes(x = lat, y = SA_int, group = LME, col = LME), method = "lm") +
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
##      Min       1Q   Median       3Q      Max 
## -19.1496  -5.8948   0.7596   5.0316  29.3945 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  44.4657    15.7599   2.821  0.00622 **
## lat          -0.3735     0.2110  -1.770  0.08100 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 8.254 on 70 degrees of freedom
## Multiple R-squared:  0.04286,	Adjusted R-squared:  0.02919 
## F-statistic: 3.135 on 1 and 70 DF,  p-value: 0.081
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
  filter(LME != "NAR") %>%
  nest_by(LME) %>%
  # Fit linear model for each arctic region
  mutate(mod = list(lm(SA_int ~ lat, data = data))) %>%
  # Find coefficients, summary (AIC, r2, etc) and fit
  summarise(coefs = list(tidy(mod)),
            summary_mod = list(glance(mod)),
            fit = list(augment(mod, interval = "confidence"))) %>%
  ungroup()

# Extract coefficients
LM_area_coefs <- LM_area %>%
  dplyr::select(LME, coefs) %>% 
  unnest(cols = c(coefs)) %>%
  mutate(term = if_else(term == "(Intercept)", "itcpt", "lat")) %>%
  rename(est = estimate, 
         sd = std.error,
         p_val = p.value) %>%
  pivot_wider(names_from = term, values_from = c(est, sd, statistic, p_val))

# Extract fit
LM_area_fit <- LM_area %>%
  dplyr::select(LME, fit) %>% 
  unnest(cols = c(fit))

# Combine data
LM_area_res <- left_join(LM_area_coefs, LM_area_fit)
```

```
## Joining, by = "LME"
```

Plot regressions. We can see that mesopelagic backscatter decreases with increasing latitude in all regions.


```r
LM_area_res %>%
  ggplot() + 
  geom_point(data = SA_lm, aes(x = lat, y = SA_int, col = LME)) +
  geom_line(aes(x = lat, y = .fitted, col = LME), size = 0.8) + 
  geom_ribbon(aes(x = lat, ymin = .lower, ymax = .upper, group = LME), 
              alpha = 0.1)
```

![](PanArctic_DSL_statistics_files/figure-html/lm-plot-area-1.png)<!-- -->

## Baffin Bay only

There was no poleward decrease in Baffin Bay. Thus, I want to make sure that this pattern is still valid if I remove data from eastern Baffin Bay where myctophid dominated the mesopelagic layer.


```r
load("data/nets/trawl_BB.RData")

SA_grid_laea %>%
  mutate(exclude = if_else(xc == -2025 & yc == -1275, T, 
                   if_else(xc == -2025 & yc == -1125, T, 
                   if_else(xc == -1875 & yc == -1125, T, 
                   if_else(xc == -1725 & yc == -825, T, F))))) %>%
  ggplot() +
  # Coastlines
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), 
               fill = "grey25") + # Coast
  geom_tile(aes(x = xc, y = yc, fill = exclude), col = "white") +
  scale_fill_viridis_d() + 
  # Stations with myctophid
  geom_point(data = trawl_BB, aes(x = xc, y = yc, col = mycto)) +
  # Exclusion line
  geom_segment(aes(x = -2200, xend = -1700, 
                   y = -1400, yend = -600)) + 
  scale_x_continuous(breaks = seq(-10000, 0, 200)) +
  scale_y_continuous(breaks = seq(-10000, 0, 200)) +
  coord_fixed(xlim = c(-2450, -1000), ylim = c(-1500, 0), expand = F)
```

![](PanArctic_DSL_statistics_files/figure-html/map-BB-station-1.png)<!-- -->



```r
# Remove data where myctophid occur
SA_BB <- SA_lm %>%
  filter(LME == "BB") %>%
  mutate(exclude = if_else(xc == -2025 & yc == -1275, T, 
                   if_else(xc == -2025 & yc == -1125, T, 
                   if_else(xc == -1875 & yc == -1125, T, 
                   if_else(xc == -1725 & yc == -825, T, F))))) %>%
  filter(exclude == F) 
# Fit LM
LM_BB <- lm(SA_int ~ lat, data = SA_BB)
summary(LM_BB)
```

```
## 
## Call:
## lm(formula = SA_int ~ lat, data = SA_BB)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -16.6413  -1.6970  -0.0396   2.9139   8.7292 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept) 21.01544   23.79007   0.883    0.385
## lat         -0.01895    0.32695  -0.058    0.954
## 
## Residual standard error: 5.281 on 26 degrees of freedom
## Multiple R-squared:  0.0001292,	Adjusted R-squared:  -0.03833 
## F-statistic: 0.00336 on 1 and 26 DF,  p-value: 0.9542
```
There are still no decrease of mesopelagic backscatter when we remove cells with myctophids.


```r
SA_BB %>%
  ggplot(aes(x = lat, y = SA_int)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_abline(intercept = 21.01544, slope = -0.01895)
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](PanArctic_DSL_statistics_files/figure-html/plot-LM-BB-1.png)<!-- -->


## Save data


```r
# Save data
save(LM_area_coefs, LM_area_fit, LM_area_res, SA_lm, 
     file = "data/statistics/LM_results.RData")
```

# HGAM - Gamma NASC

## Data preparation

I select data at 318 m depth because it fits well with the DSL diurnal centre of mass. It also prevent from removing too much data like the deeper depth levels would do.


```r
SA_df <- stat_laea %>%
  filter(depth == 318) %>%
  dplyr::select(year, xc, yc, LME, NASC_int, SA_int, v, depth) %>%
  mutate(year = factor(year))
```

Plot data.


```r
SA_df %>%
  ggplot(aes(x = v, y = SA_int, col = LME)) +
  geom_point() +
  facet_grid(depth ~ LME, scales = "free") +
  theme(legend.position = "none")
```

![](PanArctic_DSL_statistics_files/figure-html/plot-data-area-1.png)<!-- -->

I do not expect regional functional responses to derive from a global functional responses because I expect differences in species assemblages and environmental conditions from one region of the Arctic to the other. So I fit models S and I following the nomenclature proposed by Pedersen et al. 2019.

## Model fitting

In model S the random intercept is already included in the `s(bs = "fs")` term.


```r
# Model G
GAM_G <- gam(NASC_int ~ 
               s(v, k = 5, bs = "tp") + # Global smoothers
               s(LME, bs = "re") + # Random effect
               s(year, bs = "re"), # Random effect
            data = SA_df, family = Gamma(link = "log"), method = "ML")
# Model S
GAM_S <- gam(NASC_int ~
                s(year, bs = "re") +
                s(LME, bs = "re") +
                s(v, LME, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "ML")
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
# Model I
GAM_I <- gam(NASC_int ~
                s(year, bs = "re") +
                s(LME, bs = "re") +
                s(v, by = LME, k = 5, bs = "tp"),
             data = SA_df, family = Gamma(link = "log"), method = "ML")
```

Check model metrics


```r
# Summary metrics
GAM_AIC <- AIC(GAM_G, 
               GAM_S,
               GAM_I) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAM_G", 
                     "GAM_S",
                     "GAM_I"),
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
<div id="htmlwidget-edca92fbb7deaff0f976" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-edca92fbb7deaff0f976">{"x":{"filter":"none","vertical":false,"data":[["GAM_I","GAM_S","GAM_G"],[11.684,14.254,8.66],[64.99,64.68,54.76],[0.11,0.1,0.07],[409.69,414.38,417.69],[813.593,819.913,828.905],[0,6.32,15.31],[0.95887,0.04068,0.00045]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## NASC_int ~ s(v, k = 5, bs = "tp") + s(LME, bs = "re") + s(year, 
##     bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    5.064      1.051    4.82 9.65e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##           edf Ref.df      F  p-value    
## s(v)    2.219  2.687  2.873 0.050018 .  
## s(LME)  1.856  2.000 10.816 0.000525 ***
## s(year) 1.842  2.000 12.658 9.93e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0676   Deviance explained = 54.8%
## -ML = 417.69  Scale est. = 2.5171    n = 69
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
## NASC_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, LME, 
##     bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.8565     0.6576   7.386 6.74e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(year)  1.771      2 9.462 0.000189 ***
## s(LME)   1.713      2 6.622 0.001502 ** 
## s(v,LME) 6.758     14 4.614 0.000161 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0993   Deviance explained = 64.7%
## -ML = 414.38  Scale est. = 1.5572    n = 69
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
## NASC_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, by = LME, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.8858     0.6788   7.197 1.24e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df      F  p-value    
## s(year)    1.797  2.000 15.316 3.11e-05 ***
## s(LME)     1.807  2.000 10.777 0.000595 ***
## s(v):LMEBF 1.000  1.000  9.644 0.002918 ** 
## s(v):LMEBB 1.000  1.000  0.701 0.405804    
## s(v):LMESV 3.290  3.729 17.232  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.113   Deviance explained =   65%
## -ML = 409.69  Scale est. = 1.3595    n = 69
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
##          k'      edf   k-index p-value
## s(year)   3 1.771413        NA      NA
## s(LME)    3 1.713144        NA      NA
## s(v,LME) 15 6.757841 0.8987912  0.5575
```


```r
appraise(GAM_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_I)
```

```
##            k'      edf   k-index p-value
## s(year)     3 1.796624        NA      NA
## s(LME)      3 1.806629        NA      NA
## s(v):LMEBF  4 1.000014 0.8722884  0.4350
## s(v):LMEBB  4 1.000018 0.8722884  0.4725
## s(v):LMESV  4 3.290189 0.8722884  0.4550
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
            ggplot(aes(x = LME, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_S_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA),
          nrow = 1)
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
            ggplot(aes(x = LME, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_I_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA), 
          nrow = 1)
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
                s(year, bs = "re") +
                s(LME, bs = "re") +
                s(v, LME, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "REML", 
             select = T)
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
summary(GAM_S_p)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, LME, 
##     bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.8603     0.7466    6.51 1.98e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df      F  p-value    
## s(year)  1.824      2 10.286 0.000137 ***
## s(LME)   1.777      2  7.129 0.001504 ** 
## s(v,LME) 6.660     14  4.695 0.000449 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.103   Deviance explained = 64.6%
## -REML = 413.84  Scale est. = 1.5358    n = 69
```

No terms can be dropped, so I refit model with REML


```r
GAM_S2 <- gam(NASC_int ~ 
                s(year, bs = "re") +
                s(LME, bs = "re") +
                s(v, LME, bs = "fs", k = 5),
              data = SA_df, family = Gamma(link = "log"), method = "REML")
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
summary(GAM_S2)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, LME, 
##     bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.8603     0.7466    6.51 1.98e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df      F  p-value    
## s(year)  1.824      2 10.286 0.000137 ***
## s(LME)   1.777      2  7.129 0.001504 ** 
## s(v,LME) 6.660     14  4.695 0.000449 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.103   Deviance explained = 64.6%
## -REML = 413.84  Scale est. = 1.5358    n = 69
```

Check residuals.


```r
appraise(GAM_S2)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS3-residuals-covariates-1.png)<!-- -->

```r
GAM_S2_resid <- bind_cols(SA_df, residuals.gam(GAM_S2)) %>%
  rename(resid = "...9")
plot_grid(GAM_S2_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_S2_resid %>%
            ggplot(aes(x = LME, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_S2_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA),
          nrow = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS3-residuals-covariates-2.png)<!-- -->

### Model I


```r
# Model I
GAM_I_p <- gam(NASC_int ~
                 s(year, bs = "re") +
                 s(LME, bs = "re") +
                 s(v, by = LME, k = 5, bs = "tp"),
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
## NASC_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, by = LME, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.8762     0.7596   6.419 2.43e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                  edf Ref.df      F  p-value    
## s(year)    1.8423830      2 18.384 1.85e-05 ***
## s(LME)     1.8394519      2 10.577  0.00073 ***
## s(v):LMEBF 1.2945419      4  2.907  0.00601 ** 
## s(v):LMEBB 0.0008912      4  0.000  0.50650    
## s(v):LMESV 3.0691543      4 20.505 1.28e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.157   Deviance explained = 64.6%
## -REML = 412.21  Scale est. = 1.3573    n = 69
```

No terms can be dropped, so I refit model with REML.


```r
GAM_I2 <- gam(NASC_int ~ 
                 s(year, bs = "re") +
                 s(LME, bs = "re") +
                 s(v, by = LME, k = 5, bs = "tp"),
              data = SA_df, family = Gamma(link = "log"), method = "REML")
summary(GAM_I2)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, by = LME, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.8864     0.7507   6.509 1.82e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df      F  p-value    
## s(year)    1.830  2.000 16.701 2.67e-05 ***
## s(LME)     1.846  2.000 11.247 0.000734 ***
## s(v):LMEBF 1.000  1.001  9.609 0.002964 ** 
## s(v):LMEBB 1.000  1.000  0.707 0.403663    
## s(v):LMESV 3.298  3.735 17.443  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.114   Deviance explained = 65.1%
## -REML = 408.41  Scale est. = 1.3392    n = 69
```

Check residuals.


```r
appraise(GAM_I2)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI3-residuals-covariates-1.png)<!-- -->

```r
GAM_I2_resid <- bind_cols(SA_df, residuals.gam(GAM_I2)) %>%
  rename(resid = "...9")
plot_grid(GAM_I2_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_I2_resid %>%
            ggplot(aes(x = LME, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_I2_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA),
          nrow = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI3-residuals-covariates-2.png)<!-- -->

## Model predictions

I predict the model to get the smooths in the response scale (NASC). I need to create a specific data frame for this with the factors of interest `LME` and `year`, and covariate `v`.


```r
# Find median, min, and max of SA_df
val <- SA_df %>%
  dplyr::select(LME, year, xc, yc, v) %>%
  group_by(LME) %>%
  summarise(min_v = min(v),
            max_v = max(v),
            median_v = median(v))
# Resolution for predictions
res = 50
# LME areas
LME_area <- levels(SA_df$LME)
# Empty data frame that we populate with new values
SA_new <- data.frame()

for (i in LME_area) {
  # Select data from which we build the 
  val_tmp <- val %>% 
    filter(LME == i)
  # Temporary data frame for each region
  SA_tmp <- data.frame(
    LME = rep(i, res), 
    year = 2015,
    var = rep("v", res), 
    # Velocity
    v = seq(subset(val, LME == i)$min_v, 
            subset(val, LME == i)$max_v, 
            length.out = res))
  # Append data
  SA_new <- bind_rows(SA_new, SA_tmp)
}
# Remove unused variables
rm(val_tmp, SA_tmp, LME_area, res)
```

Get GAM predictions (`GAM_S` and `GAM_I`) for the new data. I then calculate the average functional response for all years.


```r
ilink_S <- family(GAM_S2)$linkinv # Get link function
pred_S <- predict(GAM_S2, SA_new, type = "link", se.fit = TRUE,
                  exclude = "s(year)") %>%
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilink_S(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilink_S(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilink_S(fit)) %>% # Calculate fit
  mutate(SA_fit = 10 * log10(fitted), # Convert to SA
         SA_lwr_ci = 10 * log10(lwr_ci),
         SA_upr_ci = 10 * log10(upr_ci))
```


```r
ilink_I <- family(GAM_I2)$linkinv # Get link function
pred_I <- predict(GAM_I2, SA_new, type = "link", se.fit = TRUE,
                  exclude = "s(year,LME)") %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilink_I(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilink_I(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilink_I(fit)) %>% # Calculate fit
  mutate(SA_fit = 10 * log10(fitted), # Convert to SA
         SA_lwr_ci = 10 * log10(lwr_ci),
         SA_upr_ci = 10 * log10(upr_ci))
```


```r
# Save data
save(pred_S, pred_I, GAM_S2, GAM_I2, file = "data/statistics/GAM_results.RData")
```

## Model visualization

Plot to see if predictions worked well.


```r
plot_grid(
  pred_S %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = fitted)) +
    geom_ribbon(aes(x = v, ymin = lwr_ci, ymax = upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = NASC_int, group = LME)) +
    facet_wrap(~ LME, scales = "free") +
    ggtitle("Model S"),
  pred_I %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = fitted)) +
    geom_ribbon(aes(x = v, ymin = lwr_ci, ymax = upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = NASC_int, group = LME)) + 
    facet_wrap(~ LME, scales = "free") +
    ggtitle("Model I"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/velo-ggplot-1.png)<!-- -->


```r
GAMS_velo <- pred_S %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = LME), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = LME),
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 40), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "NASC Model S",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAMI_velo <- pred_I %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = LME), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = LME), 
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 40), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "NASC Model I",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
# Combine plots
plot_grid(GAMS_velo, GAMI_velo,
          ncol = 1, rel_heights = c(1, 1))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS-pretty-ggplot-1.png)<!-- -->

```r
# Remove unused variables
rm(GAMS_velo, GAMI_velo)
```

# HGAM - Gaussian SA

## Data preparation

I select data at 318 m depth because it fits well with the DSL diurnal centre of mass. It also prevent from removing too much data like the deeper depth levels would do.


```r
SA_df <- stat_laea %>%
  filter(depth == 318) %>%
  dplyr::select(year, xc, yc, LME, NASC_int, SA_int, v, depth) %>%
  mutate(year = factor(year))
```


## Model fitting

In model S the random intercept is already included in the `s(bs = "fs")` term.


```r
# Model G
GAMSA_G <- gam(SA_int ~ 
                 s(v, k = 5, bs = "tp") + # Global smoothers
                 s(LME, bs = "re") + # Random effect
                 s(year, bs = "re"),# + # Random effect
               data = SA_df, family = "gaussian", method = "ML")
# Model S
GAMSA_S <- gam(SA_int ~
                s(year, bs = "re") +
                s(LME, bs = "re") +
                s(v, LME, bs = "fs", k = 5),
               data = SA_df, family = "gaussian", method = "ML")
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
# Model I
GAMSA_I <- gam(SA_int ~
                s(year, bs = "re") +
                s(LME, bs = "re") +
                s(v, by = LME, k = 5, bs = "tp"),
             data = SA_df, family = "gaussian", method = "ML")
```

Check model metrics


```r
# Summary metrics
GAMSA_AIC <- AIC(GAMSA_G, 
                 GAMSA_S,
                 GAMSA_I) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAMSA_G", 
                     "GAMSA_S",
                     "GAMSA_I"),
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
<div id="htmlwidget-c8e95c78ea54e43f1276" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-c8e95c78ea54e43f1276">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_S","GAMSA_G"],[9.987,10.734,7.428],[40.25,40.57,30.78],[0.34,0.34,0.26],[232.48,234.65,236.54],[467.81,468.933,472.842],[0,1.12,5.03],[0.60564,0.34543,0.04893]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## SA_int ~ s(v, k = 5, bs = "tp") + s(LME, bs = "re") + s(year, 
##     bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.952      2.126   7.974 3.71e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(v)    2.2344  2.696 2.690 0.142507    
## s(LME)  0.4572  2.000 0.364 0.219560    
## s(year) 1.6495  2.000 6.747 0.000728 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.261   Deviance explained = 30.8%
## -ML = 236.54  Scale est. = 48.432    n = 69
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
## SA_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, LME, bs = "fs", 
##     k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    16.59       2.10   7.896 6.34e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##             edf Ref.df     F  p-value    
## s(year)  1.6665      2 7.031 0.000559 ***
## s(LME)   0.4664      2 0.375 0.201885    
## s(v,LME) 4.4367     14 1.133 0.007712 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.342   Deviance explained = 40.6%
## -ML = 234.65  Scale est. = 43.092    n = 69
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
## SA_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, by = LME, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.857      2.402   7.017 2.14e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df     F p-value    
## s(year)    1.693  2.000 7.946  0.0005 ***
## s(LME)     1.218  2.000 2.530  0.0355 *  
## s(v):LMEBF 1.000  1.000 3.440  0.0685 .  
## s(v):LMEBB 1.000  1.000 2.759  0.1018    
## s(v):LMESV 1.964  2.357 3.561  0.0382 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.335   Deviance explained = 40.2%
## -ML = 232.48  Scale est. = 43.54     n = 69
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
##          k'       edf  k-index p-value
## s(year)   3 1.6665136       NA      NA
## s(LME)    3 0.4664448       NA      NA
## s(v,LME) 15 4.4366928 1.043272  0.5875
```


```r
appraise(GAMSA_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAMSA_I)
```

```
##            k'      edf k-index p-value
## s(year)     3 1.693057      NA      NA
## s(LME)      3 1.218332      NA      NA
## s(v):LMEBF  4 1.000028 1.01646  0.5150
## s(v):LMEBB  4 1.000005 1.01646  0.4900
## s(v):LMESV  4 1.963541 1.01646  0.5125
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
            ggplot(aes(x = LME, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAMSA_S_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA), 
          nrow = 1)
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
            ggplot(aes(x = LME, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAMSA_I_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA),
          nrow = 1)
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
                   s(year, bs = "re") +
                   s(LME, bs = "re") +
                   s(v, LME, bs = "fs", k = 5),
             data = SA_df, family = "gaussian", method = "REML", 
             select = T)
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
summary(GAMSA_S_p)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, LME, bs = "fs", 
##     k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.574      2.527   6.558  1.3e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(year)  1.773      2 7.525 0.000555 ***
## s(LME)   0.550      2 0.429 0.218692    
## s(v,LME) 4.387     14 1.136 0.009979 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.343   Deviance explained = 40.8%
## -REML = 232.91  Scale est. = 43.017    n = 69
```

`s(LME)` can be dropped. Refit the model without those terms.


```r
# Model S
GAMSA_S2 <- gam(SA_int ~ 
                  s(year, bs = "re") +
                  # s(LME, bs = "re") +
                  s(v, LME, bs = "fs", k = 5),
                data = SA_df, family = "gaussian", method = "ML")
summary(GAMSA_S2)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(v, LME, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.590      2.079   7.979 4.51e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(year)  1.668      2 6.979 0.000504 ***
## s(v,LME) 4.794     14 1.221 0.005164 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.34   Deviance explained = 40.3%
## -ML = 234.69  Scale est. = 43.208    n = 69
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
<div id="htmlwidget-7598e9be962aa164be3b" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-7598e9be962aa164be3b">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_S2","GAMSA_S"],[10.544,10.734],[40.3,40.57],[0.34,0.34],[234.69,234.65],[468.861,468.933],[0,0.07],[0.509,0.491]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

#### Final model 
Refit model with REML


```r
GAMSA_S3 <- gam(SA_int ~ 
                  s(year, bs = "re") +
                  # s(LME, bs = "re") +
                  s(v, LME, bs = "fs", k = 5),
                data = SA_df, family = "gaussian", method = "REML")
summary(GAMSA_S3)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(v, LME, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.577      2.517   6.585 1.16e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(year)  1.775      2 7.443 0.000501 ***
## s(v,LME) 4.821     14 1.255 0.006282 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.342   Deviance explained = 40.5%
## -REML = 232.95  Scale est. = 43.126    n = 69
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
            ggplot(aes(x = LME, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAMSA_S3_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA),
          nrow = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS3-residuals-covariates-2.png)<!-- -->

### Model I


```r
# Model I
GAMSA_I_p <- gam(SA_int ~
                   s(year, bs = "re") +
                   s(LME, bs = "re") +
                   s(v, by = LME, k = 5, bs = "tp"),
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
## SA_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, by = LME, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    16.51       2.61   6.325 3.17e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F  p-value    
## s(year)    1.7747      2 7.898 0.000483 ***
## s(LME)     0.9015      2 1.005 0.128378    
## s(v):LMEBF 0.7035      4 0.758 0.071297 .  
## s(v):LMEBB 1.1521      4 0.915 0.065444 .  
## s(v):LMESV 1.7363      4 1.970 0.028429 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.338   Deviance explained = 39.9%
## -REML = 233.11  Scale est. = 43.334    n = 69
```

I remove `s(LME)` from the model.


```r
# Model I
GAMSA_I2 <- gam(SA_int ~
                  s(year, bs = "re") +
                  # s(LME, bs = "re") +
                  s(v, by = LME, k = 5, bs = "tp"),
                data = SA_df, family = "gaussian", method = "ML")
summary(GAMSA_I2)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(v, by = LME, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.269      2.093   7.771 1.02e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df     F p-value    
## s(year)    1.668  2.000 6.945 0.00058 ***
## s(v):LMEBF 1.000  1.000 3.114 0.08254 .  
## s(v):LMEBB 1.811  2.258 2.486 0.10281    
## s(v):LMESV 1.897  2.285 3.504 0.03799 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.329   Deviance explained = 39.1%
## -ML = 233.39  Scale est. = 43.981    n = 69
```

Compare models


```r
# # Summary metrics
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
<div id="htmlwidget-f2241ef97423f77e75ea" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f2241ef97423f77e75ea">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_I2"],[9.987,9.484],[40.25,39.15],[0.34,0.33],[232.48,233.39],[467.81,468.061],[0,0.25],[0.53133,0.46867]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

#### Final model 

Refit model with REML


```r
GAMSA_I3 <- gam(SA_int ~
                  s(year, bs = "re") +
                  # s(LME, bs = "re") +
                  s(v, by = LME, k = 5, bs = "tp"),
                data = SA_df, family = "gaussian", method = "REML")
summary(GAMSA_I3)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(v, by = LME, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.073      2.473     6.5 1.68e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df     F  p-value    
## s(year)    1.762  2.000 7.029 0.000772 ***
## s(v):LMEBF 1.000  1.000 3.122 0.082239 .  
## s(v):LMEBB 2.189  2.681 2.013 0.097877 .  
## s(v):LMESV 2.176  2.605 3.572 0.033094 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.336   Deviance explained = 40.6%
## -REML = 225.68  Scale est. = 43.484    n = 69
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
            ggplot(aes(x = LME, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAMSA_I3_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA), 
          nrow = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI3-residuals-covariates-2.png)<!-- -->

## Model predictions

I predict the model to get the smooths in the response scale (NASC). I need to create a specific data frame for this with the factors of interest `LME` and `year`, and covariate `v`.


```r
# Find median, min, and max of SA_df
val <- SA_df %>%
  dplyr::select(LME, year, xc, yc, v) %>%
  group_by(LME) %>%
  summarise(min_v = min(v),
            max_v = max(v),
            median_v = median(v))
# Resolution for predictions
res = 50
# LME areas
LME_area <- levels(SA_df$LME)
# Empty data frame that we populate with new values
SA_new <- data.frame()

for (i in LME_area) {
  # Select data from which we build the 
  val_tmp <- val %>% 
    filter(LME == i)
  # Temporary data frame for each region
  SA_tmp <- data.frame(
    LME = rep(i, res), 
    year = 2015,
    var = rep("v", res), 
    # Velocity
    v = seq(subset(val, LME == i)$min_v, 
            subset(val, LME == i)$max_v, 
            length.out = res))
  # Append data
  SA_new <- bind_rows(SA_new, SA_tmp)
}
# Remove unused variables
rm(val_tmp, SA_tmp, LME_area, res)
```

Get GAM predictions (`GAM_S` and `GAM_I`) for the new data. I then calculate the average functional response for all years.


```r
ilinkSA_S <- family(GAMSA_S3)$linkinv # Get link function
predSA_S <- predict(GAMSA_S3, SA_new, type = "link", se.fit = TRUE,
                    exclude = "s(year)") %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilinkSA_S(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilinkSA_S(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilinkSA_S(fit)) # Calculate fit
```


```r
ilinkSA_I <- family(GAMSA_I3)$linkinv # Get link function
predSA_I <- predict(GAMSA_I3, SA_new, type = "link", se.fit = TRUE,
                    exclude = "s(year)") %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilinkSA_I(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilinkSA_I(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilinkSA_I(fit)) # Calculate fit
```


```r
# Save data
save(predSA_S, predSA_I, GAMSA_S3, GAMSA_I3, SA_df, 
     file = "data/statistics/GAMSA_results.RData")
```

## Model visualization

Plot to see if predictions worked well.


```r
plot_grid(
  predSA_S %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = fitted)) +
    geom_ribbon(aes(x = v, ymin = lwr_ci, ymax = upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = SA_int, group = LME)) +
    facet_wrap(~ LME, scales = "free") +
    ggtitle("SA Model S"),
  predSA_I %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = fitted)) +
    geom_ribbon(aes(x = v, ymin = lwr_ci, ymax = upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = SA_int, group = LME)) + 
    facet_wrap(~ LME, scales = "free") +
    ggtitle("SA Model I"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/veloSA-ggplot-1.png)<!-- -->


```r
GAMSAS_velo <- predSA_S %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = fitted, col = LME), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = lwr_ci, ymax = upr_ci, fill = LME),
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 30), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "SA Model S",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAMSAI_velo <- predSA_I %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = fitted, col = LME), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = lwr_ci, ymax = upr_ci, fill = LME), 
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 30), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(title = "SA Model I",
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
# Combine plots
plot_grid(GAMSAS_velo, GAMSAI_velo, ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS-pretty-ggplot-1.png)<!-- -->
