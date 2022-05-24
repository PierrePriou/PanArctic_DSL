---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/24 at 19:13"
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
cell_res <- 100
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

I combine backscatter and environmental data which were gridded on the Lambert Azimuthal Equal Area grid (EPSG:6931) with a cell resolution of 100 x 100 km. Since I do not have a lot of samples from the the Beaufort Sea and West Arctic Ocean, I combine those two regions due to their geographical proximity.


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
    filter(depth == 318) %>%
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

This is because the coverage of acoustic data and remote sensing products slightly differ close to the coasts. For sea ice concentration (and open water duration) the NAs are close to land while for ocean colour they are in ice covered areas. I use the average value of neighbouring cells (8 cells around the pixel of interest) to fill missing values.


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
                       xc >= xc_v - cell_res & xc <= xc_v + cell_res & 
                         yc >= yc_v - cell_res & yc <= yc_v + cell_res & 
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
    filter(depth == 318) %>%
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
    filter(depth == 318) %>%
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
## 1 WAO_BF      14
## 2 CAA         13
## 3 BB          32
## 4 DS          16
## 5 EAO         21
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
##   .y.          n statistic    df           p method         var     
##   <chr>    <int>     <dbl> <int>       <dbl> <chr>          <chr>   
## 1 NASC_int    97     21.4      4 0.000265    Kruskal-Wallis IHO_area
## 2 NASC_int    97      8.97     2 0.0113      Kruskal-Wallis year    
## 3 CM          97     36.9      4 0.000000189 Kruskal-Wallis IHO_area
## 4 CM          97      1.03     2 0.596       Kruskal-Wallis year
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
## 1 WAO_BF   NASC_int    14     5.90      2 0.0523 Kruskal-Wallis
## 2 CAA      NASC_int    13     7.14      2 0.0281 Kruskal-Wallis
## 3 BB       NASC_int    32     1.38      2 0.502  Kruskal-Wallis
## 4 DS       NASC_int    16     0.176     2 0.916  Kruskal-Wallis
## 5 EAO      NASC_int    22     3.96      2 0.138  Kruskal-Wallis
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
## -18.895  -4.900  -0.493   4.221  33.162 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  45.9856    13.7259   3.350  0.00116 **
## lat          -0.3861     0.1834  -2.105  0.03797 * 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 8.806 on 95 degrees of freedom
## Multiple R-squared:  0.04454,	Adjusted R-squared:  0.03449 
## F-statistic: 4.429 on 1 and 95 DF,  p-value: 0.03797
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
SA_df %>%
  ggplot(aes(x = v, y = NASC_int, col = IHO)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, col = "grey20") +
  facet_wrap(~ IHO, scales = "free") +
  theme(legend.position = "none")
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

![](PanArctic_DSL_statistics_files/figure-html/plot-data-area-2.png)<!-- -->

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
               # Second order effects
               s(IHO, year, bs = "re"), # Random slope
            data = SA_df, family = Gamma(link = "log"), method = "REML")
# Model S
GAM_S <- gam(NASC_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "REML")
# Model I
GAM_I <- gam(NASC_int ~
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
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
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-3e399fe16efdabcd80a9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-3e399fe16efdabcd80a9">{"x":{"filter":"none","vertical":false,"data":[["GAM_I","GAM_S","GAM_G"],[27.173,30.313,18.42],[71.67,72.47,62.93],[0.42,0.47,0.28],[569.74,583.46,589.01],[1161.442,1164.578,1174.239],[0,3.14,12.8],[0.82636,0.17226,0.00138]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   4.8060     0.5546   8.665 3.52e-13 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df     F p-value  
## s(v)        1.0022  1.004 2.540  0.1148  
## s(chl)      2.1976  2.676 3.338  0.0196 *
## s(IHO)      2.0616  4.000 3.544  0.4802  
## s(year)     0.8597  2.000 4.195  0.3525  
## s(IHO,year) 7.4555 14.000 3.136  0.1254  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.284   Deviance explained = 62.9%
## -REML = 589.01  Scale est. = 2.2942    n = 96
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
## NASC_int ~ s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, 
##     bs = "re") + s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", 
##     k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.4361     0.4854    9.14 9.77e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                  edf Ref.df     F  p-value    
## s(year)     1.141854      2 5.739  0.22603    
## s(IHO)      0.001956      4 0.000  0.27082    
## s(year,IHO) 7.287030     14 1.698  0.03469 *  
## s(chl,IHO)  5.978133     24 1.795 3.85e-05 ***
## s(v,IHO)    7.234360     24 1.048  0.00941 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.468   Deviance explained = 72.5%
## -REML = 583.46  Scale est. = 1.7834    n = 96
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
## NASC_int ~ s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, 
##     bs = "re") + s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.5882     0.4057   11.31   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df      F  p-value    
## s(year)          1.211e+00  2.000  4.627  0.17965    
## s(IHO)           2.082e-05  4.000  0.000  0.71145    
## s(year,IHO)      5.470e+00 14.000  0.967  0.10059    
## s(chl):IHOWAO_BF 1.000e+00  1.001  9.581  0.00277 ** 
## s(chl):IHOCAA    1.000e+00  1.000  0.025  0.87374    
## s(chl):IHOBB     2.240e+00  2.673  4.027  0.01620 *  
## s(chl):IHODS     1.000e+00  1.000  0.015  0.90281    
## s(chl):IHOEAO    1.000e+00  1.000 33.197 3.09e-07 ***
## s(v):IHOWAO_BF   1.002e+00  1.003  2.500  0.11816    
## s(v):IHOCAA      1.000e+00  1.000  0.001  0.97013    
## s(v):IHOBB       2.572e+00  3.039  3.529  0.01810 *  
## s(v):IHODS       1.000e+00  1.000  0.774  0.38176    
## s(v):IHOEAO      2.540e+00  3.035  5.121  0.00272 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.416   Deviance explained = 71.7%
## -REML = 569.74  Scale est. = 1.6401    n = 96
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
##             k'         edf   k-index p-value
## s(year)      3 1.141853948        NA      NA
## s(IHO)       5 0.001956388        NA      NA
## s(year,IHO) 15 7.287029664        NA      NA
## s(chl,IHO)  25 5.978132860 0.8048709  0.2125
## s(v,IHO)    25 7.234359596 0.9153094  0.6525
```


```r
appraise(GAM_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_I)
```

```
##                  k'          edf   k-index p-value
## s(year)           3 1.211474e+00        NA      NA
## s(IHO)            5 2.081778e-05        NA      NA
## s(year,IHO)      15 5.469597e+00        NA      NA
## s(chl):IHOWAO_BF  4 1.000293e+00 0.7655921  0.0900
## s(chl):IHOCAA     4 1.000044e+00 0.7655921  0.1025
## s(chl):IHOBB      4 2.240052e+00 0.7655921  0.0850
## s(chl):IHODS      4 1.000036e+00 0.7655921  0.0900
## s(chl):IHOEAO     4 1.000068e+00 0.7655921  0.0775
## s(v):IHOWAO_BF    4 1.001700e+00 0.9262131  0.6475
## s(v):IHOCAA       4 1.000064e+00 0.9262131  0.7175
## s(v):IHOBB        4 2.572119e+00 0.9262131  0.6600
## s(v):IHODS        4 1.000082e+00 0.9262131  0.6775
## s(v):IHOEAO       4 2.540037e+00 0.9262131  0.7100
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
acf(resid(GAM_S), lag.max = 36, main = "ACF")
pacf(resid(GAM_S), lag.max = 36, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS-acf-pacf-1.png)<!-- -->


```r
par(mfrow = c(1, 2))
acf(resid(GAM_I), lag.max = 36, main = "ACF")
pacf(resid(GAM_I), lag.max = 36, main = "pACF")
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
## NASC_int ~ s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, 
##     bs = "re") + s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", 
##     k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.4361     0.4854    9.14 9.77e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                  edf Ref.df     F  p-value    
## s(year)     1.141854      2 5.739  0.22603    
## s(IHO)      0.001956      4 0.000  0.27082    
## s(year,IHO) 7.287030     14 1.698  0.03469 *  
## s(chl,IHO)  5.978133     24 1.795 3.85e-05 ***
## s(v,IHO)    7.234360     24 1.048  0.00941 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.468   Deviance explained = 72.5%
## -REML = 583.46  Scale est. = 1.7834    n = 96
```

The random effects `s(IHO)` and `s(year)` also have low edf but I decide to keep them in the model as their interaction `s(year,IHO)` is significant.


```r
# Model S
GAM_S2 <- gam(NASC_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                # s(IHO, bs = "re") +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = Gamma(link = "log"), method = "REML")
summary(GAM_S2)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(year, bs = "re") + s(year, IHO, bs = "re") + s(chl, 
##     IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.4360     0.4854    9.14 9.78e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F  p-value    
## s(year)     1.142      2 5.739  0.22601    
## s(year,IHO) 7.287     14 1.698  0.03465 *  
## s(chl,IHO)  5.978     24 1.794 3.85e-05 ***
## s(v,IHO)    7.236     24 1.048  0.00941 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.468   Deviance explained = 72.5%
## -REML = 583.46  Scale est. = 1.7834    n = 96
```

Compare models.


```r
GAM_S_AIC <- AIC(GAM_S, GAM_S2) %>% 
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
  full_join(., GAM_S_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 10)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-2e4e01734d173c1301ad" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2e4e01734d173c1301ad">{"x":{"filter":"none","vertical":false,"data":[["GAM_S","GAM_S2"],[30.313,30.313],[72.47,72.47],[0.47,0.47],[583.46,583.46],[1164.578,1164.578],[0,0],[0.5,0.5]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

`GAM_S` has better metics. So I keep that model.

Check residuals.


```r
appraise(GAM_S2)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS2-residuals-covariates-1.png)<!-- -->

```r
GAM_S2_resid <- bind_cols(SA_df, residuals.gam(GAM_S2)) %>%
  rename(resid = "...9")
plot_grid(GAM_S2_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_S2_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_S2_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_S2_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMS2-residuals-covariates-2.png)<!-- -->

### Model I


```r
# Model I
GAM_I_p <- gam(NASC_int ~
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  s(year, bs = "re") +
                  s(IHO, bs = "re") +
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
## NASC_int ~ s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, 
##     bs = "re") + s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.4219     0.4768   9.274 3.45e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df      F p-value    
## s(year)          1.397e+00      2 13.718 0.11500    
## s(IHO)           5.722e-04      4  0.000 0.54837    
## s(year,IHO)      7.835e+00     14  2.157 0.05252 .  
## s(chl):IHOWAO_BF 8.523e-01      4  4.939 0.01582 *  
## s(chl):IHOCAA    8.521e-05      4  0.000 0.95020    
## s(chl):IHOBB     1.803e+00      4  4.012 0.01256 *  
## s(chl):IHODS     6.154e-05      4  0.000 0.98680    
## s(chl):IHOEAO    2.309e+00      4 13.330 5.1e-06 ***
## s(v):IHOWAO_BF   4.750e-03      4  0.001 0.47033    
## s(v):IHOCAA      1.146e-04      4  0.000 0.86653    
## s(v):IHOBB       1.915e+00      4  2.909 0.00786 ** 
## s(v):IHODS       7.346e-04      4  0.000 0.40779    
## s(v):IHOEAO      1.616e+00      4  2.208 0.01496 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.59   Deviance explained = 72.5%
## -REML = 580.88  Scale est. = 1.5869    n = 96
```

The random effects `s(IHO)` and `s(year)` also have low edf but I decide to keep them in the model as their interaction `s(year,IHO)` is significant. I refit the model without the `s(v,year)` and `s(chl,year)` terms.


```r
# Model I
GAM_I2 <- gam(NASC_int ~
                 # First order effects
                 # s(chl, bs = "tp", k = 5) +
                 # s(v, bs = "tp", k = 5) +
                 # s(year, bs = "re") +
                 # s(IHO, bs = "re") +
                 # Second order random effects
                 # s(year, IHO, bs = "re") +
                 # Second order functional effects
                 s(chl, by = IHO, k = 5, bs = "tp") + 
                 s(v, by = IHO, k = 5, bs = "tp"),
               data = SA_df, family = Gamma(link = "log"), method = "REML")
summary(GAM_I2)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.9909     0.2401   20.79   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df      F  p-value    
## s(chl):IHOWAO_BF 1.779  2.222  8.414 0.000353 ***
## s(chl):IHOCAA    1.000  1.000  0.397 0.530262    
## s(chl):IHOBB     2.258  2.719  3.431 0.026217 *  
## s(chl):IHODS     1.000  1.000  0.106 0.745617    
## s(chl):IHOEAO    1.000  1.001 35.288  < 2e-16 ***
## s(v):IHOWAO_BF   1.000  1.001  5.370 0.023044 *  
## s(v):IHOCAA      1.000  1.000  0.008 0.930078    
## s(v):IHOBB       2.362  2.835  1.849 0.167846    
## s(v):IHODS       1.000  1.000  0.476 0.492372    
## s(v):IHOEAO      2.431  2.928  3.655 0.024274 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.65   Deviance explained = 65.4%
## -REML = 571.55  Scale est. = 2.396     n = 96
```

Compare models.


```r
GAM_I_AIC <- AIC(GAM_I, GAM_I2) %>% 
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
  full_join(., GAM_I_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 10)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-ad952f3a2b53cc250ddd" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ad952f3a2b53cc250ddd">{"x":{"filter":"none","vertical":false,"data":[["GAM_I","GAM_I2"],[27.173,18.707],[71.67,65.4],[0.42,0.65],[569.74,571.55],[1161.442,1166.726],[0,5.28],[0.9335161995,0.0664838005]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

`GAM_I` has better metics. So I keep that model.

Check residuals.


```r
appraise(GAM_I2)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI2-residuals-covariates-1.png)<!-- -->

```r
GAM_I2_resid <- bind_cols(SA_df, residuals.gam(GAM_I2)) %>%
  rename(resid = "...9")
plot_grid(GAM_S2_resid %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          GAM_S2_resid %>%
            ggplot(aes(x = chl, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "chlorophyll") +
            geom_point(), 
          GAM_S2_resid %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          GAM_S2_resid %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI2-residuals-covariates-2.png)<!-- -->

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
## # A tibble: 5 × 7
##   IHO    min_v max_v median_v min_chl max_chl median_chl
##   <fct>  <dbl> <dbl>    <dbl>   <dbl>   <dbl>      <dbl>
## 1 WAO_BF 0.606  4.32     2.49  0.110    1.07       0.251
## 2 CAA    0.633  2.64     1.02  0.309    1.07       0.423
## 3 BB     0.297  3.97     1.75  0.309    1.82       0.855
## 4 DS     0.598  2.94     1.60  0.318    0.748      0.432
## 5 EAO    0.292  4.66     1.99  0.0892   2.20       1.04
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
}
# Reorder factors
SA_new <- SA_new %>%
  mutate(IHO = factor(IHO, levels = c("WAO_BF","CAA", "BB", "DS", "EAO")))
```

Get GAM predictions (`GAM_S` and `GAM_I`) for the new data. I then calculate the average functional response for all years.


```r
ilink_S <- family(GAM_S)$linkinv # Get link function
pred_S <- predict(GAM_S, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
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
ilink_I <- family(GAM_I)$linkinv # Get link function
pred_I <- predict(GAM_I, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
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
save(pred_S, pred_I, GAM_S, GAM_I, file = "data/statistics/GAM_results.RData")
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
               # Second order effects
               s(IHO, year, bs = "re"), # Random slope
            data = SA_df, family = "gaussian", method = "REML")
# Model S
GAMSA_S <- gam(SA_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, IHO, bs = "fs", k = 5) +
                s(v, IHO, bs = "fs", k = 5),
             data = SA_df, family = "gaussian", method = "REML")
# Model I
GAMSA_I <- gam(SA_int ~
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                s(year, bs = "re") +
                s(IHO, bs = "re") +
                # Second order random effects
                s(year, IHO, bs = "re") +
                # Second order functional effects
                s(chl, by = IHO, k = 5, bs = "tp") + 
                s(v, by = IHO, k = 5, bs = "tp"),
             data = SA_df, family = "gaussian", method = "REML")
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
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-db3958cb65cac0525503" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-db3958cb65cac0525503">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_S","GAMSA_G"],[25.064,23.479,18.886],[54.24,51.99,36.74],[0.43,0.43,0.27],[306.06,331.77,336.93],[666.285,667.732,685.024],[0,1.45,18.74],[0.67334,0.3266,6e-05]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
##     s(IHO, bs = "re") + s(year, bs = "re") + s(IHO, year, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   17.042      1.859   9.166 3.33e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df     F p-value
## s(v)        2.2031  2.696 1.549   0.296
## s(chl)      2.5769  3.066 1.085   0.307
## s(IHO)      1.1639  4.000 0.660   0.437
## s(year)     0.7162  2.000 1.279   0.343
## s(IHO,year) 6.2709 14.000 1.136   0.120
## 
## R-sq.(adj) =  0.268   Deviance explained = 36.7%
## -REML = 336.93  Scale est. = 58.039    n = 96
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
## SA_int ~ s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, 
##     bs = "re") + s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", 
##     k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.563      2.405   6.886  1.2e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                 edf Ref.df     F p-value   
## s(year)     1.49144      2 6.315  0.0722 . 
## s(IHO)      0.00077      4 0.000  0.0998 . 
## s(year,IHO) 4.27786     14 0.531  0.1334   
## s(chl,IHO)  7.46194     24 1.956  0.0018 **
## s(v,IHO)    2.19864     24 0.139  0.1680   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.427   Deviance explained =   52%
## -REML = 331.77  Scale est. = 45.433    n = 96
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
## SA_int ~ s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, 
##     bs = "re") + s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.693      2.215   7.537 8.58e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                       edf Ref.df     F p-value   
## s(year)          1.453625  2.000 5.030 0.10368   
## s(IHO)           0.001634  4.000 0.000 0.38556   
## s(year,IHO)      3.930004 14.000 0.564 0.13353   
## s(chl):IHOWAO_BF 1.000024  1.000 5.794 0.01851 * 
## s(chl):IHOCAA    1.000020  1.000 0.013 0.91079   
## s(chl):IHOBB     1.879422  2.277 2.436 0.09292 . 
## s(chl):IHODS     1.000031  1.000 0.315 0.57643   
## s(chl):IHOEAO    2.555437  3.039 5.786 0.00128 **
## s(v):IHOWAO_BF   1.000061  1.000 1.477 0.22806   
## s(v):IHOCAA      1.000025  1.000 0.238 0.62702   
## s(v):IHOBB       2.428626  2.901 2.949 0.03928 * 
## s(v):IHODS       1.000043  1.000 0.312 0.57807   
## s(v):IHOEAO      1.000054  1.000 0.944 0.33434   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.426   Deviance explained = 54.2%
## -REML = 306.06  Scale est. = 45.483    n = 96
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
##             k'          edf   k-index p-value
## s(year)      3 1.4914351710        NA      NA
## s(IHO)       5 0.0007699674        NA      NA
## s(year,IHO) 15 4.2778591432        NA      NA
## s(chl,IHO)  25 7.4619396229 0.9794999  0.3750
## s(v,IHO)    25 2.1986389681 1.0385734  0.6025
```


```r
appraise(GAMSA_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAMSA_I)
```

```
##                  k'         edf   k-index p-value
## s(year)           3 1.453625176        NA      NA
## s(IHO)            5 0.001633903        NA      NA
## s(year,IHO)      15 3.930004344        NA      NA
## s(chl):IHOWAO_BF  4 1.000023655 0.9507955  0.2600
## s(chl):IHOCAA     4 1.000019798 0.9507955  0.2600
## s(chl):IHOBB      4 1.879421884 0.9507955  0.3200
## s(chl):IHODS      4 1.000031144 0.9507955  0.2875
## s(chl):IHOEAO     4 2.555436938 0.9507955  0.2575
## s(v):IHOWAO_BF    4 1.000060612 1.0398818  0.6450
## s(v):IHOCAA       4 1.000024510 1.0398818  0.5975
## s(v):IHOBB        4 2.428625579 1.0398818  0.5875
## s(v):IHODS        4 1.000042979 1.0398818  0.5925
## s(v):IHOEAO       4 1.000053799 1.0398818  0.6475
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
acf(resid(GAMSA_S), lag.max = 36, main = "ACF")
pacf(resid(GAMSA_S), lag.max = 36, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS-acf-pacf-1.png)<!-- -->


```r
par(mfrow = c(1, 2))
acf(resid(GAMSA_I), lag.max = 36, main = "ACF")
pacf(resid(GAMSA_I), lag.max = 36, main = "pACF")
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
## SA_int ~ s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, 
##     bs = "re") + s(chl, IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", 
##     k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.563      2.405   6.886  1.2e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                 edf Ref.df     F p-value   
## s(year)     1.49144      2 6.315  0.0722 . 
## s(IHO)      0.00077      4 0.000  0.0998 . 
## s(year,IHO) 4.27786     14 0.531  0.1334   
## s(chl,IHO)  7.46194     24 1.956  0.0018 **
## s(v,IHO)    2.19864     24 0.139  0.1680   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.427   Deviance explained =   52%
## -REML = 331.77  Scale est. = 45.433    n = 96
```

The random effects `s(IHO)` and `s(year)` also have low edf but I decide to keep them in the model as their interaction `s(year,IHO)` is significant.

### Model I


```r
# Model I
GAMSA_I_p <- gam(SA_int ~
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  s(year, bs = "re") +
                  s(IHO, bs = "re") +
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
## SA_int ~ s(year, bs = "re") + s(IHO, bs = "re") + s(year, IHO, 
##     bs = "re") + s(chl, by = IHO, k = 5, bs = "tp") + s(v, by = IHO, 
##     k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.662      2.306   7.226 2.44e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df     F  p-value    
## s(year)          1.532e+00      2 9.234 0.065534 .  
## s(IHO)           9.567e-01      4 0.436 0.403143    
## s(year,IHO)      5.200e+00     14 0.958 0.135065    
## s(chl):IHOWAO_BF 7.400e-01      4 2.050 0.053229 .  
## s(chl):IHOCAA    2.519e-05      4 0.000 0.849089    
## s(chl):IHOBB     1.524e+00      4 2.269 0.030548 *  
## s(chl):IHODS     2.955e-05      4 0.000 0.763772    
## s(chl):IHOEAO    1.962e+00      4 6.365 0.000592 ***
## s(v):IHOWAO_BF   3.455e-04      4 0.000 0.336614    
## s(v):IHOCAA      3.693e-05      4 0.000 0.635799    
## s(v):IHOBB       1.859e+00      4 2.542 0.008396 ** 
## s(v):IHODS       5.135e-05      4 0.000 0.507434    
## s(v):IHOEAO      7.085e-02      4 0.019 0.299247    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.461   Deviance explained = 53.9%
## -REML = 329.39  Scale est. = 42.731    n = 96
```

The random effects `s(IHO)` and `s(year)` also have low edf but I decide to keep them in the model as their interaction `s(year,IHO)` is significant. I refit the model without the `s(v,year)` and `s(chl,year)` terms.

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
## # A tibble: 5 × 7
##   IHO    min_v max_v median_v min_chl max_chl median_chl
##   <fct>  <dbl> <dbl>    <dbl>   <dbl>   <dbl>      <dbl>
## 1 WAO_BF 0.606  4.32     2.49  0.110    1.07       0.251
## 2 CAA    0.633  2.64     1.02  0.309    1.07       0.423
## 3 BB     0.297  3.97     1.75  0.309    1.82       0.855
## 4 DS     0.598  2.94     1.60  0.318    0.748      0.432
## 5 EAO    0.292  4.66     1.99  0.0892   2.20       1.04
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
}
# Reorder factors
SA_new <- SA_new %>%
  mutate(IHO = factor(IHO, levels = c("WAO_BF","CAA", "BB", "DS", "EAO")))
```

Get GAM predictions (`GAM_S` and `GAM_I`) for the new data. I then calculate the average functional response for all years.


```r
ilinkSA_S <- family(GAMSA_S)$linkinv # Get link function
predSA_S <- predict(GAMSA_S, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
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
ilinkSA_I <- family(GAMSA_I)$linkinv # Get link function
predSA_I <- predict(GAMSA_I, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
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
save(predSA_S, predSA_I, GAMSA_S, GAMSA_I, file = "data/statistics/GAMSA_results.RData")
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
