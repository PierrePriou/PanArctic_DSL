---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/30 at 17:31"
output: 
  html_document:
    code_folding: hide
    keep_md: yes
    toc: true
    toc_float: true
    toc_collapsed: true
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
cell_res <- 150
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
    geom_point(aes(col = chl_mean)) +
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
    geom_point(aes(col = chl_median)) +
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

This is because the coverage of acoustic data and remote sensing products slightly differ close to the coasts. For sea ice concentration (and open water duration) the NAs are close to land while for ocean colour they are in ice covered areas. I use the average value of neighbouring cells (24 cells around the pixel of interest) to fill missing values.


```r
stat_laea <- stat_laea %>%
  rowwise() %>%
  mutate( 
    # Ocean colour
    xc_v = if_else(is.na(v_mean) == T, xc, NaN),
    yc_v = if_else(is.na(v_mean) == T, yc, NaN),
    year_v = if_else(is.na(v_mean) == T, year, NaN), 
    v = if_else(
      is.na(v_mean) == T,
      mean(pull(subset(phy_grid_laea, 
                       xc >= xc_v - 1 * cell_res & xc <= xc_v + 1 * cell_res & 
                         yc >= yc_v - 1 * cell_res & yc <= yc_v + 1 * cell_res & 
                         year == year_v,
                       select = v_mean),
                v_mean),
           na.rm = T),
      v_mean),
    # Ocean colour
    xc_chl = if_else(is.na(chl_median) == T, xc, NaN),
    yc_chl = if_else(is.na(chl_median) == T, yc, NaN),
    year_chl = if_else(is.na(chl_median) == T, year, NaN), 
    chl = if_else(
      is.na(chl_mean) == T,
      mean(pull(subset(chl_grid_laea, 
                       xc >= xc_chl - 1 * cell_res & xc <= xc_chl + 1 * cell_res & 
                         yc >= yc_chl - 1 * cell_res & yc <= yc_chl + 1 * cell_res & 
                         year == year_chl,
                       select = chl_mean),
                chl_mean),
           na.rm = T),
      chl_mean),
    # Convert veloctity to cm/s
    v = v * 100, 
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
    geom_point(aes(col = v)) +
    scale_colour_cmocean(name = "speed", na.value = "red", limits = c(0, 5)) + 
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
  filter(is.na(chl) == F & is.na(v) == F)
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
    geom_point(aes(col = v)) +
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
## 1 WAO_BF      12
## 2 CAA         14
## 3 BB          26
## 4 DS           9
## 5 EAO         14
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
  geom_tile(aes(fill = chl_mean)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("chl (mg m-3)", name = "algae", limits = c(0,3)) +
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
  geom_tile(aes(fill = v_mean * 100)) +
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
  geom_tile(aes(fill = v_mean * 100)) +
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
  dplyr::select(year, IHO_area, NASC_int, SA_int, CM) %>%
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
##   .y.          n statistic    df         p method         var     
##   <chr>    <int>     <dbl> <int>     <dbl> <chr>          <chr>   
## 1 NASC_int    74     18.0      4 0.00122   Kruskal-Wallis IHO_area
## 2 NASC_int    74      9.98     2 0.00682   Kruskal-Wallis year    
## 3 CM          74     26.0      4 0.0000314 Kruskal-Wallis IHO_area
## 4 CM          74      3.23     2 0.199     Kruskal-Wallis year
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
## 1 WAO_BF   NASC_int    12     5.11      2 0.0777 Kruskal-Wallis
## 2 CAA      NASC_int    13     6.79      2 0.0336 Kruskal-Wallis
## 3 BB       NASC_int    26     2.01      2 0.366  Kruskal-Wallis
## 4 DS       NASC_int     9     0.167     2 0.92   Kruskal-Wallis
## 5 EAO      NASC_int    14     2.33      2 0.312  Kruskal-Wallis
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
##      Min       1Q   Median       3Q      Max 
## -20.0310  -4.9443   0.0176   4.3320  29.6689 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)  35.4415    15.4249   2.298   0.0245 *
## lat          -0.2344     0.2070  -1.132   0.2613  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 8.166 on 72 degrees of freedom
## Multiple R-squared:  0.01749,	Adjusted R-squared:  0.003847 
## F-statistic: 1.282 on 1 and 72 DF,  p-value: 0.2613
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
  dplyr::select(year, xc, yc, IHO_area, NASC_int, SA_int, v, chl) %>%
  group_by(IHO_area) %>%
  mutate(year = factor(year)) %>%
  ungroup() %>%
  rename(IHO = IHO_area) 
```

Plot data.


```r
plot_grid(SA_df %>%
            ggplot(aes(x = v, y = NASC_int, col = IHO)) +
            geom_point() +
            # geom_smooth(method = "gam", se = F, col = "grey20") +
            facet_wrap(~ IHO, scales = "free") +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = chl, y = NASC_int, col = IHO)) +
            geom_point() +
            # geom_smooth(method = "gam", se = F, col = "grey20") +
            facet_wrap(~ IHO, scales = "free") +
            theme(legend.position = "none"),
          ncol = 2)
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
               # te(xc, yc, by = year, k = 5) +
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
                # te(xc, yc, by = year, k = 5) +
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
                # te(xc, yc, by = year, k = 5) +
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
<div id="htmlwidget-104693b95379fec89712" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-104693b95379fec89712">{"x":{"filter":"none","vertical":false,"data":[["GAM_S","GAM_I","GAM_G"],[31.615,27.353,20.42],[85.49,81.37,76.68],[0.33,0.27,0.13],[454.78,448.44,458.9],[891.273,901.007,906.478],[0,9.73,15.2],[0.99187,0.00763,0.0005]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   4.5408     0.6788   6.689 9.14e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df      F  p-value    
## s(v)        2.729  3.218  7.314 0.000240 ***
## s(chl)      2.862  3.336  6.481 0.000530 ***
## s(IHO)      3.403  4.000 22.293 0.000511 ***
## s(year)     1.358  2.000  9.789 0.046074 *  
## s(IHO,year) 4.764 14.000  1.426 0.048061 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.127   Deviance explained = 76.7%
## -ML =  458.9  Scale est. = 1.2771    n = 75
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
## (Intercept)   5.0505     0.4463   11.31 1.96e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                   edf Ref.df     F  p-value    
## s(year)     1.862e-04      2 0.000  0.46077    
## s(IHO)      5.793e-05      4 0.000  0.03839 *  
## s(year,IHO) 4.933e+00     14 0.955  0.00552 ** 
## s(chl,IHO)  1.120e+01     24 4.695  < 2e-16 ***
## s(v,IHO)    7.561e+00     24 1.846 1.03e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.326   Deviance explained = 85.5%
## -ML = 454.78  Scale est. = 0.82125   n = 75
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
## (Intercept)   5.0453     0.5859   8.611 1.33e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                     edf Ref.df      F  p-value    
## s(year)          0.5521  2.000  0.826  0.30182    
## s(IHO)           2.7883  4.000 14.191  0.00030 ***
## s(year,IHO)      4.9114 14.000  1.014  0.03372 *  
## s(chl):IHOWAO_BF 1.0000  1.000  5.290  0.02549 *  
## s(chl):IHOCAA    1.0000  1.000  0.123  0.72762    
## s(chl):IHOBB     1.0000  1.000  3.554  0.06498 .  
## s(chl):IHODS     1.0000  1.000  0.000  0.99263    
## s(chl):IHOEAO    1.0000  1.000 23.691 1.14e-05 ***
## s(v):IHOWAO_BF   1.0000  1.000  7.314  0.00923 ** 
## s(v):IHOCAA      1.0000  1.000  0.049  0.82492    
## s(v):IHOBB       2.5877  3.083  5.366  0.00259 ** 
## s(v):IHODS       1.0000  1.000  1.195  0.27936    
## s(v):IHOEAO      2.8479  3.266 12.922 1.81e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.273   Deviance explained = 81.4%
## -ML = 448.44  Scale est. = 1.095     n = 75
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
##             k'          edf   k-index p-value
## s(year)      3 1.862448e-04        NA      NA
## s(IHO)       5 5.792909e-05        NA      NA
## s(year,IHO) 15 4.933321e+00        NA      NA
## s(chl,IHO)  25 1.120141e+01 0.8800380  0.3425
## s(v,IHO)    25 7.561252e+00 0.7491447  0.0350
```


```r
appraise(GAM_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_I)
```

```
##                  k'       edf   k-index p-value
## s(year)           3 0.5520565        NA      NA
## s(IHO)            5 2.7882951        NA      NA
## s(year,IHO)      15 4.9113733        NA      NA
## s(chl):IHOWAO_BF  4 1.0000002 0.8992914  0.4300
## s(chl):IHOCAA     4 1.0000076 0.8992914  0.4125
## s(chl):IHOBB      4 1.0000316 0.8992914  0.3700
## s(chl):IHODS      4 1.0000113 0.8992914  0.3900
## s(chl):IHOEAO     4 1.0000015 0.8992914  0.3650
## s(v):IHOWAO_BF    4 1.0000320 0.7773980  0.0600
## s(v):IHOCAA       4 1.0000098 0.7773980  0.0575
## s(v):IHOBB        4 2.5876860 0.7773980  0.0575
## s(v):IHODS        4 1.0000037 0.7773980  0.0725
## s(v):IHOEAO       4 2.8479169 0.7773980  0.0775
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
                # te(xc, yc, by = year, k = 5) +
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
## (Intercept)   4.9974     0.4873   10.26 6.79e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                   edf Ref.df     F  p-value    
## s(year)     6.798e-04      2 0.000  0.44519    
## s(IHO)      1.592e-04      4 0.000  0.07363 .  
## s(year,IHO) 5.021e+00     14 1.007  0.00453 ** 
## s(chl,IHO)  1.149e+01     24 4.857  < 2e-16 ***
## s(v,IHO)    7.461e+00     24 1.781 1.63e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.34   Deviance explained = 85.7%
## -REML = 454.56  Scale est. = 0.80754   n = 75
```

`s(IHO)` and `s(year)` can be dropped. Refit model without those terms.


```r
# Model S
GAM_S2 <- gam(NASC_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                # s(year, bs = "re") +
                # s(IHO, bs = "re") +
                # te(xc, yc, by = year, k = 5) +
                # Second order random effects
                s(year, IHO, bs = "re") +
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
## NASC_int ~ s(year, IHO, bs = "re") + s(chl, IHO, bs = "fs", k = 5) + 
##     s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   5.0505     0.4463   11.31 1.96e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df     F  p-value    
## s(year,IHO)  4.933     14 0.955  0.00552 ** 
## s(chl,IHO)  11.202     24 4.695  < 2e-16 ***
## s(v,IHO)     7.561     24 1.846 1.03e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.326   Deviance explained = 85.5%
## -ML = 454.78  Scale est. = 0.82123   n = 75
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
<div id="htmlwidget-5abecdc8ccda2cc6c25c" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-5abecdc8ccda2cc6c25c">{"x":{"filter":"none","vertical":false,"data":[["GAM_S","GAM_S2"],[31.615,31.616],[85.49,85.49],[0.33,0.33],[454.78,454.78],[891.273,891.274],[0,0],[0.50012,0.49988]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAM_S3 <- gam(NASC_int ~ 
                # First order effects
                # s(chl, bs = "tp", k = 5) +
                # s(v, bs = "tp", k = 5) +
                # s(year, bs = "re") +
                # s(IHO, bs = "re") +
                # te(xc, yc, by = year, k = 5) +
                # Second order random effects
                s(year, IHO, bs = "re") +
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
## NASC_int ~ s(year, IHO, bs = "re") + s(chl, IHO, bs = "fs", k = 5) + 
##     s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.9976     0.4873   10.26 6.79e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df     F  p-value    
## s(year,IHO)  5.022     14 1.007  0.00452 ** 
## s(chl,IHO)  11.490     24 4.857  < 2e-16 ***
## s(v,IHO)     7.461     24 1.781 1.63e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.34   Deviance explained = 85.7%
## -REML = 454.56  Scale est. = 0.80752   n = 75
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
                 # te(xc, yc, by = year, k = 5) +
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
## (Intercept)    4.911      0.585   8.395 2.16e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df      F  p-value    
## s(year)          0.0002099      2  0.000  0.42969    
## s(IHO)           3.5654745      4 19.663  0.00741 ** 
## s(year,IHO)      5.6284931     14  2.103  0.01617 *  
## s(chl):IHOWAO_BF 0.8113599      4  3.741  0.02132 *  
## s(chl):IHOCAA    0.0001485      4  0.000  0.85730    
## s(chl):IHOBB     1.3681335      4  1.947  0.09123 .  
## s(chl):IHODS     0.0001356      4  0.000  0.68266    
## s(chl):IHOEAO    3.0644557      4 88.978  < 2e-16 ***
## s(v):IHOWAO_BF   0.8606060      4  3.646  0.00668 ** 
## s(v):IHOCAA      0.0001942      4  0.000  0.94369    
## s(v):IHOBB       2.2107359      4  8.754 4.98e-06 ***
## s(v):IHODS       0.3527712      4  0.202  0.18380    
## s(v):IHOEAO      1.7483209      4  2.759  0.01191 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.534   Deviance explained = 85.4%
## -REML = 450.74  Scale est. = 0.74979   n = 75
```

`s(year)` can be dropped. Refit model without those terms.


```r
# Model I
GAM_I2 <- gam(NASC_int ~
                 # First order effects
                 # s(chl, bs = "tp", k = 5) +
                 # s(v, bs = "tp", k = 5) +
                 # s(year, bs = "re") +
                 s(IHO, bs = "re") +
                 # te(xc, yc, by = year, k = 5) +
                 # Second order random effects
                 s(year, IHO, bs = "re") +
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
## NASC_int ~ s(IHO, bs = "re") + s(year, IHO, bs = "re") + s(chl, 
##     by = IHO, k = 5, bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   5.0451     0.5616   8.984 3.49e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df      F  p-value    
## s(IHO)           2.744  4.000 14.405 0.000357 ***
## s(year,IHO)      5.486 14.000  1.135 0.009729 ** 
## s(chl):IHOWAO_BF 1.000  1.000  5.594 0.021788 *  
## s(chl):IHOCAA    1.000  1.000  0.105 0.747683    
## s(chl):IHOBB     1.000  1.000  3.291 0.075423 .  
## s(chl):IHODS     1.000  1.000  0.002 0.960552    
## s(chl):IHOEAO    1.000  1.000 24.237 9.42e-06 ***
## s(v):IHOWAO_BF   1.000  1.000  7.951 0.006787 ** 
## s(v):IHOCAA      1.000  1.000  0.054 0.817469    
## s(v):IHOBB       2.601  3.097  5.619 0.001951 ** 
## s(v):IHODS       1.000  1.000  1.161 0.286160    
## s(v):IHOEAO      2.864  3.282 13.381 1.32e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.273   Deviance explained = 81.3%
## -ML = 448.48  Scale est. = 1.072     n = 75
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
<div id="htmlwidget-c6f6c946d0e835ed2831" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-c6f6c946d0e835ed2831">{"x":{"filter":"none","vertical":false,"data":[["GAM_I","GAM_I2"],[27.353,27.388],[81.37,81.33],[0.27,0.27],[448.44,448.48],[901.007,901.286],[0,0.28],[0.53482,0.46518]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAM_I3 <- gam(NASC_int ~ 
                 # First order effects
                 # s(chl, bs = "tp", k = 5) +
                 # s(v, bs = "tp", k = 5) +
                 # s(year, bs = "re") +
                 s(IHO, bs = "re") +
                 # te(xc, yc, by = year, k = 5) +
                 # Second order random effects
                 s(year, IHO, bs = "re") +
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
## NASC_int ~ s(IHO, bs = "re") + s(year, IHO, bs = "re") + s(chl, 
##     by = IHO, k = 5, bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   5.0654     0.6204   8.165 1.06e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df      F  p-value    
## s(IHO)           3.094  4.000 12.533 8.58e-05 ***
## s(year,IHO)      5.273 14.000  1.304 0.002209 ** 
## s(chl):IHOWAO_BF 1.704  1.966  3.840 0.042829 *  
## s(chl):IHOCAA    1.000  1.000  0.135 0.714988    
## s(chl):IHOBB     1.592  1.938  3.477 0.065530 .  
## s(chl):IHODS     1.000  1.000  0.032 0.858228    
## s(chl):IHOEAO    3.180  3.408 20.905  < 2e-16 ***
## s(v):IHOWAO_BF   1.000  1.000 10.431 0.002215 ** 
## s(v):IHOCAA      1.000  1.000  0.054 0.816448    
## s(v):IHOBB       2.762  3.223  6.966 0.000509 ***
## s(v):IHODS       1.000  1.000  1.520 0.223455    
## s(v):IHOEAO      2.301  2.660  3.326 0.026850 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.526   Deviance explained = 86.2%
## -REML =  438.7  Scale est. = 0.7512    n = 75
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
    ggtitle("NASC Model S"),
  pred_I %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = NASC_fit)) +
    geom_ribbon(aes(x = chl, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("NASC Model I"),
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
  labs(title = "NASC Model S",
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
  labs(title = "NASC Model I",
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
                 # te(xc, yc, by = year, k = 5) +
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
                # te(xc, yc, by = year, k = 5) +
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
                # te(xc, yc, by = year, k = 5) +
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
<div id="htmlwidget-ce1f4087b524af687bde" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ce1f4087b524af687bde">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_G","GAMSA_S"],[27.435,16.322,23.739],[71.61,53.85,61.92],[0.6,0.46,0.52],[243.5,252.04,251.05],[486.813,501.027,501.453],[0,14.21,14.64],[0.99852,0.00082,0.00066]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   17.453      1.943   8.984 7.01e-13 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F  p-value    
## s(v)        2.175  2.619 3.024 0.071144 .  
## s(chl)      2.986  3.437 5.844 0.000957 ***
## s(IHO)      2.709  4.000 3.119 0.019731 *  
## s(year)     1.200  2.000 2.485 0.074106 .  
## s(IHO,year) 2.029 14.000 0.213 0.277891    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.457   Deviance explained = 53.8%
## -ML = 252.04  Scale est. = 35.987    n = 75
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
## (Intercept)   19.419      1.564   12.42   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                   edf Ref.df     F  p-value    
## s(year)     1.0653251      2 1.970 0.083698 .  
## s(IHO)      0.0001275      4 0.000 0.337403    
## s(year,IHO) 2.4773644     14 0.275 0.124203    
## s(chl,IHO)  5.8562273     24 1.063 0.000392 ***
## s(v,IHO)    6.3236830     24 1.258 0.000224 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.516   Deviance explained = 61.9%
## -ML = 251.05  Scale est. = 32.053    n = 75
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
## (Intercept)   19.158      1.763   10.87 5.57e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df     F  p-value    
## s(year)          7.306e-01  2.000 3.428 0.157866    
## s(IHO)           4.949e-05  4.000 0.000 0.197844    
## s(year,IHO)      7.240e+00 14.000 2.026 0.011269 *  
## s(chl):IHOWAO_BF 1.000e+00  1.000 7.625 0.007937 ** 
## s(chl):IHOCAA    1.000e+00  1.000 0.472 0.495066    
## s(chl):IHOBB     1.000e+00  1.000 3.428 0.069767 .  
## s(chl):IHODS     1.000e+00  1.000 0.031 0.861808    
## s(chl):IHOEAO    1.000e+00  1.000 8.733 0.004688 ** 
## s(v):IHOWAO_BF   1.000e+00  1.000 3.636 0.062069 .  
## s(v):IHOCAA      1.000e+00  1.000 0.012 0.911806    
## s(v):IHOBB       2.383e+00  2.877 3.275 0.027167 *  
## s(v):IHODS       1.000e+00  1.000 0.536 0.467464    
## s(v):IHOEAO      3.703e+00  3.894 6.208 0.000436 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.596   Deviance explained = 71.6%
## -ML =  243.5  Scale est. = 26.807    n = 75
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
## s(year)      3 1.0653251124        NA      NA
## s(IHO)       5 0.0001274537        NA      NA
## s(year,IHO) 15 2.4773644245        NA      NA
## s(chl,IHO)  25 5.8562272649 1.0488638  0.6275
## s(v,IHO)    25 6.3236830476 0.8341003  0.0475
```


```r
appraise(GAMSA_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAMSA_I)
```

```
##                  k'          edf   k-index p-value
## s(year)           3 7.305644e-01        NA      NA
## s(IHO)            5 4.949212e-05        NA      NA
## s(year,IHO)      15 7.239932e+00        NA      NA
## s(chl):IHOWAO_BF  4 1.000001e+00 0.9771774  0.3425
## s(chl):IHOCAA     4 1.000002e+00 0.9771774  0.3500
## s(chl):IHOBB      4 1.000002e+00 0.9771774  0.3450
## s(chl):IHODS      4 1.000002e+00 0.9771774  0.3825
## s(chl):IHOEAO     4 1.000000e+00 0.9771774  0.3600
## s(v):IHOWAO_BF    4 1.000005e+00 0.8523421  0.0675
## s(v):IHOCAA       4 1.000002e+00 0.8523421  0.0950
## s(v):IHOBB        4 2.382724e+00 0.8523421  0.0875
## s(v):IHODS        4 1.000001e+00 0.8523421  0.0825
## s(v):IHOEAO       4 3.702662e+00 0.8523421  0.0975
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
                # te(xc, yc, by = year, k = 5) +
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
## (Intercept)   19.389      1.772   10.95 9.32e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                  edf Ref.df     F  p-value    
## s(year)     1.348089      2 2.882 0.071413 .  
## s(IHO)      0.000171      4 0.000 0.309684    
## s(year,IHO) 2.381668     14 0.258 0.142743    
## s(chl,IHO)  5.671724     24 1.021 0.000484 ***
## s(v,IHO)    6.362679     24 1.293 0.000208 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.516   Deviance explained = 61.9%
## -REML = 249.63  Scale est. = 32.096    n = 75
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
                  # te(xc, yc, by = year, k = 5) +
                  # Second order random effects
                  s(year, IHO, bs = "re") +
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
## SA_int ~ s(year, bs = "re") + s(year, IHO, bs = "re") + s(chl, 
##     IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   19.419      1.564   12.42   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F  p-value    
## s(year)     1.065      2 1.970 0.083698 .  
## s(year,IHO) 2.477     14 0.275 0.124201    
## s(chl,IHO)  5.856     24 1.063 0.000392 ***
## s(v,IHO)    6.324     24 1.258 0.000224 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.516   Deviance explained = 61.9%
## -ML = 251.05  Scale est. = 32.053    n = 75
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
<div id="htmlwidget-a79eaadb426eb0dfb5ff" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-a79eaadb426eb0dfb5ff">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_S","GAMSA_S2"],[23.739,23.739],[61.92,61.92],[0.52,0.52],[251.05,251.05],[501.453,501.453],[0,0],[0.5,0.5]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAMSA_S3 <- gam(SA_int ~ 
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  s(year, bs = "re") +
                  # s(IHO, bs = "re") +
                  # te(xc, yc, by = year, k = 5) +
                  # Second order random effects
                  s(year, IHO, bs = "re") +
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
## SA_int ~ s(year, bs = "re") + s(year, IHO, bs = "re") + s(chl, 
##     IHO, bs = "fs", k = 5) + s(v, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   19.389      1.772   10.95 9.32e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F  p-value    
## s(year)     1.348      2 2.882 0.071412 .  
## s(year,IHO) 2.382     14 0.258 0.142735    
## s(chl,IHO)  5.672     24 1.021 0.000484 ***
## s(v,IHO)    6.363     24 1.293 0.000208 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.516   Deviance explained = 61.9%
## -REML = 249.63  Scale est. = 32.096    n = 75
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
                   # te(xc, yc, by = year, k = 5) +
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
## (Intercept)   18.836      2.098   8.979 1.19e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df     F p-value   
## s(year)          1.380e+00      2 4.754 0.07834 . 
## s(IHO)           2.641e+00      4 2.718 0.11265   
## s(year,IHO)      3.331e+00     14 0.514 0.17467   
## s(chl):IHOWAO_BF 8.165e-01      4 2.644 0.02226 * 
## s(chl):IHOCAA    1.744e-05      4 0.000 0.66686   
## s(chl):IHOBB     7.795e-01      4 1.921 0.03363 * 
## s(chl):IHODS     1.700e-05      4 0.000 0.93978   
## s(chl):IHOEAO    1.335e+00      4 4.194 0.00238 **
## s(v):IHOWAO_BF   6.348e-01      4 0.679 0.10202   
## s(v):IHOCAA      1.466e-05      4 0.000 0.99890   
## s(v):IHOBB       1.900e+00      4 3.106 0.00406 **
## s(v):IHODS       3.237e-05      4 0.000 0.47637   
## s(v):IHOEAO      1.898e+00      4 3.607 0.00264 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.552   Deviance explained = 64.1%
## -REML = 247.96  Scale est. = 29.717    n = 75
```

I remove `s(year)` from the model.


```r
# Model I
GAMSA_I2 <- gam(SA_int ~
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  # s(year, bs = "re") +
                  s(IHO, bs = "re") +
                  # te(xc, yc, by = year, k = 5) +
                  # Second order random effects
                  s(year, IHO, bs = "re") +
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
## SA_int ~ s(IHO, bs = "re") + s(year, IHO, bs = "re") + s(chl, 
##     by = IHO, k = 5, bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   19.056      1.577   12.08   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df     F  p-value    
## s(IHO)           1.417e-05  4.000 0.000 0.259204    
## s(year,IHO)      8.324e+00 14.000 2.285 0.000286 ***
## s(chl):IHOWAO_BF 1.000e+00  1.000 7.281 0.009374 ** 
## s(chl):IHOCAA    1.000e+00  1.000 0.412 0.523607    
## s(chl):IHOBB     1.000e+00  1.000 3.128 0.082821 .  
## s(chl):IHODS     1.000e+00  1.000 0.093 0.761421    
## s(chl):IHOEAO    1.000e+00  1.000 8.947 0.004241 ** 
## s(v):IHOWAO_BF   1.000e+00  1.000 3.558 0.064857 .  
## s(v):IHOCAA      1.000e+00  1.000 0.013 0.909243    
## s(v):IHOBB       2.412e+00  2.908 3.424 0.022698 *  
## s(v):IHODS       1.000e+00  1.000 0.478 0.492481    
## s(v):IHOEAO      3.719e+00  3.900 6.349 0.000358 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.599   Deviance explained = 72.1%
## -ML = 243.63  Scale est. = 26.561    n = 75
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
<div id="htmlwidget-03f12ba122db251569af" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-03f12ba122db251569af">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I2","GAMSA_I"],[26.87,27.435],[72.09,71.61],[0.6,0.6],[243.63,243.5],[484.412,486.813],[0,2.4],[0.76861,0.23139]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAMSA_I3 <- gam(SA_int ~
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  # s(year, bs = "re") +
                  s(IHO, bs = "re") +
                  # te(xc, yc, by = year, k = 5) +
                  # Second order random effects
                  s(year, IHO, bs = "re") +
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
## SA_int ~ s(IHO, bs = "re") + s(year, IHO, bs = "re") + s(chl, 
##     by = IHO, k = 5, bs = "tp") + s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   20.358      2.168   9.388 1.54e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df     F  p-value    
## s(IHO)           1.883  4.000 2.290 0.058709 .  
## s(year,IHO)      5.946 14.000 1.245 0.007802 ** 
## s(chl):IHOWAO_BF 1.915  2.193 4.190 0.020166 *  
## s(chl):IHOCAA    1.000  1.000 0.469 0.496795    
## s(chl):IHOBB     1.000  1.001 3.688 0.060636 .  
## s(chl):IHODS     1.000  1.000 0.005 0.946767    
## s(chl):IHOEAO    2.884  3.122 6.572 0.000681 ***
## s(v):IHOWAO_BF   1.000  1.000 4.453 0.039977 *  
## s(v):IHOCAA      1.000  1.000 0.004 0.951627    
## s(v):IHOBB       2.581  3.077 3.919 0.013577 *  
## s(v):IHODS       1.000  1.000 0.635 0.429303    
## s(v):IHOEAO      3.692  3.889 5.533 0.001476 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.659   Deviance explained = 77.4%
## -REML = 220.87  Scale est. = 22.575    n = 75
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
save(predSA_S, predSA_I, GAMSA_S3, GAMSA_I3, SA_df, file = "data/statistics/GAMSA_results.RData")
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
    ggtitle("SA Model S"),
  predSA_I %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = SA_fit)) +
    geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = SA_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("SA Model I"),
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
    ggtitle("SA Model S"),
  predSA_I %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = SA_fit)) +
    geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = SA_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free") +
    ggtitle("SA Model I"),
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
  labs(title = "SA Model S",
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
  labs(title = "SA Model I",
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
