---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/06/01 at 15:11"
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

This is because the coverage of acoustic data and remote sensing products slightly differ close to the coasts. For sea ice concentration (and open water duration) the NAs are close to land while for ocean colour they are in ice covered areas. I use the average value of neighbouring cells (8 cells around the pixel of interest) to fill missing values.


```r
# stat_laea <- stat_laea %>%
#   rowwise() %>%
#   mutate( 
#     # Ocean colour
#     xc_v = if_else(is.na(v_mean) == T, xc, NaN),
#     yc_v = if_else(is.na(v_mean) == T, yc, NaN),
#     year_v = if_else(is.na(v_mean) == T, year, NaN), 
#     v = if_else(
#       is.na(v_mean) == T,
#       mean(pull(subset(phy_grid_laea, 
#                        xc >= xc_v - 1 * cell_res & xc <= xc_v + 1 * cell_res & 
#                          yc >= yc_v - 1 * cell_res & yc <= yc_v + 1 * cell_res & 
#                          year == year_v,
#                        select = v_mean),
#                 v_mean),
#            na.rm = T),
#       v_mean),
#     # Ocean colour
#     xc_chl = if_else(is.na(chl_median) == T, xc, NaN),
#     yc_chl = if_else(is.na(chl_median) == T, yc, NaN),
#     year_chl = if_else(is.na(chl_median) == T, year, NaN), 
#     chl = if_else(
#       is.na(chl_mean) == T,
#       mean(pull(subset(chl_grid_laea, 
#                        xc >= xc_chl - 1 * cell_res & xc <= xc_chl + 1 * cell_res & 
#                          yc >= yc_chl - 1 * cell_res & yc <= yc_chl + 1 * cell_res & 
#                          year == year_chl,
#                        select = chl_mean),
#                 chl_mean),
#            na.rm = T),
#       chl_mean),
#     # Convert veloctity to cm/s
#     v = v * 100, 
#     IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
#                                 IHO_area == "West Arctic Ocean" ~ "WAO_BF",
#                                 IHO_area == "Beaufort Sea" ~ "WAO_BF",
#                                 IHO_area == "The Northwestern Passages" ~ "CAA",
#                                 IHO_area == "Baffin Bay" ~ "BB",
#                                 IHO_area == "Davis Strait" ~ "DS"),
#                       levels = c("WAO_BF", "CAA", "BB", "DS", "EAO"))) %>%
#   ungroup() %>%
#   dplyr::select(-xc_chl, -yc_chl, -year_chl)
```

Check if I replace NAs correctly.


```r
# plot_grid(
#   stat_laea %>% 
#     ggplot(aes(x = xc,  y = yc)) +
#     geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
#                  fill = "grey80") +
#     geom_point(aes(col = chl)) +
#     scale_colour_cmocean(name = "algae", na.value = "red") + 
#     facet_wrap(~ year) +
#     coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
#     theme(axis.text = element_blank(), 
#           axis.ticks = element_blank(), 
#           axis.title = element_blank()),
#   stat_laea %>%
#     filter(depth == 318) %>%
#     ggplot(aes(x = xc,  y = yc)) +
#     geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
#                  fill = "grey80") +
#     geom_point(aes(col = v)) +
#     scale_colour_cmocean(name = "speed", na.value = "red", limits = c(0, 5)) + 
#     facet_wrap(~ year) +
#     coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
#     theme(axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.title = element_blank()),
#   ncol = 1, align = "hv", axis = "tblr")
```

There are still one missing chl values, so this cell is excluded from further analyses.


```r
stat_laea <- stat_laea %>%
  mutate(chl = chl_mean, 
         v = v_mean * 100,
         IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                     IHO_area == "West Arctic Ocean" ~ "WAO_BF",
                                     IHO_area == "Beaufort Sea" ~ "WAO_BF",
                                     IHO_area == "The Northwestern Passages" ~ "CAA",
                                     IHO_area == "Baffin Bay" ~ "BB",
                                     IHO_area == "Davis Strait" ~ "DS"),
                           levels = c("WAO_BF", "CAA", "BB", "DS", "EAO"))) %>%
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
## 3 BB          25
## 4 DS           9
## 5 EAO         12
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
# # Model S
# GAM_S <- gam(NASC_int ~ 
#                 # First order effects
#                 # s(chl, bs = "tp", k = 5) +
#                 # s(v, bs = "tp", k = 5) +
#                 s(year, bs = "re") +
#                 s(IHO, bs = "re") +
#                 # te(xc, yc, by = year, k = 5) +
#                 # Second order random effects
#                 s(year, IHO, bs = "re") +
#                 # Second order functional effects
#                 s(chl, IHO, bs = "fs", k = 5) +
#                 s(v, IHO, bs = "fs", k = 5),
#              data = SA_df, family = Gamma(link = "log"), method = "ML")
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
GAM_AIC <- AIC(GAM_G, 
               # GAM_S,
               GAM_I) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAM_G", 
                     # "GAM_S", 
                     "GAM_I"),
           reml = round(c(GAM_G$gcv.ubre,
                          # GAM_S$gcv.ubre, 
                          GAM_I$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAM_G$deviance / GAM_G$null.deviance)) * 100,
             # (1 - (GAM_S$deviance / GAM_S$null.deviance)) * 100,
             (1 - (GAM_I$deviance / GAM_I$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAM_G)$r.sq,
                        # summary(GAM_S)$r.sq, 
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
<div id="htmlwidget-b716e09b8c8f1cab09a5" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-b716e09b8c8f1cab09a5">{"x":{"filter":"none","vertical":false,"data":[["GAM_G","GAM_I"],[20.423,28.454],[78.8,82.33],[0.15,0.39],[441.22,435.04],[870.794,872.156],[0,1.36],[0.66396,0.33604]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   4.6069     0.6896   6.681 1.18e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df      F  p-value    
## s(v)        2.574  3.053  7.028 0.000407 ***
## s(chl)      2.762  3.254  7.034 0.000351 ***
## s(IHO)      3.427  4.000 27.730 0.000377 ***
## s(year)     1.362  2.000 11.087 0.053071 .  
## s(IHO,year) 5.030 14.000  1.727 0.038599 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.15   Deviance explained = 78.8%
## -ML = 441.22  Scale est. = 1.1591    n = 72
```

Model S summary.


```r
summary(GAM_S)
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
## (Intercept)   4.9284     0.5628   8.758 1.62e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df      F  p-value    
## s(year)          0.0001687  2.000  0.000 0.222470    
## s(IHO)           2.4076505  4.000 11.882 0.009368 ** 
## s(year,IHO)      7.4165714 14.000  2.231 0.001937 ** 
## s(chl):IHOWAO_BF 1.0000041  1.000  4.791 0.033499 *  
## s(chl):IHOCAA    1.0000035  1.000  0.217 0.643708    
## s(chl):IHOBB     1.0000127  1.000  0.987 0.325481    
## s(chl):IHODS     1.0000063  1.000  0.019 0.892162    
## s(chl):IHOEAO    1.0000126  1.000 21.454  2.8e-05 ***
## s(v):IHOWAO_BF   1.0000157  1.000  6.466 0.014276 *  
## s(v):IHOCAA      1.0000053  1.000  0.029 0.866136    
## s(v):IHOBB       2.6591574  3.161  4.850 0.004548 ** 
## s(v):IHODS       1.0000020  1.000  1.075 0.304904    
## s(v):IHOEAO      2.5256034  2.914  6.782 0.000489 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.393   Deviance explained = 82.3%
## -ML = 435.04  Scale est. = 1.1133    n = 72
```

## Model checking

### Basis size and residual distribution

First I check the basis size k. `k-indexes` are \> 1 or close to 1 so the basis size is large enough. The residual plot look good too.


```r
appraise(GAM_S, method = "simulate")
k.check(GAM_S)
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
## s(year)           3 0.0001686593        NA      NA
## s(IHO)            5 2.4076505437        NA      NA
## s(year,IHO)      15 7.4165714240        NA      NA
## s(chl):IHOWAO_BF  4 1.0000040953 0.8254481  0.1950
## s(chl):IHOCAA     4 1.0000034978 0.8254481  0.1400
## s(chl):IHOBB      4 1.0000126627 0.8254481  0.1925
## s(chl):IHODS      4 1.0000063020 0.8254481  0.1225
## s(chl):IHOEAO     4 1.0000126173 0.8254481  0.1375
## s(v):IHOWAO_BF    4 1.0000157092 0.8258448  0.1725
## s(v):IHOCAA       4 1.0000052920 0.8258448  0.1700
## s(v):IHOBB        4 2.6591573911 0.8258448  0.1700
## s(v):IHODS        4 1.0000019600 0.8258448  0.1675
## s(v):IHOEAO       4 2.5256034072 0.8258448  0.1925
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

Check residuals.


```r
appraise(GAM_S3)
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
## (Intercept)   4.8209     0.5765   8.362 3.55e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df      F  p-value    
## s(year)          5.615e-01      2  1.362  0.30921    
## s(IHO)           3.439e+00      4 26.838  0.00895 ** 
## s(year,IHO)      5.896e+00     14  2.548  0.03981 *  
## s(chl):IHOWAO_BF 7.860e-01      4  2.834  0.03309 *  
## s(chl):IHOCAA    1.397e-04      4  0.000  0.76736    
## s(chl):IHOBB     5.709e-01      4  0.411  0.25092    
## s(chl):IHODS     1.765e-04      4  0.000  0.73482    
## s(chl):IHOEAO    2.626e+00      4 32.752  < 2e-16 ***
## s(v):IHOWAO_BF   8.346e-01      4  2.971  0.01480 *  
## s(v):IHOCAA      9.846e-05      4  0.000  0.96274    
## s(v):IHOBB       2.313e+00      4 13.206 4.48e-07 ***
## s(v):IHODS       3.467e-01      4  0.214  0.19453    
## s(v):IHOEAO      1.868e+00      4  5.021  0.00239 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.51   Deviance explained = 84.8%
## -REML = 437.89  Scale est. = 0.80125   n = 72
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
## (Intercept)   4.9284     0.5627   8.758 1.62e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df      F  p-value    
## s(IHO)           2.408  4.000 11.882 0.009369 ** 
## s(year,IHO)      7.417 14.000  2.231 0.001936 ** 
## s(chl):IHOWAO_BF 1.000  1.000  4.791 0.033498 *  
## s(chl):IHOCAA    1.000  1.000  0.217 0.643711    
## s(chl):IHOBB     1.000  1.000  0.987 0.325489    
## s(chl):IHODS     1.000  1.000  0.019 0.892146    
## s(chl):IHOEAO    1.000  1.000 21.454  2.8e-05 ***
## s(v):IHOWAO_BF   1.000  1.000  6.466 0.014275 *  
## s(v):IHOCAA      1.000  1.000  0.029 0.866136    
## s(v):IHOBB       2.659  3.161  4.850 0.004548 ** 
## s(v):IHODS       1.000  1.000  1.075 0.304908    
## s(v):IHOEAO      2.526  2.914  6.782 0.000489 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.393   Deviance explained = 82.3%
## -ML = 435.04  Scale est. = 1.1133    n = 72
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
<div id="htmlwidget-bd2b9f167bb955f27421" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bd2b9f167bb955f27421">{"x":{"filter":"none","vertical":false,"data":[["GAM_I","GAM_I2"],[28.454,28.454],[82.33,82.33],[0.39,0.39],[435.04,435.04],[872.156,872.156],[0,0],[0.5,0.5]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   5.0527     0.6143   8.225 1.26e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df      F  p-value    
## s(IHO)           3.154  4.000 14.976 4.64e-05 ***
## s(year,IHO)      5.077 14.000  1.240  0.00283 ** 
## s(chl):IHOWAO_BF 1.696  1.955  3.766  0.03974 *  
## s(chl):IHOCAA    1.000  1.000  0.128  0.72167    
## s(chl):IHOBB     1.699  2.060  1.588  0.23112    
## s(chl):IHODS     1.000  1.000  0.052  0.82100    
## s(chl):IHOEAO    3.194  3.415 21.171  < 2e-16 ***
## s(v):IHOWAO_BF   1.000  1.000 10.495  0.00220 ** 
## s(v):IHOCAA      1.000  1.000  0.056  0.81367    
## s(v):IHOBB       2.673  3.137  4.433  0.00831 ** 
## s(v):IHODS       1.000  1.000  1.505  0.22608    
## s(v):IHOEAO      1.946  2.261  1.551  0.19389    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.476   Deviance explained = 85.9%
## -REML = 425.16  Scale est. = 0.74513   n = 72
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
save(#pred_S,
     pred_I, 
     #GAM_S3, 
     GAM_I3, file = "data/statistics/GAM_results.RData")
```

## Model visualization

Plot to see if predictions worked well.


```r
plot_grid(
  # pred_S %>%
  #   filter(var == "v") %>%
  #   ggplot() +
  #   geom_line(aes(x = v, y = NASC_fit)) +
  #   geom_ribbon(aes(x = v, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
  #   geom_point(data = SA_df, aes(x = v, y = NASC_int, group = IHO)) + 
  #   facet_wrap(~ IHO, scales = "free") +
  #   ggtitle("Model S"),
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
  # pred_S %>%
  #   filter(var == "chl") %>%
  #   ggplot() +
  #   geom_line(aes(x = chl, y = NASC_fit)) +
  #   geom_ribbon(aes(x = chl, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
  #   geom_point(data = SA_df, aes(x = chl, y = NASC_int, group = IHO)) +
  #   facet_wrap(~ IHO, scales = "free") +
  #   ggtitle("NASC Model S"),
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
# # Model S
# GAMSA_S <- gam(SA_int ~ 
#                 # First order effects
#                 # s(chl, bs = "tp", k = 5) +
#                 # s(v, bs = "tp", k = 5) +
#                 s(year, bs = "re") +
#                 s(IHO, bs = "re") +
#                 # te(xc, yc, by = year, k = 5) +
#                 # Second order random effects
#                 s(year, IHO, bs = "re") +
#                 # Second order functional effects
#                 s(chl, IHO, bs = "fs", k = 5) +
#                 s(v, IHO, bs = "fs", k = 5),
#                data = SA_df, family = "gaussian", method = "ML")
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
GAMSA_AIC <- AIC(GAMSA_G, 
                 # GAMSA_S,
                 GAMSA_I) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAMSA_G", 
                     # "GAMSA_S", 
                     "GAMSA_I"),
           reml = round(c(GAMSA_G$gcv.ubre,
                          # GAMSA_S$gcv.ubre, 
                          GAMSA_I$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAMSA_G$deviance / GAMSA_G$null.deviance)) * 100,
             # (1 - (GAMSA_S$deviance / GAMSA_S$null.deviance)) * 100,
             (1 - (GAMSA_I$deviance / GAMSA_I$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAMSA_G)$r.sq,
                        # summary(GAMSA_S)$r.sq, 
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
<div id="htmlwidget-6010a8e3a3b228e3f1dc" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-6010a8e3a3b228e3f1dc">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_G"],[27.177,15.186],[70.77,56.02],[0.58,0.48],[234.3,239.68],[470.856,476.3],[0,5.44],[0.93831,0.06169]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   17.721      2.183   8.117 2.89e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F  p-value    
## s(v)        1.686  2.052 1.512 0.231594    
## s(chl)      2.847  3.329 6.669 0.000425 ***
## s(IHO)      2.932  4.000 3.904 0.007285 ** 
## s(year)     1.455  2.000 4.206 0.027220 *  
## s(IHO,year) 1.450 14.000 0.142 0.314281    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.485   Deviance explained =   56%
## -ML = 239.68  Scale est. = 34.037    n = 72
```

Model S summary.


```r
summary(GAMSA_S)
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
## (Intercept)   19.167      2.139   8.963 6.22e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                     edf Ref.df      F p-value   
## s(year)          1.1335  2.000  6.732 0.10375   
## s(IHO)           0.8571  4.000  0.713 0.14866   
## s(year,IHO)      5.9099 14.000  1.510 0.03737 * 
## s(chl):IHOWAO_BF 1.0000  1.000  6.052 0.01746 * 
## s(chl):IHOCAA    1.0000  1.000  0.479 0.49224   
## s(chl):IHOBB     1.0000  1.000  2.864 0.09694 . 
## s(chl):IHODS     1.0000  1.000  0.007 0.93572   
## s(chl):IHOEAO    1.0000  1.000 11.211 0.00157 **
## s(v):IHOWAO_BF   1.0000  1.000  2.935 0.09302 . 
## s(v):IHOCAA      1.0000  1.000  0.011 0.91576   
## s(v):IHOBB       2.2403  2.740  2.290 0.08709 . 
## s(v):IHODS       1.0000  1.000  0.563 0.45676   
## s(v):IHOEAO      3.4445  3.758  4.520 0.00543 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.58   Deviance explained = 70.8%
## -ML =  234.3  Scale est. = 27.752    n = 72
```

## Model checking

### Basis size and residual distribution

First I check the basis size k. `k-indexes` are \> 1 or close to 1 so the basis size is large enough. The residual plot look good too.


```r
appraise(GAMSA_S, method = "simulate")
k.check(GAMSA_S)
```


```r
appraise(GAMSA_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAMSA_I)
```

```
##                  k'       edf   k-index p-value
## s(year)           3 1.1335299        NA      NA
## s(IHO)            5 0.8570963        NA      NA
## s(year,IHO)      15 5.9099077        NA      NA
## s(chl):IHOWAO_BF  4 1.0000002 0.9541349  0.3050
## s(chl):IHOCAA     4 1.0000013 0.9541349  0.3050
## s(chl):IHOBB      4 1.0000162 0.9541349  0.3000
## s(chl):IHODS      4 1.0000023 0.9541349  0.2900
## s(chl):IHOEAO     4 1.0000018 0.9541349  0.3225
## s(v):IHOWAO_BF    4 1.0000035 0.9397667  0.2950
## s(v):IHOCAA       4 1.0000020 0.9397667  0.2850
## s(v):IHOBB        4 2.2402550 0.9397667  0.2800
## s(v):IHODS        4 1.0000010 0.9397667  0.2750
## s(v):IHOEAO       4 3.4445105 0.9397667  0.2725
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

Check residuals.


```r
appraise(GAMSA_S3)
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
## (Intercept)    18.72       2.54   7.372 6.65e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df     F p-value   
## s(year)          1.678e+00      2 7.805 0.00831 **
## s(IHO)           2.983e+00      4 3.394 0.02478 * 
## s(year,IHO)      1.697e+00     14 0.188 0.26911   
## s(chl):IHOWAO_BF 7.697e-01      4 1.858 0.04034 * 
## s(chl):IHOCAA    1.544e-05      4 0.000 0.56392   
## s(chl):IHOBB     1.906e+00      4 6.237 0.00265 **
## s(chl):IHODS     6.982e-06      4 0.000 0.99542   
## s(chl):IHOEAO    1.263e+00      4 2.945 0.00463 **
## s(v):IHOWAO_BF   5.440e-01      4 0.451 0.14214   
## s(v):IHOCAA      1.011e-05      4 0.000 0.98972   
## s(v):IHOBB       2.314e-04      4 0.000 0.53039   
## s(v):IHODS       2.057e-05      4 0.000 0.45989   
## s(v):IHOEAO      1.669e+00      4 1.783 0.03077 * 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.531   Deviance explained = 61.4%
## -REML = 237.91  Scale est. = 30.995    n = 72
```

I remove `s(year)` and `s(IHO)` from the model.


```r
# Model I
GAMSA_I2 <- gam(SA_int ~
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  # s(year, bs = "re") +
                  # s(IHO, bs = "re") +
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
## SA_int ~ s(year, IHO, bs = "re") + s(chl, by = IHO, k = 5, bs = "tp") + 
##     s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   18.937      1.654   11.45 2.29e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df      F p-value    
## s(year,IHO)      8.738 14.000  2.613 0.00013 ***
## s(chl):IHOWAO_BF 1.000  1.000  6.323 0.01532 *  
## s(chl):IHOCAA    1.000  1.000  0.415 0.52257    
## s(chl):IHOBB     1.000  1.000  1.927 0.17148    
## s(chl):IHODS     1.000  1.000  0.123 0.72710    
## s(chl):IHOEAO    1.000  1.000 11.714 0.00128 ** 
## s(v):IHOWAO_BF   1.000  1.000  3.137 0.08288 .  
## s(v):IHOCAA      1.000  1.000  0.014 0.90670    
## s(v):IHOBB       2.371  2.878  2.750 0.05117 .  
## s(v):IHODS       1.000  1.000  0.442 0.50956    
## s(v):IHOEAO      3.578  3.829  4.673 0.00301 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.595   Deviance explained = 72.5%
## -ML = 234.53  Scale est. = 26.746    n = 72
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
<div id="htmlwidget-a559d13e43398e0bf731" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-a559d13e43398e0bf731">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I2","GAMSA_I"],[27.073,27.177],[72.46,70.77],[0.6,0.58],[234.53,234.3],[466.367,470.856],[0,4.49],[0.90418,0.09582]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAMSA_I3 <- gam(SA_int ~
                  # First order effects
                  # s(chl, bs = "tp", k = 5) +
                  # s(v, bs = "tp", k = 5) +
                  # s(year, bs = "re") +
                  # s(IHO, bs = "re") +
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
## SA_int ~ s(year, IHO, bs = "re") + s(chl, by = IHO, k = 5, bs = "tp") + 
##     s(v, by = IHO, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   20.322      1.653   12.29 3.99e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df     F  p-value    
## s(year,IHO)      7.907 14.000 1.871 0.000988 ***
## s(chl):IHOWAO_BF 2.145  2.441 4.358 0.013596 *  
## s(chl):IHOCAA    1.000  1.000 0.679 0.414014    
## s(chl):IHOBB     1.479  1.770 2.088 0.188561    
## s(chl):IHODS     1.000  1.000 0.001 0.979688    
## s(chl):IHOEAO    2.573  2.863 7.201 0.003677 ** 
## s(v):IHOWAO_BF   1.000  1.000 4.919 0.031537 *  
## s(v):IHOCAA      1.000  1.000 0.019 0.892066    
## s(v):IHOBB       2.381  2.864 1.872 0.131340    
## s(v):IHODS       1.000  1.000 0.672 0.416697    
## s(v):IHOEAO      3.586  3.827 3.971 0.006817 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.64   Deviance explained = 76.7%
## -REML = 212.66  Scale est. = 23.808    n = 72
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
save(#predSA_S,
     predSA_I,
     # GAMSA_S3,
     GAMSA_I3, 
     SA_df, file = "data/statistics/GAMSA_results.RData")
```

## Model visualization

Plot to see if predictions worked well.


```r
plot_grid(
  # predSA_S %>%
  #   filter(var == "v") %>%
  #   ggplot() +
  #   geom_line(aes(x = v, y = SA_fit)) +
  #   geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci), alpha = 0.1) +
  #   geom_point(data = SA_df, aes(x = v, y = SA_int, group = IHO)) + 
  #   facet_wrap(~ IHO, scales = "free") +
  #   ggtitle("SA Model S"),
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
  # predSA_S %>%
  #   filter(var == "chl") %>%
  #   ggplot() +
  #   geom_line(aes(x = chl, y = SA_fit)) +
  #   geom_ribbon(aes(x = chl, ymin = SA_lwr_ci, ymax = SA_upr_ci), alpha = 0.1) +
  #   geom_point(data = SA_df, aes(x = chl, y = SA_int, group = IHO)) + 
  #   facet_wrap(~ IHO, scales = "free") +
  #   ggtitle("SA Model S"),
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
