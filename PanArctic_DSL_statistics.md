---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/08/17 at 17:35"
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
    filter(depth == 380) %>%
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
  filter(depth == 380) %>%
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
## 4 BB       31
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
  filter(depth == 380) %>%
  group_by(LME) %>% 
  summarise(n = n())
```

```
## # A tibble: 3 × 2
##   LME       n
##   <fct> <int>
## 1 BF       21
## 2 BB       31
## 3 SV       15
```


```r
plot_grid(
  stat_laea %>%
    filter(depth == 380) %>%
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
    filter(depth == 380) %>%
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
  filter(depth == 380) %>%
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = v_mean * 100)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("velo (cm/s)", name = "speed", limits = c(0, 8)) +
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

Nares Strait data is excluded because of the low number of data points from that region.


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
  filter(LME != "NAR") %>%
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
  mutate(var = c("LME", "year", "LME", "year")) %>%
  kbl() %>%
  kable_classic()
```

<table class=" lightable-classic" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> .y. </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:right;"> statistic </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:left;"> method </th>
   <th style="text-align:left;"> var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 7.002972 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.030200 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
   <td style="text-align:left;"> LME </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 17.629804 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.000149 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
   <td style="text-align:left;"> year </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CM </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 7.945393 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.018800 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
   <td style="text-align:left;"> LME </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CM </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 1.541519 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.463000 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
   <td style="text-align:left;"> year </td>
  </tr>
</tbody>
</table>

I use Dunn's test to know which group are significantly different.


```r
dunn_test(SA_diff, NASC_int ~ LME) %>%
  kbl() %>%
  kable_classic()
```

<table class=" lightable-classic" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> .y. </th>
   <th style="text-align:left;"> group1 </th>
   <th style="text-align:left;"> group2 </th>
   <th style="text-align:right;"> n1 </th>
   <th style="text-align:right;"> n2 </th>
   <th style="text-align:right;"> statistic </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:right;"> p.adj </th>
   <th style="text-align:left;"> p.adj.signif </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:left;"> BF </td>
   <td style="text-align:left;"> BB </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 2.555379 </td>
   <td style="text-align:right;"> 0.0106072 </td>
   <td style="text-align:right;"> 0.0318217 </td>
   <td style="text-align:left;"> * </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:left;"> BF </td>
   <td style="text-align:left;"> SV </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.686710 </td>
   <td style="text-align:right;"> 0.4922655 </td>
   <td style="text-align:right;"> 0.4922655 </td>
   <td style="text-align:left;"> ns </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:left;"> BB </td>
   <td style="text-align:left;"> SV </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> -1.551510 </td>
   <td style="text-align:right;"> 0.1207795 </td>
   <td style="text-align:right;"> 0.2415591 </td>
   <td style="text-align:left;"> ns </td>
  </tr>
</tbody>
</table>

## Within areas

I check whether there are inter-annual variability in S\~A\~ across years within areas. This will be used to decide whether a random effect is needed in the GAM. I used a non parametric Kruskal-Wallis test to test for interannual variability within areas.

### NASC


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(LME) %>%
  kruskal_test(NASC_int ~ year) %>%
  kbl() %>%
  kable_classic()
```

<table class=" lightable-classic" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> LME </th>
   <th style="text-align:left;"> .y. </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:right;"> statistic </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:left;"> method </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> BF </td>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 14.308637 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.000781 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BB </td>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 4.514738 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.105000 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SV </td>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 3.617143 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.164000 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
  </tr>
</tbody>
</table>
Dunn test for BF.


```r
SA_diff %>%
  filter(LME == "BF") %>%
  dunn_test(NASC_int ~ year) %>%
  kbl() %>%
  kable_classic()
```

<table class=" lightable-classic" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> .y. </th>
   <th style="text-align:left;"> group1 </th>
   <th style="text-align:left;"> group2 </th>
   <th style="text-align:right;"> n1 </th>
   <th style="text-align:right;"> n2 </th>
   <th style="text-align:right;"> statistic </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:right;"> p.adj </th>
   <th style="text-align:left;"> p.adj.signif </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:left;"> 2015 </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1.419030 </td>
   <td style="text-align:right;"> 0.1558902 </td>
   <td style="text-align:right;"> 0.1558902 </td>
   <td style="text-align:left;"> ns </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:left;"> 2015 </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 3.781775 </td>
   <td style="text-align:right;"> 0.0001557 </td>
   <td style="text-align:right;"> 0.0004671 </td>
   <td style="text-align:left;"> *** </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NASC_int </td>
   <td style="text-align:left;"> 2016 </td>
   <td style="text-align:left;"> 2017 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1.903094 </td>
   <td style="text-align:right;"> 0.0570282 </td>
   <td style="text-align:right;"> 0.1140564 </td>
   <td style="text-align:left;"> ns </td>
  </tr>
</tbody>
</table>

### CM


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(LME) %>%
  kruskal_test(CM ~ year) %>%
  kbl() %>%
  kable_classic()
```

<table class=" lightable-classic" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> LME </th>
   <th style="text-align:left;"> .y. </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:right;"> statistic </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:left;"> method </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> BF </td>
   <td style="text-align:left;"> CM </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 4.384086 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.112 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BB </td>
   <td style="text-align:left;"> CM </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 1.252238 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.535 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SV </td>
   <td style="text-align:left;"> CM </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 2.935238 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.230 </td>
   <td style="text-align:left;"> Kruskal-Wallis </td>
  </tr>
</tbody>
</table>

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
##     Min      1Q  Median      3Q     Max 
## -26.990  -5.764   1.233   6.155  31.556 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  50.2048    18.9078   2.655  0.00983 **
## lat          -0.4696     0.2532  -1.855  0.06790 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 9.937 on 69 degrees of freedom
## Multiple R-squared:  0.04749,	Adjusted R-squared:  0.03369 
## F-statistic:  3.44 on 1 and 69 DF,  p-value: 0.0679
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
  mutate(mycto_TF = if_else(xc == -2025 & yc == -1275, T, 
                 if_else(xc == -2025 & yc == -1125, T, 
                 if_else(xc == -1875 & yc == -1125, T, 
                 if_else(xc == -1875 & yc == -975, T, 
                 if_else(xc == -1725 & yc == -825, T, F)))))) %>%
  ggplot() +
  # Coastlines
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), 
               fill = "grey25") + # Coast
  geom_tile(aes(x = xc, y = yc, fill = mycto_TF), col = "white") +
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
  mutate(mycto_TF = if_else(xc == -2025 & yc == -1275, T, 
                    if_else(xc == -2025 & yc == -1125, T, 
                    if_else(xc == -1875 & yc == -1125, T, 
                    if_else(xc == -1875 & yc == -975, T, 
                    if_else(xc == -1725 & yc == -825, T, F)))))) 
# Fit LM
LM_BB_all <- SA_BB %>%
  # filter(mycto_TF == T) %>%
  lm(SA_int ~ lat, data = .)
summary(LM_BB_all)
```

```
## 
## Call:
## lm(formula = SA_int ~ lat, data = .)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -16.7948  -2.2964  -0.4967   2.6743   9.2521 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept) 12.60417   19.97055   0.631    0.533
## lat          0.08967    0.27646   0.324    0.748
## 
## Residual standard error: 4.822 on 30 degrees of freedom
## Multiple R-squared:  0.003495,	Adjusted R-squared:  -0.02972 
## F-statistic: 0.1052 on 1 and 30 DF,  p-value: 0.7479
```

```r
# Fit LM
LM_BB_myctoT <- SA_BB %>%
  filter(mycto_TF == T) %>%
  lm(SA_int ~ lat, data = .)
summary(LM_BB_myctoT)
```

```
## 
## Call:
## lm(formula = SA_int ~ lat, data = .)
## 
## Residuals:
##        1        2        3        4        5        6 
## -0.08945  2.12569  1.46837 -2.35473 -0.98407 -0.16581 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)  -9.8014    37.0480  -0.265    0.804
## lat           0.3816     0.5290   0.721    0.511
## 
## Residual standard error: 1.818 on 4 degrees of freedom
## Multiple R-squared:  0.1151,	Adjusted R-squared:  -0.1061 
## F-statistic: 0.5203 on 1 and 4 DF,  p-value: 0.5106
```

```r
# Fit LM
LM_BB_myctoF <- SA_BB %>%
  filter(mycto_TF == F) %>%
  lm(SA_int ~ lat, data = .)
summary(LM_BB_myctoF)
```

```
## 
## Call:
## lm(formula = SA_int ~ lat, data = .)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -16.8409  -2.6286  -0.0738   3.9164   8.8008 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept) 22.80138   23.55328   0.968    0.343
## lat         -0.04441    0.32382  -0.137    0.892
## 
## Residual standard error: 5.205 on 24 degrees of freedom
## Multiple R-squared:  0.0007832,	Adjusted R-squared:  -0.04085 
## F-statistic: 0.01881 on 1 and 24 DF,  p-value: 0.8921
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
save(SA_BB, LM_BB_all, LM_BB_myctoT, LM_BB_myctoF,
     file = "data/statistics/LM_BB_results.RData")
```

# HGAM - Gamma NASC

## Data preparation

I select data at 380 m depth because it fits well with the DSL diurnal centre of mass. It also prevent from removing too much data like the deeper depth levels would do.


```r
SA_df <- stat_laea %>%
  filter(depth == 380) %>%
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
<div id="htmlwidget-21347e10553fdef5dd17" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-21347e10553fdef5dd17">{"x":{"filter":"none","vertical":false,"data":[["GAM_I","GAM_S","GAM_G"],[11.722,13.832,8.517],[62.39,60.29,52.94],[-0.05,0.06,0.06],[384.19,389.07,390.62],[762.594,771.64,774.818],[0,9.05,12.22],[0.9871,0.01072,0.00219]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)    4.872      1.163   4.188 9.33e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##           edf Ref.df      F  p-value    
## s(v)    2.074  2.543  2.084  0.16419    
## s(LME)  1.831  2.000 10.779  0.00112 ** 
## s(year) 1.884  2.000 19.482 3.92e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0597   Deviance explained = 52.9%
## -ML = 390.62  Scale est. = 2.8104    n = 67
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
## (Intercept)   4.7751     0.8792   5.431 1.24e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df      F  p-value    
## s(year)  1.837      2 14.148 1.23e-05 ***
## s(LME)   1.784      2  9.838 0.000248 ***
## s(v,LME) 6.132     14  3.181 0.004441 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0633   Deviance explained = 60.3%
## -ML = 389.07  Scale est. = 1.9951    n = 67
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
## (Intercept)   4.8410     0.8332    5.81 2.94e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df      F  p-value    
## s(year)    1.807  2.000 16.171 5.93e-05 ***
## s(LME)     1.837  2.000 13.489 0.000359 ***
## s(v):LMEBF 1.000  1.000 12.021 0.001008 ** 
## s(v):LMEBB 1.000  1.000  1.851 0.179062    
## s(v):LMESV 3.353  3.759 17.426  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  -0.0465   Deviance explained = 62.4%
## -ML = 384.19  Scale est. = 1.7195    n = 67
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
##          k'      edf  k-index p-value
## s(year)   3 1.836964       NA      NA
## s(LME)    3 1.784141       NA      NA
## s(v,LME) 15 6.131791 0.933916    0.82
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
## s(year)     3 1.806705        NA      NA
## s(LME)      3 1.836556        NA      NA
## s(v):LMEBF  4 1.000173 0.9365695  0.7500
## s(v):LMEBB  4 1.000015 0.9365695  0.7325
## s(v):LMESV  4 3.353482 0.9365695  0.7575
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
## (Intercept)   4.7744     0.9976   4.786 1.27e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df      F  p-value    
## s(year)  1.875      2 15.281 8.27e-06 ***
## s(LME)   1.829      2 10.358 0.000275 ***
## s(v,LME) 5.924     14  3.148 0.008624 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0682   Deviance explained =   60%
## -REML = 388.25  Scale est. = 1.9836    n = 67
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
## (Intercept)   4.7744     0.9976   4.786 1.27e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df      F  p-value    
## s(year)  1.875      2 15.281 8.27e-06 ***
## s(LME)   1.829      2 10.358 0.000275 ***
## s(v,LME) 5.924     14  3.148 0.008624 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0682   Deviance explained =   60%
## -REML = 388.25  Scale est. = 1.9836    n = 67
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
## (Intercept)   4.7996     0.9383   5.115 3.76e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df      F  p-value    
## s(year)    1.8614      2 20.489 1.98e-05 ***
## s(LME)     1.8646      2 13.552  0.00047 ***
## s(v):LMEBF 1.2853      4  3.535  0.00245 ** 
## s(v):LMEBB 0.3666      4  0.178  0.21027    
## s(v):LMESV 3.1171      4 25.955 2.29e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0585   Deviance explained =   62%
## -REML = 386.94  Scale est. = 1.6921    n = 67
```

No terms can be dropped, so I refit model with REML.


```r
GAM_I2 <- gam(NASC_int ~ 
                 s(year, bs = "re") +
                 s(LME, bs = "re") +
                 s(v, by = LME, k =5, bs = "tp"),
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
## (Intercept)   4.8414     0.9107   5.316 1.87e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df      F  p-value    
## s(year)    1.834  2.000 16.405 5.82e-05 ***
## s(LME)     1.873  2.000 14.745 0.000212 ***
## s(v):LMEBF 1.410  1.699  9.225 0.002519 ** 
## s(v):LMEBB 1.001  1.001  1.942 0.168811    
## s(v):LMESV 3.354  3.759 18.501  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  -0.0726   Deviance explained = 63.1%
## -REML = 382.35  Scale est. = 1.6341    n = 67
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

```r
k.check(GAM_I2)
```

```
##            k'      edf   k-index p-value
## s(year)     3 1.833930        NA      NA
## s(LME)      3 1.873210        NA      NA
## s(v):LMEBF  4 1.410153 0.9507057  0.8125
## s(v):LMEBB  4 1.000510 0.9507057  0.7725
## s(v):LMESV  4 3.353545 0.9507057  0.8050
```

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
       x = expression("Current velocity at 380 m (cm s"^-1*")"),
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
       x = expression("Current velocity at 380 m (cm s"^-1*")"),
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

I select data at 380 m depth because it fits well with the DSL diurnal centre of mass. It also prevent from removing too much data like the deeper depth levels would do.


```r
SA_df <- stat_laea %>%
  filter(depth == 380) %>%
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
<div id="htmlwidget-ded66a9eb4c07efa27dc" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ded66a9eb4c07efa27dc">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_S","GAMSA_G"],[10.176,10.961,6.708],[46.29,44.54,31.36],[0.4,0.38,0.27],[236.36,239.61,241.52],[474.55,478.269,484.056],[0,3.72,9.51],[0.85883,0.13376,0.00741]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   15.143      3.162   4.789 1.08e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##           edf Ref.df     F  p-value    
## s(v)    1.000      1 0.610 0.437515    
## s(LME)  1.309      2 3.331 0.022335 *  
## s(year) 1.706      2 8.649 0.000379 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.269   Deviance explained = 31.4%
## -ML = 241.52  Scale est. = 71.123    n = 67
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
## (Intercept)   14.982      2.965   5.053 4.47e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(year)  1.696      2 7.997 0.000459 ***
## s(LME)   1.322      2 3.214 0.016549 *  
## s(v,LME) 3.783     14 0.899 0.021942 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.382   Deviance explained = 44.5%
## -ML = 239.61  Scale est. = 60.163    n = 67
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
## (Intercept)   15.229      3.183   4.785 1.19e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df     F  p-value    
## s(year)    1.710  2.000 9.014 0.000435 ***
## s(LME)     1.513  2.000 6.047 0.003822 ** 
## s(v):LMEBF 1.000  1.000 4.238 0.043943 *  
## s(v):LMEBB 1.000  1.000 2.615 0.111159    
## s(v):LMESV 1.946  2.353 3.002 0.038986 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.397   Deviance explained = 46.3%
## -ML = 236.36  Scale est. = 58.629    n = 67
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
##          k'      edf k-index p-value
## s(year)   3 1.695677      NA      NA
## s(LME)    3 1.322275      NA      NA
## s(v,LME) 15 3.782634 0.97081  0.3625
```


```r
appraise(GAMSA_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAMSA_I)
```

```
##            k'      edf   k-index p-value
## s(year)     3 1.710307        NA      NA
## s(LME)      3 1.512509        NA      NA
## s(v):LMEBF  4 1.000031 0.9508999    0.30
## s(v):LMEBB  4 1.000015 0.9508999    0.31
## s(v):LMESV  4 1.946118 0.9508999    0.30
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
          GAMSA_S_resid %>%
            ggplot(aes(x = xc, y = yc, col = sign(resid), size = abs(resid))) +
            geom_point() +
            metR::scale_color_divergent() + 
            coord_fixed() +
            scale_size(range = c(0.5, 4)),
          nrow = 2)
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
          GAMSA_I_resid %>%
            ggplot(aes(x = xc, y = yc, col = sign(resid), size = abs(resid))) +
            geom_point() +
            metR::scale_color_divergent() + 
            coord_fixed() +
            scale_size(range = c(0.5, 4)),
          nrow = 2)
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
GAMSA_S_resid %>%
  ggplot(aes(x = xc, y = yc, col = resid)) +
  geom_point() +
  metR::scale_color_divergent() + 
  coord_fixed() +
  facet_wrap(~ year)
```

![](PanArctic_DSL_statistics_files/figure-html/unnamed-chunk-1-1.png)<!-- -->



```r
par(mfrow = c(1, 2))
acf(resid(GAMSA_I), lag.max = 20, main = "ACF")
pacf(resid(GAMSA_I), lag.max = 20, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-acf-pacf-1.png)<!-- -->


```r
GAMSA_I_resid %>%
  ggplot(aes(x = xc, y = yc, col = resid)) +
  geom_point() +
  metR::scale_color_divergent() + 
  coord_fixed() +
  facet_wrap(~ year)
```

![](PanArctic_DSL_statistics_files/figure-html/unnamed-chunk-2-1.png)<!-- -->
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
## (Intercept)   14.969      3.402     4.4  4.6e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(year)  1.775      2 8.455 0.000483 ***
## s(LME)   1.407      2 3.492 0.018526 *  
## s(v,LME) 3.754     14 0.907 0.026471 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.383   Deviance explained = 44.7%
## -REML = 237.54  Scale est. = 60.077    n = 67
```

No terms can be dropped.

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
## (Intercept)    15.15       3.55   4.268 7.18e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F  p-value    
## s(year)    1.7835      2 9.965 0.000327 ***
## s(LME)     1.5349      2 5.914 0.006368 ** 
## s(v):LMEBF 0.7578      4 1.021 0.046595 *  
## s(v):LMEBB 0.6110      4 0.450 0.114438    
## s(v):LMESV 1.7746      4 2.076 0.042884 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.393   Deviance explained = 45.3%
## -REML = 237.24  Scale est. = 59.049    n = 67
```

No terms can be dropped.

## Final models

### Model S

Refit model with REML


```r
GAMSA_S3 <- gam(SA_int ~ 
                  s(year, bs = "re") +
                  s(LME, bs = "re") +
                  s(v, LME, bs = "fs", k = 5),
                data = SA_df, family = "gaussian", method = "REML")
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
summary(GAMSA_S3)
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
## (Intercept)   14.969      3.402     4.4  4.6e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(year)  1.775      2 8.455 0.000483 ***
## s(LME)   1.407      2 3.492 0.018526 *  
## s(v,LME) 3.754     14 0.907 0.026471 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.383   Deviance explained = 44.7%
## -REML = 237.54  Scale est. = 60.077    n = 67
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

Refit model with REML


```r
GAMSA_I3 <- gam(SA_int ~
                  s(year, bs = "re") +
                  s(LME, bs = "re") +
                  s(v, by = LME, k = 6, bs = "tp"),
                data = SA_df, family = "gaussian", method = "REML")
summary(GAMSA_I3)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(LME, bs = "re") + s(v, by = LME, 
##     k = 6, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   15.374      3.785   4.061 0.000153 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df      F  p-value    
## s(year)    1.847  2.000 15.466 2.44e-05 ***
## s(LME)     1.622  2.000  8.967 0.002569 ** 
## s(v):LMEBF 1.000  1.000  4.700 0.034416 *  
## s(v):LMEBB 1.000  1.000  3.183 0.079805 .  
## s(v):LMESV 4.462  4.867  5.720 0.000957 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.522   Deviance explained = 59.4%
## -REML = 227.05  Scale est. = 46.512    n = 67
```

Check residuals.


```r
appraise(GAMSA_I3)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI3-residuals-covariates-1.png)<!-- -->

```r
k.check(GAMSA_I3)
```

```
##            k'      edf  k-index p-value
## s(year)     3 1.846634       NA      NA
## s(LME)      3 1.622055       NA      NA
## s(v):LMEBF  5 1.000176 1.071799   0.710
## s(v):LMEBB  5 1.000166 1.071799   0.685
## s(v):LMESV  5 4.462338 1.071799   0.650
```

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

```r
GAMSA_I3_resid %>%
  ggplot(aes(x = xc, y = yc, col = resid)) +
  geom_point() +
  metR::scale_color_divergent() + 
  coord_fixed() +
  facet_wrap(~ year)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI3-residuals-covariates-3.png)<!-- -->

```r
par(mfrow = c(1, 2))
acf(resid(GAMSA_I3), lag.max = 20, main = "ACF")
pacf(resid(GAMSA_I3), lag.max = 20, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI3-residuals-covariates-4.png)<!-- -->


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
res = 100
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
       x = expression("Current velocity at 380 m (cm s"^-1*")"),
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
       x = expression("Current velocity at 380 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
# Combine plots
plot_grid(GAMSAS_velo, GAMSAI_velo, ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS-pretty-ggplot-1.png)<!-- -->


```r
library(plotly)
```

```
## 
## Attaching package: 'plotly'
```

```
## The following object is masked from 'package:plotfunctions':
## 
##     add_bars
```

```
## The following object is masked from 'package:raster':
## 
##     select
```

```
## The following object is masked from 'package:ggplot2':
## 
##     last_plot
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```
## The following object is masked from 'package:graphics':
## 
##     layout
```

```r
predSA_I %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = fitted, col = LME), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = lwr_ci, ymax = upr_ci, fill = LME), 
              alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 30), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(breaks = seq(0, 50, 10)) +
  theme(panel.grid = element_line(color = "grey80"),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
```

![](PanArctic_DSL_statistics_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

