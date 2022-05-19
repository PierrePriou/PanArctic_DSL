---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/19 at 19:15"
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
## # A tibble: 5 × 2
##   IHO_area     n
##   <fct>    <int>
## 1 WAO_BF     106
## 2 CAA        101
## 3 BB         258
## 4 DS         167
## 5 EAO        157
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
<div id="htmlwidget-92f8c3bb0e8f5579d556" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-92f8c3bb0e8f5579d556">{"x":{"filter":"none","vertical":false,"data":[["GAM_GI","GAM_GS","GAM_G"],[32.757,32.25,20.63],[75.48,73.56,63.39],[0.14,0,0.02],[773.83,791.47,804.43],[1560.654,1570.577,1597.349],[0,9.92,36.69],[0.99305,0.00695,0]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)    4.790      0.692   6.922 2.76e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df      F p-value  
## s(v)        2.0798  2.541  1.929  0.1484  
## s(chl)      2.3076  2.803  1.973  0.1382  
## s(IHO)      2.4457  4.000 16.737  0.2802  
## s(year)     0.7288  2.000  5.750  0.4689  
## s(IHO,year) 8.3573 14.000  7.612  0.0911 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0235   Deviance explained = 63.4%
## -REML = 804.43  Scale est. = 2.5392    n = 131
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
## (Intercept)   4.7370     0.4091   11.58   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                   edf Ref.df     F  p-value    
## s(chl)      1.0000999      1 2.891  0.09200 .  
## s(v)        1.0001924      1 0.732  0.39441    
## s(year)     0.5361772      2 1.113  0.45265    
## s(IHO)      0.0006894      4 0.000  0.43770    
## s(year,IHO) 8.4731208     14 2.554  0.00451 ** 
## s(chl,IHO)  7.4905385     23 4.751 2.01e-05 ***
## s(v,IHO)    6.0464317     23 1.197  0.00190 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.00331   Deviance explained = 73.6%
## -REML = 791.47  Scale est. = 1.5614    n = 131
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
## (Intercept)   4.7643     0.2905    16.4   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf    Ref.df     F  p-value    
## s(chl)           1.000e+00 1.000e+00 0.013 0.909883    
## s(v)             1.878e+00 2.282e+00 1.182 0.306491    
## s(year)          5.668e-02 2.000e+00 0.036 0.381708    
## s(IHO)           1.614e-04 4.000e+00 0.000 0.767014    
## s(year,IHO)      7.312e+00 1.400e+01 1.551 0.001286 ** 
## s(chl):IHOWAO_BF 1.152e+00 1.240e+00 4.247 0.041117 *  
## s(chl):IHOCAA    1.000e+00 1.000e+00 0.028 0.868563    
## s(chl):IHOBB     2.833e+00 3.264e+00 1.916 0.114304    
## s(chl):IHODS     6.231e-05 1.093e-04 0.000 0.999858    
## s(chl):IHOEAO    2.844e+00 3.250e+00 6.273 0.000741 ***
## s(v):IHOWAO_BF   1.000e+00 1.000e+00 0.758 0.385990    
## s(v):IHOCAA      1.000e+00 1.000e+00 0.437 0.510196    
## s(v):IHOBB       1.859e+00 2.225e+00 1.180 0.313733    
## s(v):IHODS       4.399e-01 7.195e-01 0.016 0.914963    
## s(v):IHOEAO      3.571e+00 3.861e+00 4.382 0.002258 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Rank: 70/72
## R-sq.(adj) =  0.143   Deviance explained = 75.5%
## -REML = 773.83  Scale est. = 1.3028    n = 131
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
##             k'          edf   k-index p-value
## s(chl)       4 1.0000998632 0.9269839  0.6750
## s(v)         4 1.0001924303 0.8362879  0.2225
## s(year)      3 0.5361772200        NA      NA
## s(IHO)       5 0.0006893627        NA      NA
## s(year,IHO) 15 8.4731208365        NA      NA
## s(chl,IHO)  25 7.4905384548 0.9269839  0.6825
## s(v,IHO)    25 6.0464317442 0.8362879  0.2250
```


```r
appraise(GAM_GI, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_GI)
```

```
##                  k'          edf   k-index p-value
## s(chl)            4 1.000117e+00 0.9276536  0.6800
## s(v)              4 1.877884e+00 0.8373218  0.2200
## s(year)           3 5.667702e-02        NA      NA
## s(IHO)            5 1.613875e-04        NA      NA
## s(year,IHO)      15 7.311715e+00        NA      NA
## s(chl):IHOWAO_BF  4 1.151523e+00 0.9276536  0.6600
## s(chl):IHOCAA     4 1.000227e+00 0.9276536  0.6450
## s(chl):IHOBB      4 2.833476e+00 0.9276536  0.6025
## s(chl):IHODS      4 6.231082e-05 0.9276536  0.6675
## s(chl):IHOEAO     4 2.843581e+00 0.9276536  0.6400
## s(v):IHOWAO_BF    4 1.000208e+00 0.8373218  0.2050
## s(v):IHOCAA       4 1.000101e+00 0.8373218  0.1800
## s(v):IHOBB        4 1.858875e+00 0.8373218  0.2000
## s(v):IHODS        4 4.398652e-01 0.8373218  0.2100
## s(v):IHOEAO       4 3.570807e+00 0.8373218  0.2350
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
## (Intercept)   4.6486     0.3996   11.63   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                  edf Ref.df     F  p-value    
## s(chl)      0.722110      4 0.570 1.83e-06 ***
## s(v)        0.001407      4 0.000 0.137060    
## s(year)     0.641707      2 1.678 0.424289    
## s(IHO)      0.002308      4 0.001 0.373199    
## s(year,IHO) 8.369545     14 2.559 0.007021 ** 
## s(chl,IHO)  8.021388     24 4.447 3.27e-05 ***
## s(v,IHO)    6.022511     24 1.443 0.000244 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0246   Deviance explained = 73.3%
## -REML = 793.12  Scale est. = 1.5604    n = 131
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
## (Intercept)   4.7633     0.4116   11.57   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                  edf Ref.df     F  p-value    
## s(chl)      1.000805  1.001 3.079 0.082089 .  
## s(year)     0.641933  2.000 1.702 0.422750    
## s(IHO)      0.001804  4.000 0.000 0.398950    
## s(year,IHO) 8.406296 14.000 2.612 0.006718 ** 
## s(chl,IHO)  7.736855 23.000 4.757 3.11e-05 ***
## s(v,IHO)    5.958314 24.000 1.470 0.000207 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0271   Deviance explained = 73.3%
## -REML = 791.65  Scale est. = 1.5458    n = 131
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
<div id="htmlwidget-0543d1e599f8d9ba10ef" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-0543d1e599f8d9ba10ef">{"x":{"filter":"none","vertical":false,"data":[["GAM_GS2","GAM_GS"],[31.313,32.25],[73.3,73.56],[0.03,0],[791.65,791.47],[1570.078,1570.577],[0,0.5],[0.56205343,0.43794657]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   4.7565     0.2306   20.62   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df      F  p-value    
## s(chl)           0.0002656      4  0.000 0.642995    
## s(v)             1.7308203      4  2.089 0.022039 *  
## s(year)          0.0002642      2  0.000 0.507723    
## s(IHO)           0.0006196      4  0.000 0.514574    
## s(year,IHO)      8.4822702     14  2.061 0.000145 ***
## s(chl):IHOWAO_BF 0.9449256      4 17.885 5.93e-05 ***
## s(chl):IHOCAA    0.0005015      4  0.000 0.769832    
## s(chl):IHOBB     0.5743704      4  0.400 0.194823    
## s(chl):IHODS     0.0002877      4  0.000 0.722881    
## s(chl):IHOEAO    2.5568529      4 20.606  < 2e-16 ***
## s(v):IHOWAO_BF   0.0003997      4  0.000 0.736161    
## s(v):IHOCAA      0.0001820      4  0.000 0.881129    
## s(v):IHOBB       0.5845854      4  0.414 0.109454    
## s(v):IHODS       1.1519118      4  1.664 0.034366 *  
## s(v):IHOEAO      2.6377838      4  5.870 5.12e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.161   Deviance explained = 73.3%
## -REML = 786.61  Scale est. = 1.3989    n = 131
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
## (Intercept)   4.9413     0.2952   16.74   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf   Ref.df      F  p-value    
## s(v)             1.5447126  1.86716  2.088 0.082281 .  
## s(year)          0.0429052  2.00000  0.025 0.385091    
## s(IHO)           0.0001016  4.00000  0.000 0.805284    
## s(year,IHO)      7.1907430 14.00000  1.813 0.000286 ***
## s(chl):IHOWAO_BF 1.3789902  1.55482  9.260 0.000530 ***
## s(chl):IHOCAA    1.0001017  1.00018  0.312 0.577765    
## s(chl):IHOBB     1.0002494  1.00049  1.599 0.208871    
## s(chl):IHODS     1.0000841  1.00015  0.282 0.596409    
## s(chl):IHOEAO    2.8286254  3.23478 14.913  < 2e-16 ***
## s(v):IHOWAO_BF   1.0029975  1.00511  4.017 0.047488 *  
## s(v):IHOCAA      1.0001697  1.00032  1.332 0.251115    
## s(v):IHOBB       2.2574940  2.67927  3.758 0.016040 *  
## s(v):IHODS       0.0026291  0.00507  0.011 0.994191    
## s(v):IHOEAO      3.5646486  3.86067  6.565 9.34e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Rank: 67/68
## R-sq.(adj) =  0.164   Deviance explained =   74%
## -REML = 775.36  Scale est. = 1.4058    n = 131
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
<div id="htmlwidget-23b98ed4725e08c64b4d" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-23b98ed4725e08c64b4d">{"x":{"filter":"none","vertical":false,"data":[["GAM_GI","GAM_GI2"],[32.757,29.961],[75.48,74.01],[0.14,0.16],[773.83,775.36],[1560.654,1563.415],[0,2.76],[0.7990712904,0.2009287096]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
SA_df %>%
  dplyr::select(IHO, year, v, chl) %>%
  # group_by(IHO, year) %>%
  group_by(IHO) %>%
  get_summary_stats(type = "common")
```

```
## # A tibble: 10 × 11
##    IHO    variable     n   min   max median   iqr  mean    sd    se    ci
##    <fct>  <chr>    <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
##  1 WAO_BF chl         16 0.111 0.857  0.245 0.331 0.371 0.245 0.061 0.131
##  2 WAO_BF v           16 0.575 5.40   2.42  2.05  2.41  1.41  0.352 0.749
##  3 CAA    chl         21 0.234 0.806  0.405 0.156 0.448 0.15  0.033 0.068
##  4 CAA    v           21 0.473 3.39   0.874 0.654 1.12  0.712 0.155 0.324
##  5 BB     chl         46 0.283 2.54   1.01  0.829 1.08  0.525 0.077 0.156
##  6 BB     v           46 0.119 4.07   1.81  1.38  1.74  0.99  0.146 0.294
##  7 DS     chl         25 0.319 0.685  0.413 0.111 0.435 0.094 0.019 0.039
##  8 DS     v           25 0.37  2.74   1.50  1.39  1.57  0.782 0.156 0.323
##  9 EAO    chl         23 0.089 4.06   1.36  1.06  1.22  0.898 0.187 0.388
## 10 EAO    v           23 0.45  5.20   2.07  2.15  2.42  1.50  0.314 0.65
```


```r
res = 20 # resolution for predictions
# WAO and BF
SA_WAO_BF <- data.frame(
  IHO = factor(rep(rep(rep("WAO_BF", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.575, 5.395, length.out = res), 3), # Velocity
        rep(2.415, res * 3)), # Dummy
  chl = c(rep(0.245, res * 3), # Dummy
          rep(seq(0.111, 0.857, length.out = res), 3))) # Chlorophyll
# CAA
SA_CAA <- data.frame(
  IHO = factor(rep(rep(rep("CAA", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.473, 3.387, length.out = res), 3), # Velocity
        rep(0.874, res * 3)), # Dummy
  chl = c(rep(0.405, res * 3), # Dummy
          rep(seq(0.234, 0.806, length.out = res), 3))) # Chlorophyll
# BB
SA_BB <- data.frame(
  IHO = factor(rep(rep(rep("BB", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.119, 4.069, length.out = res), 3), # Velocity
        rep(1.814, res * 3)), # Dummy
  chl = c(rep(1.005, res * 3), # Dummy
          rep(seq(0.283, 2.539, length.out = res), 3))) # Chlorophyll
# DS
SA_DS <- data.frame(
  IHO = factor(rep(rep(rep("DS", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.370, 2.743, length.out = res), 3), # Velocity
        rep(1.495, res * 3)), # Dummy
  chl = c(rep(0.413, res * 3), # Dummy
          rep(seq(0.319, 0.685, length.out = res), 3))) # Chlorophyll
# EAO
SA_EAO <- data.frame(
  IHO = factor(rep(rep(rep("EAO", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.450, 5.204, length.out = res), 3), # Velocity
        rep(2.070, res * 3)), # Dummy
  chl = c(rep(1.355, res * 3), # Dummy
          rep(seq(0.089, 4.063, length.out = res), 3))) # Chlorophyll
# Combine datasets
SA_new <- bind_rows(SA_WAO_BF, SA_CAA, SA_BB, SA_DS, SA_EAO)
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
