---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/09 at 15:30"
output: 
  html_document:
    keep_md: yes
  github_document:
    always_allow_html: true
---

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
library(tidymv)       # Predict GAM
library(MuMIn)        # AIC weights
library(DT)           # Interactive table
library(rstatix)      # Pipe-friendly stats
# Custom figure theme
theme_set(theme_bw())
theme_update(axis.text = element_text(size = 9),
             axis.title = element_text(size = 9),
             strip.text.x = element_text(size = 9, face = "plain", hjust = 0.5),
             strip.background = element_rect(colour = "transparent", fill = "transparent"),
             legend.title = element_text(size = 9),
             legend.margin = margin(0, 0, 0, 0),
             legend.box.margin = margin(0, 0, -8, 0),
             panel.grid = element_blank(), 
             plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
options(dplyr.summarise.inform = F) # Suppress summarise() warning
```

I want to test whether temperature and salinity at mesopelagic depth, sea-ice concentration, open-water duration (a proxy for productivity) have an effect on the backscatter anomalies observed per year. I therefore combined gridded acoustic data—integrated mesopelagic NASC—with gridded CTD, and remote sensing data projected on either the WGS84 or the EASE-Grid 2.0 North.


```r
# Map projections
cell_res <- 50 # Cell resolution in km
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931") # Seaice projection
projection(arctic_laea) <- gsub("units=m", "units=km", projection(arctic_laea)) # Convert proj unit from m to km
res(arctic_laea) <- c(cell_res, cell_res) # Define the 100 km cell resolution

arctic_latlon <- raster(extent(-155, 35, 66, 85), # Base projection for acoustic and CTD data
                        crs = "EPSG:4326", 
                        res = c(2, 1)) # cells of 2 degree longitude per 1 degree latitude

# Coastline shapefiles
coast_10m_laea <- readOGR("data/bathy/ne_10m_land.shp", verbose = F) %>% # Coastline in laea
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  crop(extent(-180, 180, 0, 90)) %>% # Crop shapefile
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to a dataframe for ggplot
  rename(xc = long, yc = lat)

# Gridded acoustic, CTD, and sea ice data
load("data/acoustics/SA_grids.RData") # Acoustic data
load("data/remote_sensing/physics_grids.RData") # Modelled physics data 
load("data/remote_sensing/seaice_grids.RData") # Remote sensing sea ice data
```

# Data preparation

I combine data using the EASE-Grid 2.0 North (EPSG:6931) and normalize the covariates per IHO area. I combined the Beaufort Sea and West Arctic Ocean data to have more points for the regression due to their geographical proximity.


```r
SA_laea <- SA_grid_laea %>%  # Tidy anomaly dataset for joining
  dplyr::select(-lat, -lon) %>% 
   filter(empty_clean == F) # Select cells with DSL
phy_laea <- phy_grid_laea %>% # Tidy remote sensing dataset for joining
  dplyr::select(year, area, xc, yc, cell_res, depth, thetao, velocity)
seaice_laea <- seaice_grid_laea %>%
  dplyr::select(year, area, xc, yc, cell_res, openwater_duration, ice_week)

stat_laea <- left_join(SA_laea, phy_laea, by = c("year", "area", "xc", "yc", "cell_res")) %>% # Join acoustic and physics
  left_join(., seaice_laea, by = c("year", "area", "xc", "yc", "cell_res")) %>%
  # For cells that have NaN, get the mean of the surrounding cells
  rowwise() %>%
  mutate(xc_na = if_else(is.na(openwater_duration) == T, xc, NaN),
         yc_na = if_else(is.na(openwater_duration) == T, yc, NaN),
         year_na = if_else(is.na(openwater_duration) == T, year, NaN), 
         openwater_duration = if_else(is.na(openwater_duration) == T, mean(pull(subset(seaice_grid_laea,
                                                                                       xc >= xc_na - 50 &
                                                                                         xc <= xc_na + 50 & 
                                                                                         yc >= yc_na - 50 & 
                                                                                         yc <= yc_na + 50 & 
                                                                                         year == year_na,
                                                                                       select = openwater_duration),
                                                                                openwater_duration),
                                                                           na.rm = T),
                                      openwater_duration),
         ice_week = if_else(is.na(ice_week) == T, mean(pull(subset(seaice_grid_laea,
                                                                   xc >= xc_na - 50 &
                                                                     xc <= xc_na + 50 & 
                                                                     yc >= yc_na - 50 & 
                                                                     yc <= yc_na + 50 & 
                                                                     year == year_na,
                                                                   select = ice_week),
                                                            ice_week),
                                                       na.rm = T),
                            ice_week)) %>%
  filter(depth == 380) %>% # Select data at 380 m depth and remove cells with no DSL
  mutate(velocity = velocity * 100, # Convert to cm/s
         IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                     IHO_area == "West Arctic Ocean" ~ "WAO_BF_CAA",
                                     IHO_area == "Beaufort Sea" ~ "WAO_BF_CAA",
                                     IHO_area == "The Northwestern Passages" ~ "WAO_BF_CAA",
                                     IHO_area == "Baffin Bay" ~ "BB",
                                     IHO_area == "Davis Strait" ~ "DS"),
                            levels = c("WAO_BF_CAA", "BB", "DS", "EAO"))) %>%
  ungroup() %>%
  dplyr::select(-xc_na, -yc_na, -year_na)
```

# Data exploration 

Maps of all variables.


```r
stat_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_point(aes(col = IHO_area)) +
  ggtitle("regions") +
  facet_wrap(~ year) +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F)
```

![](PanArctic_DSL_statistics_files/figure-html/map-SA-int-1.png)<!-- -->


```r
stat_laea %>% # Openwater duration
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_tile(aes(fill = openwater_duration)) +
  scale_fill_viridis_c("Day", option = "plasma", na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Openwater duration") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-ow-1.png)<!-- -->


```r
stat_laea %>% # Openwater duration
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_tile(aes(fill = ice_week)) +
  scale_fill_viridis_c("Week", option = "viridis", na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Week ice breakup") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-ice-week-1.png)<!-- -->


```r
stat_laea %>% # Temperature
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_tile(aes(fill = thetao)) +
  scale_fill_cmocean("Temp (dC)", name = "thermal", na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Temperature at 380 m depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-thetao-1.png)<!-- -->


```r
stat_laea %>% # Ice concentration
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_tile(aes(fill = velocity)) +
  scale_fill_cmocean("velo (cm/s)", name = "speed", na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Current velocity at 380 m depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-velocity-1.png)<!-- -->

# Spatial and interannual variability


```r
SA_diff <- stat_laea %>%
  dplyr::select(year, IHO_area, SA_int_clean, NASC_int, CM) %>%
  mutate(year = factor(year))
```

## All areas

I check whether there are inter-annual and spatial variability in S~A~ and the centre of mass across years and areas.


```r
# Visualize data
plot_grid(SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = SA_int_clean)),
          SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = CM)))
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-1.png)<!-- -->

```r
# Check assumptions
## Outliers
SA_diff %>%
  group_by(year) %>%
  identify_outliers(SA_int_clean) 
```

```
## # A tibble: 2 x 7
##   year  IHO_area SA_int_clean NASC_int    CM is.outlier is.extreme
##   <fct> <fct>           <dbl>    <dbl> <dbl> <lgl>      <lgl>     
## 1 2017  EAO              41.6   14314.  412. TRUE       FALSE     
## 2 2017  EAO              48.3   67316.  392. TRUE       FALSE
```

```r
SA_diff %>%
  group_by(year) %>%
  identify_outliers(CM) 
```

```
## # A tibble: 2 x 7
##   year  IHO_area   SA_int_clean NASC_int    CM is.outlier is.extreme
##   <fct> <fct>             <dbl>    <dbl> <dbl> <lgl>      <lgl>     
## 1 2015  WAO_BF_CAA         7.69     7.72  591. TRUE       TRUE      
## 2 2015  WAO_BF_CAA         3.81     5.78  572. TRUE       TRUE
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(SA_int_clean) 
```

```
## # A tibble: 5 x 7
##   IHO_area year  SA_int_clean NASC_int    CM is.outlier is.extreme
##   <fct>    <fct>        <dbl>    <dbl> <dbl> <lgl>      <lgl>     
## 1 BB       2017          3.10     5.08  285. TRUE       FALSE     
## 2 BB       2017          2.83     3.83  284. TRUE       FALSE     
## 3 DS       2016         30.2   1036.    364. TRUE       FALSE     
## 4 EAO      2017         41.6  14314.    412. TRUE       FALSE     
## 5 EAO      2017         48.3  67316.    392. TRUE       FALSE
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(CM) 
```

```
## # A tibble: 7 x 7
##   IHO_area   year  SA_int_clean NASC_int    CM is.outlier is.extreme
##   <fct>      <fct>        <dbl>    <dbl> <dbl> <lgl>      <lgl>     
## 1 WAO_BF_CAA 2015          7.69     7.72  591. TRUE       TRUE      
## 2 WAO_BF_CAA 2015          3.81     5.78  572. TRUE       TRUE      
## 3 BB         2015         14.3     28.8   415. TRUE       FALSE     
## 4 BB         2016         13.5     24.7   393. TRUE       FALSE     
## 5 BB         2016         16.6     48.0   427. TRUE       FALSE     
## 6 BB         2016         16.6     47.4   401. TRUE       FALSE     
## 7 BB         2017         13.0     20.8   265. TRUE       FALSE
```

```r
## Normality: Checked by looking at the anova residuals with QQ plot
ggqqplot(SA_diff, "SA_int_clean", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-2.png)<!-- -->

```r
ggqqplot(SA_diff, "CM", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-3.png)<!-- -->

```r
ggqqplot(SA_diff, "SA_int_clean", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-4.png)<!-- -->

```r
ggqqplot(SA_diff, "CM", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-5.png)<!-- -->

```r
## Homogeneity of variance
plot(lm(SA_int_clean ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-6.png)<!-- -->

```r
plot(lm(CM ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-7.png)<!-- -->

```r
plot(lm(SA_int_clean ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-8.png)<!-- -->

```r
plot(lm(CM ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-9.png)<!-- -->

```r
SA_diff %>%
  levene_test(SA_int_clean ~ year)
```

```
## # A tibble: 1 x 4
##     df1   df2 statistic     p
##   <int> <int>     <dbl> <dbl>
## 1     2   117     0.465 0.629
```

I used a non parametric Kruskal-Wallis test to test for interannual and spatial variability because the variance is not homogenous and not very normal.


```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(SA_int_clean ~ year)
```

```
## # A tibble: 1 x 6
##   .y.              n statistic    df      p method        
## * <chr>        <int>     <dbl> <int>  <dbl> <chr>         
## 1 SA_int_clean   120      5.17     2 0.0754 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ year)
```

```
## # A tibble: 1 x 6
##   .y.       n statistic    df     p method        
## * <chr> <int>     <dbl> <int> <dbl> <chr>         
## 1 CM      120      2.94     2  0.23 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(SA_int_clean ~ IHO_area)
```

```
## # A tibble: 1 x 6
##   .y.              n statistic    df      p method        
## * <chr>        <int>     <dbl> <int>  <dbl> <chr>         
## 1 SA_int_clean   120      8.37     3 0.0389 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ IHO_area)
```

```
## # A tibble: 1 x 6
##   .y.       n statistic    df            p method        
## * <chr> <int>     <dbl> <int>        <dbl> <chr>         
## 1 CM      120      38.1     3 0.0000000266 Kruskal-Wallis
```

The center of mass did not vary significantly between years (H = 1.652899, p = 0.438) but varied significantly among areas (*H* = 57.8495, *p* < 0.001). Mesopelagic backscatter S~A~ varied significantly between areas (*H* = 26.19394, *p* < 0.001) but not years (*H* = 8.6464, *p* = 0.013). Thus, in the S~A~ HGAM a random effect for change in S~A~ per area is likely needed.

## Within areas

I check whether there are inter-annual variability in S~A~ across years within areas. This will be used to decide whether a random effect is needed in the GAM. I used a non parametric Kruskal-Wallis test to test for interannual variability within areas.


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(IHO_area) %>%
  kruskal_test(SA_int_clean ~ year)
```

```
## # A tibble: 4 x 7
##   IHO_area   .y.              n statistic    df        p method        
## * <fct>      <chr>        <int>     <dbl> <int>    <dbl> <chr>         
## 1 WAO_BF_CAA SA_int_clean    29   17.6        2 0.000149 Kruskal-Wallis
## 2 BB         SA_int_clean    42    4.83       2 0.0893   Kruskal-Wallis
## 3 DS         SA_int_clean    24    0.0857     2 0.958    Kruskal-Wallis
## 4 EAO        SA_int_clean    25    4.87       2 0.0878   Kruskal-Wallis
```

There was no interannual differences in SA within each region (WAO_BF - H = 4.8395483, *p* = 0.0889; CAA - H = 6.9278752, *p* = 0.0313; BB - H = 0.0968, *p* = 0.0889; DS - H = 0.9470, *p* = 0.0889; EAO - H = 5.0263899, *p* = 0.0810). Therefore, there seems to be no need for a random effect for interannual variability within years.

# HGAM - Gaussian SA

## Data preparation


```r
SA_df <- stat_laea %>%
  dplyr::select(year, xc, yc, IHO_area, SA_int_clean, velocity, thetao, openwater_duration, ice_week) %>%
  group_by(IHO_area) %>%
  mutate(year = factor(year),
         SA_int_n = (SA_int_clean - min(SA_int_clean)) / (max(SA_int_clean) - min(SA_int_clean)),
         v_n = (velocity - min(velocity)) / (max(velocity) - min(velocity)),
         t_n = (thetao - min(thetao)) / (max(thetao) - min(thetao)),
         o_n = (openwater_duration - min(openwater_duration)) / (max(openwater_duration) - min(openwater_duration)),
         s_n = (ice_week - min(ice_week)) / (max(ice_week) - min(ice_week))) %>%
  ungroup() %>%
  rename(IHO = IHO_area, 
         v = velocity,
         t = thetao,
         o = openwater_duration,
         s = ice_week)
```

Plot data.


```r
plot_grid(SA_df %>%
            ggplot(aes(x = v, y = SA_int_clean, col = IHO)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = t, y = SA_int_clean, col = IHO)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = o, y = SA_int_clean, col = IHO)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = s, y = SA_int_clean, col = IHO)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/plot-data-area-1.png)<!-- -->


```r
plot_grid(SA_df %>%
            ggplot(aes(x = v, fill = IHO)) +
            geom_histogram() +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = t, fill = IHO)) +
            geom_histogram() +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = o, fill = IHO)) +
            geom_histogram() +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = s, fill = IHO)) +
            geom_histogram() +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/plot-distribution-1.png)<!-- -->

Check correlations.


```r
corr_EAO <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "EAO") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_EAO, type = "lower", lab = T, title = "corr EAO")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-1.png)<!-- -->

```r
corr_BB <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "BB") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_BB, type = "lower", lab = T, title = "corr BB")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-2.png)<!-- -->

```r
corr_DS <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "DS") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_DS, type = "lower", lab = T, title = "corr BB")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-3.png)<!-- -->

```r
corr_WAO_BF_CAA <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "WAO_BF_CAA") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_WAO_BF_CAA, type = "lower", lab = T, title = "corr WAO")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-4.png)<!-- -->

## Model fitting

I decided to fit hierarchical generalized additive model due to their flexibility in modelling non linear relationships. I fit several models with different structures (random intercept, random "slope") following Pedersen et al. 2019.


```r
# Model G
GAM1 <- gam(SA_int_n ~ s(v, k = 5, bs = "tp") + s(t, k = 5, bs = "tp") + s(o, k = 5, bs = "tp"),
            data = SA_df, family = "gaussian", method = "REML")
# Model S
GAM2 <- gam(SA_int_n ~ s(v, IHO, bs = "fs", k = 5) + s(t, IHO, bs = "fs", k = 5) + s(o, IHO, bs = "fs", k = 5), 
            data = SA_df, family = "gaussian", method = "REML")
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
# Model I
GAM3 <- gam(SA_int_n ~ IHO + s(v, by = IHO, k = 5, bs = "tp") + s(t, by = IHO, k = 5, bs = "tp") + 
              s(o, by = IHO, k = 5, bs = "tp"),
            data = SA_df, family = "gaussian", method = "REML")
# Model GS
GAM4 <- gam(SA_int_n ~ IHO + s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + s(v, IHO, bs = "fs", k = 5) +
              s(t, IHO, bs = "fs", k = 5) + s(o, IHO, bs = "fs", k = 5),
            data = SA_df, family = "gaussian", method = "REML")
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
# Model GI
GAM5 <- gam(SA_int_n ~ IHO + s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + s(v, by = IHO, k = 5, bs = "tp") + 
              s(t, by = IHO, k = 5, bs = "tp") + s(o, by = IHO, k = 5, bs = "tp"),
            data = SA_df, family = "gaussian", method = "REML")

# Summary metrics
GAM_AIC <- AIC(GAM1, GAM2, GAM3, GAM4, GAM5) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
summ_GAM <- data.frame(model = c("GAM1", "GAM2", "GAM3", "GAM4", "GAM5"),
                       reml = round(c(GAM1$gcv.ubre, GAM2$gcv.ubre, GAM3$gcv.ubre, GAM4$gcv.ubre, GAM5$gcv.ubre), 2), 
                       dev_expl = round(c((1 - (GAM1$deviance / GAM1$null.deviance)) * 100,
                                          (1 - (GAM2$deviance / GAM2$null.deviance)) * 100,
                                          (1 - (GAM3$deviance / GAM3$null.deviance)) * 100,
                                          (1 - (GAM4$deviance / GAM4$null.deviance)) * 100,
                                          (1 - (GAM5$deviance / GAM5$null.deviance)) * 100), 2),
                       r2 = round(c(summary(GAM1)$r.sq, summary(GAM2)$r.sq, summary(GAM3)$r.sq, summary(GAM4)$r.sq, 
                                    summary(GAM5)$r.sq), 2)) %>%
  full_join(., GAM_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 10)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
summ_GAM
```

```{=html}
<div id="htmlwidget-15bc36f43093385136d4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-15bc36f43093385136d4">{"x":{"filter":"none","vertical":false,"data":[["GAM5","GAM3","GAM4","GAM1","GAM2"],[27.088,25.201,24.297,10.836,25.413],[62.97,59.8,54.74,41.98,54.06],[0.54,0.51,0.47,0.38,0.46],[0.8,-0.83,-0.55,-2.41,-4.49],[-44.6,-38.5,-26.084,-23.211,-22.058],[0,6.1,18.52,21.39,22.54],[0.9546633465,0.0452118292,9.10232e-05,2.16415e-05,1.21595e-05]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

I select `GAM4 model GS` because I expect to have a similar functional response within each region but that response can vary between region.


```r
summary(GAM4)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int_n ~ IHO + s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + s(v, 
##     IHO, bs = "fs", k = 5) + s(t, IHO, bs = "fs", k = 5) + s(o, 
##     IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)   0.3480     0.1069   3.256  0.00154 **
## IHOBB         0.1595     0.1114   1.431  0.15554   
## IHODS         0.1829     0.1901   0.962  0.33828   
## IHOEAO        0.1326     0.1393   0.952  0.34337   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##             edf Ref.df     F  p-value    
## s(v)     1.2233  1.300 3.136 0.076955 .  
## s(t)     3.0861  3.543 3.026 0.037283 *  
## s(o)     2.3794  2.726 8.530 0.000357 ***
## s(v,IHO) 4.6301 15.000 0.620 0.023297 *  
## s(t,IHO) 0.8142 14.000 0.107 0.118222    
## s(o,IHO) 2.8830 15.000 0.366 0.027976 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.467   Deviance explained = 54.7%
## -REML = -0.54669  Scale est. = 0.037341  n = 120
```

## Model checking

First I check the basis size k. `k-indexes` are > 1 or close to 1 so the basis size is large enough. The residual plot look good too.


```r
par(mfrow = c(2, 2))
gam.check(GAM4, rep = 500)
```

![](PanArctic_DSL_statistics_files/figure-html/basis-size-residuals-1.png)<!-- -->

```
## 
## Method: REML   Optimizer: outer newton
## full convergence after 15 iterations.
## Gradient range [-1.38077e-06,5.605702e-06]
## (score -0.5466875 & scale 0.03734113).
## eigenvalue range [-6.525242e-13,56.66375].
## Model rank =  76 / 76 
## 
## Basis dimension (k) checking results. Low p-value (k-index<1) may
## indicate that k is too low, especially if edf is close to k'.
## 
##              k'    edf k-index p-value
## s(v)      4.000  1.223    1.08    0.76
## s(t)      4.000  3.086    1.04    0.64
## s(o)      4.000  2.379    0.90    0.11
## s(v,IHO) 20.000  4.630    1.08    0.80
## s(t,IHO) 20.000  0.814    1.04    0.61
## s(o,IHO) 20.000  2.883    0.90    0.11
```

```r
worm_plot(GAM4, method = "simulate", point_col = "steelblue", point_alpha = 0.4)
```

![](PanArctic_DSL_statistics_files/figure-html/basis-size-residuals-2.png)<!-- -->

Check residuals versus covariates.


```r
resid_GAM <- bind_cols(SA_df, residuals.gam(GAM4)) %>%
  rename(resid = "...15")
plot_grid(resid_GAM %>%
            ggplot(aes(x = v, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "velocity") +
            geom_point(), 
          resid_GAM %>%
            ggplot(aes(x = t, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "temperature") +
            geom_point(), 
          resid_GAM %>%
            ggplot(aes(x = o, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "open water") +
            geom_point(), 
          resid_GAM %>%
            ggplot(aes(x = IHO, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "area") +
            geom_boxplot(fill = NA), 
          resid_GAM %>%
            ggplot(aes(x = year, y = resid)) + 
            geom_hline(yintercept = 0, col = "red") +
            labs(x = "year") +
            geom_boxplot(fill = NA))
```

![](PanArctic_DSL_statistics_files/figure-html/residuals-covariates-1.png)<!-- -->

Check auto-correlation wihtin the model. All seems in order.


```r
par(mfrow = c(1, 2), mar = c(4, 4, 4, 0.5))
acf(resid(GAM4), main = "ACF")
pacf(resid(GAM4), main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/acf-pacf-1.png)<!-- -->

## Model selection

I turn on the double penalty (`select = TRUE`) and check the covariates that has been shrunk. 


```r
GAM4_p <- gam(SA_int_clean ~ year + s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + s(v, IHO, bs = "fs", k = 5) + 
              s(t, IHO, bs = "fs", k = 5) + s(o, IHO, bs = "fs", k = 5),
            data = SA_df, family = "gaussian", method = "REML", select=TRUE)
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
summary(GAM4_p)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int_clean ~ year + s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + 
##     s(v, IHO, bs = "fs", k = 5) + s(t, IHO, bs = "fs", k = 5) + 
##     s(o, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.150      1.617   9.989   <2e-16 ***
## year2016       2.996      1.848   1.621    0.108    
## year2017       2.033      1.731   1.174    0.243    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df     F  p-value    
## s(v)     0.0002123      4 0.000 0.153662    
## s(t)     0.2200515      4 0.071 0.008110 ** 
## s(o)     1.1096904      4 1.018 0.004604 ** 
## s(v,IHO) 6.7540544     19 0.783 0.011258 *  
## s(t,IHO) 2.7589777     18 0.462 0.007070 ** 
## s(o,IHO) 4.0161195     19 0.798 0.000177 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.35   Deviance explained = 44.2%
## -REML = 402.12  Scale est. = 39.768    n = 120
```

## Model visualization

Quick partial effects plot with `gratia::draw`.


```r
draw(GAM4, residuals = TRUE, scales = "fixed")
```

![](PanArctic_DSL_statistics_files/figure-html/draw-GAM-1.png)<!-- -->
I predict the data and then use that dataframe to plot it in a "clean" way with `ggplot`.


```r
# Fitted GAMs
pred_v <- fitted_values(GAM4, terms = c("s(v)", "s(v,IHO)")) %>% # Fit velocity
  rename(v_fit = fitted, 
         v_se = se,
         v_lower = lower,
         v_upper = upper) 
pred_t <- fitted_values(GAM4, terms = c("s(t)", "s(t,IHO)")) %>% # Fit temperature
  rename(t_fit = fitted, 
         t_se = se,
         t_lower = lower,
         t_upper = upper)
pred_o <- fitted_values(GAM4, terms = c("s(o)", "s(o,IHO)")) %>% # Fit open water
  rename(o_fit = fitted, 
         o_se = se,
         o_lower = lower,
         o_upper = upper)
pred_GAM <- left_join(SA_df, pred_v, by = c("IHO", "v", "t", "o")) %>%
  left_join(., pred_t, by = c("IHO", "v", "t", "o")) %>%
  left_join(., pred_o, by = c("IHO", "v", "t", "o")) 

# Partial effects
part_eff_GAM <- smooth_estimates(GAM4) %>% # Evaluate smooths
  add_confint()
part_res_SA <- SA_df %>% # add partial residuals to data
  add_partial_residuals(GAM4)

# Save predictions
save(pred_GAM, GAM4, file = "data/statistics/GAM_results.RData")

rm(pred_v, pred_t, pred_o)
```

Plot fitted GAM with `ggplot`.


```r
# col_pal <- c("#5BBCD6", "#00A08A", "#F2AD00", "#F98400", "#FF0000")# wesanderson::wes_palette("Darjeeling1", type = "discrete") 
col_pal <- c("#5BBCD6",  "#00A08A", "#F98400", "#FF0000")# wesanderson::wes_palette("Darjeeling1", type = "discrete") 

GAM_plot <- plot_grid(get_legend(ggplot(data = pred_GAM, aes(x = v, y = v_fit, fill = IHO, col = IHO)) + 
                                   geom_point() +
                                   scale_colour_manual(values = col_pal) + 
                                   scale_fill_manual(values = col_pal) + 
                                   theme(legend.position = "top", legend.title = element_blank())),
          plot_grid(
            pred_GAM %>% # velocity
              ggplot(aes(x = v)) +
              geom_ribbon(aes(ymin = v_lower, ymax = v_upper, fill = IHO), alpha = 0.1) +
              geom_line(aes(y = v_fit, col = IHO), size = 1) +
              scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
              scale_colour_manual(values = col_pal) + 
              scale_fill_manual(values = col_pal) + 
              labs(x = expression("Velocity at 380 m (cm s"^-1*")"), 
                   y = expression("S"[A]*"")) + 
              theme(legend.position = "none"),
            pred_GAM %>% # temperature
              ggplot(aes(x = t)) +
              geom_ribbon(aes(ymin = t_lower, ymax = t_upper, fill = IHO), alpha = 0.1) +
              geom_line(aes(y = t_fit, col = IHO), size = 1) +
              scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
              scale_colour_manual(values = col_pal) + 
              scale_fill_manual(values = col_pal) + 
              labs(x = "Temp. at 380 m (°C)", 
                   y = expression("S"[A]*"")) + 
              theme(legend.position = "none"),
            pred_GAM %>% # open water
              ggplot(aes(x = o)) +
              geom_ribbon(aes(ymin = o_lower, ymax = o_upper, fill = IHO), alpha = 0.1) +
              geom_line(aes(y = o_fit, col = IHO), size = 1) +
              scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
              scale_colour_manual(values = col_pal) + 
              scale_fill_manual(values = col_pal) + 
              labs(x = "Open water (days)", 
                   y = expression("S"[A]*"")) + 
              theme(legend.position = "none")), 
          nrow = 2, rel_heights = c(0.1, 1))
GAM_plot
```

![](PanArctic_DSL_statistics_files/figure-html/ggplot-GAM-1.png)<!-- -->

Plot partial effects.


```r
plot_grid(part_eff_GAM %>%
            filter(smooth == "s(v)") %>%
            ggplot() +
            geom_rug(data = part_res_SA, aes(x = v), sides = "b", length = grid::unit(0.02, "npc")) +
            geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = v), alpha = 0.2) +
            geom_point(data = part_res_SA, aes(x = v, y = `s(v)`), col = "steelblue3", cex = 1.5) +
            geom_line(aes(x = v, y = est), lwd = 1.2) +
            labs(y = "Partial effect"),
          part_eff_GAM %>%
            filter(smooth == "s(t)") %>%
            ggplot() +
            geom_rug(data = part_res_SA, aes(x = t), sides = "b", length = grid::unit(0.02, "npc")) +
            geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = t), alpha = 0.2) +
            geom_point(data = part_res_SA, aes(x = t, y = `s(t)`), col = "steelblue3", cex = 1.5) +
            geom_line(aes(x = t, y = est), lwd = 1.2) +
            labs(y = "Partial effect"),
          part_eff_GAM %>%
            filter(smooth == "s(o)") %>%
            ggplot() +
            geom_rug(data = part_res_SA, aes(x = o), sides = "b", length = grid::unit(0.02, "npc")) +
            geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = o), alpha = 0.2) +
            geom_point(data = part_res_SA, aes(x = o, y = `s(o)`), col = "steelblue3", cex = 1.5) +
            geom_line(aes(x = o, y = est), lwd = 1.2) +
            labs(y = "Partial effect"))
```

![](PanArctic_DSL_statistics_files/figure-html/partial-effect-plots-1.png)<!-- -->
