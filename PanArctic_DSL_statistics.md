---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/10 at 16:18"
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
  filter(depth == 222) %>% # Select data at 222 m depth and remove cells with no DSL
  mutate(velocity = velocity * 100, # Convert to cm/s
         IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                     IHO_area == "West Arctic Ocean" ~ "WAO_BF",
                                     IHO_area == "Beaufort Sea" ~ "WAO_BF",
                                     IHO_area == "The Northwestern Passages" ~ "CAA",
                                     IHO_area == "Baffin Bay" ~ "BB",
                                     IHO_area == "Davis Strait" ~ "DS"),
                            levels = c("WAO_BF", "CAA", "BB", "DS", "EAO"))) %>%
  ungroup() %>%
  dplyr::select(-xc_na, -yc_na, -year_na)

stat_laea %>% group_by(IHO_area) %>% summarise(n = n())
```

```
## # A tibble: 5 x 2
##   IHO_area     n
##   <fct>    <int>
## 1 WAO_BF      12
## 2 CAA         21
## 3 BB          44
## 4 DS          24
## 5 EAO         25
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
  ggtitle("Temperature at 222 m depth") +
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
  ggtitle("Current velocity at 222 m depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-velocity-1.png)<!-- -->

# Spatial and interannual variability


```r
SA_diff <- SA_grid_laea %>%
  filter(empty_clean == F) %>%
  dplyr::select(year, IHO_area, NASC_int_clean, SA_int_clean, NASC_int, CM) %>%
  mutate(IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                     IHO_area == "West Arctic Ocean" ~ "WAO_BF",
                                     IHO_area == "Beaufort Sea" ~ "WAO_BF",
                                     IHO_area == "The Northwestern Passages" ~ "CAA",
                                     IHO_area == "Baffin Bay" ~ "BB",
                                     IHO_area == "Davis Strait" ~ "DS"),
                           levels = c("WAO_BF", "CAA", "BB", "DS", "EAO")),
         year = factor(year))
```

## All areas

I check whether there are inter-annual and spatial variability in S~A~ and the centre of mass across years and areas.


```r
# Visualize data
plot_grid(SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = NASC_int_clean)),
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
  identify_outliers(NASC_int_clean) 
```

```
## # A tibble: 12 x 8
##    year  IHO_area NASC_int_clean SA_int_clean NASC_int    CM is.outlier
##    <fct> <fct>             <dbl>        <dbl>    <dbl> <dbl> <lgl>     
##  1 2015  BB                 946.         29.8     947.  313. TRUE      
##  2 2015  BB                1528.         31.8    1529.  292. TRUE      
##  3 2016  BB                 765.         28.8     765.  353. TRUE      
##  4 2016  BB                 548.         27.4     548.  330. TRUE      
##  5 2016  BB                1029.         30.1    1030.  338. TRUE      
##  6 2016  EAO               2434.         33.9    2435.  471. TRUE      
##  7 2016  DS                1036.         30.2    1036.  364. TRUE      
##  8 2017  CAA                650.         28.1     650.  393. TRUE      
##  9 2017  CAA                543.         27.3     543.  303. TRUE      
## 10 2017  EAO              14313.         41.6   14314.  412. TRUE      
## 11 2017  EAO              67315.         48.3   67316.  392. TRUE      
## 12 2017  EAO                649.         28.1     650.  310. TRUE      
## # ... with 1 more variable: is.extreme <lgl>
```

```r
SA_diff %>%
  group_by(year) %>%
  identify_outliers(CM) 
```

```
## # A tibble: 2 x 8
##   year  IHO_area NASC_int_clean SA_int_clean NASC_int    CM is.outlier
##   <fct> <fct>             <dbl>        <dbl>    <dbl> <dbl> <lgl>     
## 1 2015  WAO_BF             5.87         7.69     7.72  591. TRUE      
## 2 2015  WAO_BF             2.40         3.81     5.78  572. TRUE      
## # ... with 1 more variable: is.extreme <lgl>
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(NASC_int_clean) 
```

```
## # A tibble: 11 x 8
##    IHO_area year  NASC_int_clean SA_int_clean NASC_int    CM is.outlier
##    <fct>    <fct>          <dbl>        <dbl>    <dbl> <dbl> <lgl>     
##  1 BB       2015            946.         29.8     947.  313. TRUE      
##  2 BB       2015           1528.         31.8    1529.  292. TRUE      
##  3 BB       2016            765.         28.8     765.  353. TRUE      
##  4 BB       2016            548.         27.4     548.  330. TRUE      
##  5 BB       2016           1029.         30.1    1030.  338. TRUE      
##  6 DS       2016           1036.         30.2    1036.  364. TRUE      
##  7 DS       2017            228.         23.6     229.  341. TRUE      
##  8 EAO      2016           2434.         33.9    2435.  471. TRUE      
##  9 EAO      2017          14313.         41.6   14314.  412. TRUE      
## 10 EAO      2017          67315.         48.3   67316.  392. TRUE      
## 11 EAO      2017            649.         28.1     650.  310. TRUE      
## # ... with 1 more variable: is.extreme <lgl>
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(CM) 
```

```
## # A tibble: 9 x 8
##   IHO_area year  NASC_int_clean SA_int_clean NASC_int    CM is.outlier
##   <fct>    <fct>          <dbl>        <dbl>    <dbl> <dbl> <lgl>     
## 1 WAO_BF   2015            5.87         7.69     7.72  591. TRUE      
## 2 WAO_BF   2015            2.40         3.81     5.78  572. TRUE      
## 3 WAO_BF   2017           49.9         17.0     50.6   300. TRUE      
## 4 BB       2015           26.8         14.3     28.8   415. TRUE      
## 5 BB       2016           84.7         19.3     86.5   261. TRUE      
## 6 BB       2016           22.2         13.5     24.7   393. TRUE      
## 7 BB       2016           45.6         16.6     48.0   427. TRUE      
## 8 BB       2016           45.9         16.6     47.4   401. TRUE      
## 9 BB       2017           19.8         13.0     20.8   265. TRUE      
## # ... with 1 more variable: is.extreme <lgl>
```

```r
## Normality: Checked by looking at the anova residuals with QQ plot
ggqqplot(SA_diff, "NASC_int_clean", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-2.png)<!-- -->

```r
ggqqplot(SA_diff, "CM", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-3.png)<!-- -->

```r
ggqqplot(SA_diff, "NASC_int_clean", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-4.png)<!-- -->

```r
ggqqplot(SA_diff, "CM", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-5.png)<!-- -->

```r
## Homogeneity of variance
plot(lm(NASC_int_clean ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-6.png)<!-- -->

```r
plot(lm(CM ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-7.png)<!-- -->

```r
plot(lm(NASC_int_clean ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-8.png)<!-- -->

```r
plot(lm(CM ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-9.png)<!-- -->

I used a non parametric Kruskal-Wallis test to test for interannual and spatial variability because the variance is not homogenous and not very normal.


```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(NASC_int_clean ~ year)
```

```
## # A tibble: 1 x 6
##   .y.                n statistic    df      p method        
## * <chr>          <int>     <dbl> <int>  <dbl> <chr>         
## 1 NASC_int_clean   127      5.74     2 0.0566 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ year)
```

```
## # A tibble: 1 x 6
##   .y.       n statistic    df     p method        
## * <chr> <int>     <dbl> <int> <dbl> <chr>         
## 1 CM      127      3.40     2 0.183 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(NASC_int_clean ~ IHO_area)
```

```
## # A tibble: 1 x 6
##   .y.                n statistic    df        p method        
## * <chr>          <int>     <dbl> <int>    <dbl> <chr>         
## 1 NASC_int_clean   127      21.3     4 0.000278 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ IHO_area)
```

```
## # A tibble: 1 x 6
##   .y.       n statistic    df        p method        
## * <chr> <int>     <dbl> <int>    <dbl> <chr>         
## 1 CM      127      48.9     4 6.22e-10 Kruskal-Wallis
```

The center of mass did not vary significantly between years (H = 1.652899, p = 0.438) but varied significantly among areas (*H* = 57.8495, *p* < 0.001). Mesopelagic backscatter S~A~ varied significantly between areas (*H* = 26.19394, *p* < 0.001) but not years (*H* = 8.6464, *p* = 0.013). Thus, in the S~A~ HGAM a random effect for change in S~A~ per area is likely needed.

## Within areas

I check whether there are inter-annual variability in S~A~ across years within areas. This will be used to decide whether a random effect is needed in the GAM. I used a non parametric Kruskal-Wallis test to test for interannual variability within areas.


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(IHO_area) %>%
  kruskal_test(NASC_int_clean ~ year)
```

```
## # A tibble: 5 x 7
##   IHO_area .y.                n statistic    df      p method        
## * <fct>    <chr>          <int>     <dbl> <int>  <dbl> <chr>         
## 1 WAO_BF   NASC_int_clean    12    3.69       2 0.158  Kruskal-Wallis
## 2 CAA      NASC_int_clean    21   11.8        2 0.0027 Kruskal-Wallis
## 3 BB       NASC_int_clean    45    5.14       2 0.0764 Kruskal-Wallis
## 4 DS       NASC_int_clean    24    0.0857     2 0.958  Kruskal-Wallis
## 5 EAO      NASC_int_clean    25    4.87       2 0.0878 Kruskal-Wallis
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
corr_WAO_BF <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "WAO_BF") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_WAO_BF, type = "lower", lab = T, title = "corr WAO BF")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-4.png)<!-- -->

```r
corr_CAA <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "CAA") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_CAA, type = "lower", lab = T, title = "corr CAA")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-5.png)<!-- -->

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
GAM4 <- gam(SA_int_n ~ IHO + s(v, k = 5, bs = "tp") + s(t, k = 5, bs = "tp") + s(o, k = 5, bs = "tp") +
              s(v, IHO, bs = "fs", k = 5) + s(t, IHO, bs = "fs", k = 5) + s(o, IHO, bs = "fs", k = 5),
            data = SA_df, family = "gaussian", method = "REML")
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
# Model GI
GAM5 <- gam(SA_int_n ~ IHO + s(v, k = 5, bs = "tp") + s(t, k = 5, bs = "tp") + s(o, k = 5, bs = "tp") +
              s(v, by = IHO, k = 5, bs = "tp") + s(t, by = IHO, k = 5, bs = "tp") + s(o, by = IHO, k = 5, bs = "tp"),
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
<div id="htmlwidget-1521e1629accc8c2a8e4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1521e1629accc8c2a8e4">{"x":{"filter":"none","vertical":false,"data":[["GAM3","GAM5","GAM2","GAM4","GAM1"],[35.642,35.409,27.112,22.458,7.959],[68.76,68.54,56.13,52.25,35.27],[0.59,0.58,0.48,0.45,0.33],[1.45,3.9,-3.36,2.2,1.41],[-46.918,-46.483,-21.188,-19.822,-10.48],[0,0.44,25.73,27.1,36.44],[0.5541604533,0.4458373822,1.4336e-06,7.241e-07,6.8e-09]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## SA_int_n ~ IHO + s(v, k = 5, bs = "tp") + s(t, k = 5, bs = "tp") + 
##     s(o, k = 5, bs = "tp") + s(v, IHO, bs = "fs", k = 5) + s(t, 
##     IHO, bs = "fs", k = 5) + s(o, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.507007   0.077452   6.546 2.04e-09 ***
## IHOCAA      -0.032185   0.089917  -0.358    0.721    
## IHOBB        0.013235   0.077978   0.170    0.866    
## IHODS       -0.127884   0.134543  -0.951    0.344    
## IHOEAO       0.002248   0.124390   0.018    0.986    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df     F  p-value    
## s(v)     2.694e+00  3.182 3.091 0.027741 *  
## s(t)     1.000e+00  1.000 4.637 0.033515 *  
## s(o)     2.133e+00  2.405 5.668 0.007431 ** 
## s(v,IHO) 4.182e-05 19.000 0.000 0.319376    
## s(t,IHO) 1.229e+00 17.000 0.098 0.214854    
## s(o,IHO) 5.718e+00 19.000 0.909 0.000784 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.448   Deviance explained = 52.3%
## -REML = 2.2014  Scale est. = 0.040778  n = 126
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
## Gradient range [-5.06926e-07,1.630614e-05]
## (score 2.201414 & scale 0.04077809).
## Hessian positive definite, eigenvalue range [2.070206e-07,59.16527].
## Model rank =  92 / 92 
## 
## Basis dimension (k) checking results. Low p-value (k-index<1) may
## indicate that k is too low, especially if edf is close to k'.
## 
##                k'      edf k-index p-value
## s(v)     4.00e+00 2.69e+00    0.98    0.33
## s(t)     4.00e+00 1.00e+00    1.12    0.89
## s(o)     4.00e+00 2.13e+00    0.92    0.18
## s(v,IHO) 2.50e+01 4.18e-05    0.98    0.32
## s(t,IHO) 2.50e+01 1.23e+00    1.12    0.90
## s(o,IHO) 2.50e+01 5.72e+00    0.92    0.16
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
GAM4_p <- gam(SA_int_n ~ IHO + s(v, k = 5, bs = "tp") + s(t, k = 5, bs = "tp") + s(o, k = 5, bs = "tp") +
                s(v, IHO, bs = "fs", k = 5) + s(t, IHO, bs = "fs", k = 5) + s(o, IHO, bs = "fs", k = 5),
              data = SA_df, family = "gaussian", method = "REML", select = TRUE)
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
## SA_int_n ~ IHO + s(v, k = 5, bs = "tp") + s(t, k = 5, bs = "tp") + 
##     s(o, k = 5, bs = "tp") + s(v, IHO, bs = "fs", k = 5) + s(t, 
##     IHO, bs = "fs", k = 5) + s(o, IHO, bs = "fs", k = 5)
## 
## Parametric coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.495570   0.076987   6.437 3.57e-09 ***
## IHOCAA      -0.009073   0.085689  -0.106    0.916    
## IHOBB        0.023893   0.077135   0.310    0.757    
## IHODS       -0.116148   0.131770  -0.881    0.380    
## IHOEAO      -0.011798   0.119167  -0.099    0.921    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                edf Ref.df     F  p-value    
## s(v)     2.155e-07      4 0.000  0.66518    
## s(t)     8.062e-01      4 1.040  0.00866 ** 
## s(o)     1.754e+00      4 3.288 7.23e-06 ***
## s(v,IHO) 4.497e+00     19 0.689  0.00361 ** 
## s(t,IHO) 1.800e+00     18 0.154  0.13321    
## s(o,IHO) 5.540e+00     20 1.245 1.82e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.466   Deviance explained = 54.5%
## -REML = -1.6833  Scale est. = 0.039476  n = 126
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
col_pal <- c("#5BBCD6", "#00A08A", "#F2AD00", "#F98400", "#FF0000")# wesanderson::wes_palette("Darjeeling1", type = "discrete")
# col_pal <- c("#5BBCD6",  "#00A08A", "#F98400", "#FF0000")# wesanderson::wes_palette("Darjeeling1", type = "discrete") 

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
              # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
              scale_colour_manual(values = col_pal) + 
              scale_fill_manual(values = col_pal) + 
              labs(x = expression("Velocity at 380 m (cm s"^-1*")"), 
                   y = expression("S"[A]*"")) + 
              theme(legend.position = "none"),
            pred_GAM %>% # temperature
              ggplot(aes(x = t)) +
              geom_ribbon(aes(ymin = t_lower, ymax = t_upper, fill = IHO), alpha = 0.1) +
              geom_line(aes(y = t_fit, col = IHO), size = 1) +
              # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
              scale_colour_manual(values = col_pal) + 
              scale_fill_manual(values = col_pal) + 
              labs(x = "Temp. at 380 m (°C)", 
                   y = expression("S"[A]*"")) + 
              theme(legend.position = "none"),
            pred_GAM %>% # open water
              ggplot(aes(x = o)) +
              geom_ribbon(aes(ymin = o_lower, ymax = o_upper, fill = IHO), alpha = 0.1) +
              geom_line(aes(y = o_fit, col = IHO), size = 1) +
              # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
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
