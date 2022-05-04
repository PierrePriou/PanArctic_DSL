---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/04 at 20:10"
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
  filter(empty == F)
phy_laea <- phy_grid_laea %>% # Tidy remote sensing dataset for joining
  dplyr::select(-lat, -lon)
seaice_laea <- seaice_grid_laea %>%
  dplyr::select(-lat, -lon)

stat_laea <- left_join(SA_laea, phy_laea, by = c("year", "area", "xc", "yc", "cell_res")) %>% # Join acoustic and physics
  left_join(., seaice_laea, by = c("year", "area", "xc", "yc", "cell_res")) %>%
  # For cells that have NaN, get the mean of the surrounding cells
  rowwise() %>%
  mutate(xc_na = if_else(is.na(mean_ice_conc) == T, xc, NaN),
         yc_na = if_else(is.na(mean_ice_conc) == T, yc, NaN),
         year_na = if_else(is.na(mean_ice_conc) == T, year, NaN), 
         mean_ice_conc = if_else(is.na(mean_ice_conc) == T, mean(pull(subset(seaice_grid_laea,
                                                                             xc >= xc_na - 50 &
                                                                               xc <= xc_na + 50 & 
                                                                               yc >= yc_na - 50 & 
                                                                               yc <= yc_na + 50 & 
                                                                               year == year_na,
                                                                             select = mean_ice_conc),
                                                                      mean_ice_conc),
                                                                 na.rm = T),
                                 mean_ice_conc),
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
         ice_break = if_else(is.na(ice_break) == T, mean(pull(subset(seaice_grid_laea,
                                                                     xc >= xc_na - 50 &
                                                                       xc <= xc_na + 50 & 
                                                                       yc >= yc_na - 50 & 
                                                                       yc <= yc_na + 50 & 
                                                                       year == year_na,
                                                                     select = ice_break),
                                                              ice_break),
                                                         na.rm = T),
                             ice_break),
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
  filter(depth == 380) %>% # Select data at 380 m depth
  mutate(SA_int = 10 * log10(NASC_int),
         velocity = velocity * 100, # Convert to cm/s
         IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                     IHO_area == "West Arctic Ocean" ~ "WAO_BF",
                                     IHO_area == "Beaufort Sea" ~ "WAO_BF",
                                     IHO_area == "The Northwestern Passages" ~ "CAA",
                                     IHO_area == "Baffin Bay" ~ "BB",
                                     IHO_area == "Davis Strait" ~ "DS"),
                            levels = c("WAO_BF", "CAA", "BB", "DS", "EAO"))) %>%
  ungroup() %>%
  dplyr::select(-vxo, -vyo, -mean_temp_area_depth, -mean_velo_area_depth, -xc_na, -yc_na, -year_na, -siconc, -sithick,
                -temp_anomaly, -velocity_anomaly)
```

# Data exploration 

Maps of all variables.


```r
stat_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_point(aes(col = IHO_area)) +
  ggtitle("regions") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-SA-int-1.png)<!-- -->


```r
stat_laea %>% # Ice concentration
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_tile(aes(fill = mean_ice_conc)) +
  scale_fill_cmocean("Ice (%)", name = "ice", na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Sea ice concentration") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-siconc-1.png)<!-- -->


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
  geom_tile(aes(fill = ice_break)) +
  scale_fill_viridis_c("Day", option = "viridis", na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Day ice breakup") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-ice-day-1.png)<!-- -->


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
stat_laea %>% # Salinity
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_tile(aes(fill = so)) +
  scale_fill_cmocean("Sal (psu)", name = "haline", na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Salinity at 380 m depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-so-1.png)<!-- -->


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


```r
stat_laea %>% # Mixed layer depth
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_tile(aes(fill = mlotst)) +
  scale_fill_viridis_c("depth (m)", option = "plasma", direction = -1, na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Mixed layer depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-mld-1.png)<!-- -->

# Spatial and interannual variability


```r
SA_diff <- stat_laea %>%
  dplyr::select(year, IHO_area, SA_int, NASC_int, CM) %>%
  mutate(year = factor(year))
```

## All areas

I check whether there are inter-annual and spatial variability in S~A~ and the centre of mass across years and areas.


```r
# Visualize data
plot_grid(SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = NASC_int)),
          SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = SA_int)),
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
  identify_outliers(NASC_int) 
```

```
## # A tibble: 11 x 7
##    year  IHO_area SA_int NASC_int    CM is.outlier is.extreme
##    <fct> <fct>     <dbl>    <dbl> <dbl> <lgl>      <lgl>     
##  1 2015  BB         29.8     947.  313. TRUE       TRUE      
##  2 2015  BB         31.8    1529.  292. TRUE       TRUE      
##  3 2016  BB         28.8     765.  353. TRUE       FALSE     
##  4 2016  BB         30.1    1030.  338. TRUE       FALSE     
##  5 2016  EAO        33.9    2435.  471. TRUE       TRUE      
##  6 2016  DS         30.2    1036.  364. TRUE       FALSE     
##  7 2017  CAA        28.1     650.  393. TRUE       FALSE     
##  8 2017  CAA        27.4     543.  303. TRUE       FALSE     
##  9 2017  EAO        41.6   14314.  412. TRUE       TRUE      
## 10 2017  EAO        48.3   67316.  392. TRUE       TRUE      
## 11 2017  EAO        28.1     650.  310. TRUE       FALSE
```

```r
SA_diff %>%
  group_by(year) %>%
  identify_outliers(SA_int) 
```

```
## # A tibble: 5 x 7
##   year  IHO_area SA_int NASC_int    CM is.outlier is.extreme
##   <fct> <fct>     <dbl>    <dbl> <dbl> <lgl>      <lgl>     
## 1 2015  CAA        11.5     14.2  332. TRUE       FALSE     
## 2 2015  BB         29.8    947.   313. TRUE       TRUE      
## 3 2015  BB         31.8   1529.   292. TRUE       TRUE      
## 4 2017  EAO        41.6  14314.   412. TRUE       TRUE      
## 5 2017  EAO        48.3  67316.   392. TRUE       TRUE
```

```r
SA_diff %>%
  group_by(year) %>%
  identify_outliers(CM) 
```

```
## [1] year       IHO_area   SA_int     NASC_int   CM         is.outlier is.extreme
## <0 rows> (or 0-length row.names)
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(NASC_int) 
```

```
## # A tibble: 9 x 7
##   IHO_area year  SA_int NASC_int    CM is.outlier is.extreme
##   <fct>    <fct>  <dbl>    <dbl> <dbl> <lgl>      <lgl>     
## 1 WAO_BF   2016    12.5     17.8  364. TRUE       FALSE     
## 2 BB       2015    29.8    947.   313. TRUE       FALSE     
## 3 BB       2015    31.8   1529.   292. TRUE       TRUE      
## 4 BB       2016    28.8    765.   353. TRUE       FALSE     
## 5 BB       2016    30.1   1030.   338. TRUE       FALSE     
## 6 DS       2016    30.2   1036.   364. TRUE       TRUE      
## 7 DS       2017    23.6    229.   341. TRUE       FALSE     
## 8 EAO      2017    41.6  14314.   412. TRUE       TRUE      
## 9 EAO      2017    48.3  67316.   392. TRUE       TRUE
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(SA_int) 
```

```
## # A tibble: 3 x 7
##   IHO_area year  SA_int NASC_int    CM is.outlier is.extreme
##   <fct>    <fct>  <dbl>    <dbl> <dbl> <lgl>      <lgl>     
## 1 WAO_BF   2016    12.5     17.8  364. TRUE       FALSE     
## 2 DS       2016    12.5     17.9  369. TRUE       FALSE     
## 3 DS       2016    30.2   1036.   364. TRUE       FALSE
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(CM) 
```

```
## # A tibble: 4 x 7
##   IHO_area year  SA_int NASC_int    CM is.outlier is.extreme
##   <fct>    <fct>  <dbl>    <dbl> <dbl> <lgl>      <lgl>     
## 1 WAO_BF   2017    17.0     50.6  300. TRUE       FALSE     
## 2 BB       2016    13.9     24.7  393. TRUE       FALSE     
## 3 BB       2016    16.8     48.0  427. TRUE       TRUE      
## 4 BB       2016    16.8     47.4  401. TRUE       FALSE
```

```r
## Normality: Checked by looking at the anova residuals with QQ plot
ggqqplot(SA_diff, "NASC_int", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-2.png)<!-- -->

```r
ggqqplot(SA_diff, "SA_int", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-3.png)<!-- -->

```r
ggqqplot(SA_diff, "CM", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-4.png)<!-- -->

```r
ggqqplot(SA_diff, "NASC_int", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-5.png)<!-- -->

```r
ggqqplot(SA_diff, "SA_int", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-6.png)<!-- -->

```r
ggqqplot(SA_diff, "CM", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-7.png)<!-- -->

```r
# Data is fully normal

## Homogeneity of variance
plot(lm(NASC_int ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-8.png)<!-- -->

```r
plot(lm(SA_int ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-9.png)<!-- -->

```r
plot(lm(CM ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-10.png)<!-- -->

```r
plot(lm(NASC_int ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-11.png)<!-- -->

```r
plot(lm(SA_int ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-12.png)<!-- -->

```r
plot(lm(CM ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-13.png)<!-- -->

```r
SA_diff %>%
  levene_test(SA_int ~ year)
```

```
## # A tibble: 1 x 4
##     df1   df2 statistic     p
##   <int> <int>     <dbl> <dbl>
## 1     2    86     0.197 0.822
```

```r
# Not homogenous
```

I used a non parametric Kruskal-Wallis test to test for interannual and spatial variability because the variance is not homogenous and not very normal.


```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(NASC_int ~ year)
```

```
## # A tibble: 1 x 6
##   .y.          n statistic    df     p method        
## * <chr>    <int>     <dbl> <int> <dbl> <chr>         
## 1 NASC_int    89      2.49     2 0.289 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(SA_int ~ year)
```

```
## # A tibble: 1 x 6
##   .y.        n statistic    df     p method        
## * <chr>  <int>     <dbl> <int> <dbl> <chr>         
## 1 SA_int    89      2.49     2 0.289 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ year)
```

```
## # A tibble: 1 x 6
##   .y.       n statistic    df       p method        
## * <chr> <int>     <dbl> <int>   <dbl> <chr>         
## 1 CM       89      9.44     2 0.00892 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(NASC_int ~ IHO_area)
```

```
## # A tibble: 1 x 6
##   .y.          n statistic    df      p method        
## * <chr>    <int>     <dbl> <int>  <dbl> <chr>         
## 1 NASC_int    89      9.55     4 0.0488 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(SA_int ~ IHO_area)
```

```
## # A tibble: 1 x 6
##   .y.        n statistic    df      p method        
## * <chr>  <int>     <dbl> <int>  <dbl> <chr>         
## 1 SA_int    89      9.55     4 0.0488 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ IHO_area)
```

```
## # A tibble: 1 x 6
##   .y.       n statistic    df         p method        
## * <chr> <int>     <dbl> <int>     <dbl> <chr>         
## 1 CM       89      25.0     4 0.0000515 Kruskal-Wallis
```

The center of mass did not vary significantly between years (H = 1.652899, p = 0.438) but varied significantly among areas (*H* = 57.8495, *p* < 0.001). Mesopelagic backscatter S~A~ varied significantly between areas (*H* = 26.19394, *p* < 0.001) but not years (*H* = 8.6464, *p* = 0.013). Thus, in the S~A~ HGAM a random effect for change in S~A~ per area is likely needed.

## Within areas

I check whether there are inter-annual variability in S~A~ across years within areas. This will be used to decide whether a random effect is needed in the GAM. I used a non parametric Kruskal-Wallis test to test for interannual variability within areas.


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(IHO_area) %>%
  kruskal_test(NASC_int ~ year)
```

```
## # A tibble: 5 x 7
##   IHO_area .y.          n statistic    df      p method        
## * <fct>    <chr>    <int>     <dbl> <int>  <dbl> <chr>         
## 1 WAO_BF   NASC_int     4    0          1 1      Kruskal-Wallis
## 2 CAA      NASC_int    15    8.13       2 0.0171 Kruskal-Wallis
## 3 BB       NASC_int    36    2.39       2 0.303  Kruskal-Wallis
## 4 DS       NASC_int    23    0.0242     2 0.988  Kruskal-Wallis
## 5 EAO      NASC_int    11    2.18       2 0.336  Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(IHO_area) %>%
  kruskal_test(SA_int ~ year)
```

```
## # A tibble: 5 x 7
##   IHO_area .y.        n statistic    df      p method        
## * <fct>    <chr>  <int>     <dbl> <int>  <dbl> <chr>         
## 1 WAO_BF   SA_int     4    0          1 1      Kruskal-Wallis
## 2 CAA      SA_int    15    8.13       2 0.0171 Kruskal-Wallis
## 3 BB       SA_int    36    2.39       2 0.303  Kruskal-Wallis
## 4 DS       SA_int    23    0.0242     2 0.988  Kruskal-Wallis
## 5 EAO      SA_int    11    2.18       2 0.336  Kruskal-Wallis
```

There was no interannual differences in SA within each region (WAO_BF - H = 4.8395483, *p* = 0.0889; CAA - H = 6.9278752, *p* = 0.0313; BB - H = 0.0968, *p* = 0.0889; DS - H = 0.9470, *p* = 0.0889; EAO - H = 5.0263899, *p* = 0.0810). Therefore, there seems to be no need for a random effect for interannual variability within years.

# HGAM - Gaussian SA

## Data preparation


```r
SA_df <- stat_laea %>%
  dplyr::select(year, xc, yc, area, IHO_area, NASC_int, SA_int, velocity, thetao, openwater_duration, ice_week) %>%
  group_by(IHO_area) %>%
  mutate(year = factor(year),
         SA_int_n = (SA_int - min(SA_int)) / (max(SA_int) - min(SA_int)),
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
            ggplot(aes(x = v, y = SA_int, col = IHO)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = t, y = SA_int, col = IHO)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = o, y = SA_int, col = IHO)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = s, y = SA_int, col = IHO)) +
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
  dplyr::select(-year, -xc, -yc, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_EAO, type = "lower", lab = T, title = "corr EAO")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-1.png)<!-- -->

```r
corr_BB <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "BB") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_BB, type = "lower", lab = T, title = "corr BB")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-2.png)<!-- -->

```r
corr_DS <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "DS") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_DS, type = "lower", lab = T, title = "corr BB")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-3.png)<!-- -->

```r
corr_WAO <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "WAO_BF") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_WAO, type = "lower", lab = T, title = "corr WAO")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-4.png)<!-- -->

```r
corr_CAA <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "CAA") %>%
  dplyr::select(-year, -xc, -yc, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
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
GAM3 <- gam(SA_int_n ~ s(v, by = IHO, k = 5, bs = "tp") + s(t, by = IHO, k = 5, bs = "tp") + 
              s(o, by = IHO, k = 5, bs = "tp") + s(IHO, bs = "re"),
            data = SA_df, family = "gaussian", method = "REML")
# Model GS
GAM4 <- gam(SA_int_n ~ s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + s(v, IHO, bs = "fs", k = 5) + 
              s(t_n, IHO, bs = "fs", k = 5) + s(o_n, IHO, bs = "fs", k = 5),
            data = SA_df, family = "gaussian", method = "REML")
```

```
## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
## d smooths of same variable.
```

```r
# Model GI
GAM5 <- gam(SA_int_n ~ s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + s(v, by = IHO, k = 5, bs = "tp") + 
              s(t, by = IHO, k = 5, bs = "tp") + s(o, by = IHO, k = 5, bs = "tp") + s(IHO, bs = "re"),
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
<div id="htmlwidget-7a1a43c1329627d4ffb2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-7a1a43c1329627d4ffb2">{"x":{"filter":"none","vertical":false,"data":[["GAM5","GAM3","GAM4","GAM2","GAM1"],[22.651,22.76,22.65,21.386,5],[63.41,63.49,62.78,56.04,17.55],[0.53,0.53,0.54,0.47,0.15],[11.04,9.05,8.61,10.02,19.51],[-8.93,-8.908,-7.404,4.877,28.076],[0,0.02,1.53,13.81,37.01],[0.4071114087,0.4026577233,0.1898220168,0.0004088475,3.7e-09]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

`GAM3` and `GAM5` are the two models with similar performances, although `GAM3` seems to perform slightly better.


```r
summary(GAM3)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int_n ~ s(v, by = IHO, k = 5, bs = "tp") + s(t, by = IHO, 
##     k = 5, bs = "tp") + s(o, by = IHO, k = 5, bs = "tp") + s(IHO, 
##     bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.50409    0.09481   5.317 1.24e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                  edf Ref.df      F  p-value    
## s(v):IHOWAO_BF 1.000  1.000  2.296 0.134304    
## s(v):IHOCAA    1.000  1.000  1.669 0.200704    
## s(v):IHOBB     1.456  1.750  1.536 0.152803    
## s(v):IHODS     1.000  1.000  3.634 0.060788 .  
## s(v):IHOEAO    2.396  2.823  5.834 0.001450 ** 
## s(t):IHOWAO_BF 1.000  1.000  0.147 0.702793    
## s(t):IHOCAA    1.000  1.000  0.240 0.626024    
## s(t):IHOBB     1.000  1.000 16.695 0.000117 ***
## s(t):IHODS     1.000  1.000  2.164 0.145804    
## s(t):IHOEAO    1.000  1.000  0.058 0.809751    
## s(o):IHOWAO_BF 1.000  1.000  0.745 0.391013    
## s(o):IHOCAA    2.272  2.527  5.955 0.001259 ** 
## s(o):IHOBB     1.000  1.000  0.087 0.768454    
## s(o):IHODS     1.000  1.000  0.393 0.532771    
## s(o):IHOEAO    1.000  1.000  0.174 0.677848    
## s(IHO)         1.216  4.000  1.028 0.035992 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.532   Deviance explained = 63.5%
## -REML = 9.0497  Scale est. = 0.041174  n = 89
```

## Model checking

First I check the basis size k. `k-indexes` are > 1 or close to 1 so the basis size is large enough. The residual plot look good too.


```r
par(mfrow = c(2, 2))
gam.check(GAM3, rep = 500)
```

![](PanArctic_DSL_statistics_files/figure-html/basis-size-residuals-1.png)<!-- -->

```
## 
## Method: REML   Optimizer: outer newton
## full convergence after 15 iterations.
## Gradient range [-2.395962e-06,1.638428e-06]
## (score 9.049679 & scale 0.04117382).
## Hessian positive definite, eigenvalue range [2.318526e-08,36.53652].
## Model rank =  66 / 66 
## 
## Basis dimension (k) checking results. Low p-value (k-index<1) may
## indicate that k is too low, especially if edf is close to k'.
## 
##                  k'  edf k-index p-value
## s(v):IHOWAO_BF 4.00 1.00    1.07    0.66
## s(v):IHOCAA    4.00 1.00    1.07    0.71
## s(v):IHOBB     4.00 1.46    1.07    0.73
## s(v):IHODS     4.00 1.00    1.07    0.72
## s(v):IHOEAO    4.00 2.40    1.07    0.69
## s(t):IHOWAO_BF 4.00 1.00    1.21    0.97
## s(t):IHOCAA    4.00 1.00    1.21    0.96
## s(t):IHOBB     4.00 1.00    1.21    0.99
## s(t):IHODS     4.00 1.00    1.21    0.99
## s(t):IHOEAO    4.00 1.00    1.21    0.97
## s(o):IHOWAO_BF 4.00 1.00    1.00    0.47
## s(o):IHOCAA    4.00 2.27    1.00    0.54
## s(o):IHOBB     4.00 1.00    1.00    0.46
## s(o):IHODS     4.00 1.00    1.00    0.49
## s(o):IHOEAO    4.00 1.00    1.00    0.44
## s(IHO)         5.00 1.22      NA      NA
```

Check residuals versus covariates.


```r
resid_GAM <- bind_cols(SA_df, residuals.gam(GAM3)) %>%
  rename(resid = "...17")
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
acf(resid(GAM3), lag.max = 36, main = "ACF")
pacf(resid(GAM3), lag.max = 36, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/acf-pacf-1.png)<!-- -->

## Model selection

I turn on the double penalty (`select = TRUE`) and check the covariates that has been shrunk. All covariates are imporant in the model. 


```r
GAM3_p <- gam(SA_int ~ s(v, by = IHO, k = 5, bs = "tp") + s(t, by = IHO, k = 5, bs = "tp") + 
                s(o, by = IHO, k = 5, bs = "tp") + s(IHO, bs = "re"),
              data = SA_df, family = "gaussian", method = "REML", select = TRUE)
summary(GAM3_p)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(v, by = IHO, k = 5, bs = "tp") + s(t, by = IHO, k = 5, 
##     bs = "tp") + s(o, by = IHO, k = 5, bs = "tp") + s(IHO, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  20.9736     0.5751   36.47   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                      edf Ref.df     F  p-value    
## s(v):IHOWAO_BF 8.784e-05      4 0.000  0.38982    
## s(v):IHOCAA    9.298e-05      4 0.000  0.44688    
## s(v):IHOBB     6.155e-01      4 0.302  0.14862    
## s(v):IHODS     7.583e-01      4 0.784  0.04117 *  
## s(v):IHOEAO    1.916e+00      4 2.538  0.00199 ** 
## s(t):IHOWAO_BF 8.083e-01      3 1.404  0.02430 *  
## s(t):IHOCAA    3.376e-05      4 0.000  0.79533    
## s(t):IHOBB     9.226e-01      4 2.977  0.00036 ***
## s(t):IHODS     5.911e-01      4 0.361  0.11049    
## s(t):IHOEAO    3.280e-01      4 0.099  0.19020    
## s(o):IHOWAO_BF 1.987e-04      4 0.000  0.36215    
## s(o):IHOCAA    2.196e+00      4 5.252 6.48e-05 ***
## s(o):IHOBB     2.939e-05      4 0.000  0.86580    
## s(o):IHODS     3.954e-05      4 0.000  0.97374    
## s(o):IHOEAO    2.072e+00      4 7.526  < 2e-16 ***
## s(IHO)         1.045e-05      4 0.000  0.84915    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.53   Deviance explained = 58.5%
## -REML = 264.39  Scale est. = 16.777    n = 89
```


```r
draw(GAM3_p, scales = "fixed")
```

![](PanArctic_DSL_statistics_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

## Model visualization

Quick plots using `visreg`.


```r
visreg(GAM3, "v", "IHO", overlay = T, scale = "response", band = F, partial = T)
```

![](PanArctic_DSL_statistics_files/figure-html/plot-GAM3-1.png)<!-- -->

```r
visreg(GAM3, "t", "IHO", overlay = T, scale = "response", band = F, partial = T)
```

![](PanArctic_DSL_statistics_files/figure-html/plot-GAM3-2.png)<!-- -->

```r
visreg(GAM3, "o", "IHO", overlay = T, scale = "response", band = F, partial = T)
```

![](PanArctic_DSL_statistics_files/figure-html/plot-GAM3-3.png)<!-- -->

I predict the data and then use that dataframe to plot it in a "clean" way with `ggplot`.


```r
pred_v <- get_gam_predictions(GAM3, series = v, .comparison = IHO, series_length = 50) %>% 
  mutate(var = "v",
         signif = case_when(IHO == "WAO_BF" ~  0.05, # Significance level
                            IHO == "CAA" ~ 1,
                            IHO == "BB" ~ 0.01,
                            IHO == "DS" ~ 0.05,
                            IHO == "EAO" ~ 1),
         signif_group = factor(if_else(signif > 0.1, F, T))) %>%
  rename(idx = ".idx",
         fit = SA_int_n,
         value = v)
pred_t <- get_gam_predictions(GAM3, series = t, .comparison = IHO, series_length = 50) %>%
  mutate(var = "t",
         signif = case_when(IHO == "WAO_BF" ~ 0.001, # significance level
                            IHO == "CAA" ~ 1,
                            IHO == "BB" ~ 0.001,
                            IHO == "DS" ~ 1,
                            IHO == "EAO" ~ 0.01),
         signif_group = factor(if_else(signif > 0.1, F, T))) %>%
    rename(idx = ".idx",
           fit = SA_int_n,
           value = t)
pred_o <- get_gam_predictions(GAM3, series = o, .comparison = IHO, series_length = 50) %>%
  mutate(var = "o", 
         signif = case_when(IHO == "WAO_BF" ~ 0.01, # significance level
                            IHO == "CAA" ~ 0.001,
                            IHO == "BB" ~ 0.001,
                            IHO == "DS" ~ 1,
                            IHO == "EAO" ~ 0.05),
         signif_group = factor(if_else(signif > 0.1, F, T))) %>%
  rename(idx = ".idx",
         fit = SA_int_n,
         value = o)
# Combine data
pred_GAM <- bind_rows(pred_v, pred_t, pred_o)

# Save data
save(pred_GAM, GAM3, file = "data/statistics/GAM_results.RData")
```

`ggplot` plot.


```r
col_pal <- c("#5BBCD6", "#00A08A", "#F2AD00", "#F98400", "#FF0000")# wesanderson::wes_palette("Darjeeling1", type = "discrete") 

plot_grid(pred_GAM %>%
            filter(var == "v" & SE < 4) %>%
            ggplot() +
            geom_line(aes(x = value, y = fit, col = IHO, linetype = IHO), size = 0.75, alpha = 0.8) +
            # geom_ribbon(aes(x = value, ymin = CI_lower, ymax = CI_upper, fill = IHO), alpha = 0.1) +
            scale_colour_manual(name = "IHO", values = col_pal) + 
            scale_fill_manual(name = "IHO", values = col_pal) + 
            scale_y_continuous(limits = c(0,1)) +
            scale_linetype_manual(name = "IHO", values = c(2, 2, 2, 1, 1)) +
            labs(x = "Velo. 380 m", y = expression("S"[A]*"")),
          pred_GAM %>%
            filter(var == "t" & SE < 4) %>%
            ggplot() +
            geom_line(aes(x = value, y = fit, col = IHO, linetype = IHO), size = 0.75, alpha = 0.8) +
            # geom_ribbon(aes(x = value, ymin = CI_lower, ymax = CI_upper, fill = IHO), alpha = 0.1) +
            scale_colour_manual(name = "IHO", values = col_pal) + 
            scale_fill_manual(name = "IHO", values = col_pal) + 
            scale_y_continuous(limits = c(0,1)) +
            scale_linetype_manual(name = "IHO", values = c(2, 2, 1, 2, 2)) +
            labs(x = "Temp. 380 m", y = expression("S"[A]*"")),
          pred_GAM %>%
            filter(var == "o" & SE < 4) %>%
            ggplot() +
            geom_line(aes(x = value, y = fit, col = IHO, linetype = IHO), size = 0.75, alpha = 0.8) +
            # geom_ribbon(aes(x = value, ymin = CI_lower, ymax = CI_upper, fill = IHO), alpha = 0.1) +
            scale_colour_manual(name = "IHO", values = col_pal) + 
            scale_fill_manual(name = "IHO", values = col_pal) + 
            scale_y_continuous(limits = c(0,1)) +
            scale_linetype_manual(name = "IHO", values = c(2, 1, 2, 2, 2)) +
            labs(x = "Open water days", y = expression("S"[A]*"")))
```

```
## Warning: Removed 49 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 3 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 31 row(s) containing missing values (geom_path).
```

![](PanArctic_DSL_statistics_files/figure-html/ggplot-gam-1.png)<!-- -->
