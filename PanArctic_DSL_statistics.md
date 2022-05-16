---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/16 at 18:43"
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
options(dplyr.summarise.inform = F) # Suppress summarise() warning
```


```r
# Laea projection
cell_res <- 25 # Cell resolution in km
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931")
projection(arctic_laea) <- gsub("units=m", "units=km", 
                                projection(arctic_laea)) # Convert from m to km
res(arctic_laea) <- c(cell_res, cell_res) # Define the 100 km cell resolution

# Coastline shapefiles reprojected to laea
coast_10m_laea <- readOGR("data/bathy/ne_10m_land.shp", verbose = F) %>% 
  spTransform(CRSobj = crs("EPSG:4326")) %>% # Verify that the shapefile proj is the right one
  crop(extent(-180, 180, 0, 90)) %>% # Crop shapefile
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to a dataframe for ggplot
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
# Prepare datasets
SA_laea <- SA_grid_laea %>%  # Gridded SA
  dplyr::select(-lat, -lon) 
phy_laea <- phy_grid_laea %>% # Remote sensing: physics reanalysis
  dplyr::select(year, area, xc, yc, cell_res, depth, thetao, velocity, so)
seaice_laea <- seaice_grid_laea %>% # Remote sensing: seaice concentration
  dplyr::select(year, area, xc, yc, cell_res, openwater_duration, mean_ice_conc)

stat_laea <- SA_laea %>%
  left_join(., phy_laea, by = c("year", "xc", "yc", "area", "cell_res")) %>% 
  left_join(., seaice_laea, by = c("year", "xc", "yc", "area", "cell_res")) %>%
  left_join(., chl_grid_laea, by = c("year", "xc", "yc", "area", "cell_res")) %>%
  filter(depth == 222) # Select physics at 222 m depth
```

There are NAs in the remote sensing dataset.


```r
plot_grid(stat_laea %>% 
            ggplot(aes(x = xc,  y = yc)) +
            geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                         fill = "grey80") +
            geom_point(aes(col = chl)) +
            scale_colour_cmocean(name = "algae", na.value = "red") + 
            facet_wrap(~ year) +
            coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
            theme(axis.text = element_blank(), axis.ticks = element_blank(), 
                  axis.title = element_blank()),
          stat_laea %>% 
            ggplot(aes(x = xc,  y = yc)) +
            geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                         fill = "grey80") +
            geom_point(aes(col = openwater_duration)) +
            scale_colour_viridis_c("ow", na.value = "red") + 
            facet_wrap(~ year) +
            coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
            theme(axis.text = element_blank(), axis.ticks = element_blank(), 
                  axis.title = element_blank()),
          ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_statistics_files/figure-html/map-NA-1.png)<!-- -->

This is because the coverage of acoustic data and remote sensing products slightly differ close to the coasts. For sea ice concentration (and open water duration) the NAs are close to land while for ocean colour they are in ice covered areas. I use the average value of neighbouring cells (8 cells around the pixel of interest) to fill missing values.


```r
stat_laea <- stat_laea %>%
  rowwise() %>%
  mutate( 
    # Physics
    xc_phy = if_else(is.na(openwater_duration) == T, xc, NaN),
    yc_phy = if_else(is.na(openwater_duration) == T, yc, NaN),
    year_phy = if_else(is.na(openwater_duration) == T, year, NaN), 
    openwater_duration = if_else(
      is.na(openwater_duration) == T,
      mean(pull(subset(seaice_grid_laea, 
                       xc >= xc_phy - cell_res & xc <= xc_phy + cell_res & 
                         yc >= yc_phy - cell_res & yc <= yc_phy + cell_res & 
                         year == year_phy,
                       select = openwater_duration),
                openwater_duration),
           na.rm = T),
      openwater_duration),
    mean_ice_conc = if_else(
      is.na(mean_ice_conc) == T, 
      mean(pull(subset(seaice_grid_laea,
                       xc >= xc_phy - cell_res & xc <= xc_phy + cell_res & 
                         yc >= yc_phy - cell_res & yc <= yc_phy + cell_res & 
                         year == year_phy,
                       select = mean_ice_conc),
                mean_ice_conc),
           na.rm = T),
      mean_ice_conc),
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
  filter(depth == 222) %>% # Select physics at 222 m depth
  ungroup() %>%
  dplyr::select(-xc_phy, -yc_phy, -year_phy, -xc_chl, -yc_chl, -year_chl)
```

Check if I replace NAs correctly.


```r
plot_grid(stat_laea %>% 
            ggplot(aes(x = xc,  y = yc)) +
            geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                         fill = "grey80") +
            geom_point(aes(col = chl)) +
            scale_colour_cmocean(name = "algae", na.value = "red") + 
            facet_wrap(~ year) +
            coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
            theme(axis.text = element_blank(), axis.ticks = element_blank(), 
                  axis.title = element_blank()),
          stat_laea %>%
            ggplot(aes(x = xc,  y = yc)) +
            geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                         fill = "grey80") +
            geom_point(aes(col = openwater_duration)) +
            scale_colour_viridis_c("ow", na.value = "red") +
            facet_wrap(~ year) +
            coord_fixed(xlim = c(-2450, 600), ylim = c(-1500, 1800), expand = F) +
            theme(axis.text = element_blank(), axis.ticks = element_blank(),
                  axis.title = element_blank()),
          ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_statistics_files/figure-html/map-NA-fix-1.png)<!-- -->

NA replacement worked well for ice concentration. There are still missing chl values, so these cells are excluded from further analyses.


```r
stat_laea <- stat_laea %>%
  filter(is.na(chl) == F) 
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
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
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
## 1 WAO_BF      16
## 2 CAA         21
## 3 BB          47
## 4 DS          25
## 5 EAO         23
```

# Data exploration

Maps of all variables.


```r
stat_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  geom_point(aes(col = IHO_area)) +
  ggtitle("regions") +
  facet_wrap(~ year) +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-SA-int-1.png)<!-- -->


```r
seaice_grid_laea %>% # Openwater duration
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = openwater_duration)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_viridis_c("Day", option = "plasma", na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Openwater duration") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-ow-1.png)<!-- -->


```r
chl_grid_laea %>% # Temperature
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = chl)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("chl (mg m-3)", name = "algae", limits = c(0,5),
                     na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("chlorophyll") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-html/map-chl-1.png)<!-- -->


```r
phy_grid_laea %>% # Ice concentration
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = velocity*100)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("velo (cm/s)", name = "speed", limits = c(0,8),
                     na.value = "red") +
  facet_wrap(~ year, ncol = 3) +
  ggtitle("Current velocity at 222 m depth") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
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

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-1.png)<!-- -->

```r
# Check assumptions
## Outliers
SA_diff %>%
  group_by(year) %>%
  identify_outliers(NASC_int) 
```

```
## # A tibble: 14 × 7
##    year  IHO_area NASC_int SA_int    CM is.outlier is.extreme
##    <fct> <fct>       <dbl>  <dbl> <dbl> <lgl>      <lgl>     
##  1 2015  BB           947.   29.8  313. TRUE       TRUE      
##  2 2015  BB          1529.   31.8  292. TRUE       TRUE      
##  3 2016  BB           765.   28.8  353. TRUE       FALSE     
##  4 2016  BB           548.   27.4  330. TRUE       FALSE     
##  5 2016  BB          1030.   30.1  338. TRUE       TRUE      
##  6 2016  BB           512.   27.1  340. TRUE       FALSE     
##  7 2016  EAO         2435.   33.9  471. TRUE       TRUE      
##  8 2016  DS          1036.   30.2  364. TRUE       TRUE      
##  9 2017  CAA          481.   26.8  289. TRUE       FALSE     
## 10 2017  CAA          650.   28.1  393. TRUE       FALSE     
## 11 2017  CAA          543.   27.4  303. TRUE       FALSE     
## 12 2017  EAO        14314.   41.6  412. TRUE       TRUE      
## 13 2017  EAO        67316.   48.3  392. TRUE       TRUE      
## 14 2017  EAO          650.   28.1  310. TRUE       FALSE
```

```r
SA_diff %>%
  group_by(year) %>%
  identify_outliers(CM) 
```

```
## # A tibble: 5 × 7
##   year  IHO_area NASC_int SA_int    CM is.outlier is.extreme
##   <fct> <fct>       <dbl>  <dbl> <dbl> <lgl>      <lgl>     
## 1 2015  WAO_BF      0.412  -3.85  602. TRUE       FALSE     
## 2 2015  WAO_BF      7.72    8.88  591. TRUE       FALSE     
## 3 2015  WAO_BF      0.545  -2.63  564. TRUE       FALSE     
## 4 2015  WAO_BF      5.78    7.62  572. TRUE       FALSE     
## 5 2015  WAO_BF      0.770  -1.13  581. TRUE       FALSE
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(NASC_int) 
```

```
## # A tibble: 14 × 7
##    IHO_area year  NASC_int SA_int    CM is.outlier is.extreme
##    <fct>    <fct>    <dbl>  <dbl> <dbl> <lgl>      <lgl>     
##  1 WAO_BF   2016      59.0   17.7  382. TRUE       FALSE     
##  2 WAO_BF   2017      50.6   17.0  300. TRUE       FALSE     
##  3 WAO_BF   2017      56.7   17.5  365. TRUE       FALSE     
##  4 BB       2015     947.    29.8  313. TRUE       TRUE      
##  5 BB       2015    1529.    31.8  292. TRUE       TRUE      
##  6 BB       2016     765.    28.8  353. TRUE       FALSE     
##  7 BB       2016     548.    27.4  330. TRUE       FALSE     
##  8 BB       2016    1030.    30.1  338. TRUE       TRUE      
##  9 DS       2016    1036.    30.2  364. TRUE       TRUE      
## 10 DS       2017     229.    23.6  341. TRUE       FALSE     
## 11 EAO      2016    2435.    33.9  471. TRUE       TRUE      
## 12 EAO      2017   14314.    41.6  412. TRUE       TRUE      
## 13 EAO      2017   67316.    48.3  392. TRUE       TRUE      
## 14 EAO      2017     650.    28.1  310. TRUE       TRUE
```

```r
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(CM) 
```

```
## # A tibble: 3 × 7
##   IHO_area year  NASC_int SA_int    CM is.outlier is.extreme
##   <fct>    <fct>    <dbl>  <dbl> <dbl> <lgl>      <lgl>     
## 1 BB       2015      28.8   14.6  415. TRUE       FALSE     
## 2 BB       2016      48.0   16.8  427. TRUE       FALSE     
## 3 BB       2016      47.4   16.8  401. TRUE       FALSE
```

```r
## Normality: Checked by looking at the anova residuals with QQ plot
ggqqplot(SA_diff, "NASC_int", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-2.png)<!-- -->

```r
ggqqplot(SA_diff, "CM", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-3.png)<!-- -->

```r
ggqqplot(SA_diff, "NASC_int", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-4.png)<!-- -->

```r
ggqqplot(SA_diff, "CM", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-5.png)<!-- -->

```r
## Homogeneity of variance
plot(lm(NASC_int ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-6.png)<!-- -->

```r
plot(lm(CM ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-7.png)<!-- -->

```r
plot(lm(NASC_int ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-8.png)<!-- -->

```r
plot(lm(CM ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-html/all-interannual-spatial-prep-9.png)<!-- -->

I used a non parametric Kruskal-Wallis test to test for interannual and spatial variability because the variance is not homogenous and not very normal.


```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(NASC_int ~ year)
```

```
## # A tibble: 1 × 6
##   .y.          n statistic    df      p method        
## * <chr>    <int>     <dbl> <int>  <dbl> <chr>         
## 1 NASC_int   139      9.06     2 0.0108 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ year)
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df     p method        
## * <chr> <int>     <dbl> <int> <dbl> <chr>         
## 1 CM      139      4.16     2 0.125 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(NASC_int ~ IHO_area)
```

```
## # A tibble: 1 × 6
##   .y.          n statistic    df          p method        
## * <chr>    <int>     <dbl> <int>      <dbl> <chr>         
## 1 NASC_int   139      29.8     4 0.00000539 Kruskal-Wallis
```

```r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ IHO_area)
```

```
## # A tibble: 1 × 6
##   .y.       n statistic    df        p method        
## * <chr> <int>     <dbl> <int>    <dbl> <chr>         
## 1 CM      139      61.2     4 1.63e-12 Kruskal-Wallis
```

The center of mass did not vary significantly between years (H = 1.652899, p = 0.438) but varied significantly among areas (*H* = 57.8495, *p* \< 0.001). Mesopelagic backscatter S\~A\~ varied significantly between areas (*H* = 26.19394, *p* \< 0.001) but not years (*H* = 8.6464, *p* = 0.013). Thus, in the S\~A\~ HGAM a random effect for change in S\~A\~ per area is likely needed.

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

There was no interannual differences in SA within each region (WAO_BF - H = 4.8395483, *p* = 0.0889; CAA - H = 6.9278752, *p* = 0.0313; BB - H = 0.0968, *p* = 0.0889; DS - H = 0.9470, *p* = 0.0889; EAO - H = 5.0263899, *p* = 0.0810). Therefore, there seems to be no need for a random effect for interannual variability within years.

# HGAM - Gamma NASC

## Data preparation


```r
SA_df <- stat_laea %>%
  dplyr::select(year, xc, yc, IHO_area, NASC_int, SA_int, velocity, thetao, so, openwater_duration, mean_ice_conc, chl) %>%
  group_by(IHO_area) %>%
  mutate(year = factor(year),
         NASC_int_n = (NASC_int - min(NASC_int)) / (max(NASC_int) - min(NASC_int)),
         SA_int_n = (SA_int - min(SA_int)) / (max(SA_int) - min(SA_int))) %>%
  ungroup() %>%
  rename(IHO = IHO_area, 
         v = velocity,
         t = thetao,
         o = openwater_duration,
         si = mean_ice_conc,
         s = so)
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
  ggplot(aes(x = o, y = NASC_int, col = IHO)) +
  geom_point() +
  geom_smooth(method = "gam", se = F, col = "grey20") +
  facet_wrap(~ IHO, scales = "free") +
  theme(legend.position = "none")
```

```
## `geom_smooth()` using formula 'y ~ s(x, bs = "cs")'
```

![](PanArctic_DSL_statistics_files/figure-html/plot-data-area-2.png)<!-- -->

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

![](PanArctic_DSL_statistics_files/figure-html/plot-data-area-3.png)<!-- -->

Check correlations.


```r
corr <- SA_df %>% # Compute Spearman correlation matrix
  dplyr::select(-year, -xc, -yc, -IHO, -NASC_int_n, -SA_int_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr, type = "lower", lab = T, title = "corr")
```

![](PanArctic_DSL_statistics_files/figure-html/corr-area-1.png)<!-- -->

# Velocity and chlorophyll

## Model fitting

In model GS the random intercept is already included in the `s(bs = "fs")` term.


```r
# Model G
GAM_G <- gam(NASC_int ~ 
               # First order effects
               s(v, k = 5, bs = "tp") + # Global smoothers
               s(chl, k = 5, bs = "tp") + # Global smoothers
               s(IHO, bs = "re") + # Random intercept
               # Second order effects
               s(IHO, year, bs = "re"), # Random effect
            data = SA_df, family = Gamma(link = "log"), method = "REML")
# Model GS
GAM_GS <- gam(NASC_int ~ 
                # First order effects
                s(v, k = 5, bs = "tp") + # Global smoothers
                s(chl, k = 5, bs = "tp") + # Global smoothers
                s(IHO, bs = "re") + # Random intercept
                # Second order effects
                s(v, IHO, bs = "fs", k = 5) + # Group-level smoothers
                s(chl, IHO, bs = "fs", k = 5) + # Group-level smoothers
                s(IHO, year, bs = "re"), # Random effect
              data = SA_df, family = Gamma(link = "log"), method = "REML")
# Model GI
GAM_GI <- gam(NASC_int ~
                # First order effects
                s(v, k = 5, bs = "tp") + # Global smoothers
                s(chl, k = 5, bs = "tp") + # Global smoothers
                s(IHO, bs = "re") + # Random intercept
                # Secnd order effects
                s(v, by = IHO, k = 5, bs = "tp") + # Group-level smoother
                s(chl, by = IHO, k = 5, bs = "tp") + # Group-level smoother
                s(year, IHO, bs = "re"), # Random effect
               data = SA_df, family = Gamma(link = "log"), method = "REML")

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
         w_AIC = round(Weights(AIC), 10)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(dAIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-dedba048c03d5f738938" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-dedba048c03d5f738938">{"x":{"filter":"none","vertical":false,"data":[["GAM_GI","GAM_GS","GAM_G"],[31.229,31.595,20.698],[75.99,74.33,63.3],[0.15,0,0.02],[778.45,795.98,811.65],[1567.222,1577.743,1611.396],[0,10.52,44.17],[0.994834117,0.0051658828,3e-10]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


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
##     s(IHO, bs = "re") + s(IHO, year, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.7769     0.6119   7.807 2.98e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df      F p-value  
## s(v)        1.958  2.433  1.615   0.201  
## s(chl)      2.421  2.921  2.010   0.129  
## s(IHO)      2.232  4.000 15.564   0.312  
## s(IHO,year) 9.383 14.000  9.706   0.027 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0171   Deviance explained = 63.3%
## -REML = 811.65  Scale est. = 2.4581    n = 132
```


```r
summary(GAM_GS)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(v, k = 5, bs = "tp") + s(chl, k = 5, bs = "tp") + 
##     s(IHO, bs = "re") + s(v, IHO, bs = "fs", k = 5) + s(chl, 
##     IHO, bs = "fs", k = 5) + s(IHO, year, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    4.722      0.367   12.87   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                   edf Ref.df     F  p-value    
## s(v)        1.0065499   1.01 0.908 0.345074    
## s(chl)      1.0001102   1.00 3.143 0.079087 .  
## s(IHO)      0.0004368   4.00 0.000 0.540801    
## s(v,IHO)    5.4681916  23.00 1.629 0.000155 ***
## s(chl,IHO)  7.9875412  23.00 6.270 2.06e-06 ***
## s(IHO,year) 8.5840549  14.00 2.454 1.20e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.00459   Deviance explained = 74.3%
## -REML = 795.98  Scale est. = 1.3962    n = 132
```


```r
summary(GAM_GI)
```

```
## 
## Family: Gamma 
## Link function: log 
## 
## Formula:
## NASC_int ~ s(v, k = 5, bs = "tp") + s(chl, k = 5, bs = "tp") + 
##     s(IHO, bs = "re") + s(v, by = IHO, k = 5, bs = "tp") + s(chl, 
##     by = IHO, k = 5, bs = "tp") + s(year, IHO, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.7646     0.2631   18.11   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf    Ref.df     F  p-value    
## s(v)             1.000e+00 1.000e+00 0.597 0.441594    
## s(chl)           1.001e+00 1.001e+00 0.199 0.656911    
## s(IHO)           1.549e-04 4.000e+00 0.000 0.957234    
## s(v):IHOWAO_BF   1.000e+00 1.001e+00 2.481 0.118212    
## s(v):IHOCAA      8.359e-05 1.518e-04 0.000 0.500000    
## s(v):IHOBB       2.367e+00 2.837e+00 3.408 0.031009 *  
## s(v):IHODS       1.000e+00 1.000e+00 0.037 0.847099    
## s(v):IHOEAO      3.636e+00 3.897e+00 4.491 0.001933 ** 
## s(chl):IHOWAO_BF 1.000e+00 1.000e+00 5.175 0.024914 *  
## s(chl):IHOCAA    1.000e+00 1.000e+00 0.426 0.515294    
## s(chl):IHOBB     3.068e+00 3.471e+00 3.631 0.008215 ** 
## s(chl):IHODS     8.770e-05 1.597e-04 0.000 0.500000    
## s(chl):IHOEAO    3.133e+00 3.527e+00 9.469 8.24e-05 ***
## s(year,IHO)      6.852e+00 1.400e+01 1.805 0.000135 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Rank: 67/69
## R-sq.(adj) =  0.149   Deviance explained =   76%
## -REML = 778.45  Scale est. = 1.1674    n = 132
```

## Model checking

### Basis size and residual distribution

First I check the basis size k. `k-indexes` are \> 1 or close to 1 so the basis size is large enough. The residual plot look good too.


```r
appraise(GAM_GS)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGS-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_GS)
```

```
##             k'          edf   k-index p-value
## s(v)         4 1.0065499359 0.9204424  0.6325
## s(chl)       4 1.0001101984 0.9260238  0.6075
## s(IHO)       5 0.0004368472        NA      NA
## s(v,IHO)    25 5.4681915561 0.9204424  0.5900
## s(chl,IHO)  25 7.9875412221 0.9260238  0.6400
## s(IHO,year) 15 8.5840548767        NA      NA
```


```r
appraise(GAM_GI)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_GI)
```

```
##                  k'          edf   k-index p-value
## s(v)              4 1.000086e+00 0.9836785  0.8325
## s(chl)            4 1.000594e+00 0.9135958  0.5250
## s(IHO)            5 1.549404e-04        NA      NA
## s(v):IHOWAO_BF    4 1.000352e+00 0.9836785  0.8325
## s(v):IHOCAA       4 8.358556e-05 0.9836785  0.8725
## s(v):IHOBB        4 2.367328e+00 0.9836785  0.8625
## s(v):IHODS        4 1.000051e+00 0.9836785  0.8650
## s(v):IHOEAO       4 3.636028e+00 0.9836785  0.8725
## s(chl):IHOWAO_BF  4 1.000201e+00 0.9135958  0.5550
## s(chl):IHOCAA     4 1.000092e+00 0.9135958  0.5675
## s(chl):IHOBB      4 3.067986e+00 0.9135958  0.4725
## s(chl):IHODS      4 8.770214e-05 0.9135958  0.5675
## s(chl):IHOEAO     4 3.132766e+00 0.9135958  0.5250
## s(year,IHO)      15 6.851543e+00        NA      NA
```

### Resiudals against covariates


```r
GAM_GS_resid <- bind_cols(SA_df, residuals.gam(GAM_GS)) %>%
  rename(resid = "...15")
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
  rename(resid = "...15")
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

I select `model GI` because it explains the most deviance, has the best metrics and the better residuals.

## Term selection

I turn on the double penalty (`select = TRUE`) and check the covariates that has been shrunk.


```r
# Model GI
GAM_GI_p <- gam(NASC_int ~
                # First order effects
                s(v, k = 5, bs = "tp") + # Global smoothers
                s(chl, k = 5, bs = "tp") + # Global smoothers
                s(IHO, bs = "re") + # Random intercept
                # Second order effects
                s(v, by = IHO, k = 5, bs = "tp") + # Group-level smoother
                s(chl, by = IHO, k = 5, bs = "tp") + # Group-level smoother
                s(year, IHO, bs = "re"), # Random effect
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
## NASC_int ~ s(v, k = 5, bs = "tp") + s(chl, k = 5, bs = "tp") + 
##     s(IHO, bs = "re") + s(v, by = IHO, k = 5, bs = "tp") + s(chl, 
##     by = IHO, k = 5, bs = "tp") + s(year, IHO, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   4.6507     0.2288   20.33   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                        edf Ref.df      F  p-value    
## s(v)             1.961e+00      4  3.894  0.00314 ** 
## s(chl)           5.766e-05      4  0.000  0.51589    
## s(IHO)           9.007e-04      4  0.000  0.38618    
## s(v):IHOWAO_BF   1.536e-04      4  0.000  0.73977    
## s(v):IHOCAA      1.275e-04      4  0.000  0.75710    
## s(v):IHOBB       5.102e-05      4  0.000  0.69049    
## s(v):IHODS       1.309e+00      4  3.051  0.00623 ** 
## s(v):IHOEAO      1.655e+00      4  2.189  0.00790 ** 
## s(chl):IHOWAO_BF 9.486e-01      4 26.269 2.32e-05 ***
## s(chl):IHOCAA    1.551e-04      4  0.000  0.86424    
## s(chl):IHOBB     1.012e+00      4  1.174  0.12452    
## s(chl):IHODS     7.175e-05      4  0.000  0.81508    
## s(chl):IHOEAO    2.675e+00      4 34.360  < 2e-16 ***
## s(year,IHO)      8.623e+00     14  2.312 5.36e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.126   Deviance explained = 72.8%
## -REML = 793.15  Scale est. = 1.34      n = 132
```

All covariates are important.

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
##  2 WAO_BF v           16 0.801 7.55   3.03  2.03  3.13  1.72  0.431 0.918
##  3 CAA    chl         21 0.234 0.806  0.405 0.156 0.448 0.15  0.033 0.068
##  4 CAA    v           21 0.564 2.10   1.22  0.749 1.12  0.472 0.103 0.215
##  5 BB     chl         47 0.283 2.54   1.00  0.834 1.07  0.525 0.077 0.154
##  6 BB     v           47 0.167 6.36   1.81  2.61  2.30  1.61  0.235 0.474
##  7 DS     chl         25 0.319 0.685  0.413 0.111 0.435 0.094 0.019 0.039
##  8 DS     v           25 0.513 4.20   2.16  2.37  2.33  1.23  0.246 0.507
##  9 EAO    chl         23 0.089 4.06   1.36  1.06  1.22  0.898 0.187 0.388
## 10 EAO    v           23 0.525 6.00   2.21  2.49  2.84  1.73  0.361 0.749
```


```r
res = 20 # resolution for predictions
# WAO and BF
SA_WAO_BF <- data.frame(
  IHO = factor(rep(rep(rep("WAO_BF", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.801, 7.554, length.out = res), 3), # Velocity
        rep(3.032, res * 3)), # Dummy
  chl = c(rep(0.245, res * 3), # Dummy
          rep(seq(0.111, 0.857, length.out = res), 3))) # Chlorophyll
# CAA
SA_CAA <- data.frame(
  IHO = factor(rep(rep(rep("CAA", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.564, 2.097, length.out = res), 3), # Velocity
        rep(1.218, res * 3)), # Dummy
  chl = c(rep(0.405, res * 3), # Dummy
          rep(seq(0.234, 0.806, length.out = res), 3))) # Chlorophyll
# BB
SA_BB <- data.frame(
  IHO = factor(rep(rep(rep("BB", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.167, 6.357, length.out = res), 3), # Velocity
        rep(1.807, res * 3)), # Dummy
  chl = c(rep(1.005, res * 3), # Dummy
          rep(seq(0.283, 2.539, length.out = res), 3))) # Chlorophyll
# DS
SA_DS <- data.frame(
  IHO = factor(rep(rep(rep("DS", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.513, 4.205, length.out = res), 3), # Velocity
        rep(2.159, res * 3)), # Dummy
  chl = c(rep(0.413, res * 3), # Dummy
          rep(seq(0.319, 0.685, length.out = res), 3))) # Chlorophyll
# EAO
SA_EAO <- data.frame(
  IHO = factor(rep(rep(rep("EAO", res), 3), 2)),
  year = factor(rep(c(rep(2015, res), rep(2016, res), rep(2017, res)), 2)),
  var = c(rep("v", res * 3), rep("chl", res * 3)),
  v = c(rep(seq(0.525, 5.997, length.out = res), 3), # Velocity
        rep(2.209, res * 3)), # Dummy
  chl = c(rep(1.355, res * 3), # Dummy
          rep(seq(0.089, 4.063, length.out = res), 3))) # Chlorophyll
# Combine datasets
SA_new <- bind_rows(SA_WAO_BF, SA_CAA, SA_BB, SA_DS, SA_EAO)
```

Get GAM predictions for the new data. I then calculate the average functional response for all years.


```r
ilink <- family(GAM_GI)$linkinv # Get link function
pred <- predict(GAM_GI, SA_new, type = "link", se.fit = TRUE) %>% # Predict data
  bind_cols(., SA_new) %>%
  transform(lwr_ci = ilink(fit - (2 * se.fit)), # Calculate lower CI
            upr_ci = ilink(fit + (2 * se.fit)), # Calculate upper CI
            fitted = ilink(fit)) %>% # Calculate fit
  group_by(IHO, var, v, chl) %>%
  summarise(NASC_fit = mean(fitted), # Calculate average functional response
            NASC_lwr_ci = mean(lwr_ci),
            NASC_upr_ci = mean(upr_ci)) %>%
  mutate(SA_fit = 10 * log10(NASC_fit), # Convert to SA
         SA_lwr_ci = 10 * log10(NASC_lwr_ci),
         SA_upr_ci = 10 * log10(NASC_upr_ci))
# Save data
save(pred, GAM_GI, file = "data/statistics/GAM_results.RData")
```

## Model visualization

Quick partial effects plot with `gratia::draw`.


```r
draw(GAM_GS, residuals = TRUE, scales = "free")
```

![](PanArctic_DSL_statistics_files/figure-html/draw-GAMGI-1.png)<!-- -->

Plot to see if predictions worked well.


```r
plot_grid(
  pred %>%
    filter(var == "v") %>%
    ggplot() +
    geom_line(aes(x = v, y = NASC_fit)) +
    geom_ribbon(aes(x = v, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = v, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free"),
  pred %>%
    filter(var == "chl") %>%
    ggplot() +
    geom_line(aes(x = chl, y = NASC_fit)) +
    geom_ribbon(aes(x = chl, ymin = NASC_lwr_ci, ymax = NASC_upr_ci), alpha = 0.1) +
    geom_point(data = SA_df, aes(x = chl, y = NASC_int, group = IHO)) + 
    facet_wrap(~ IHO, scales = "free"),
  ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGI-ggplot-1.png)<!-- -->

Pretty plot.


```r
pred_plot <- pred %>%
  mutate(IHO = factor(
    case_when(IHO == "EAO" ~ "East Arctic\nOcean",
              IHO == "WAO_BF" ~ "West Arctic Ocean\n&Beaufort Sea",
              IHO == "CAA" ~ "Can. Arctic\nArchipelago",
              IHO == "BB" ~ "Baffin Bay",
              IHO == "DS" ~ "Davis Strait"),
    levels = c("West Arctic Ocean\n&Beaufort Sea", "Can. Arctic\nArchipelago",
               "Baffin Bay", "Davis Strait", "East Arctic\nOcean")))
GAM_velo <- pred_plot %>%
  filter(var == "v") %>%
  ggplot() +
  geom_line(aes(x = v, y = SA_fit, col = IHO), size = 0.8) +
  geom_ribbon(aes(x = v, ymin = SA_lwr_ci, ymax = SA_upr_ci, fill = IHO), alpha = 0.1) +
  scale_colour_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  scale_fill_manual(values = met.brewer("Johnson", n = 6, direction = -1)) +
  coord_cartesian(ylim = c(-2, 50), expand = T) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(x = expression("Current velocity at 222 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
GAM_chl <- pred_plot %>%
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
plot_grid(GAM_velo, GAM_chl, ncol = 1, rel_heights = c(1,0.8))
```

![](PanArctic_DSL_statistics_files/figure-html/GAMGI-pretty-ggplot-1.png)<!-- -->
