PanArctic DSL - Statistics
================
[Pierre Priou](mailto:pierre.priou@mi.mun.ca)
2022/05/03 at 15:13

# Package loading

``` r
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

I want to test whether temperature and salinity at mesopelagic depth,
sea-ice concentration, open-water duration (a proxy for productivity)
have an effect on the backscatter anomalies observed per year. I
therefore combined gridded acoustic data—integrated mesopelagic
NASC—with gridded CTD, and remote sensing data projected on either the
WGS84 or the EASE-Grid 2.0 North.

``` r
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
SA_grid_laea <- SA_grid_laea %>%
  mutate(SA_int = 10 * log10(NASC_int))
load("data/remote_sensing/physics_grids.RData") # Modelled physics data 
load("data/remote_sensing/seaice_grids.RData") # Remote sensing sea ice data
```

# Data preparation

I combine data using the EASE-Grid 2.0 North (EPSG:6931) and normalize
the covariates per IHO area. I combined the Beaufort Sea and West Arctic
Ocean data to have more points for the regression due to their
geographical proximity.

``` r
SA_laea <- SA_grid_laea %>%  # Tidy anomaly dataset for joining
  dplyr::select(-lat, -lon)
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
  dplyr::select(-mean_NASC_area_year, -sd_NASC_area_year, -NASC_anomaly_d, -cell_res, -vxo, -vyo,
                -mean_temp_area_depth, -mean_velo_area_depth, -xc_na, -yc_na, -year_na, -siconc, -sithick,
                -temp_anomaly, -velocity_anomaly)
```

# Data exploration

Maps of all variables.

``` r
stat_laea %>% 
  ggplot(aes(x = xc,  y = yc)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_point(aes(col = IHO_area)) +
  ggtitle("regions") +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

![](PanArctic_DSL_statistics_files/figure-gfm/map-SA-int-1.png)<!-- -->

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/map-siconc-1.png)<!-- -->

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/map-ow-1.png)<!-- -->

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/map-ice-day-1.png)<!-- -->

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/map-ice-week-1.png)<!-- -->

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/map-thetao-1.png)<!-- -->

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/map-so-1.png)<!-- -->

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/map-velocity-1.png)<!-- -->

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/map-mld-1.png)<!-- -->

# Spatial and interannual variability

``` r
SA_diff <- stat_laea %>%
  dplyr::select(year, IHO_area, SA_int, NASC_int, CM) %>%
  mutate(year = factor(year))
```

## All areas

I check whether there are inter-annual and spatial variability in
S<sub>A</sub> and the centre of mass across years and areas.

``` r
# Visualize data
plot_grid(SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = SA_int)),
          SA_diff %>% 
            ggplot() +
            geom_boxplot(aes(x = year, y = CM)))
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-1.png)<!-- -->

``` r
# Check assumptions
## Outliers
SA_diff %>%
  group_by(year) %>%
  identify_outliers(SA_int) 
```

    ## # A tibble: 2 x 7
    ##   year  IHO_area SA_int NASC_int    CM is.outlier is.extreme
    ##   <fct> <fct>     <dbl>    <dbl> <dbl> <lgl>      <lgl>     
    ## 1 2017  EAO        41.6   14296.  412. TRUE       FALSE     
    ## 2 2017  EAO        48.3   67316.  392. TRUE       FALSE

``` r
# There are no extreme outlier in CAA
SA_diff %>%
  group_by(year) %>%
  identify_outliers(CM) 
```

    ## # A tibble: 5 x 7
    ##   year  IHO_area SA_int NASC_int    CM is.outlier is.extreme
    ##   <fct> <fct>     <dbl>    <dbl> <dbl> <lgl>      <lgl>     
    ## 1 2015  WAO_BF    -3.85    0.412  601. TRUE       FALSE     
    ## 2 2015  WAO_BF     8.88    7.72   591. TRUE       FALSE     
    ## 3 2015  WAO_BF    -2.63    0.545  564. TRUE       FALSE     
    ## 4 2015  WAO_BF     7.62    5.78   572. TRUE       FALSE     
    ## 5 2015  WAO_BF    -1.13    0.770  581. TRUE       FALSE

``` r
# There are no extreme outliers
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(SA_int) 
```

    ## # A tibble: 6 x 7
    ##   IHO_area year  SA_int NASC_int    CM is.outlier is.extreme
    ##   <fct>    <fct>  <dbl>    <dbl> <dbl> <lgl>      <lgl>     
    ## 1 BB       2016    1.98     1.58  295. TRUE       FALSE     
    ## 2 BB       2016    3.58     2.28  354. TRUE       FALSE     
    ## 3 BB       2017    4.26     2.67  302. TRUE       FALSE     
    ## 4 DS       2016   30.4   1101.    367. TRUE       FALSE     
    ## 5 EAO      2017   41.6  14296.    412. TRUE       FALSE     
    ## 6 EAO      2017   48.3  67316.    392. TRUE       FALSE

``` r
# There are no extreme outlier in CAA
SA_diff %>%
  group_by(IHO_area) %>%
  identify_outliers(CM) 
```

    ## # A tibble: 5 x 7
    ##   IHO_area year  SA_int NASC_int    CM is.outlier is.extreme
    ##   <fct>    <fct>  <dbl>    <dbl> <dbl> <lgl>      <lgl>     
    ## 1 BB       2015    14.6     28.8  415. TRUE       FALSE     
    ## 2 BB       2016    13.9     24.7  393. TRUE       FALSE     
    ## 3 BB       2016    16.8     48.0  427. TRUE       FALSE     
    ## 4 BB       2016    16.8     47.4  401. TRUE       FALSE     
    ## 5 BB       2017    13.2     20.8  265. TRUE       FALSE

``` r
# There are no extreme outliers

## Normality: Checked by looking at the anova residuals with QQ plot
ggqqplot(SA_diff, "SA_int", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-2.png)<!-- -->

``` r
ggqqplot(SA_diff, "CM", facet.by = "year")
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-3.png)<!-- -->

``` r
ggqqplot(SA_diff, "SA_int", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-4.png)<!-- -->

``` r
ggqqplot(SA_diff, "CM", facet.by = "IHO_area")
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-5.png)<!-- -->

``` r
# Data is fully normal

## Homogeneity of variance
plot(lm(SA_int ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-6.png)<!-- -->

``` r
plot(lm(CM ~ year, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-7.png)<!-- -->

``` r
plot(lm(SA_int ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-8.png)<!-- -->

``` r
plot(lm(CM ~ IHO_area, data = SA_diff), 1)
```

![](PanArctic_DSL_statistics_files/figure-gfm/all-interannual-spatial-prep-9.png)<!-- -->

``` r
SA_diff %>%
  levene_test(SA_int ~ year)
```

    ## # A tibble: 1 x 4
    ##     df1   df2 statistic     p
    ##   <int> <int>     <dbl> <dbl>
    ## 1     2   129     0.333 0.718

``` r
# Not homogenous
```

I used a non parametric Kruskal-Wallis test to test for interannual and
spatial variability because the variance is not homogenous and not very
normal.

``` r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(SA_int ~ year)
```

    ## # A tibble: 1 x 6
    ##   .y.        n statistic    df      p method        
    ## * <chr>  <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 SA_int   132      8.65     2 0.0133 Kruskal-Wallis

``` r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ year)
```

    ## # A tibble: 1 x 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 CM      132      1.65     2 0.438 Kruskal-Wallis

``` r
SA_diff %>% # Kruskal wallis interannual difference SA
  kruskal_test(SA_int ~ IHO_area)
```

    ## # A tibble: 1 x 6
    ##   .y.        n statistic    df         p method        
    ## * <chr>  <int>     <dbl> <int>     <dbl> <chr>         
    ## 1 SA_int   132      26.2     4 0.0000289 Kruskal-Wallis

``` r
SA_diff %>% # Kruskal wallis interannual difference CM
  kruskal_test(CM ~ IHO_area)
```

    ## # A tibble: 1 x 6
    ##   .y.       n statistic    df        p method        
    ## * <chr> <int>     <dbl> <int>    <dbl> <chr>         
    ## 1 CM      132      57.8     4 8.21e-12 Kruskal-Wallis

The center of mass did not vary significantly between years (H =
1.652899, p = 0.438) but varied significantly among areas (*H* =
57.8495, *p* &lt; 0.001). Mesopelagic backscatter S<sub>A</sub> varied
significantly between areas (*H* = 26.19394, *p* &lt; 0.001) but not
years (*H* = 8.6464, *p* = 0.013). Thus, in the S<sub>A</sub> HGAM a
random effect for change in S<sub>A</sub> per area is likely needed.

## Within areas

I check whether there are inter-annual variability in S<sub>A</sub>
across years within areas. This will be used to decide whether a random
effect is needed in the GAM. I used a non parametric Kruskal-Wallis test
to test for interannual variability within areas.

``` r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(IHO_area) %>%
  kruskal_test(SA_int ~ year)
```

    ## # A tibble: 5 x 7
    ##   IHO_area .y.        n statistic    df      p method        
    ## * <fct>    <chr>  <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 WAO_BF   SA_int    16     4.84      2 0.0889 Kruskal-Wallis
    ## 2 CAA      SA_int    18     6.93      2 0.0313 Kruskal-Wallis
    ## 3 BB       SA_int    46     4.67      2 0.0968 Kruskal-Wallis
    ## 4 DS       SA_int    24     0.108     2 0.947  Kruskal-Wallis
    ## 5 EAO      SA_int    28     5.03      2 0.081  Kruskal-Wallis

There was no interannual differences in SA within each region (WAO\_BF -
H = 4.8395483, *p* = 0.0889; CAA - H = 6.9278752, *p* = 0.0313; BB - H =
0.0968, *p* = 0.0889; DS - H = 0.9470, *p* = 0.0889; EAO - H =
5.0263899, *p* = 0.0810). Therefore, there seems to be no need for a
random effect for interannual variability within years.

# HGAM

## Data preparation

``` r
SA_df <- stat_laea %>%
  dplyr::select(year, xc, yc, area, IHO_area_large, IHO_area, SA_int, velocity, thetao, openwater_duration, ice_week) %>%
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

``` r
plot_grid(SA_df %>%
            ggplot(aes(x = v, y = SA_int, col = IHO_area_large)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = t, y = SA_int, col = IHO_area_large)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = o, y = SA_int, col = IHO_area_large)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = s, y = SA_int, col = IHO_area_large)) +
            geom_point() +
            geom_smooth(method = "lm", se = F, col = "grey20") +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-gfm/plot-data-area-1.png)<!-- -->

``` r
plot_grid(SA_df %>%
            ggplot(aes(x = v, fill = IHO_area_large)) +
            geom_histogram() +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = t, fill = IHO_area_large)) +
            geom_histogram() +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = o, fill = IHO_area_large)) +
            geom_histogram() +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          SA_df %>%
            ggplot(aes(x = s, fill = IHO_area_large)) +
            geom_histogram() +
            facet_grid(~ IHO) +
            theme(legend.position = "none"),
          ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-gfm/plot-distribution-1.png)<!-- -->

Check correlations.

``` r
corr_EAO <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "EAO") %>%
  dplyr::select(-year, -xc, -yc, -IHO_area_large, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_EAO, type = "lower", lab = T, title = "corr EAO")
```

![](PanArctic_DSL_statistics_files/figure-gfm/corr-area-1.png)<!-- -->

``` r
corr_BB <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "BB") %>%
  dplyr::select(-year, -xc, -yc, -IHO_area_large, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_BB, type = "lower", lab = T, title = "corr BB")
```

![](PanArctic_DSL_statistics_files/figure-gfm/corr-area-2.png)<!-- -->

``` r
corr_DS <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "DS") %>%
  dplyr::select(-year, -xc, -yc, -IHO_area_large, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_DS, type = "lower", lab = T, title = "corr BB")
```

![](PanArctic_DSL_statistics_files/figure-gfm/corr-area-3.png)<!-- -->

``` r
corr_WAO <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "WAO_BF") %>%
  dplyr::select(-year, -xc, -yc, -IHO_area_large, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_WAO, type = "lower", lab = T, title = "corr WAO")
```

![](PanArctic_DSL_statistics_files/figure-gfm/corr-area-4.png)<!-- -->

``` r
corr_CAA <- SA_df %>% # Compute Spearman correlation matrix
  filter(IHO == "CAA") %>%
  dplyr::select(-year, -xc, -yc, -IHO_area_large, -IHO, -area, -SA_int_n, -v_n, -t_n, -o_n, -s_n) %>%
  cor(., method = "spearman") %>%
  round(., 2)
ggcorrplot(corr_CAA, type = "lower", lab = T, title = "corr CAA")
```

![](PanArctic_DSL_statistics_files/figure-gfm/corr-area-5.png)<!-- -->

## Model fitting

I decided to fit hierarchical generalized additive model due to their
flexibility in modelling non linear relationships. I fit several models
with different structures (random intercept, random “slope”) following
Pedersen et al. 2019.

``` r
# Model G
GAM1 <- gam(SA_int ~ s(v, k = 5, bs = "tp") + s(t, k = 5, bs = "tp") + s(o, k = 5, bs = "tp"),
            data = SA_df, family = "gaussian", method = "REML")
# Model S
GAM2 <- gam(SA_int ~ s(v, IHO, bs = "fs", k = 5) + s(t, IHO, bs = "fs", k = 5) + s(o, IHO, bs = "fs", k = 5), 
            data = SA_df, family = "gaussian", method = "REML")
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
    ## d smooths of same variable.

``` r
# Model I
GAM3 <- gam(SA_int ~ s(v, by = IHO, k = 5, bs = "tp") + s(t, by = IHO, k = 5, bs = "tp") + 
              s(o, by = IHO, k = 5, bs = "tp") + s(IHO, bs = "re"),
            data = SA_df, family = "gaussian", method = "REML")
# Model GS
GAM4 <- gam(SA_int ~ s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + s(v, IHO, bs = "fs", k = 5) + 
              s(t_n, IHO, bs = "fs", k = 5) + s(o_n, IHO, bs = "fs", k = 5),
            data = SA_df, family = "gaussian", method = "REML")
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has repeated 1-
    ## d smooths of same variable.

``` r
# Model GI
GAM5 <- gam(SA_int ~ s(v, k = 5) + s(t, k = 5) + s(o, k = 5) + s(v, by = IHO, k = 5, bs = "tp") + 
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

![](PanArctic_DSL_statistics_files/figure-gfm/HGAM-fit-1.png)<!-- -->

`GAM3` and `GAM5` are the two models with similar performances, although
`GAM3` seems to perform slightly better.

``` r
summary(GAM3)
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
    ## (Intercept)   18.562      1.741   10.66   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                   edf Ref.df      F  p-value    
    ## s(v):IHOWAO_BF 1.0000  1.000  0.485  0.48770    
    ## s(v):IHOCAA    1.0000  1.000  0.243  0.62327    
    ## s(v):IHOBB     1.0004  1.001  6.195  0.01430 *  
    ## s(v):IHODS     1.0000  1.000  1.302  0.25629    
    ## s(v):IHOEAO    1.0001  1.000  0.755  0.38671    
    ## s(t):IHOWAO_BF 1.7901  1.961  5.536  0.01183 *  
    ## s(t):IHOCAA    1.0000  1.000  0.265  0.60790    
    ## s(t):IHOBB     3.0982  3.524  3.188  0.01363 *  
    ## s(t):IHODS     1.0000  1.000  0.007  0.93573    
    ## s(t):IHOEAO    1.0000  1.000  7.031  0.00919 ** 
    ## s(o):IHOWAO_BF 1.0002  1.000  3.896  0.05087 .  
    ## s(o):IHOCAA    2.5439  2.925  4.709  0.00245 ** 
    ## s(o):IHOBB     1.0000  1.000 18.966 3.04e-05 ***
    ## s(o):IHODS     1.0001  1.000  0.045  0.83216    
    ## s(o):IHOEAO    1.2179  1.369  5.655  0.04707 *  
    ## s(IHO)         0.9169  4.000  0.497  0.09668 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =   0.44   Deviance explained = 52.8%
    ## -REML = 409.83  Scale est. = 42.286    n = 132

## Model checking

First I check the basis size k. `k-indexes` are &gt; 1 or close to 1 so
the basis size is large enough. The residual plot look good too.

``` r
par(mfrow = c(2, 2))
gam.check(GAM3, rep = 500)
```

![](PanArctic_DSL_statistics_files/figure-gfm/basis-size-residuals-1.png)<!-- -->

    ## 
    ## Method: REML   Optimizer: outer newton
    ## full convergence after 13 iterations.
    ## Gradient range [-0.0001038437,0.0003489404]
    ## (score 409.8302 & scale 42.2863).
    ## Hessian positive definite, eigenvalue range [5.778968e-06,58.03567].
    ## Model rank =  66 / 66 
    ## 
    ## Basis dimension (k) checking results. Low p-value (k-index<1) may
    ## indicate that k is too low, especially if edf is close to k'.
    ## 
    ##                   k'   edf k-index p-value
    ## s(v):IHOWAO_BF 4.000 1.000    1.06    0.72
    ## s(v):IHOCAA    4.000 1.000    1.06    0.73
    ## s(v):IHOBB     4.000 1.000    1.06    0.76
    ## s(v):IHODS     4.000 1.000    1.06    0.72
    ## s(v):IHOEAO    4.000 1.000    1.06    0.75
    ## s(t):IHOWAO_BF 4.000 1.790    1.12    0.88
    ## s(t):IHOCAA    4.000 1.000    1.12    0.89
    ## s(t):IHOBB     4.000 3.098    1.12    0.94
    ## s(t):IHODS     4.000 1.000    1.12    0.92
    ## s(t):IHOEAO    4.000 1.000    1.12    0.88
    ## s(o):IHOWAO_BF 4.000 1.000    1.02    0.51
    ## s(o):IHOCAA    4.000 2.544    1.02    0.55
    ## s(o):IHOBB     4.000 1.000    1.02    0.58
    ## s(o):IHODS     4.000 1.000    1.02    0.59
    ## s(o):IHOEAO    4.000 1.218    1.02    0.54
    ## s(IHO)         5.000 0.917      NA      NA

Check residuals versus covariates.

``` r
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

![](PanArctic_DSL_statistics_files/figure-gfm/residuals-covariates-1.png)<!-- -->

Check auto-correlation wihtin the model. All seems in order.

``` r
par(mfrow = c(1, 2), mar = c(4, 4, 4, 0.5))
acf(resid(GAM3), lag.max = 36, main = "ACF")
pacf(resid(GAM3), lag.max = 36, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-gfm/acf-pacf-1.png)<!-- -->

## Model selection

I turn on the double penalty (`select = TRUE`) and check the covariates
that has been shrunk. All covariates are imporant in the model.

``` r
GAM3_p <- gam(SA_int ~ s(v, by = IHO, k = 5, bs = "tp") + s(t, by = IHO, k = 5, bs = "tp") + 
                s(o, by = IHO, k = 5, bs = "tp") + s(IHO, bs = "re"),
              data = SA_df, family = "gaussian", method = "REML", select = TRUE)
summary(GAM3_p)
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
    ## (Intercept)  17.8786     0.6966   25.67   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                      edf Ref.df     F  p-value    
    ## s(v):IHOWAO_BF 4.376e-05      4 0.000 0.895751    
    ## s(v):IHOCAA    1.086e-04      4 0.000 0.419790    
    ## s(v):IHOBB     1.700e+00      4 2.283 0.003999 ** 
    ## s(v):IHODS     3.321e-01      4 0.105 0.266644    
    ## s(v):IHOEAO    1.377e-01      4 0.040 0.265774    
    ## s(t):IHOWAO_BF 1.157e+00      4 4.612 1.96e-05 ***
    ## s(t):IHOCAA    3.546e-05      4 0.000 0.462514    
    ## s(t):IHOBB     6.828e-01      4 0.538 0.074136 .  
    ## s(t):IHODS     1.032e-05      4 0.000 0.461164    
    ## s(t):IHOEAO    7.918e-01      4 0.950 0.027320 *  
    ## s(o):IHOWAO_BF 7.902e-01      4 0.942 0.030198 *  
    ## s(o):IHOCAA    1.879e+00      4 3.465 0.000585 ***
    ## s(o):IHOBB     9.027e-01      4 2.316 0.001407 ** 
    ## s(o):IHODS     2.100e-05      4 0.000 0.589685    
    ## s(o):IHOEAO    9.237e-01      4 3.028 0.000176 ***
    ## s(IHO)         1.410e-05      4 0.000 0.490977    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.406   Deviance explained = 44.8%
    ## -REML = 448.81  Scale est. = 44.824    n = 132

``` r
draw(GAM3_p, scales = "free")
```

![](PanArctic_DSL_statistics_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Model visualization

Quick plots using `visreg`.

``` r
visreg(GAM3, "v", "IHO", overlay = F, scale = "response", band = F, partial = T)
```

![](PanArctic_DSL_statistics_files/figure-gfm/plot-GAM3-1.png)<!-- -->

``` r
visreg(GAM3, "t", "IHO", overlay = F, scale = "response", band = F, partial = T)
```

![](PanArctic_DSL_statistics_files/figure-gfm/plot-GAM3-2.png)<!-- -->

``` r
visreg(GAM3, "o", "IHO", overlay = F, scale = "response", band = F, partial = T)
```

![](PanArctic_DSL_statistics_files/figure-gfm/plot-GAM3-3.png)<!-- -->

I predict the data and then use that dataframe to plot it in a “clean”
way with `ggplot`.

``` r
pred_v <- get_gam_predictions(GAM3, series = v, .comparison = IHO, series_length = 50) %>% 
  mutate(var = "v",
         signif = case_when(IHO == "WAO_BF" ~  0.05, # Significance level
                            IHO == "CAA" ~ 1,
                            IHO == "BB" ~ 0.01,
                            IHO == "DS" ~ 0.05,
                            IHO == "EAO" ~ 1),
         signif_group = factor(if_else(signif > 0.1, F, T))) %>%
  rename(idx = ".idx",
         fit = SA_int,
         value = v)
pred_t <- get_gam_predictions(GAM3, series = t, .comparison = IHO, series_length = 50) %>%
  mutate(var = "t",
         signif = case_when(IHO == "WAO_BF" ~ 0.001, # significance level
                            IHO == "CAA" ~ 0.1,
                            IHO == "BB" ~ 0.001,
                            IHO == "DS" ~ 1,
                            IHO == "EAO" ~ 0.1),
         signif_group = factor(if_else(signif > 0.1, F, T))) %>%
    rename(idx = ".idx",
           fit = SA_int,
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
         fit = SA_int,
         value = o)
# Combine data
pred_GAM <- bind_rows(pred_v, pred_t, pred_o)

# Save data
save(pred_GAM, GAM3, file = "data/statistics/GAM_results.RData")
```

`ggplot` plot.

``` r
col_pal <- c("#5BBCD6", "#00A08A", "#F2AD00", "#F98400", "#FF0000")# wesanderson::wes_palette("Darjeeling1", type = "discrete") 

plot_grid(pred_GAM %>%
            filter(var == "v" & SE < 5) %>%
            ggplot() +
            geom_line(aes(x = value, y = fit, col = IHO, linetype = IHO), size = 0.75, alpha = 0.8) +
            # geom_ribbon(aes(x = value, ymin = CI_lower, ymax = CI_upper, fill = IHO), alpha = 0.1) +
            scale_colour_manual(name = "IHO", values = col_pal) + 
            scale_fill_manual(name = "IHO", values = col_pal) + 
            scale_linetype_manual(name = "IHO", values = c(2, 2, 1, 2, 2)) +
            labs(x = "Velo. 380 m", y = expression("S"[A]*"")),
          pred_GAM %>%
            filter(var == "t" & SE < 5) %>%
            ggplot() +
            geom_line(aes(x = value, y = fit, col = IHO, linetype = IHO), size = 0.75, alpha = 0.8) +
            # geom_ribbon(aes(x = value, ymin = CI_lower, ymax = CI_upper, fill = IHO), alpha = 0.1) +
            scale_colour_manual(name = "IHO", values = col_pal) + 
            scale_fill_manual(name = "IHO", values = col_pal) + 
            scale_linetype_manual(name = "IHO", values = c(1, 2, 1, 2, 1)) +
            labs(x = "Temp. 380 m", y = expression("S"[A]*"")),
          pred_GAM %>%
            filter(var == "o" & SE < 5) %>%
            ggplot() +
            geom_line(aes(x = value, y = fit, col = IHO, linetype = IHO), size = 0.75, alpha = 0.8) +
            # geom_ribbon(aes(x = value, ymin = CI_lower, ymax = CI_upper, fill = IHO), alpha = 0.1) +
            scale_colour_manual(name = "IHO", values = col_pal) + 
            scale_fill_manual(name = "IHO", values = col_pal) + 
            scale_linetype_manual(name = "IHO", values = c(1, 1, 1, 2, 1)) +
            labs(x = "Open water days", y = expression("S"[A]*"")))
```

![](PanArctic_DSL_statistics_files/figure-gfm/ggplot-gam-1.png)<!-- -->
