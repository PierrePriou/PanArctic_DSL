---
title: "PanArctic DSL - Remote sensing"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/24 at 15:40"
output: 
  html_document:
    keep_md: yes
  github_document:
    always_allow_html: true
---

# Package loading


```r
# Load packages
library(tidyverse)  # Tidy code
library(dtplyr)     # Speed up code
library(lubridate)  # Deal with dates
library(tidync)     # Read NetCDF
library(ncmeta)     # Metadata NetCDF
library(sf)         # Spatial data
library(raster)     # Rasterize data
library(rgdal)      # Read shapefiles
library(ecmwfr)     # Download Copernicus data
library(cowplot)    # Plots on a grid
library(cmocean)    # Oceanographic colour palettes
library(metR)
# Custom figure theme
theme_set(theme_bw())
theme_update(axis.text = element_text(size = 9),
             axis.title = element_text(size = 9),
             strip.text.x = element_text(size = 9, face = "plain", hjust = 0.5),
             strip.background = element_rect(colour = "transparent",
                                             fill = "transparent"),
             legend.title = element_text(size = 9, vjust = 1),
             legend.margin = margin(0, 0, 0, 0),
             legend.box.margin = margin(0, 0, -8, 0),
             panel.grid = element_blank(), 
             plot.margin = unit(c(0.02, 0.02, 0.02, 0.02), "in"),
             plot.title = element_text(size = 9, face = "bold"))
# Suppress summarise() warning
options(dplyr.summarise.inform = F) 
```


```r
# Laea projection
cell_res <- 12.5
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931")
projection(arctic_laea) <- gsub("units=m", "units=km", projection(arctic_laea))
res(arctic_laea) <- c(cell_res, cell_res)
```

# Sea ice data

Code that downloads sea ice data from the [Copernicus Climate website](https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-sea-ice-concentration?tab=overview). We use the daily sea ice concentration derived from satellite observations from ECMWF. Data is downloaded in two batches because ECMWF changed their Climate Data Record (CDR) on 2016-01-01. I changed the sea ice data projection from Lambert's azimuthal Equal Area (EPSG:6931) to WGS84 (EPSG:4326) and match the resolution to that of backscatter anomalies (2 deg longitude \* 1 deg latitude).

## CDR data 2015-01-01 - 2015-12-31

Download "raw" data that I only crop to match the area of interest to reduce the file size.


```r
# Date range of interest
dates <- seq(ymd("2015-01-01"), ymd("2015-12-31"), 1)
# Download data
for (d in dates) {
  # Format years, month, and day to be correctly read in the request
  year = as.character(year(as.Date(d, origin = "1970-01-01")))
  month = 
    if_else(as.numeric(month(as.Date(d, origin =" 1970-01-01"))) < 10,
            paste0("0", as.character(month(as.Date(d, origin = "1970-01-01")))),
            as.character(month(as.Date(d, origin = "1970-01-01"))))
  day =
    if_else(as.numeric(day(as.Date(d, origin = "1970-01-01"))) < 10,
            paste0("0", as.character(day(as.Date(d, origin = "1970-01-01")))),
            as.character(day(as.Date(d, origin = "1970-01-01"))))
  # Create request
  request <- list(version = "v2", 
                  variable = "all",
                  format = "zip",
                  origin = "eumetsat_osi_saf",
                  region = "northern_hemisphere",
                  cdr_type = "cdr",
                  year = year,
                  month = month,
                  day = day,
                  dataset_short_name = "satellite-sea-ice-concentration",
                  target = "download.zip")
  print(paste0("Downloading ", request$dataset_short_name, 
               " from ", year, "-", month, "-", day))
  # Download requested file
  tmp_raw <- wf_request(user = "113650",
                        request = request, 
                        transfer = T, 
                        verbose = F,
                        path = "data/remote_sensing/sea_ice/tmp") %>%
    # Unzip file
    unzip(exdir = "data/remote_sensing/sea_ice/tmp") %>%  
    tidync() %>% 
    activate("D2,D3") %>%
    # Crop data to fit area of interest
    hyper_filter(xc = xc >= -2700 & xc <= 2700, yc = yc >= -2700 & yc <= 2700) 
  # Extract latlon and laea coordinates
  tmp_coord <- tmp_raw %>% 
    activate("D2,D3") %>%
    hyper_tibble()
  tmp_seaice <- tmp_raw %>%
    activate("D2,D3,D0") %>%
    hyper_tibble() %>%
    # Join latlon coordinates to sea ice data
    left_join(., tmp_coord, by = c("xc", "yc")) %>%
    mutate(date = ymd(format(as_datetime(time, origin = "1978-01-01", tz = "UTC"),
                             format = "%Y%m%d"))) 
  # Save as a csv file
  write_csv(tmp_seaice, file = paste0("data/remote_sensing/sea_ice/",
                                      year, month, day, "_laea_ice_conc.csv")) 
  # Delete temporary netcdf and zip files
  file.remove(list.files("data/remote_sensing/sea_ice/tmp", full.names = T)) 
}
# Remove unused variables
rm(request, tmp_raw, tmp_coord, tmp_seaice)
```

## ICDR data 2016-01-01 - 2017-12-31


```r
# Date range of interest
dates <- seq(ymd("2016-01-10"), ymd("2017-12-31"), 1)
# Download data
for (d in dates) {
  # Format years, month, and day to be correctly read in the request
  year = as.character(year(as.Date(d, origin = "1970-01-01")))
  month = 
    if_else(as.numeric(month(as.Date(d, origin =" 1970-01-01"))) < 10,
            paste0("0", as.character(month(as.Date(d, origin = "1970-01-01")))),
            as.character(month(as.Date(d, origin = "1970-01-01"))))
  day = 
    if_else(as.numeric(day(as.Date(d, origin = "1970-01-01"))) < 10,
            paste0("0", as.character(day(as.Date(d, origin = "1970-01-01")))),
            as.character(day(as.Date(d, origin = "1970-01-01"))))
  # Create request
  request <- list(version = "v2",
                  variable = "all",
                  format = "zip",
                  origin = "eumetsat_osi_saf",
                  region = "northern_hemisphere",
                  cdr_type = "icdr",
                  year = year,
                  month = month,
                  day = day,
                  dataset_short_name = "satellite-sea-ice-concentration",
                  target = "download.zip")
  print(paste0("Downloading ", request$dataset_short_name, 
               " from ", year, "-", month, "-", day))
  # Download requested file
  tmp_raw <- wf_request(user = "113650", 
                        request = request, 
                        transfer = T, 
                        verbose = F,
                        path = "data/remote_sensing/sea_ice/tmp") %>%
    # Unzip file
    unzip(exdir = "data/remote_sensing/sea_ice/tmp") %>%
    tidync() %>% 
    activate("D2,D3") %>%
    # Crop data to fit area of interest
    hyper_filter(xc = xc >= -2700 & xc <= 2700, yc = yc >= -2700 & yc <= 2700) 
  # Extract latlon and laea coordinates
  tmp_coord <- tmp_raw %>% 
    activate("D2,D3") %>%
    hyper_tibble()
  tmp_seaice <- tmp_raw %>%
    activate("D2,D3,D0") %>%
    hyper_tibble() %>%
    # Join latlon coordinates to sea ice data
    left_join(., tmp_coord, by = c("xc", "yc")) %>% 
    mutate(date = ymd(format(as_datetime(time, origin = "1978-01-01", tz = "UTC"),
                             format = "%Y%m%d"))) 
  # Save as a csv file
  write_csv(tmp_seaice, file = paste0("data/remote_sensing/sea_ice/", 
                                      year, month, day, "_laea_ice_conc.csv")) 
  # Delete temporary netcdf and zip files
  file.remove(list.files("data/remote_sensing/sea_ice/tmp", full.names = T))
}
# Remove unused variables
rm(request, tmp_raw, tmp_coord, tmp_seaice)
```

## Data tidying and gridding

First, I plot data to see if download worked.


```r
plot_grid(
  read_csv("data/remote_sensing/sea_ice/20151231_laea_ice_conc.csv",
           show_col_types = F) %>%
    ggplot(aes(x = xc, y = yc, fill = ice_conc)) + 
    geom_tile() + 
    scale_fill_cmocean("Ice concentration (%)", name = "ice") +
    ggtitle("2015-02-12 - CDR") + 
    coord_fixed(expand = F) +
    theme(legend.position = "top", 
          legend.key.height = unit(0.1, "in"),
          legend.key.width = unit(0.3, "in")),
  read_csv("data/remote_sensing/sea_ice/20170212_laea_ice_conc.csv",
           show_col_types = F) %>%
    ggplot(aes(x = xc, y = yc, fill = ice_conc)) + 
    geom_tile() + 
    scale_fill_cmocean("Ice concentration (%)", name = "ice") +
    ggtitle("2017-02-12 - ICDR") + 
    coord_fixed(expand = F) +
    theme(legend.position = "top",
          legend.key.height = unit(0.1, "in"), 
          legend.key.width = unit(0.3, "in")),
  ncol = 2, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_remote_sensing_files/figure-html/plot-seaice-latlon-cdr-1.png)<!-- -->

I calculate sea ice and open-water duration, and mean sea ice concentration for each cell per year. We define a cell being ice-covered when the mean sea ice concentration is above 15 % [(as per NSIDC Sea Ice Index)](https://nsidc.org/cryosphere/glossary/term/ice-extent).


```r
# Empty dataframe 
seaice_year <- data.frame() 

# Loop through data per year
for (i in seq(2015, 2017, 1)) { 
  # List files in folder
  seaice_tmp <- list.files("data/remote_sensing/sea_ice", 
                            pattern = paste0(i),
                            full.names = T) %>%
    # Read sea ice files
    map_dfr(.f = ~ read_csv(., show_col_types = F)) %>% 
    # Call dtplyr which considerably speed up calculations
    lazy_dt() %>% 
    mutate(year = i,
           # calculate day of the year
           yday = yday(date),
           # calculate week of the year
           week = week(date), 
           # > 15% ice cover = ice covered day
           ice_covered_day = if_else(ice_conc > 15, 1, 0)) %>% 
    # Remove date to speed up calculations
    dplyr::select(-date) %>%
    # Remove data from land (1), lakes (2), and possible false ice (16)
    filter(status_flag != c(1, 2, 16) & is.na(ice_conc) == F) %>% 
    # Calculate metrics for each cell per year
    group_by(year, xc, yc, lat, lon) %>% 
    summarise(
      # Total days per year (2016 was a leap year)
      total_day_year = max(yday), 
      # Duration of sea ice cover
      seaice_duration = sum(ice_covered_day), 
      # Mean ice concentration
      mean_ice_conc = round(mean(ice_conc, na.rm = T), 2),
      # Ice breakup day of the year
      ice_break = first(subset(., ice_conc < 50)$yday), 
      # Ice breakup week of the year
      ice_week = first(subset(., ice_conc < 50)$week)) %>% 
    # Duration of open water and ice break day and week
    mutate(openwater_duration = total_day_year - seaice_duration, 
           ice_break = if_else(is.na(ice_break) == T, 365, ice_break), 
           ice_week = if_else(is.na(ice_week) == T, 52, ice_week)) %>% 
    as_tibble()
  seaice_year <- bind_rows(seaice_year, seaice_tmp)
  print(paste0("Processing seaice data of ", i, "  finished."))
}
```

```
## [1] "Processing seaice data of 2015  finished."
## [1] "Processing seaice data of 2016  finished."
## [1] "Processing seaice data of 2017  finished."
```

```r
# Save data
save(seaice_year, file = "data/remote_sensing/remote_sensing_seaice_year.RData") 
# Remove temporary data
rm(seaice_tmp, i) 
```

To match sea ice records with acoustic data, I rasterized sea ice data on the same grid as the acoustic data. I tried two different grids; the WGS84 projection (EPSG:4326) with grid cells of 2°lon \* 1°lat, and the EASE-Grid 2.0 North (EPSG:6931) which is the default grid for sea-ice data. For each cell I calculated the mean ice concentration, open water, sea ice duration, day of ice breakup, and week of ice breakup. The ice breakup day is defined as the first day of the year when sea ice concentration fell below 50 %. Similarly, ice breakup week is the first week of the year when sea ice concentration fell below 50 %.

### EPSG:6931 - EASE-Grid 2.0 North (Lambert's equal-area, azimuthal)

More info on this projection can be found on the [NSIDC website](https://nsidc.org/data/ease/). Because "raw" sea ice data has a 25 km \* 25 km resolution, I match that resolution to the final projection (150 km \* 150 km).


```r
# Empty dataframe
seaice_grid_laea <- data.frame() 

# Data gridding
for (i in seq(2015, 2017, 1)) { 
  seaice_tmp <- seaice_year %>%
    filter(year == i)
  # Rasterize data in latlon
  seaice_tmp_laea <- SpatialPointsDataFrame(
    SpatialPoints(cbind(seaice_tmp$xc, seaice_tmp$yc), 
                  proj4string = CRS("EPSG:6931")),
    data.frame(lat = seaice_tmp$lat,
               lon = seaice_tmp$lon,
               seaice_duration = seaice_tmp$seaice_duration,
               openwater_duration = seaice_tmp$openwater_duration,
               mean_ice_conc = seaice_tmp$mean_ice_conc,
               ice_break = seaice_tmp$ice_break,
               ice_week = seaice_tmp$ice_week)) %>%
    # Rasterize data in laea
    rasterize(., arctic_laea, fun = mean, na.rm = T) %>%
    # Remove ID layer
    dropLayer(1) %>%
    # Convert raster to data frame
    rasterToPoints() %>% 
    as.data.frame() %>%
    # Rename variables
    rename(xc = x, yc = y) %>% 
    mutate(year = i, 
           area = factor(
             case_when(lon > -155 & lon <= -95 & lat > 65 & lat <= 82 ~ "BF_CAA",
                       lon > -95 & lon <= -50 & lat > 66 & lat <= 82 ~ "BB",
                       lon >= -25 & lon <= 145 & lat > 77 & lat <= 90 ~ "SV"),
             levels = c("BF_CAA", "BB", "SV"))) %>%
    dplyr::select(year, area, lat, lon, xc, yc, seaice_duration, 
                  openwater_duration, mean_ice_conc, ice_break, ice_week)
  seaice_grid_laea <- bind_rows(seaice_grid_laea, seaice_tmp_laea)
}
# Add cell resolution to dataframe
seaice_grid_laea <- seaice_grid_laea %>%
  mutate(cell_res = cell_res) 
# Remove temporary data
rm(seaice_tmp, seaice_tmp_laea, i, cell_res)
```

Plot mean ice concentration and open water duration per year with the EASE-Grid 2.0 - North projection.


```r
# Coastline in laea
coast_10m_laea <- readOGR("data/bathy/ne_10m_land.shp", verbose = F) %>% 
  spTransform(CRSobj = crs("EPSG:4326")) %>% 
  crop(extent(-180, 180, 0, 90)) %>% 
  spTransform(CRSobj = crs(arctic_laea)) %>%
  fortify() %>%
  rename(xc = long, yc = lat)

# Plot
plot_grid(
  # Map mean ice concentration
  seaice_grid_laea %>% 
    # Remove land
    filter(is.na(mean_ice_conc) == F) %>% 
    ggplot() +
    geom_tile(aes(x = xc, y = yc, fill = mean_ice_conc)) +
    scale_fill_cmocean("Mean ice concentration (%)", name = "ice") +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
                 fill = "grey80") +
    facet_wrap(~ year, ncol = 3) +
    coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) +
    theme(legend.position = "top", 
          legend.key.height = unit(0.1, "in"), 
          legend.key.width = unit(0.3, "in"),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title = element_blank()),
  # Map openwater duration
  seaice_grid_laea %>% 
    # Remove land
    filter(is.na(mean_ice_conc) == F) %>%
    ggplot() +
    geom_tile(aes(x = xc, y = yc, fill = openwater_duration)) +
    scale_fill_viridis_c("Open water duration (days)", direction = -1) +
    geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill =
                   "grey80") +
    facet_wrap(~ year, ncol = 3) +
    coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
    theme(legend.position = "top",
          legend.key.height = unit(0.1, "in"),
          legend.key.width = unit(0.3, "in"),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title = element_blank()),
  ncol = 1, align = "hv", axis = "tblr")
```

![](PanArctic_DSL_remote_sensing_files/figure-html/EPSG-6931-map-seaice-1.png)<!-- -->

# Arctic Ocean Physics Reanalysis

The [Arctic Ocean Physics Reanalysis](https://resources.marine.copernicus.eu/product-detail/ARCTIC_MULTIYEAR_PHY_002_003/INFORMATION) dataset was downloaded using a Jupyter notebook. The dataset has a 12.5 km^2^ grid resolution and a specific polar stereographic projection. Because PROJ4 strings are not have been replaced by WKT2 strings in `rgdal`, I need to convert the PROJ4 custom projection of Copernicus into a WKT2 string. I need to verify if the `st_crs` function converts the proj4 string into a WKT2 string correctly. Data is downloaded at 0, 222, 266, 318, 380, 644, 1062 m depth.

The variables are:

-   thetao: Temperature degree Celsius
-   so: Salinity psu
-   vxo: Zonal velocity m/s (x)
-   vyo: Meridional velocity m/s (y)
-   x: x-coordinate (not longitude, need to be reprojected)
-   y: y-coordinate (not latitude, need to be reprojected)
-   mlotst: mixed layer thickness (m)
-   siconc: sea ice concentration (1)
-   sithick: sea ice thickness (m)
-   depth: depth in meters
-   time: time in hours since 1950-01-01 00:00:00
-   latitude: latitude
-   longitude: longitude

## Exploratory code


```r
# Read sea ice data and mixed layer depth
test_ice <- tidync("data/remote_sensing/physics/20150115_222m_phy_CMEMS.nc") %>% 
  activate("D0,D1,D2") %>%
  hyper_tibble() 
# Read temperature, salinity, and current velocity
test_oce <- tidync("data/remote_sensing/physics/20150115_222m_phy_CMEMS.nc") %>% 
  activate("D0,D1,D3,D2") %>% # D0 = x, D1 = y, D3 = depth, D2 = time
  hyper_tibble()
# Read sea ice data and mixed layer depth
test_meta <- tidync("data/remote_sensing/physics/20150115_222m_phy_CMEMS.nc") %>%
  activate("D0,D1") %>%
  hyper_tibble()

# Combine data
test <- left_join(test_ice, test_oce, by = c("x", "y", "time")) %>% 
  left_join(., test_meta, by = c("x", "y")) %>%
  # Remove data from land or below seafloor
  filter(is.na(depth) == F) %>% 
  # Calculate date
  mutate(date = as.POSIXct(time * 3600, 
                           origin = "1950-01-01 00:00:00",
                           tz = "UTC")) %>% 
  rename(lat = latitude, lon = longitude) %>%
  dplyr::select(-time)

# Define projection
proj_original <- "+proj=stere +a=6378273.0 +b=6378273.0 +lon_0=-45.0 +lat_0=90.0 +lat_ts=90.0 +ellps=sphere +units=m" 

# Convert dataframe to sf
test_sf <- st_as_sf(test, coords = c("x", "y"), crs = proj_original)
plot(test_sf["vxo"], axes = TRUE)
```

![](PanArctic_DSL_remote_sensing_files/figure-html/processing-trial-sf-1.png)<!-- -->

```r
# Convert dataframe to sf
test_sf2 <- st_as_sf(test, coords = c("x", "y"), crs = proj_original) %>% 
  # Transform to EPSG:4326
  st_transform(., crs = "EPSG:4326") 
# Not working
plot(test_sf2["vxo"], axes = TRUE) 
```

![](PanArctic_DSL_remote_sensing_files/figure-html/processing-trial-sf-2.png)<!-- -->

```r
# Convert dataframe to sf with correct CRS
test_sf3 <- st_as_sf(test, coords = c("lon", "lat"), crs = "EPSG:4326") 
# Working!
plot(test_sf3["vxo"], axes = TRUE) 
```

![](PanArctic_DSL_remote_sensing_files/figure-html/processing-trial-sf-3.png)<!-- -->

```r
# Convert dataframe to sf with correct CRS
test_sf4 <- st_as_sf(test, coords = c("lon", "lat"), crs = "EPSG:4326") %>% 
  # Convert to EPSG:6931
  st_transform(., crs = "EPSG:6931") 
# Working
plot(test_sf4["vxo"], axes = TRUE) 
```

![](PanArctic_DSL_remote_sensing_files/figure-html/processing-trial-sf-4.png)<!-- -->

```r
# Remove temporary data
rm(test_ice, test_meta, test_oce, test_sf, test_sf2, test_sf3, test_sf4, test)
```

## Batch processing

Batch processing of Arctic Ocean Physics Reanalysis data. For convenience, all NetCDF data is combined per year and saved as csv. Since I downloaded physics data in two batches—0, 222, 380, 644, 1062 m depth and 266 and 318 m depth, I split the tidying into two chunks.


```r
# Loop through data per year
for (y in seq(2015, 2017, 1)) { 
  # Empty dataframe
  phy_year <-  data.frame() 
  # List files per year
  file_list <- list.files("data/remote_sensing/physics", 
                          pattern = paste0(y), 
                          full.names = T) 
  # Loop through each file
  for (i in file_list) { 
    # Read NetCDF 
    phy_ice_tmp <- tidync(i) %>%
      # sea ice data and mixed layer depth
      activate("D0,D1,D2") %>% 
      hyper_tibble()
    phy_oce_tmp <- tidync(i) %>%
      # temperature, salinity, depth, and current velocity
      activate("D0,D1,D3,D2") %>% 
      hyper_tibble()
    phy_coord_tmp <- tidync(i) %>%
      # latitude and longitude
      activate("D0,D1") %>% 
      hyper_tibble()
    # Combine datasets
    phy_tmp <- left_join(phy_ice_tmp, phy_oce_tmp, by = c("x", "y", "time")) %>%
      left_join(., phy_coord_tmp, by = c("x", "y")) %>%
      # Remove data from land or below seafloor
      filter(is.na(depth) == F) %>% 
      # Calculate date
      mutate(date = as.POSIXct(time * 3600, 
                               origin = "1950-01-01 00:00:00",
                               tz = "UTC")) %>% 
      rename(lat = latitude, lon = longitude)
    # Append data together
    phy_year <- bind_rows(phy_year, phy_tmp)
  }
  # Write one csv per year
  write_csv(phy_year, 
            file = paste0("data/remote_sensing/physics/physics_reanalysis_",
                          y, ".csv")) 
  print(paste0("Processing of ", y, " data finished."))
}
# Remove temprorary data
rm(phy_coord_tmp, phy_ice_tmp, phy_oce_tmp, phy_tmp, phy_year)
```


```r
# Loop through data per year
for (y in seq(2015, 2017, 1)) {
  # Empty dataframe
  phy_year <-  data.frame() 
  # List files per year
  file_list <- list.files("data/remote_sensing/physics/phy_266_318m", 
                          pattern = paste0(y), 
                          full.names = T) 
  # Loop through each file
  for (i in file_list) { 
    # Read NetCDF 
    phy_ice_tmp <- tidync(i) %>%
      # sea ice data and mixed layer depth
      activate("D0,D1,D2") %>% 
      hyper_tibble()
    phy_oce_tmp <- tidync(i) %>%
      # temperature, salinity, depth, and velocity
      activate("D0,D1,D3,D2") %>% 
      hyper_tibble()
    phy_coord_tmp <- tidync(i) %>%
      # latitude and longitude
      activate("D0,D1") %>% 
      hyper_tibble()
    # Combine datasets
    phy_tmp <- left_join(phy_ice_tmp, phy_oce_tmp, by = c("x", "y", "time")) %>%
      left_join(., phy_coord_tmp, by = c("x", "y")) %>%
      # Remove data from land or below seafloor
      filter(is.na(depth) == F) %>% 
      # Calculate date
      mutate(date = as.POSIXct(time * 3600, 
                               origin = "1950-01-01 00:00:00",
                               tz = "UTC")) %>% 
      rename(lat = latitude, lon = longitude)
    # Append data together
    phy_year <- bind_rows(phy_year, phy_tmp)
  }
  # Write one csv per year
  write_csv(
    phy_year, 
    file = paste0("data/remote_sensing/physics/physics_reanalysis_266_318m_",
                  y, ".csv")) 
  print(paste0("Processing of ", y, " data finished."))
}
# Remove temporary data
rm(phy_coord_tmp, phy_ice_tmp, phy_oce_tmp, phy_tmp, phy_year)
```

Tidy csv data and calculate current velocity and angle. I also calculate the annual average of each variable. The annual average dataset will be used for comparison with acoustic data.


```r
# List files
phy_year <- list.files("data/remote_sensing/physics",
                  pattern = "*.csv",
                  full.names = T) %>% 
  # Read files 
  map_dfr(.f = ~ read_csv(., show_col_types = F)) %>% 
  lazy_dt() %>%
  mutate(year = year(date),
         month = month(date),
         # Current velocity
         velocity = sqrt(vxo ^ 2 + vyo ^ 2), 
         # Current direction
         v_angle = atan(vyo / vxo) * (180 / pi)) %>% 
  arrange(year, date, depth) %>%
  dplyr::select(- date, - month) %>%
  # Group by cell by year by depth
  group_by(year, x, y, depth) %>% 
  # Calculate annual average for each variable
  summarise_all(mean) %>% 
  ungroup() %>%
  as_tibble()
```

Plot average current velocity in 2015 in Svalbard at 318 m depth.


```r
phy_year %>%
  filter(year == 2015 & depth == 318 & between(x, 7, 15) & between(y, -10, -6)) %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = velocity)) + 
  geom_vector(aes(x = x, y = y, dx = vxo, dy = vyo, mag = velocity),
              arrow.angle = 30, arrow.type = "open", arrow.length = .5,
              pivot = 0, preserve.dir = TRUE, direction = "ccw") +
  scale_fill_cmocean(name = "speed") +
  coord_cartesian()
```

![](PanArctic_DSL_remote_sensing_files/figure-html/plot-current-1.png)<!-- -->

### EPSG:6931 - EASE-Grid 2.0 North (Lambert's equal-area, azimuthal)

Data is projected on the EPSG:6931 projection. I also calculate temperature (x - mean(x over the whole time series) and current anomalies.


```r
# Projection, cell resolution in km
cell_res <- 12.5
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931") 
projection(arctic_laea) <- gsub("units=m", "units=km", projection(arctic_laea)) 
res(arctic_laea) <- c(cell_res, cell_res)

# Empty dataframe
phy_grid_laea <- data.frame() 

# Loop through each year
for (i in seq(2015, 2017, 1)) {
  # Loop through each depth bin
  for (j in c(0, 222, 266, 318, 380, 644, 1062)) { 
    phy_year_tmp <- phy_year %>%
      # Select depth bin
      filter(year == i & depth == j) 
    # Rasterize data in laea
    phy_tmp_laea <- SpatialPointsDataFrame(
      SpatialPoints(cbind(phy_year_tmp$lon, phy_year_tmp$lat),
                    proj4string = CRS("EPSG:4326")), 
      data.frame(lat = phy_year_tmp$lat,
                 lon = phy_year_tmp$lon,
                 depth = phy_year_tmp$depth,
                 mlotst = phy_year_tmp$mlotst,
                 siconc = phy_year_tmp$siconc,
                 sithick = phy_year_tmp$sithick,
                 thetao = phy_year_tmp$thetao,
                 so = phy_year_tmp$so,
                 velocity = phy_year_tmp$velocity,
                 vxo = phy_year_tmp$vxo,
                 vyo = phy_year_tmp$vyo)) %>%
      # Change projection to EPSG:6931
      spTransform(., CRSobj = crs(arctic_laea)) %>% 
      # Rasterize 
      rasterize(., arctic_laea, fun = mean, na.rm = T) %>% 
      # Remove ID layer
      dropLayer(1) %>% 
      # Convert raster to data frame
      rasterToPoints() %>% 
      as.data.frame() %>%
      # Rename variables
      rename(xc = x, yc = y) %>% 
      mutate(year = i, 
             depth = j,
             area = factor(
               case_when(lon > -155 & lon <= -95 & lat > 65 & lat <= 82 ~ "BF_CAA",
                         lon > -95 & lon <= -50 & lat > 66 & lat <= 82 ~ "BB",
                         lon >= -25 & lon <= 145 & lat > 77 & lat <= 90 ~ "SV"),
               levels = c("BF_CAA", "BB", "SV"))) %>%
      dplyr::select(year, area, lat, lon, xc, yc, depth, mlotst, siconc, sithick,
                    thetao, so, velocity, vxo, vyo) 
    # Combine data
    phy_grid_laea <- bind_rows(phy_grid_laea, phy_tmp_laea)
  }
}

# Add cell resolution
phy_grid_laea <- phy_grid_laea %>%
  mutate(cell_res = cell_res)

# Remove temporary data
rm(phy_year_tmp, phy_tmp_laea, i, j) 
```

Gridding worked correctly, the following plot shows current velocity at all depth levels.


```r
coast_10m_laea <- readOGR("data/bathy/ne_10m_land.shp", verbose = F) %>% 
  spTransform(CRSobj = crs("EPSG:4326")) %>% 
  crop(extent(-180, 180, 0, 90)) %>% # Crop shapefile
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to a dataframe for ggplot
  rename(xc = long, yc = lat)

phy_grid_laea %>% 
  ggplot() +
  geom_tile(aes(x = xc, y = yc, fill = velocity * 100)) +
  scale_fill_cmocean("Current velocity (cm/s)", name = "speed") +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  facet_grid(depth ~ year) +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) + 
  theme(legend.position = "top",
        legend.key.height = unit(0.1, "in"),
        legend.key.width = unit(0.3, "in"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
```

![](PanArctic_DSL_remote_sensing_files/figure-html/EPSG-6931-map-physics-1.png)<!-- -->

# Ocean colour data

Data downloaded from Copernicus Marine. The ocean colour data comes from the [OCEANCOLOUR_ARC_CHL_L4_REP_OBSERVATIONS_009_088](https://resources.marine.copernicus.eu/product-detail/OCEANCOLOUR_ARC_CHL_L4_REP_OBSERVATIONS_009_088/INFORMATION) dataset. Data is monthly averaged and has a cell resolution of 1 x 1 km on a WGS 84 / Plate Carree (EPSG:32662) projection.

I load the data and project it on the arctic laea grid (EPSG:6931 25 x 25 km). The time variable is seconds since 1970-01-01 00:00:00.


```r
# Sea ice projection
cell_res <- 12.5
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931") 
projection(arctic_laea) <- gsub("units=m", "units=km", projection(arctic_laea)) 
res(arctic_laea) <- c(cell_res, cell_res)

# List files
file_list <- list.files("data/remote_sensing/ocean_colour/nc",
                        pattern = ".nc", full.names = T) 
# Loop through each file per year
for (i in file_list) { 
  # Extract date
  date_file <- str_remove(i, pattern = "data/remote_sensing/ocean_colour/nc/") %>%
    str_remove(., pattern = "_chl_CMEMS.nc") 
  # Chlorophyll data
  chl_tmp <- tidync(i) %>%
    activate("D2,D1,D0") %>% 
    hyper_tibble() %>%
    rename(lat = latitude, lon = longitude, chl = CHL)
  
  # If data is not empty then rasterize
  if (nrow(chl_tmp > 0)){
    # Rasterize data in laea
    chl_tmp_laea <- SpatialPointsDataFrame(
      SpatialPoints(cbind(chl_tmp$lon, chl_tmp$lat),
                    proj4string = CRS("EPSG:4326")), 
      data.frame(lat = chl_tmp$lat,
                 lon = chl_tmp$lon,
                 time = chl_tmp$time,
                 chl = chl_tmp$chl)) %>%
      # Project to EPSG:6931
      spTransform(., CRSobj = crs(arctic_laea)) %>% 
      # Rasterize on the arctic_laea grid (with correct grid resolution)
      rasterize(., arctic_laea, fun = mean, na.rm = T) %>% 
      # Remove ID layer
      dropLayer(1) %>% 
      # Convert raster to data frame
      rasterToPoints() %>% 
      as.data.frame() %>%
      # Rename variables
      rename(xc = x, yc = y) %>% 
      # Tidy data
      mutate(date = as.POSIXct(time , origin = "1970-01-01 00:00:00", tz = "UTC"),
             year = year(date), 
             cell_res = cell_res) 
    # Write csv
    write_csv(chl_tmp_laea, 
              file = paste0("data/remote_sensing/ocean_colour/ocean_colour_",
                            date_file, "_", cell_res, "km.csv"))
    print(paste0("Processing of ", date_file, " chl data finished."))
    # Remove temporary file
    rm(chl_tmp_laea)
  }
  # Else move on to next file
  else {
    print(paste0("No chl data for ", date_file, ". Moving to next file."))
  }
  # Remove temporary variables
  rm(chl_tmp)
}
```

```
## [1] "No chl data for 20150101. Moving to next file."
## [1] "Processing of 20150201 chl data finished."
## [1] "Processing of 20150301 chl data finished."
## [1] "Processing of 20150401 chl data finished."
## [1] "Processing of 20150501 chl data finished."
## [1] "Processing of 20150601 chl data finished."
## [1] "Processing of 20150701 chl data finished."
## [1] "Processing of 20150801 chl data finished."
## [1] "Processing of 20150901 chl data finished."
## [1] "Processing of 20151001 chl data finished."
## [1] "No chl data for 20151101. Moving to next file."
## [1] "No chl data for 20151201. Moving to next file."
## [1] "No chl data for 20160101. Moving to next file."
## [1] "Processing of 20160201 chl data finished."
## [1] "Processing of 20160301 chl data finished."
## [1] "Processing of 20160401 chl data finished."
## [1] "Processing of 20160501 chl data finished."
## [1] "Processing of 20160601 chl data finished."
## [1] "Processing of 20160701 chl data finished."
## [1] "Processing of 20160801 chl data finished."
## [1] "Processing of 20160901 chl data finished."
## [1] "Processing of 20161001 chl data finished."
## [1] "No chl data for 20161101. Moving to next file."
## [1] "No chl data for 20161201. Moving to next file."
## [1] "No chl data for 20170101. Moving to next file."
## [1] "Processing of 20170201 chl data finished."
## [1] "Processing of 20170301 chl data finished."
## [1] "Processing of 20170401 chl data finished."
## [1] "Processing of 20170501 chl data finished."
## [1] "Processing of 20170601 chl data finished."
## [1] "Processing of 20170701 chl data finished."
## [1] "Processing of 20170801 chl data finished."
## [1] "Processing of 20170901 chl data finished."
## [1] "Processing of 20171001 chl data finished."
## [1] "No chl data for 20171101. Moving to next file."
## [1] "No chl data for 20171201. Moving to next file."
```

```r
# Remove temporary variables
rm(date_file, i, file_list)
```

After rasterizing data I reload the csv files and combine the monthly data into a single dataframe. Then, I calculate the average annual concentration of chl a per cell.


```r
# List files
chl_grid_laea <- list.files("data/remote_sensing/ocean_colour", 
                            pattern = "*.csv",
                            full.names = T) %>% 
  # Read files
  map_dfr(.f = ~ read_csv(., show_col_types = F)) %>% 
  mutate(area = factor(
    case_when(lon > -155 & lon <= -95 & lat > 65 & lat <= 82 ~ "BF_CAA",
              lon > -95 & lon <= -50 & lat > 66 & lat <= 82 ~ "BB",
              lon >= -25 & lon <= 145 & lat > 77 & lat <= 90 ~ "SV"),
    levels = c("BF_CAA", "BB", "SV"))) %>%
  group_by(year, area, xc, yc, cell_res) %>%
  summarise(chl = mean(chl)) %>%
  ungroup()
```

# Save data


```r
# Save data
save(seaice_grid_laea,
     file = paste0("data/remote_sensing/seaice_grids_", cell_res, "km.RData"))
save(phy_grid_laea, 
     file = paste0("data/remote_sensing/physics_grids_", cell_res, "km.RData"))
save(chl_grid_laea, 
     file = paste0("data/remote_sensing/chl_grids_", cell_res, "km.RData"))
```
