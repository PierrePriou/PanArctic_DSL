---
title: "PanArctic DSL - Acoustic gridding"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/05/05 at 14:12"
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
library(dtplyr)       # Speed up data processing
library(lubridate)    # Deal with dates
library(suncalc)      # Sun positions
library(marmap)       # Get depth
library(cowplot)      # Plots on a grid
library(raster)       # Data gridding
library(sf)           # Spatial data
library(fs)           # List files
library(rgdal)        # Read shapefiles
library(cmocean)      # Oceanographic color palettes
library(RColorBrewer) # Diverging color palette
# Custom figure theme
theme_set(theme_bw())
theme_update(axis.text = element_text(size = 9),
             axis.title = element_text(size = 9),
             strip.text.x = element_text(size = 9, face = "plain", hjust = 0.5),
             strip.background = element_rect(colour = "transparent", fill = "transparent"),
             legend.title = element_text(size = 9, vjust = 1),
             legend.margin = margin(0, 0, 0, 0),
             legend.box.margin = margin(0, 0, -8, 0),
             panel.grid = element_blank(), 
             plot.margin = unit(c(0.02, 0.02, 0.02, 0.02), "in"),
             plot.title = element_text(size = 9, face = "bold"))
options(dplyr.summarise.inform = F) # Suppress summarise() warning
```

I use the area definitions from IHO Sea Areas (International Hydrographic Organization) to combine data.


```r
# Projections
arctic_latlon <- raster(extent(-155, 35, 66, 85), # Base projection for acoustic and CTD data
                        crs = "EPSG:4326", 
                        res = c(2, 1)) # cells of 2 degree longitude per 1 degree latitude
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931") # Seaice projection
projection(arctic_laea) <- gsub("units=m", "units=km", projection(arctic_laea)) # Convert proj unit from m to km
cell_res <- 50 # Cell resolution in km
res(arctic_laea) <- c(cell_res, cell_res) # Define the 100 km cell resolution

# IHO areas
IHO_sf_laea <- dir_ls("data/arctic_regions/", glob = "*.shp") %>% 
  tibble(fname = .) %>%
  mutate(data = map(fname, read_sf)) %>%
  unnest(data) %>%
  st_as_sf() %>%
  st_transform(crs = crs(arctic_latlon)) %>%
  st_transform(crs = crs(arctic_laea)) %>%
  st_make_valid()

# Bathy data from marmap: for extracting bottom depth
bathy <- getNOAA.bathy(lon1 = -152, lon2 = 35, lat1 = 65, lat2 = 84, resolution = 2, keep = T, path = "data/bathy/")
```

```
## File already exists ; loading 'marmap_coord_-152;65;35;84_res_2.csv'
```

# Integrated s~A~ and centre of mass

Acoustic data were collected continuously at 38 kHz. We selected data when the ship were stationary (< 1 knot). To match acoustic records to CTD and remote sensing data, and to homogenize data spatio-temporally, I rasterized acoustic data per year. I used two different grids; the WGS84 projection (EPSG:4326) with grid cells of 2°lon x 1°lat, and the EASE-Grid 2.0 North (EPSG:6931)—which is the default grid for sea-ice data—with grid cells of 150 km x 150 km. For each cell I calculated the mean integrated nautical area scattering coefficient (s~A~; m^2^ nmi^-2^) and centre of mass over mesopelagic depths (200-1000 m). Then, I calculated the normalized backscatter anomaly for each grid cell per area—Beaufort Sea and Canadian Arctic Archipelago, Baffin Bay, and Svalbard—per year.


```r
# Load and tidy files
MVBS_raw <- list.files("C:/Users/cfer/PhD/Chapt 2 - Arctic Mesopelagic DSL/data/acoustics/90dB threshold", 
                       pattern = "*.csv", full.names = TRUE) %>% # list files in folder
  set_names() %>%
  map_dfr(.f = read_csv, .id = "filename") %>% # reads file
  dplyr::select(-Interval, -Layer, -Dist_S, -Sv_min, -Sv_max) %>%
  rename("layer_depth_min" = "Layer_depth_min",
         "layer_depth_max" = "Layer_depth_max",
         "lat" = "Lat_S",
         "lon" = "Lon_S",
         "frequency" = "Frequency") %>%
  unite(date, Date_S, Time_S, sep = " ", remove = FALSE) %>%
  mutate(filename = str_remove(filename, pattern = "C:/Users/cfer/PhD/Chapt 2 - Arctic Mesopelagic DSL/data/acoustics/90dB threshold"),
         date = ymd_hm(format(ymd_hms(date, tz = "UTC"), format = '%Y%m%d %H:%M'), tz = "UTC"),
         Date_S = as.character(ymd(Date_S)),
         date_num = as.integer(as.Date(as.character(ymd(Date_S)))), # Used for left_join with suncalc dataset
         area = factor(case_when(lon > -155 & lon <= -95 & lat > 65 & lat <= 82 ~ "BF_CAA",
                                 lon > -95 & lon <= -50 & lat > 66 & lat <= 82 ~ "BB",
                                 lon >= -25 & lon <= 145 & lat > 77 & lat <= 90 ~ "SV"),
                       levels = c("BF_CAA", "BB", "SV")))

# Find dawn and dusk times. Suncalc requires the variables to be named "date", "lat", and "lon".
MVBS_suncalc <- MVBS_raw %>%
  group_by(filename, area, Date_S, date_num) %>%
  summarise(lat = mean(lat, na.rm = T),
            lon = mean(lon, na.rm = T)) %>%
  ungroup() %>%
  # find dusk and dawn time
  mutate(suncalc = getSunlightTimes(data = data.frame(date = as.Date(Date_S), lat = lat, lon=  lon), 
                                      keep = c("dusk", "dawn"), tz="UTC"),
         dusk = suncalc$dusk,
         dawn = suncalc$dawn) %>%
  dplyr::select(filename, area, date_num, dawn, dusk)

# Find bottom depth for each MVBS cell
MVBS_bottom_depth <- MVBS_raw %>%
  group_by(filename, area, date, lat, lon) %>%
  summarise(bottom_depth_MVBS = max(layer_depth_min, na.rm = T)) %>%
  mutate(get_depth = get.depth(bathy, x = lon, y = lat, locator = F), 
         bottom_depth_GEBCO = as.numeric(abs(get_depth$depth)),
         bottom_depth = if_else(bottom_depth_MVBS < 995, bottom_depth_MVBS, bottom_depth_GEBCO)) %>%
  ungroup()

# Append dawn and dusk times, bottom depth to main dataframe
MVBS <- MVBS_raw %>%
  left_join(., MVBS_suncalc, by=c("filename", "area", "date_num")) %>%
  left_join(., MVBS_bottom_depth) %>%
  mutate(year = year(date), 
         month = month(date), 
         day_night = factor(if_else(is.na(dawn) & is.na(dusk) & month >= 4 & month <= 10, "day", # polar day case
                            if_else(is.na(dawn) & is.na(dusk) & month < 4 | month > 10, "night", # polar night case
                            if_else(date >= dawn & date <= dusk, "day", "night"))),
                            levels = c("day", "night"))) %>% 
  filter(layer_depth_min <= 995 & bottom_depth > 200 & lat != 999 & lon != 999 & Sv_mean < -30) %>% # Tidy data
  dplyr::select(year, area, date, day_night, lat, lon, layer_depth_min, Sv_mean, NASC, bottom_depth, frequency) 

rm(MVBS_suncalc, MVBS_bottom_depth)
```

Calculate integrated NASC over mesopelagic depth (200 - 1000 m depth) and centre of mass for each 10 min cell.


```r
SA_integrated <- MVBS %>%
  mutate(sv_lin = 10 ^ (Sv_mean / 10)) %>%  # linear volume backscattering coefficient
  filter(layer_depth_min >= 200 & layer_depth_min <= 995) %>% # Select mesopelagic depths
  group_by(year, area, date, day_night, lat, lon, bottom_depth) %>%
  summarise(sv_lin_mean = mean(sv_lin),
            depth_integration = max(layer_depth_min) - min(layer_depth_min) + 5, # +5 because each cell is 5m high
            NASC_int = sum(NASC), # integrated backscatter (linear)
            CM = sum(layer_depth_min * NASC) / sum(NASC)) %>% # centre of mass in meters
  ungroup() %>%
  mutate(SA_int = 10 * log10(NASC_int))
```

Save data.


```r
save(MVBS, SA_integrated, file = "data/acoustics/MVBS_2015_2017.RData")
```

## 2D Integrated backscatter - EPSG:6931 - EASE-Grid 2.0 North (Lambert's equal-area, azimuthal)

More info on this projection can be found on the [NSIDC website](https://nsidc.org/data/ease/).


```r
SA_grid_laea <- data.frame() # Empty dataframe that will be filled with gridded CTD data

for (i in seq(2015, 2017, 1)) { # Data gridding
  SA_tmp <- SA_integrated %>%
    filter(year == i & day_night == "day") %>%
    dplyr::select(year, lat, lon, NASC_int, CM)
  # Rasterize data in latlon
  SA_tmp_laea <- SpatialPointsDataFrame(SpatialPoints(cbind(SA_tmp$lon, SA_tmp$lat), 
                                                        proj4string = CRS("EPSG:4326")),
                                          data.frame(lat = SA_tmp$lat,
                                                     lon = SA_tmp$lon,
                                                     NASC_int = SA_tmp$NASC_int,
                                                     CM = SA_tmp$CM)) %>%
    spTransform(., CRSobj = crs(arctic_laea)) %>% # Change projection to EPSG:6931
    rasterize(., arctic_laea, fun = mean, na.rm = T) %>% # Rasterize data in latlon
    dropLayer(1) %>% # Remove ID layer
    rasterToPoints() %>% # Convert raster to data frame
    as.data.frame() %>%
    rename(xc = x, yc = y) %>% # Rename variables
    mutate(year = i, 
           area = factor(case_when(lon > -155 & lon <= -95 & lat > 65 & lat <= 82 ~ "BF_CAA",
                                   lon > -95 & lon <= -50 & lat > 66 & lat <= 82 ~ "BB",
                                   lon >= -25 & lon <= 145 & lat > 77 & lat <= 90 ~ "SV"),
                         levels = c("BF_CAA", "BB", "SV"))) %>%
    dplyr::select(year, area, lat, lon, xc, yc, NASC_int, CM) 
  SA_grid_laea <- bind_rows(SA_grid_laea, SA_tmp_laea)
}
rm(SA_tmp, SA_tmp_laea, i) # Remove temporary data

SA_grid_laea <- SA_grid_laea %>% # Calculate normalized backscatter anomalies
  mutate(cell_res = cell_res, # Add cell resolution to dataframe
         SA_int = 10 * log10(NASC_int)) %>% # Calculate integrated backscatter strength
  st_as_sf(coords = c("xc", "yc"), crs = st_crs(arctic_laea), remove = F) %>%
  st_join(., IHO_sf_laea, join = st_within) %>% # Append IHO region
  st_drop_geometry() %>%
  mutate(name = if_else(is.na(name) == T & yc < 0, "Baffin Bay",
                if_else(is.na(name) == T & yc > 0 & lat > 70, "The Northwestern Passages",
                if_else(name == "The Northwestern Passages" & yc < 0, "Baffin Bay",
                if_else(name == "Arctic Ocean" & name_3 == "West Arctic Ocean", "West Arctic Ocean",
                if_else(name == "Arctic Ocean" & name_3 == "East Arctic Ocean", "East Arctic Ocean", name))))),
         empty = factor(if_else(SA_int < 0, T, F))) %>%
    rename(IHO_area = name, area = area.x) %>%
  dplyr::select(year, xc, yc, lon, lat, IHO_area, area, NASC_int, SA_int, CM, empty, cell_res)
```

Map to check whether the IHO regions are well implemented.


```r
coast_10m_laea <- readOGR("data/bathy/ne_10m_land.shp", verbose = F) %>% # Coastline in laea
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  crop(extent(-180, 180, 0, 90)) %>% # Crop shapefile
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to a dataframe for ggplot
  rename(xc = long, yc = lat)

SA_grid_laea %>%
  ggplot() +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey80") +
  geom_point(aes(x = xc, y = yc, col = IHO_area)) +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F)
```

![](PanArctic_DSL_acoustics_files/figure-html/plot-areas-1.png)<!-- -->

## 3D S~V~ profiles - EPSG:6931 - EASE-Grid 2.0 North (Lambert's equal-area, azimuthal)

More info on this projection can be found on the [NSIDC website](https://nsidc.org/data/ease/).I calculate the median Sv for each depth interval of each georeferenced grid cell.


```r
Sv_grid_laea <- data.frame() # Empty dataframe that will be filled with gridded CTD data

for (i in seq(2015, 2017, 1)) { # Loop through every year
  for (j in seq(20, 995,5)) { # Loop through each depth bin
    Sv_tmp <- MVBS %>%
     filter(year == i & day_night == "day" & layer_depth_min == j) %>%
      dplyr::select(year, lat, lon, layer_depth_min, Sv_mean) %>%
      mutate(sv_lin = 10 ^ (Sv_mean / 10)) # Transform into backscattering coefficient (lin)
    # Calculate mean
    Sv_tmp_median <- SpatialPointsDataFrame(SpatialPoints(cbind(Sv_tmp$lon, Sv_tmp$lat), 
                                                       proj4string = CRS("EPSG:4326")),
                                         data.frame(lat = Sv_tmp$lat,
                                                    lon = Sv_tmp$lon,
                                                    sv_lin = Sv_tmp$sv_lin)) %>%
      spTransform(., CRSobj = crs(arctic_laea)) %>% # Change projection to EPSG:6931
      rasterize(., arctic_laea, fun = function(x, ...) {quantile(unique(na.omit(x)), 0.5)}, na.rm = F) %>% 
      dropLayer(1) %>% # Remove ID layer
      rasterToPoints() %>% # Convert raster to data frame
      as.data.frame() %>%
      rename(xc = x, yc = y, sv_lin_median = sv_lin) %>% # Rename variables
      mutate(year = i, depth = j)
    Sv_grid_laea <- bind_rows(Sv_grid_laea, Sv_tmp_median) # Combine data
  }
}
rm(Sv_tmp, Sv_tmp_median) # Remove temporary data

coord_equivalences <- SA_grid_laea %>% # Calculate coord equivalences between 4326 and 6931
  group_by(xc, yc, area, IHO_area) %>%
  summarise()

Sv_grid_laea <- Sv_grid_laea %>% # Tidy data
    mutate(Sv_median = 10 * log10(sv_lin_median), 
           area = factor(case_when(lon > -155 & lon <= -95 & lat > 65 & lat <= 82 ~ "BF_CAA",
                                   lon > -95 & lon <= -50 & lat > 66 & lat <= 82 ~ "BB",
                                   lon >= -25 & lon <= 145 & lat > 77 & lat <= 90 ~ "SV"),
                         levels = c("BF_CAA", "BB", "SV")),
           cell_res = cell_res,
           empty = factor(if_else(Sv_median < -90, T, F))) %>%
  ungroup() %>%
  left_join(., coord_equivalences, by = c("xc", "yc", "area")) %>%
  dplyr::select(year, xc, yc, lon, lat, IHO_area, area, depth, sv_lin_median, Sv_median, empty, cell_res)
```

# Save data 


```r
save(SA_grid_laea, file = "data/acoustics/SA_grids.RData") # Save data
save(Sv_grid_laea, file = "data/acoustics/Sv_grids.RData") # Save data
```
