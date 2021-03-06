PanArctic DSL - Figures
================
[Pierre Priou](mailto:pierre.priou@mi.mun.ca)
2022/05/02 at 18:12

# Package and data loading

``` r
# Load packages
library(tidyverse)              # Tidy code
library(lubridate)              # Deal with dates
library(kableExtra)             # Pretty tables
library(cowplot)                # Plots on a grid
library(raster)                 # Data gridding
library(rgdal)                  # Read shapefiles
library(rworldxtra)             # Higher resolution coastline data
library(marmap)                 # Bathymetry data of the Arctic Ocean
library(ggforce)                # Draw polygons
library(cmocean)                # Pretty color palettes
library(ggpubr)                 # Display stats with ggplot
library(RColorBrewer)           # Diverging color palette
library(ggnewscale)             # Multiple color scales on a single plot
library(lemon)                  # Repeated axis on facets
source("R/getNOAA.ice.bathy.R") # Load bathy data from NOAA
source("R/rounder.R")
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
             plot.margin = unit(c(0, 0, 0, 0), "in"))
options(dplyr.summarise.inform = F) # Suppress summarise() warning
```

``` r
# Map projection
arctic_latlon <- raster(extent(-155, 35, 66, 85), # Base projection for acoustic and CTD data
                        crs = "EPSG:4326", 
                        res = c(2, 1)) # cells of 2 degree longitude per 1 degree latitude

# Coastline shapefiles
coast_10m_latlon <- readOGR("data/bathy/ne_10m_land.shp", verbose = F) %>% # Coastline in laea
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  crop(extent(-180, 180, 0, 90)) %>% # Crop shapefile
  fortify() %>% # Convert to a dataframe for ggplot
  rename(lon = long)

# Bathy data
bathy_df <- getNOAA.ice.bathy(lon1 = -180, lon2 = 180, lat1 = 50, lat2 = 90,
                              resolution = 2, keep = T, path = "data/bathy/") %>%
  fortify.bathy() %>%
  rename(lon = x, lat = y, depth = z) %>%
  mutate(depth_d = factor(case_when(between(depth, -50, Inf) ~ "0-50",
                                  between(depth, -100, -50) ~ "50-100",
                                  between(depth, -200, -100) ~ "100-200",
                                  between(depth, -500, -200) ~ "200-500",
                                  between(depth, -1000, -500) ~ "500-1000",
                                  between(depth, -2000, -1000) ~ "1000-2000",
                                  between(depth, -4000, -2000) ~ "2000-4000",
                                  between(depth, -Inf, -4000) ~ "4000-8000"),
                        levels=c(">0","0-50","50-100","100-200","200-500","500-1000","1000-2000","2000-4000","4000-8000")))

# Arctic circle and areas of interest
arctic_circle_latlon <- data.frame(lon = seq(-180, 180, 0.1)) %>%
  mutate(lat = 66.5)
area_BF_CAA_latlon <- data.frame(lon = c(-158, -158, -93, -93, -158), lat = c(79, 68.5, 68.5, 79, 79))
area_BB_latlon <- data.frame(lon = c(-86, -65, -51.5, -51.5, -86), lat = c(72, 65.5, 65.5, 83, 83))
area_SV_latlon <- data.frame(lon = c(46, 46, 0, 0, 46), lat = c(85, 76, 76, 85, 85))

# Acoustic data
load("data/acoustics/MVBS_2015_2017.RData")
load("data/acoustics/SA_grids.RData")
load("data/acoustics/Sv_grids.RData")
MVBS_latlon <- MVBS %>% # Location of acoustic data
  mutate(lat = round(lat, 1),
         lon = round(lon, 1)) %>%
  group_by(lat,lon) %>%
  summarise(lat = mean(lat), 
            lon = mean(lon)) %>%
  ungroup() 

# CTD data
load("data/CTD/CTD_2015_2019.RData")
load("data/CTD/CTD_grids.RData")

# Trawl data
load("data/nets/trawl_mwt_2015_2017.RData")
trawl_latlon <- trawl_station %>% # Location of trawls
  group_by(date) %>%
  summarise(lat = mean(lat),
            lon = mean(lon))

# Remote sensing data
load("data/remote_sensing/remote_sensing_chl.RData")
load("data/remote_sensing/seaice_grids.RData")
```

I re-project data from the WGS84 projection to the EASE-Grid 2.0 North
projection.

``` r
# Laea projection for bathy and all other files
arctic_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931") # Seaice projection
projection(arctic_laea) <- gsub("units=m", "units=km", projection(arctic_laea)) # Convert proj unit from m to km
cell_res <- 50 # Cell resolution in km
res(arctic_laea) <- c(cell_res, cell_res) # Define the 100 km cell resolution

proj_bathy_laea <- raster(extent(-2700, 2700, -2700, 2700), crs = "EPSG:6931", res = c(5, 5)) # Seaice projection
projection(proj_bathy_laea) <- gsub("units=m", "units=km", projection(arctic_laea)) # Convert proj unit from m to km

# Project bathymetric data 
bathy_laea <- SpatialPointsDataFrame(SpatialPoints(cbind(bathy_df$lon, bathy_df$lat), 
                                                   proj4string = CRS("EPSG:4326")),
                                     data.frame(depth = bathy_df$depth)) %>%
  spTransform(., CRSobj = crs(proj_bathy_laea)) %>% # Change projection to EPSG:6931
  rasterize(., proj_bathy_laea, fun = mean, na.rm = T) %>% # Rasterize data in latlon
  dropLayer(1) %>% # Remove ID layer
  rasterToPoints() %>% # Convert raster to data frame
  as.data.frame() %>%
  rename(xc = x, yc = y) %>%
  ungroup() %>%
  mutate(depth_d = factor(case_when(between(depth, -100, Inf) ~ "0-100",
                                    between(depth, -200, -100) ~ "100-200",
                                    between(depth, -500, -200) ~ "200-500",
                                    between(depth, -1000, -500) ~ "500-1000",
                                    between(depth, -2000, -1000) ~ "1000-2000",
                                    between(depth, -3000, -2000) ~ "2000-3000",
                                    between(depth, -4000, -3000) ~ "3000-4000",
                                    between(depth, -Inf, -4000) ~ "4000-5500"),
                          levels = c(">0", "0-100", "100-200", "200-500", "500-1000",
                                     "1000-2000", "2000-3000", "3000-4000", "4000-5500")))

# Coastline shapefile
coast_10m_laea <- readOGR("data/bathy/ne_10m_land.shp", verbose = F) %>% # Coastline in laea
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  crop(extent(-180, 180, 0, 90)) %>% # Crop shapefile
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to a dataframe for ggplot
  rename(xc = long, yc = lat)
```

    ## Regions defined for each Polygons

``` r
# Location of trawls
trawl_laea <- SpatialPoints(cbind(trawl_latlon$lon, trawl_latlon$lat), 
                            proj4string = CRS("EPSG:4326")) %>%
  spTransform(., CRSobj = crs(arctic_laea)) %>% # Change projection to EPSG:6931
  as.data.frame() %>%
  rename(xc = coords.x1, yc = coords.x2)

# Location of acoustic data
MVBS_laea <- SpatialPoints(cbind(MVBS_latlon$lon, MVBS_latlon$lat), 
                           proj4string = CRS("EPSG:4326")) %>%
  spTransform(., CRSobj = crs(arctic_laea)) %>% # Change projection to EPSG:6931
  as.data.frame() %>%
  rename(xc = coords.x1, yc = coords.x2) # Rename variables

# Arctic circle and areas of interest
arctic_circle_laea <- SpatialPoints(cbind(arctic_circle_latlon$lon, arctic_circle_latlon$lat), 
                                    proj4string = CRS("EPSG:4326")) %>%
  spTransform(., CRSobj = crs(arctic_laea)) %>% # Change projection to EPSG:6931
  as.data.frame() %>%
  rename(xc = coords.x1, yc = coords.x2) # Rename variables
area_BF_CAA_laea <- SpatialPoints(cbind(c(-152, -152, -140, -122.5, -93, -93), c(78, 73, 68, 68, 68, 78)), 
                                  proj4string = CRS("EPSG:4326")) %>%
  spTransform(., CRSobj = crs(arctic_laea)) %>% # Change projection to EPSG:6931
  as.data.frame() %>%
  rename(xc = coords.x1, yc = coords.x2)
area_BB_laea <- SpatialPoints(cbind(area_BB_latlon$lon, area_BB_latlon$lat), 
                              proj4string = CRS("EPSG:4326")) %>%
  spTransform(., CRSobj = crs(arctic_laea)) %>% # Change projection to EPSG:6931
  as.data.frame() %>%
  rename(xc = coords.x1, yc = coords.x2)
area_SV_laea <- SpatialPoints(cbind(c(40, 40, 0, 0), c(85, 78, 78, 85)), 
                              proj4string = CRS("EPSG:4326")) %>%
  spTransform(., CRSobj = crs(arctic_laea)) %>% # Change projection to EPSG:6931
  as.data.frame() %>%
  rename(xc = coords.x1, yc = coords.x2)

# IHO regions
IHO_EAO <- readOGR("data/arctic_regions/iho_eastern_arctic_ocean.shp", verbose = F) %>% # Eastern Arctic Ocean
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to dataframe for ggplot
  rename(xc = long, yc = lat)
```

    ## Regions defined for each Polygons

``` r
IHO_WAO <- readOGR("data/arctic_regions/iho_western_arctic_ocean.shp", verbose = F) %>% # Western Arctic Ocean
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to dataframe for ggplot
  rename(xc = long, yc = lat)
```

    ## Regions defined for each Polygons

``` r
IHO_BF <- readOGR("data/arctic_regions/iho_beaufort_sea.shp", verbose = F) %>% # Eastern Arctic Ocean
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to dataframe for ggplot
  rename(xc = long, yc = lat)
```

    ## Regions defined for each Polygons

``` r
IHO_CAA <- readOGR("data/arctic_regions/iho_northwestern_passages.shp", verbose = F) %>% # Eastern Arctic Ocean
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to dataframe for ggplot
  rename(xc = long, yc = lat)
```

    ## Regions defined for each Polygons

``` r
IHO_BB <- readOGR("data/arctic_regions/iho_baffin_bay.shp", verbose = F) %>% # Eastern Arctic Ocean
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to dataframe for ggplot
  rename(xc = long, yc = lat)
```

    ## Regions defined for each Polygons

``` r
IHO_DS <- readOGR("data/arctic_regions/iho_davis_strait.shp", verbose = F) %>% # Eastern Arctic Ocean
  spTransform(CRSobj = crs(arctic_latlon)) %>% # Make sure that the shapefile is in the right projection
  spTransform(CRSobj = crs(arctic_laea)) %>% # Project shapefile in laea
  fortify() %>% # Convert to dataframe for ggplot
  rename(xc = long, yc = lat)
```

    ## Regions defined for each Polygons

``` r
IHO_regions <- bind_rows(IHO_EAO, IHO_WAO, IHO_BF, IHO_CAA, IHO_BB, IHO_DS) # Combine IHO definitions

# Remove unused variable
rm(cell_res)
```

# Plots

## Figure 1. Bathy map with sampling locations

Careful this map takes forever to compute (&gt; 18 h) when bathy has to
be drawn.

``` r
MVBS %>%
  mutate(lat = round(lat, 1),
         lon = round(lon, 1)) %>%
  group_by(lat,lon) %>%
  summarise(lat = mean(lat), 
            lon = mean(lon)) %>%
  ungroup() %>%
  ggplot(aes(x = lon, y = lat)) +
  # Plot bathy
  # geom_tile(data = bathy_df, aes(x = lon, y = lat, fill = depth_d)) +
  # scale_fill_cmocean("Depth (m)", name = "deep", discrete = T, na.value = NA, alpha = 0.6) +
  # Plot coastlines
  geom_polygon(data = coast_10m_latlon, aes(x = lon, y = lat, group = group), fill = "grey70") +
  # Arctic circle
  geom_segment(aes(x = -180, xend = 0, y = 66.5, yend = 66.5), col = "grey40", lty = 2, size = 0.3) +
  geom_segment(aes(x = 0, xend = 180, y = 66.5, yend = 66.5), col = "grey40", lty = 2, size = 0.3) +
  # Areas of interest
  geom_shape(data = area_BF_CAA_latlon, aes(x = lon, y = lat), fill = NA, col = "black", lwd = 0.2, lty = 1) +
  geom_shape(data = area_BB_latlon, aes(x = lon, y = lat), fill = NA, col = "black", lwd = 0.2, lty = 1) +
  geom_shape(data = area_SV_latlon, aes(x = lon, y = lat), fill = NA, col = "black", lwd = 0.2, lty = 1) +
  # Acoustic data
  geom_point(size = 0.7, shape = 4, color = "black") +
  # Trawl locations
  geom_point(data = trawl_latlon, aes(x = lon, y = lat), shape = 23, size = 1.4, color = "black", fill = "orange") +
  # Change projection
  coord_map("azequalarea", ylim = c(66, 90), xlim = c(-160, 70), orientation = c(90,0,-60)) +
  guides(fill = guide_legend(keywidth = 0.2, keyheight = 0.2, default.unit = "in"), color = "none") +
  theme(legend.position = "right", panel.border = element_rect(fill = NA), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
```

To overcome that, I project the bathy data on the EASE-Grid 2.0 North,
prior to plotting. This results in a much quicker plotting time.

``` r
# Bathymetric map of the Arctic with location of acoustic and trawl
p_MVBS_laea <- MVBS_laea %>%
  ggplot(aes(x = xc, y = yc)) +
  # Plot bathy
  geom_raster(data = bathy_laea, aes(x = xc, y = yc, fill = depth_d)) +
  scale_fill_cmocean("Depth (m)", name = "deep", discrete = T, na.value = NA, alpha = 0.6) +
  # Areas of interest
  geom_polygon(data = IHO_WAO, aes(x = xc, y = yc, group = group), fill = NA, col = "black", lwd = 0.3, lty = 1) +
  geom_shape(data = IHO_BF, aes(x = xc, y = yc, group = group), fill = NA, col = "black", lwd = 0.3, lty = 1) +
  geom_shape(data = IHO_CAA, aes(x = xc, y = yc, group = group), fill = NA, col = "black", lwd = 0.3, lty = 1) +
  geom_shape(data = IHO_BB, aes(x = xc, y = yc, group = group), fill = NA, col = "black", lwd = 0.3, lty = 1) +
  geom_shape(data = IHO_DS, aes(x = xc, y = yc, group = group), fill = NA, col = "black", lwd = 0.3, lty = 1) +
  geom_shape(data = IHO_EAO, aes(x = xc, y = yc, group = group), fill = NA, col = "black", lwd = 0.3, lty = 1) +
  # Plot coastlines
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey30") +
  # Arctic circle
  geom_path(data = arctic_circle_laea, aes(x = xc, y = yc), col = "grey10", lty = 2, size = 0.3) +
  # Acoustic data
  geom_point(size = 0.7, shape = 4, color = "black") +
  # Trawl locations
  geom_point(data = trawl_laea, aes(x = xc, y = yc), shape = 23, size = 1.4, color = "black", fill = "orange") +
  coord_fixed(xlim = c(-2700, 1500), ylim = c(-1900, 2000), expand = F) +
  guides(fill = guide_legend(keywidth = 0.2, keyheight = 0.2, default.unit = "in"), color = "none") +
  theme(legend.position = "right", panel.border = element_rect(fill = NA), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
# ggsave("plots/Fig1_map_stations.png", p_MVBS_laea, height = 4, width = 6.5, dpi = 600, units = "in")
p_MVBS_laea
```

<img src="PanArctic_DSL_figures_files/figure-gfm/fig1-map-station-laea-1.png" style="display: block; margin: auto;" />

## Figure 2. Vertical S<sub>V</sub> profiles

### Old plot

``` r
col_pal <- c("#80CBB1", "#347BA5", "#302346") # Custom colour palette

Sv_prof_mean <- Sv_grid_laea %>% # Plot median profiles
  group_by(depth, IHO_area) %>% # Calculate mean and median Sv profiles per area per year
  mutate(IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                     IHO_area == "West Arctic Ocean" ~ "WAO_BF",
                                     IHO_area == "Beaufort Sea" ~ "WAO_BF",
                                     IHO_area == "The Northwestern Passages" ~ "CAA",
                                     IHO_area == "Baffin Bay" ~ "BB",
                                     IHO_area == "Davis Strait" ~ "DS"),
                           levels = c("WAO_BF", "CAA", "BB", "DS", "EAO"))) %>%
  summarise(sv_lin_q0.05_areayear = median(sv_lin_q0.05),
            sv_lin_q0.50_areayear = median(sv_lin_q0.50),
            sv_lin_q0.95_areayear = median(sv_lin_q0.95)) %>%
  mutate(Sv_q0.05_median = 10 * log10(sv_lin_q0.05_areayear), 
         Sv_q0.50_median = 10 * log10(sv_lin_q0.50_areayear),
         Sv_q0.95_median = 10 * log10(sv_lin_q0.95_areayear))

p_Sv_grid_laea <- Sv_grid_laea %>% # Plot median profiles
  group_by(depth, IHO_area) %>% # Calculate mean and median Sv profiles per area per year
  mutate(IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "EAO",
                                     IHO_area == "West Arctic Ocean" ~ "WAO_BF",
                                     IHO_area == "Beaufort Sea" ~ "WAO_BF",
                                     IHO_area == "The Northwestern Passages" ~ "CAA",
                                     IHO_area == "Baffin Bay" ~ "BB",
                                     IHO_area == "Davis Strait" ~ "DS"),
                           levels = c("WAO_BF", "CAA", "BB", "DS", "EAO"))) %>%
  summarise(sv_lin_q0.05_areayear = median(sv_lin_q0.05),
            sv_lin_q0.25_areayear = median(sv_lin_q0.25),
            sv_lin_q0.50_areayear = median(sv_lin_q0.50),
            sv_lin_q0.75_areayear = median(sv_lin_q0.75),
            sv_lin_q0.95_areayear = median(sv_lin_q0.95)) %>%
  mutate(Sv_q0.05_median = 10 * log10(sv_lin_q0.05_areayear), 
         Sv_q0.25_median = 10 * log10(sv_lin_q0.25_areayear),
         Sv_q0.50_median = 10 * log10(sv_lin_q0.50_areayear),
         Sv_q0.25_median = 10 * log10(sv_lin_q0.25_areayear),
         Sv_q0.75_median = 10 * log10(sv_lin_q0.75_areayear),
         Sv_q0.95_median = 10 * log10(sv_lin_q0.95_areayear)) %>%
  ungroup() %>%
  # Plotting details for allowing geom_ribbon to work
  mutate(Sv_q0.05_median = case_when(Sv_q0.05_median == -999 ~ NaN,
                                     Sv_q0.05_median > -999 & Sv_q0.05_median <= -105 ~ -105,
                                     Sv_q0.05_median > -105 ~ Sv_q0.05_median),
         Sv_q0.25_median = case_when(Sv_q0.25_median == -999 ~ NaN,
                                     Sv_q0.25_median > -999 & Sv_q0.25_median <= -105 ~ -105,
                                     Sv_q0.25_median > -105 ~ Sv_q0.25_median),
         Sv_q0.50_median = case_when(Sv_q0.50_median == -999 ~ NaN,
                                     Sv_q0.50_median > -999 & Sv_q0.50_median <= -105 ~ -105,
                                     Sv_q0.50_median > -105 ~ Sv_q0.50_median),
         Sv_q0.75_median = case_when(Sv_q0.75_median == -999 ~ NaN,
                                     Sv_q0.75_median > -999 & Sv_q0.75_median <= -105 ~ -105,
                                     Sv_q0.75_median > -105 ~ Sv_q0.75_median),
         Sv_q0.95_median = case_when(Sv_q0.95_median == -999 ~ NaN,
                                     Sv_q0.95_median > -999 & Sv_q0.95_median <= -105 ~ -105,
                                     Sv_q0.95_median > -105 ~ Sv_q0.95_median)) %>%
  arrange(IHO_area, depth) %>%
  # Plot
  ggplot() +
  geom_vline(xintercept = 200, lty = 2, col = "grey20") +
  # geom_ribbon(aes(x = depth, ymin = Sv_q0.05_median, ymax = Sv_q0.95_median, fill = as.factor(year), group = year),
  #             alpha = 0.2, na.rm = T) +
  geom_ribbon(aes(x = depth, ymin = Sv_q0.05_median, ymax = Sv_q0.95_median), alpha = 0.2, na.rm = T) +
  geom_line(aes(x = depth, y = Sv_q0.50_median), na.rm = T) +
  # geom_line(aes(x = depth, y = Sv_q0.50_median, col = as.factor(year), group = year), na.rm = T) +
  scale_x_reverse("Depth (m)", breaks = seq(0, 1000, 200)) +
  geom_hline(yintercept = -105, col = "white") +
  scale_y_continuous(expression("S"[V]*" (dB re 1 m"^-1*")"), breaks = seq(-200, 0, 10), expand = c(0, 0)) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  coord_flip(ylim = c(-105, -65)) +
  facet_grid(~ IHO_area) +
  theme(legend.title = element_blank(),
        legend.position = c(0.94, 0.275),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line())
# ggsave("plots/Fig2_Sv_profiles_area_year.png", p_Sv_grid_laea, height = 3, width = 6.5, dpi = 600, units = "in")
p_Sv_grid_laea
```

<img src="PanArctic_DSL_figures_files/figure-gfm/fig2-Sv-profiles-area-year-1.png" style="display: block; margin: auto;" />

### Lines

``` r
# wesanderson::wes_palette("Darjeeling1", type = "discrete") 1st and 5th colours are inverted
col_pal <- c("#5BBCD6", "#00A08A", "#F2AD00", "#F98400", "#FF0000")

# Calculate median Sv profile for each grid cell
Sv_prof <- Sv_grid_laea %>% 
    mutate(IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "East Arctic Ocean",
                                     IHO_area == "West Arctic Ocean" ~ "West Arctic Ocean\n& Beaufort Sea",
                                     IHO_area == "Beaufort Sea" ~ "West Arctic Ocean\n& Beaufort Sea",
                                     IHO_area == "The Northwestern Passages" ~ "Canadian Arctic\nArchipelago",
                                     IHO_area == "Baffin Bay" ~ "Baffin Bay",
                                     IHO_area == "Davis Strait" ~ "Davis Strait"),
                           levels = c("West Arctic Ocean\n& Beaufort Sea", "Canadian Arctic\nArchipelago",
                                      "Baffin Bay", "Davis Strait", "East Arctic Ocean"))) %>%
  mutate(depth_round = rounder(depth, -10)) %>% # round down to nearest 10
  group_by(depth_round, IHO_area, xc, yc, year) %>%
  summarise(sv_lin50 = median(sv_lin_q0.50)) %>% 
  ungroup() %>%
  mutate(Sv_median = 10 * log10(sv_lin50),
         year = as.factor(year)) %>%
  unite(ID, xc, yc, year)

# Calculate median Sv profile 
Sv_prof_median <- Sv_grid_laea %>%
    mutate(IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "East Arctic Ocean",
                                     IHO_area == "West Arctic Ocean" ~ "West Arctic Ocean\n& Beaufort Sea",
                                     IHO_area == "Beaufort Sea" ~ "West Arctic Ocean\n& Beaufort Sea",
                                     IHO_area == "The Northwestern Passages" ~ "Canadian Arctic\nArchipelago",
                                     IHO_area == "Baffin Bay" ~ "Baffin Bay",
                                     IHO_area == "Davis Strait" ~ "Davis Strait"),
                           levels = c("West Arctic Ocean\n& Beaufort Sea", "Canadian Arctic\nArchipelago",
                                      "Baffin Bay", "Davis Strait", "East Arctic Ocean"))) %>%
  mutate(depth_round = rounder(depth, -10)) %>% # round down to nearest 20
  group_by(depth_round, IHO_area) %>%
  summarise(sv_lin50 = median(sv_lin_q0.50)) %>% 
  ungroup() %>%
  mutate(Sv_median = 10 * log10(sv_lin50))

# Plot
p_Sv_prof_lines <- Sv_prof_median %>%
  ggplot() +
  geom_line(data = Sv_prof, aes(x = depth_round, y = Sv_median, group = ID), col = "grey80", na.rm = T) +
  geom_vline(xintercept = 200, lty = 2) +
  geom_line(aes(x = depth_round, y = Sv_median, col = IHO_area), na.rm = T, size = 1) +
  scale_x_reverse("Depth (m)", breaks = seq(0, 1000, 100)) + 
  scale_y_continuous(expression("S"[V]*" (dB re 1 m"^-1*")"), breaks = seq(-200, 0, 10)) +
  coord_flip(ylim = c(-105, -70), xlim = c(800, 20)) +
  scale_colour_manual(values = col_pal) +
  facet_rep_grid(~ IHO_area) +
  guides(colour = "none") +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        panel.spacing = unit(0, "in"),
        plot.margin = unit(c(0, 0.1, 0, 0), "in"))
ggsave("plots/Fig2_Sv_profiles_lines.png", p_Sv_prof_lines, height = 4.5, width = 7, dpi = 600, units = "in")
p_Sv_prof_lines
```

![](PanArctic_DSL_figures_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### Points

``` r
p_Sv_prof_points <- Sv_grid_laea %>%
  mutate(IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "East Arctic Ocean",
                                     IHO_area == "West Arctic Ocean" ~ "West Arctic Ocean\n& Beaufort Sea",
                                     IHO_area == "Beaufort Sea" ~ "West Arctic Ocean\n& Beaufort Sea",
                                     IHO_area == "The Northwestern Passages" ~ "Canadian Arctic\nArchipelago",
                                     IHO_area == "Baffin Bay" ~ "Baffin Bay",
                                     IHO_area == "Davis Strait" ~ "Davis Strait"),
                           levels = c("West Arctic Ocean\n& Beaufort Sea", "Canadian Arctic\nArchipelago",
                                      "Baffin Bay", "Davis Strait", "East Arctic Ocean"))) %>%
  mutate(depth_round = rounder(depth, -20)) %>% # round down to nearest 20
  group_by(depth_round, IHO_area, year, xc, yc) %>%
  summarise(sv_lin = median(sv_lin_q0.50)) %>%
  ungroup() %>%
  mutate(Sv_median = 10*log10(sv_lin)) %>%
  filter(Sv_median != -999) %>%
  ggplot() +
  geom_vline(xintercept = 200, col = "red", lty = 2) +
  geom_pointrange(mapping = aes(y = Sv_median, x = depth_round, col = depth_round < 200),
                  stat = "summary", fun = median, fun.min = function(z) {quantile(z, 0.05)}, 
                  fun.max = function(z) {quantile(z, 0.95)}, size = 0.2) + 
  scale_x_reverse("Depth (m)", breaks = seq(0, 1000, 100)) + 
  scale_y_continuous(expression("S"[V]*" (dB re 1 m"^-1*")"), breaks = seq(-200, 0, 20)) +
  scale_colour_manual(values = c("black", "grey60")) + 
  coord_flip(ylim = c(-110, -60), xlim = c(800, 20)) +
  facet_rep_grid(~ IHO_area) +
  guides(colour = "none") +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        panel.spacing = unit(0, "in"),
        plot.margin = unit(c(0, 0.1, 0, 0), "in"))
ggsave("plots/Fig2_Sv_profiles_points.png", p_Sv_prof_points, height = 4.5, width = 7, dpi = 600, units = "in")
p_Sv_prof_points
```

![](PanArctic_DSL_figures_files/figure-gfm/Fig2-plot-1.png)<!-- -->

### Boxplots

``` r
p_Sv_prof_boxplot <- Sv_grid_laea %>%
  filter(depth <= 800) %>%
  mutate(IHO_area = factor(case_when(IHO_area == "East Arctic Ocean" ~ "East Arctic Ocean",
                                     IHO_area == "West Arctic Ocean" ~ "West Arctic Ocean\n& Beaufort Sea",
                                     IHO_area == "Beaufort Sea" ~ "West Arctic Ocean\n& Beaufort Sea",
                                     IHO_area == "The Northwestern Passages" ~ "Canadian Arctic\nArchipelago",
                                     IHO_area == "Baffin Bay" ~ "Baffin Bay",
                                     IHO_area == "Davis Strait" ~ "Davis Strait"),
                           levels = c("West Arctic Ocean\n& Beaufort Sea", "Canadian Arctic\nArchipelago",
                                      "Baffin Bay", "Davis Strait", "East Arctic Ocean"))) %>%
  mutate(depth_round = rounder(depth, -25)) %>% # round down to nearest 20
  group_by(depth_round, IHO_area, year, xc, yc) %>%
  summarise(sv_lin = median(sv_lin_q0.50)) %>%
  ungroup() %>%
  mutate(Sv_median = 10*log10(sv_lin), 
         depth = factor(depth_round)) %>%
  filter(Sv_median != -999) %>%
  ggplot() +
  geom_vline(xintercept = 200, col = "red", lty = 2) +
  geom_boxplot(aes(x = depth, y = Sv_median, col = depth_round < 200), 
               outlier.size = 0.5, size = 0.5) + 
  scale_x_discrete("Depth (m)", breaks = seq(0, 1000, 100), limits = rev) +
  scale_y_continuous(expression("S"[V]*" (dB re 1 m"^-1*")"), breaks = seq(-200, 0, 20)) +
  scale_colour_manual(values = c("black", "grey60")) +
  coord_flip(ylim = c(-110, -60)) +
  facet_rep_wrap(~ IHO_area, nrow = 1) +
  guides(colour = "none") +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        panel.spacing = unit(0, "in"),
        plot.margin = unit(c(0, 0.1, 0, 0), "in"))
ggsave("plots/Fig2_Sv_profiles_boxplots.png", p_Sv_prof_boxplot, height = 4.5, width = 6.5, dpi = 600, units = "in")
p_Sv_prof_boxplot
```

![](PanArctic_DSL_figures_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Figure 3. Map of mesopelagic backscatter hotspots

``` r
col_pal <- c("#80cdc1",  "#f5f5f5", "#dfc27d", "#a6611a")
# Calculate mean sa anomaly over three years of data per area
map_sa_anomaly <- SA_anomaly_2deg %>%
  group_by(lat, lon, area, region) %>%
  summarise(year_sampled = n(),
            mean_ano_sa_int = mean(anomaly_NASC_int, na.rm = T)) %>%
  ungroup() %>%
  mutate(mean_ano_sa_int_d = factor(case_when(mean_ano_sa_int < -1 ~ "<-1",
                                              mean_ano_sa_int >= -1 & mean_ano_sa_int < -0.5 ~ "[-1; -0.5[",
                                              mean_ano_sa_int >= -0.5 & mean_ano_sa_int <= 0.5 ~ "[-0.5; 0.5]",
                                              mean_ano_sa_int > 0.5 & mean_ano_sa_int <= 1 ~ "]0.5; 1]",
                                              mean_ano_sa_int > 1 ~ ">1"),
                                    levels = c("<-1", "[-1; -0.5[", "[-0.5; 0.5]", "]0.5; 1]", ">1"))) %>%
  ggplot(aes(x = lon, y = lat)) +
  # Longitude lines
  geom_segment(aes(x = -180, xend = -180, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = -150, xend = -150, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = -120, xend = -120, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = -90, xend = -90, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = -60, xend = -60, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = -30, xend = -30, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 0, xend = 0, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 180, xend = 180, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 150, xend = 150, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 120, xend = 120, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 90, xend = 90, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 60, xend = 60, y = 60, yend = 90), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 30, xend = 30, y = 60, yend = 90), col = "grey90", size = 0.3) +
  # Latitude lines
  geom_segment(aes(x = -180, xend = 0, y = 60, yend = 60), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 0, xend = 180, y = 60, yend = 60), col = "grey90", size = 0.3) +
  geom_segment(aes(x = -180, xend = 0, y = 70, yend = 70), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 0, xend = 180, y = 70, yend = 70), col = "grey90", size = 0.3) +
  geom_segment(aes(x = -180, xend = 0, y = 80, yend = 80), col = "grey90", size = 0.3) +
  geom_segment(aes(x = 0, xend = 180, y = 80, yend = 80), col = "grey90", size = 0.3) +
  # Plot 500m isobath
  geom_contour(data = bathy_df, aes(x = lon, y = lat, z = depth), breaks = -500, size = 0.25, colour = "black", lty = 3) +
  # Plot coastlines
  geom_polygon(data = coast_10m_latlon, aes(x = lon, y = lat, group = group), fill = "grey70") +
  # Areas of interest
  geom_shape(data = area_BF_CAA_latlon, aes(x = lon, y = lat), fill = NA, col = "black", lwd = 0.2, lty = 1) +
  geom_shape(data = area_BB_latlon, aes(x = lon, y = lat), fill = NA, col = "black", lwd = 0.2, lty = 1) +
  geom_shape(data = area_SV_latlon, aes(x = lon, y = lat), fill = NA, col = "black", lwd = 0.2, lty = 1) +
  # Arctic circle
  geom_segment(aes(x = -180, xend = 100, y = 66.5, yend = 66.5), col = "#666766", lty = 2, size = 0.3) +
  # Plot data
  geom_tile(aes(fill = mean_ano_sa_int_d), color = "grey30", alpha = 0.8, size = 0.15) +
  scale_fill_manual("2015-2017\nAverage\nbackscatter\nanomaly", values = col_pal) +
  coord_map("azequalarea", ylim = c(66, 90), xlim = c(-160, 30), orientation = c(90, 0, -60)) +
  # Theme
  guides(fill = guide_legend(keywidth = 0.2, keyheight = 0.2, default.unit = "in"), color = "none") +
  theme(legend.position = "right", 
        legend.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA),
        strip.text.x = element_blank(), 
        strip.background = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
map_sa_anomaly
```

    ## Warning: Removed 7197600 rows containing non-finite values (stat_contour).

<img src="PanArctic_DSL_figures_files/figure-gfm/fig3-map-sa-anomalies-1.png" style="display: block; margin: auto;" />

``` r
SA_grid_laea %>% 
  ggplot() +
  # Plot bathy
  geom_raster(data = bathy_laea, aes(x = xc, y = yc, fill = depth_d)) +
  scale_fill_cmocean("Depth (m)", name = "deep", discrete = T, na.value = NA, alpha = 0.3) +
  # scale_fill_brewer("Depth (m)", type = "seq", palette = "Greys") +
  guides(fill = "none") + # Remove legend
  new_scale_fill() +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group), fill = "grey30") +
  # Plot backscatter anomalies
  geom_tile(aes(x = xc, y = yc, fill = NASC_anomaly_d), color = "grey30", size = 0.15) +
  scale_fill_manual("Normalized sA anomaly", values = rev(brewer.pal(n = 5, name = "BrBG"))) +
  # scale_fill_manual("Normalized sA anomaly", values = rev(brewer.pal(n = 5, name = "RdBu"))) +
  facet_wrap(~ year, ncol = 3) +
  coord_fixed(xlim = c(-2600, 1100), ylim = c(-1800, 1900), expand = F) +
  guides(fill = guide_legend(keywidth = 0.2, keyheight = 0.2, default.unit = "in"), color = "none") +
  theme(legend.position = "right", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
```

<img src="PanArctic_DSL_figures_files/figure-gfm/fig3-map-sa-anomalies-laea-1.png" style="display: block; margin: auto;" />

## Table 1. Abundance mesopelagic fish

``` r
relative_abundance_area <- trawl_station %>%
  group_by(area) %>%
  summarise(total_fish = sum(number_fish)) %>% # Calculate total number of fish caught per area
  ungroup() %>%
  right_join(., trawl_area, by = c("area")) %>% # Join with relative abundance dataset
  mutate(area = factor(case_when(area == "BF_CAA" ~ "Beaufort Sea &\nCan. Arct. Archipelago",
                                 area == "BB" ~ "Baffin Bay",
                                 area == "SV" ~ "Svalbard"),
                       levels = c("Beaufort Sea &\nCan. Arct. Archipelago", "Baffin Bay", "Svalbard"))) %>%
  unite(area_label, area, total_fish, sep = "\nn = ", remove = F) %>% # Create labels for x axis
  # Format species name appropriately
  mutate(fish_species = str_replace(fish_species, pattern = "_", replacement = " "),
         fish_species = str_replace(fish_species, pattern = " sp", replacement = " sp."),
         fish_species = str_to_sentence(fish_species),
         fish_species = if_else(fish_species == "Notolepis rissoi", "Arctozenus risso", fish_species),
         fish_species = factor(fish_species, levels = c("Icelus bicornis",
                                                        "Anarhichas lupus",
                                                        "Reinhardtius hippoglossoides",
                                                        "Arctozenus risso",
                                                        "Melanogrammus aeglefinus",
                                                        "Gadus morhua",
                                                        "Leptoclinus maculatus",
                                                        "Sebastes sp.",
                                                        "Myctophidae",
                                                        "Liparidae",
                                                        "Boreogadus saida")),
         color_text = if_else(RA_area >= 40, "white", "black"),
         RA_area_round = 10 * floor(RA_area / 10),
         area_label = factor(area_label, levels = c("Beaufort Sea &\nCan. Arct. Archipelago\nn = 923",
                                                    "Baffin Bay\nn = 232",
                                                    "Svalbard\nn = 464"))) %>%
  # Plot
  ggplot() +
  geom_tile(aes(x = area_label, y = fish_species, fill = RA_area_round), color = "white", lwd = 0.6) +
  geom_text(aes(x = area_label, y = fish_species, label = sprintf("%0.1f", round(RA_area, 1)), color = color_text), size = 3) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(labels = c("Twohorn sculpin",
                              "Atlantic wolffish",
                              "Greenland halibut",
                              "Spotted barracudina",
                              "Haddock",
                              "Atlantic cod",
                              "Daubed shanny",
                              "Redfish",
                              "Lanternfish",
                              "Snailfish",
                              "Polar cod"),
                   expand = c(0, 0)) +
  scale_colour_identity() +
  scale_fill_viridis_c("Relative\nabundance (%)  ", option = "mako", breaks = seq(0, 100, 20), limits = c(0, 100), direction = -1) +
  guides(fill = guide_colorbar(label.position = "top")) +
  theme(axis.title = element_blank(),
        legend.position = "top",
        legend.title = element_text(vjust = 0),
        legend.title.align = 0.5,
        legend.key.height = unit(0.1, "in"),
        legend.key.width = unit(0.3, "in"),
        legend.box.margin = margin(0,0,-5,0))
relative_abundance_area
```

<img src="PanArctic_DSL_figures_files/figure-gfm/table1-heatmap-abundance-1.png" style="display: block; margin: auto;" />

# Supplementary data

## Table S3. Acoustic settings

<table class=" lightable-classic" style="font-family: Arial; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Year
</th>
<th style="text-align:left;">
Area
</th>
<th style="text-align:left;">
Vessel
</th>
<th style="text-align:center;">
Echosounder
</th>
<th style="text-align:center;">
Transducers depth (m)
</th>
<th style="text-align:center;">
Pulse length (msec)
</th>
<th style="text-align:center;">
Power (W)
</th>
<th style="text-align:center;">
Acoustic data duration (h)
</th>
</tr>
</thead>
<tbody>
<tr grouplength="3">
<td colspan="8" style="border-bottom: 0;">
<strong>2015</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Beaufort Sea
</td>
<td style="text-align:left;">
CCGV Amundsen
</td>
<td style="text-align:center;">
EK60
</td>
<td style="text-align:center;">
7
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
2000
</td>
<td style="text-align:center;">
117
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Baffin Bay
</td>
<td style="text-align:left;">
CCGV Amundsen
</td>
<td style="text-align:center;">
EK60
</td>
<td style="text-align:center;">
7
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
2000
</td>
<td style="text-align:center;">
70
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Svalbard
</td>
<td style="text-align:left;">
RV Polarstern
</td>
<td style="text-align:center;">
EK60
</td>
<td style="text-align:center;">
11
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
1000
</td>
<td style="text-align:center;">
230
</td>
</tr>
<tr grouplength="3">
<td colspan="8" style="border-bottom: 0;">
<strong>2016</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Beaufort Sea
</td>
<td style="text-align:left;">
CCGV Amundsen
</td>
<td style="text-align:center;">
EK60
</td>
<td style="text-align:center;">
7
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
2000
</td>
<td style="text-align:center;">
48
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Baffin Bay
</td>
<td style="text-align:left;">
CCGV Amundsen
</td>
<td style="text-align:center;">
EK60
</td>
<td style="text-align:center;">
7
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
2000
</td>
<td style="text-align:center;">
145
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Svalbard
</td>
<td style="text-align:left;">
RV Helmer Hanssen
</td>
<td style="text-align:center;">
EK60
</td>
<td style="text-align:center;">
8
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
2000
</td>
<td style="text-align:center;">
18
</td>
</tr>
<tr grouplength="3">
<td colspan="8" style="border-bottom: 0;">
<strong>2017</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Beaufort Sea
</td>
<td style="text-align:left;">
FV Frosti
</td>
<td style="text-align:center;">
EK80
</td>
<td style="text-align:center;">
6
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
2000
</td>
<td style="text-align:center;">
119
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Baffin Bay
</td>
<td style="text-align:left;">
CCGV Amundsen
</td>
<td style="text-align:center;">
EK60
</td>
<td style="text-align:center;">
7
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
2000
</td>
<td style="text-align:center;">
115
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
</td>
<td style="text-align:left;">
Svalbard
</td>
<td style="text-align:left;">
RV Polarstern
</td>
<td style="text-align:center;">
EK60
</td>
<td style="text-align:center;">
11
</td>
<td style="text-align:center;">
1.024
</td>
<td style="text-align:center;">
1000
</td>
<td style="text-align:center;">
463
</td>
</tr>
</tbody>
</table>
