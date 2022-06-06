---
title: "PanArctic DSL - Statistics"
author: "[Pierre Priou](mailto:pierre.priou@mi.mun.ca)"
date: "2022/06/06 at 19:11"
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
There are two missing chl values, so this cell is excluded from further analyses.


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
## 2 BF       24
## 3 NAR       3
## 4 BB       32
## 5 SV       12
```
Because we sampled at the border of the CAO and there is only a limited amount of data in this region we combine CAO data based on their proximity to other LME (Barents Sea and Beaufort Sea). 


```r
stat_laea <- stat_laea %>%
  mutate(LME = as.character(LME),
         LME = factor(if_else(LME == "CAO" & yc <= 0, "SV",
                      if_else(LME == "CAO" & yc > 0, "BF",
                      if_else(LME == "NAR", "BB", LME))),
                      levels = c("BF", "BB", "SV", "NAR"))) %>%
  filter(is.na(v) == F) %>%
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
## 1 BF       24
## 2 BB       34
## 3 SV       14
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

![](PanArctic_DSL_statistics_files/figure-html/map-chl-NA-1.png)<!-- -->

Check the total number of observations per group.



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
  filter(depth == 318) %>%
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = v_mean * 100)) +
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

![](PanArctic_DSL_statistics_files/figure-html/map-velocity318-1.png)<!-- -->


```r
# Velocity
phy_grid_laea %>% 
  filter(depth == 380) %>%
  ggplot(aes(x = xc,  y = yc)) +
  geom_tile(aes(fill = v_mean * 100)) +
  geom_polygon(data = coast_10m_laea, aes(x = xc, y = yc, group = group),
               fill = "grey80") +
  scale_fill_cmocean("velo (cm/s)", name = "speed", limits = c(0, 10)) +
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


```r
SA_diff <- SA_grid_laea %>%
  mutate(LME = case_when(LME == "Beaufort Sea LME" ~ "BF",
                         LME == "Central Arctic LME" ~ "CAO",
                         LME == "Baffin Bay LME" ~ "BB",
                         LME == "Barents Sea LME" ~ "SV",
                         LME == "Northern Canadian Archipelago LME" ~
                           "NAR"),
         LME = factor(if_else(LME == "CAO" & yc <= 0, "SV",
                      if_else(LME == "CAO" & yc > 0, "BF",
                      if_else(LME == "NAR", "BB", LME))),
                      levels = c("BF", "NAR", "BB", "SV")),
         year = factor(year)) %>%
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
  mutate(var = c("LME", "year", "LME", "year"))
```

```
## # A tibble: 4 × 7
##   .y.          n statistic    df       p method         var  
##   <chr>    <int>     <dbl> <int>   <dbl> <chr>          <chr>
## 1 NASC_int    74      5.54     2 0.0627  Kruskal-Wallis LME  
## 2 NASC_int    74     12.2      2 0.00224 Kruskal-Wallis year 
## 3 CM          74      8.90     2 0.0117  Kruskal-Wallis LME  
## 4 CM          74      3.74     2 0.154   Kruskal-Wallis year
```

Mesopelagic backscatter S\~A\~ varied significantly between areas (*H* = 29.792952, *p* \< 0.001) and years (*H* = 29.792952, *p* = 0.011). The center of mass varied significantly between areas (*H* = 61.185994, *p* \< 0.001) but not between years (*H* = 4.164743, *p* = 0.125).

## Within areas

I check whether there are inter-annual variability in S\~A\~ across years within areas. This will be used to decide whether a random effect is needed in the GAM. I used a non parametric Kruskal-Wallis test to test for interannual variability within areas.


```r
SA_diff %>% # Kruskal wallis interannual difference SA within group
  group_by(LME) %>%
  kruskal_test(NASC_int ~ year)
```

```
## # A tibble: 3 × 7
##   LME   .y.          n statistic    df       p method        
## * <fct> <chr>    <int>     <dbl> <int>   <dbl> <chr>         
## 1 BF    NASC_int    25     14.6      2 0.00066 Kruskal-Wallis
## 2 BB    NASC_int    35      2.13     2 0.344   Kruskal-Wallis
## 3 SV    NASC_int    14      3.14     2 0.208   Kruskal-Wallis
```

There were some interannual differences in SA within some region (WAO_BF - H = 6.4670868, *p* = 0.0394; CAA - H = 11.8260406, *p* = 0.0027; BB - H = 5.0981791, *p* = 0.0782; DS - H = 0.3282051, *p* = 0.8490; EAO - H = 4.9779498, *p* = 0.0830).

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
                      if_else(LME == "CAO" & yc > 0, "BF",
                      if_else(LME == "NAR", "BB", LME))),
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
##      Min       1Q   Median       3Q      Max 
## -18.9777  -4.2336   0.2723   4.2158  29.8301 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  48.5998    14.7680   3.291  0.00155 **
## lat          -0.4294     0.1982  -2.167  0.03355 * 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 7.83 on 72 degrees of freedom
## Multiple R-squared:  0.06123,	Adjusted R-squared:  0.04819 
## F-statistic: 4.696 on 1 and 72 DF,  p-value: 0.03355
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
  geom_point(aes(x = lat, y = SA_int, col = LME)) +
  geom_line(aes(x = lat, y = .fitted, col = LME), size = 0.8) + 
  geom_ribbon(aes(x = lat, ymin = .lower, ymax = .upper, group = LME), 
              alpha = 0.1)
```

![](PanArctic_DSL_statistics_files/figure-html/lm-plot-area-1.png)<!-- -->

Save data


```r
# Save data
save(LM_area_coefs, LM_area_fit, LM_area_res, 
     file = "data/statistics/LM_results.RData")
```

# HGAM - Gamma NASC

## Data preparation

I select data at 318 m depth because it fits well with the DSL diurnal centre of mass. It also prevent from removing too much data like the deeper depth levels would do.


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
<div id="htmlwidget-f783d5cb8c98f9c22d29" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f783d5cb8c98f9c22d29">{"x":{"filter":"none","vertical":false,"data":[["GAM_I","GAM_S","GAM_G"],[13.606,15.622,9.143],[68.64,68.41,57.97],[0.14,0.12,0.09],[417.44,420.74,423.97],[824.296,829.389,840.511],[0,5.09,16.21],[0.92708,0.07264,0.00028]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   4.9944     0.9945   5.022 4.28e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##           edf Ref.df      F  p-value    
## s(v)    2.653  3.175  6.752 0.000425 ***
## s(LME)  1.889  2.000 16.345 1.34e-05 ***
## s(year) 1.794  2.000 10.069 0.000662 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0878   Deviance explained =   58%
## -ML = 423.97  Scale est. = 2.1958    n = 72
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
## (Intercept)   4.7756     0.7384   6.467 2.09e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df      F  p-value    
## s(year)  1.739      2  8.041  0.00108 ** 
## s(LME)   1.845      2 13.010 4.76e-05 ***
## s(v,LME) 8.008     14 12.036  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.116   Deviance explained = 68.4%
## -ML = 420.74  Scale est. = 1.5144    n = 72
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
## (Intercept)   4.7222     0.8078   5.846 2.16e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df      F  p-value    
## s(year)    1.758  2.000  9.285 0.000663 ***
## s(LME)     1.878  2.000 16.070 1.76e-05 ***
## s(v):LMEBF 1.000  1.000 10.254 0.002168 ** 
## s(v):LMEBB 2.477  2.915  5.272 0.003215 ** 
## s(v):LMESV 3.305  3.736 16.132  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.137   Deviance explained = 68.6%
## -ML = 417.44  Scale est. = 1.512     n = 72
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
## s(year)   3 1.738795       NA      NA
## s(LME)    3 1.845495       NA      NA
## s(v,LME) 15 8.007890 0.939657    0.68
```


```r
appraise(GAM_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAM_I)
```

```
##            k'      edf  k-index p-value
## s(year)     3 1.757820       NA      NA
## s(LME)      3 1.877991       NA      NA
## s(v):LMEBF  4 1.000089 0.953206  0.6750
## s(v):LMEBB  4 2.476970 0.953206  0.7225
## s(v):LMESV  4 3.305120 0.953206  0.7025
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
## (Intercept)   4.7780     0.8383     5.7 3.99e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df      F  p-value    
## s(year)  1.777      2  8.445 0.000998 ***
## s(LME)   1.886      2 13.620 4.51e-05 ***
## s(v,LME) 8.008     14 14.591  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.122   Deviance explained = 68.5%
## -REML = 420.09  Scale est. = 1.4856    n = 72
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
## (Intercept)   4.7780     0.8383     5.7 3.99e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df      F  p-value    
## s(year)  1.777      2  8.445 0.000998 ***
## s(LME)   1.886      2 13.620 4.51e-05 ***
## s(v,LME) 8.008     14 14.591  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.122   Deviance explained = 68.5%
## -REML = 420.09  Scale est. = 1.4856    n = 72
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
## (Intercept)   4.8254     0.8605   5.608 5.12e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df      F  p-value    
## s(year)    1.7863      2 11.433 0.000386 ***
## s(LME)     1.8988      2 15.203 3.59e-05 ***
## s(v):LMEBF 0.8941      4  3.290 0.003118 ** 
## s(v):LMEBB 1.5791      4  1.982 0.037814 *  
## s(v):LMESV 3.0665      4 38.772  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.174   Deviance explained = 66.2%
## -REML = 420.17  Scale est. = 1.4814    n = 72
```

No terms can be dropped, so I refit model with REML.


```r
GAM_I2 <- gam(NASC_int ~ 
                 s(year, bs = "re") +
                 s(LME, bs = "re") +
                 s(v, by = LME, k = 5, bs = "tp"),
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
## (Intercept)   4.7017     0.9233   5.092 3.74e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df      F  p-value    
## s(year)    1.790  2.000  9.732 0.000557 ***
## s(LME)     1.911  2.000 17.452 1.28e-05 ***
## s(v):LMEBF 1.001  1.001 10.361 0.002075 ** 
## s(v):LMEBB 2.725  3.157  5.986 0.001056 ** 
## s(v):LMESV 3.331  3.754 16.481  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.137   Deviance explained = 69.1%
## -REML = 414.76  Scale est. = 1.4779    n = 72
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

## Model predictions

I predict the model to get the smooths in the response scale (NASC). I need to create a specific data frame for this with the factors of interest `LME` and `year`, and covariates `v` and `chl`.


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
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
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
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
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

I select data at 318 m depth because it fits well with the DSL diurnal centre of mass. It also prevent from removing too much data like the deeper depth levels would do.


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
<div id="htmlwidget-e9f4f862c84a8c0a6300" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-e9f4f862c84a8c0a6300">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_S","GAMSA_G"],[11.225,12.768,8.403],[45.68,46.21,32.36],[0.39,0.39,0.27],[241.51,244.29,247.2],[483.682,486.071,493.834],[0,2.39,10.15],[0.76388,0.23135,0.00477]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
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
## (Intercept)   16.377      1.731   9.462 6.72e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F p-value   
## s(v)    3.0009  3.499 3.822 0.00776 **
## s(LME)  0.6445  2.000 0.593 0.17302   
## s(year) 1.4500  2.000 3.684 0.00934 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.271   Deviance explained = 32.4%
## -ML =  247.2  Scale est. = 48.23     n = 72
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
## (Intercept)   16.099      1.737   9.269 2.34e-13 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##             edf Ref.df     F  p-value    
## s(year)  1.5420      2 4.728 0.003714 ** 
## s(LME)   0.5607      2 0.474 0.193777    
## s(v,LME) 6.2653     14 2.318 0.000156 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.39   Deviance explained = 46.2%
## -ML = 244.29  Scale est. = 40.36     n = 72
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
## (Intercept)   15.647      1.779   8.797 1.39e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F p-value   
## s(year)    1.5578  2.000 4.884 0.00373 **
## s(LME)     0.6841  2.000 0.623 0.18150   
## s(v):LMEBF 1.0000  1.000 5.557 0.02152 * 
## s(v):LMEBB 2.3860  2.821 6.203 0.00129 **
## s(v):LMESV 2.0124  2.400 3.597 0.03671 * 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.391   Deviance explained = 45.7%
## -ML = 241.51  Scale est. = 40.285    n = 72
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
##          k'       edf  k-index p-value
## s(year)   3 1.5420335       NA      NA
## s(LME)    3 0.5606755       NA      NA
## s(v,LME) 15 6.2653111 1.106975    0.79
```


```r
appraise(GAMSA_I, method = "simulate")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-basis-size-residuals-1.png)<!-- -->

```r
k.check(GAMSA_I)
```

```
##            k'       edf k-index p-value
## s(year)     3 1.5577669      NA      NA
## s(LME)      3 0.6840837      NA      NA
## s(v):LMEBF  4 1.0000310 1.10725  0.7850
## s(v):LMEBB  4 2.3859926 1.10725  0.7625
## s(v):LMESV  4 2.0123517 1.10725  0.8250
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
          nrow = 1)
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
          nrow = 1)
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
par(mfrow = c(1, 2))
acf(resid(GAMSA_I), lag.max = 20, main = "ACF")
pacf(resid(GAMSA_I), lag.max = 20, main = "pACF")
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI-acf-pacf-1.png)<!-- -->

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
## (Intercept)   16.094      2.092   7.691 1.32e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##             edf Ref.df     F  p-value    
## s(year)  1.6844      2 5.221 0.003764 ** 
## s(LME)   0.7297      2 0.604 0.210499    
## s(v,LME) 6.2535     14 2.383 0.000208 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.394   Deviance explained = 46.8%
## -REML = 242.73  Scale est. = 40.135    n = 72
```

`s(LME)` can be dropped. Refit the model without those terms.


```r
# Model S
GAMSA_S2 <- gam(SA_int ~ 
                  s(year, bs = "re") +
                  # s(LME, bs = "re") +
                  s(v, LME, bs = "fs", k = 5),
                data = SA_df, family = "gaussian", method = "ML")
summary(GAMSA_S2)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(v, LME, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.079      1.725   9.321  1.9e-13 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F  p-value    
## s(year)  1.541      2 4.634 0.003497 ** 
## s(v,LME) 6.778     14 2.387 0.000133 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.389   Deviance explained = 46.1%
## -ML = 244.31  Scale est. = 40.416    n = 72
```

Compare models.


```r
# Summary metrics
GAMSAS_AIC <- AIC(GAMSA_S, GAMSA_S2) %>% 
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAMSA_S", "GAMSA_S2"),
           reml = round(c(GAMSA_S$gcv.ubre,
                          GAMSA_S2$gcv.ubre), 
                        2), 
           dev_expl = round(c(
             (1 - (GAMSA_S$deviance / GAMSA_S$null.deviance)) * 100,
             (1 - (GAMSA_S2$deviance / GAMSA_S2$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAMSA_S)$r.sq,
                        summary(GAMSA_S2)$r.sq),
                      2)) %>%
  full_join(., GAMSAS_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(AIC) %>% 
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-5a7ccb0d1ac917c01fc5" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-5a7ccb0d1ac917c01fc5">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_S","GAMSA_S2"],[12.768,12.696],[46.21,46.09],[0.39,0.39],[244.29,244.31],[486.071,486.083],[0,0.01],[0.5015,0.4985]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAMSA_S3 <- gam(SA_int ~ 
                  s(year, bs = "re") +
                  # s(LME, bs = "re") +
                  s(v, LME, bs = "fs", k = 5),
                data = SA_df, family = "gaussian", method = "REML")
summary(GAMSA_S3)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(v, LME, bs = "fs", k = 5)
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   16.066      2.082   7.717 1.18e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##            edf Ref.df     F p-value    
## s(year)  1.685      2 5.063 0.00350 ** 
## s(v,LME) 6.929     14 2.474 0.00018 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.393   Deviance explained = 46.6%
## -REML = 242.75  Scale est. = 40.192    n = 72
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
## (Intercept)   15.788      2.088   7.561 1.94e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##               edf Ref.df     F  p-value    
## s(year)    1.6760      2 5.170 0.004267 ** 
## s(LME)     0.8685      2 0.816 0.187609    
## s(v):LMEBF 0.8152      4 1.368 0.023201 *  
## s(v):LMEBB 2.0852      4 4.616 0.000353 ***
## s(v):LMESV 1.7142      4 2.088 0.014377 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.385   Deviance explained = 44.7%
## -REML = 242.27  Scale est. = 40.675    n = 72
```

I remove `s(LME)` from the model.


```r
# Model I
GAMSA_I2 <- gam(SA_int ~
                  s(year, bs = "re") +
                  # s(LME, bs = "re") +
                  s(v, by = LME, k = 5, bs = "tp"),
                data = SA_df, family = "gaussian", method = "ML")
summary(GAMSA_I2)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(v, by = LME, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   15.547      1.668   9.322 1.54e-13 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df     F  p-value    
## s(year)    1.556  2.000 4.882 0.003108 ** 
## s(v):LMEBF 1.000  1.000 5.187 0.026105 *  
## s(v):LMEBB 2.414  2.860 6.363 0.000944 ***
## s(v):LMESV 1.947  2.331 3.355 0.044310 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.378   Deviance explained = 43.9%
## -ML = 241.62  Scale est. = 41.16     n = 72
```

Compare models


```r
# # Summary metrics
GAMSAI_AIC <- AIC(GAMSA_I, GAMSA_I2) %>%
  rownames_to_column() %>%
  rename(model = rowname)
# Metrics data frame
data.frame(model = c("GAMSA_I", "GAMSA_I2"),
           reml = round(c(GAMSA_I$gcv.ubre,
                          GAMSA_I2$gcv.ubre),
                        2),
           dev_expl = round(c(
             (1 - (GAMSA_I$deviance / GAMSA_I$null.deviance)) * 100,
             (1 - (GAMSA_I2$deviance / GAMSA_I2$null.deviance)) * 100),
             2),
           r2 = round(c(summary(GAMSA_I)$r.sq,
                        summary(GAMSA_I2)$r.sq),
                      2)) %>%
  full_join(., GAMSAI_AIC, by = "model") %>%
  mutate(df = round(df, 3),
         AIC = round(AIC, 3),
         dAIC = round(AIC - min(AIC), 2),
         w_AIC = round(Weights(AIC), 5)) %>%
  dplyr::select(model, df, dev_expl, r2, reml, AIC, dAIC, w_AIC) %>%
  arrange(AIC) %>%
  datatable(class = "cell-border stribe", rownames = F)
```

```{=html}
<div id="htmlwidget-5a1e2a857f6b265bc474" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-5a1e2a857f6b265bc474">{"x":{"filter":"none","vertical":false,"data":[["GAMSA_I","GAMSA_I2"],[11.225,10.089],[45.68,43.87],[0.39,0.38],[241.51,241.62],[483.682,483.776],[0,0.09],[0.51175,0.48825]],"container":"<table class=\"cell-border stribe\">\n  <thead>\n    <tr>\n      <th>model<\/th>\n      <th>df<\/th>\n      <th>dev_expl<\/th>\n      <th>r2<\/th>\n      <th>reml<\/th>\n      <th>AIC<\/th>\n      <th>dAIC<\/th>\n      <th>w_AIC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Refit model with REML


```r
GAMSA_I3 <- gam(SA_int ~
                  s(year, bs = "re") +
                  # s(LME, bs = "re") +
                  s(v, by = LME, k = 5, bs = "tp"),
                data = SA_df, family = "gaussian", method = "REML")
summary(GAMSA_I3)
```

```
## 
## Family: gaussian 
## Link function: identity 
## 
## Formula:
## SA_int ~ s(year, bs = "re") + s(v, by = LME, k = 5, bs = "tp")
## 
## Parametric coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   15.530      2.018   7.695 1.16e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##              edf Ref.df     F  p-value    
## s(year)    1.698  2.000 5.349 0.003020 ** 
## s(v):LMEBF 1.000  1.000 5.047 0.028110 *  
## s(v):LMEBB 2.571  3.021 6.288 0.000815 ***
## s(v):LMESV 2.208  2.630 3.398 0.040791 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.383   Deviance explained = 44.8%
## -REML = 234.18  Scale est. = 40.854    n = 72
```

Check residuals.


```r
appraise(GAMSA_I3)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAI3-residuals-covariates-1.png)<!-- -->

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

## Model predictions

I predict the model to get the smooths in the response scale (NASC). I need to create a specific data frame for this with the factors of interest `LME` and `year`, and covariates `v` and `chl`.


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
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
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
       x = expression("Current velocity at 318 m (cm s"^-1*")"),
       y = expression("S"[A]*" (dB re 1 m"^2*" nmi"^-2*")")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        legend.title = element_blank(), 
        legend.position = "top")
# Combine plots
plot_grid(GAMSAS_velo, GAMSAI_velo, ncol = 1)
```

![](PanArctic_DSL_statistics_files/figure-html/GAMSAS-pretty-ggplot-1.png)<!-- -->
