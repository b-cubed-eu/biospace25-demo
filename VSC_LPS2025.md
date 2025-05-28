# Virtual Suitability Data Cube Tutorial

Species suitability refers to **how favorable** an environment is for a species to survive, reproduce, and grow in a specific area and time. It considers factors like climate, landscape, and resource availability.

**Species Distribution Models** (SDMs) are tools that use environmental and species occurrence data to study and predict the distribution of species across time and space. SDMs help identify **suitable habitats**, forecast the movements of invasive species, and illustrate how species distributions might change due to climate change. They are essential for conservation planning and biodiversity protection.

Studying species suitability under different environmental conditions is crucial for understanding population dynamics, planning conservation actions, and monitoring the effects of climate change and human activities.

To observe species suitability over multiple species and time periods, we developed a framework using **data cubes** — multidimensional arrays that organize data in a structured way. In **R**, the [`stars`](https://r-spatial.github.io/stars/) package provides a convenient way to represent these cubes.

### What are `stars` objects?
`stars` objects are multidimensional arrays with spatial (and optionally temporal) structure. They allow for:
- **Slicing** across dimensions (e.g., extracting a time step or species)
- **Aggregation** (e.g., averaging suitability over time or space)
- **Visualization** using base R or `ggplot2`

### Tutorial Overview
In this tutorial, we outline the steps to create a `stars` object with three dimensions: time, space (represented as grid cells), and species, with suitability as the main attribute.

We:
1. Download **bioclimatic variables** for South Africa from WorldClim (current) and CMIP6 (future).
2. Select three relevant variables: annual mean temperature (`bio1`), annual precipitation (`bio12`), and precipitation seasonality (`bio15`).
3. Define response curves for **three virtual species**, each with distinct ecological preferences.
4. Generate **suitability maps** for each species under both present (2020) and future (2070) climate conditions.
5. Create a **spatial grid** (as polygons) over the study area.
6. Aggregate the suitability values over this grid.
7. Combine all layers into a **`stars` data cube** with dimensions:
   - **cell**: spatial polygons (grid)
   - **species**: species_1, species_2, species_3
   - **time**: 2020, 2070

This structured data format makes it easy to observe changes in species suitability over space and time, and lays the foundation for more advanced modeling workflows (e.g., validation, ensemble modeling, thresholding).

The code implementation follows this exact logic, using `terra`, `virtualspecies`, `sf`, and `stars`. Virtual species are particularly useful for this demonstration because their ecological preferences are defined explicitly, making them ideal for controlled experiments or teaching. Each virtual species is associated with a response to climate variables, and we use these responses to compute suitability scores across the landscape.

The final output is a unified, three-dimensional `stars` object that is well suited for exploratory analysis, visualization, and integration into biodiversity decision-support tools.

---
## 0. Load packages
```r
# load required packages
library(terra)           # for raster operations
library(virtualspecies)  # for generating virtual species
library(geodata)         # for downloading bioclimatic data
library(sf)              # for spatial vector operations
library(stars)           # for creating data cubes
library(ggplot2)         # for plotting
```

## 1. Import Climate Variables

Thanks to the functions `worldclim_country` and `cmip6_world` from the `geodata` R package, we can download bioclimatic variables for both the present and the future

```r
# download current and future climate data for South Africa
BIO_current <- worldclim_country("ZAF", var = "bio", res = 10, path = tempdir())
BIO_future_global <- cmip6_world(model = "CNRM-CM6-1", ssp = "585", time = "2061-2080", 
                                 var = "bio", res = 10, version = "2.1", path = tempdir())
# crop future to current extent
BIO_future <- crop(BIO_future_global, BIO_current)
```

## 2. Select Bioclimatic Variables

We select just a subset of the variables, as an example. 
```r
# select variables of interest: bio1 (temperature), bio12 (precipitation), bio15 (seasonality)
BIO_current <- BIO_current[[c("wc2.1_30s_bio_1", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15")]]
BIO_future  <- BIO_future[[c("bio01", "bio12", "bio15")]]

# rename layers for consistency
names(BIO_current) <- names(BIO_future) <- c("bio1", "bio12", "bio15")

# plot one of the climate variables for inspection
plot(BIO_current[["bio1"]])
```

## 3. Define Response Functions and Generate Virtual Species
Virtual species are... and their response functions to the environment can be customized to create a suitability map that we will use... 
Each species respond differently to each environmental factor, according to their requirements. Here, we simulate 3 types of species... 
```r
# define Gaussian response curves for 3 virtual species and for each variable
curve_list <- list(
  formatFunctions(bio1 = c(fun = 'dnorm', mean = 15, sd = 2),  
                  bio12 = c(fun = 'dnorm', mean = 800, sd = 200),  
                  bio15 = c(fun = 'dnorm', mean = 60, sd = 10)),
  formatFunctions(bio1 = c(fun = 'dnorm', mean = 22, sd = 3),  
                  bio12 = c(fun = 'dnorm', mean = 1200, sd = 300), 
                  bio15 = c(fun = 'dnorm', mean = 80, sd = 10)),
  formatFunctions(bio1 = c(fun = 'dnorm', mean = 10, sd = 2),  
                  bio12 = c(fun = 'dnorm', mean = 400, sd = 100),  
                  bio15 = c(fun = 'dnorm', mean = 40, sd = 8))
)

# generate virtual species under current and future climate
species_current <- list()
species_future  <- list()
for (i in 1:3) {
  species_current[[paste0("species_", i)]] <- generateSpFromFun(raster.stack = BIO_current, parameters = curve_list[[i]], rescale = FALSE)
  species_future[[paste0("species_", i)]]  <- generateSpFromFun(raster.stack = BIO_future,  parameters = curve_list[[i]], rescale = FALSE)
}
```

## 4. Create Spatial Grid
The aim is to merge all the infos in one object that contains everything we need to know about the suitability of our species in space. A stars object can contain n dimensions, but for simplicity our datacube will have 3 dimensions: space as cells, time and species, and one attribute: the suitability. 
First, grid. 

```r
# create a regular polygon grid over the area
bbox <- st_bbox(BIO_current[[1]])
sf_bbox <- st_as_sfc(bbox) |> st_sf()
lux_grid <- st_make_grid(sf_bbox, cellsize = 0.9, what = "polygons") |> st_as_sf() |> mutate(cell_id = 1:n())
```

## 5. Aggregate Suitability over Grid
Takes each suitability raster and calculate the mean over polygons

```r
# aggregate suitability values per grid cell for each species and scenario
dataset <- list()
for (i in 1:3) {
  sp <- paste0("species_", i)
  r_2020 <- st_as_stars(species_current[[sp]]$suitab.raster)
  r_2070 <- st_as_stars(species_future[[sp]]$suitab.raster)
  dataset[[paste0(sp, "_2020")]] <- aggregate(r_2020, lux_grid, FUN = mean, na.rm = TRUE)
  dataset[[paste0(sp, "_2070")]] <- aggregate(r_2070, lux_grid, FUN = mean, na.rm = TRUE)
}
```

## 6. Combine into a Stars Data Cube
Merge into one object
I need 
```r
# Build stars cube with dimensions: cell x species x time
current_cube <- c(dataset$species_1_2020, dataset$species_2_2020, dataset$species_3_2020)
future_cube  <- c(dataset$species_1_2070, dataset$species_2_2070, dataset$species_3_2070)

current_cube <- st_redimension(current_cube)
current_cube <- st_set_dimensions(current_cube, 2, values = paste0("species_", 1:3), names = "species")

future_cube <- st_redimension(future_cube)
future_cube <- st_set_dimensions(future_cube, 2, values = paste0("species_", 1:3), names = "species")

suitability_cube <- c(current_cube, future_cube, along = list(time = c(2020, 2070)))
suitability_cube <- st_set_dimensions(suitability_cube, 1, names = "cell")
names(suitability_cube) <- "suitability"
```

## 7. Visualize the Suitability Cube

```r
ggplot() +
  geom_stars(data = suitability_cube) +
  facet_grid(time ~ species) +
  scale_fill_viridis_c() +
  theme_minimal()
```

This concludes the step-by-step tutorial for building a suitability cube with virtual species using climate data and spatial aggregation.


