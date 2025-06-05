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
1. Download **bioclimatic variables** for Madagascar from WorldClim (current) and CMIP6 (future).
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
### 0. Load Libraries
``` r
library(terra)          # for raster data handling
library(virtualspecies) # for generating virtual species with response curves
library(geodata)        # for downloading climate and administrative data
library(sf)             # for vector data and spatial features
library(stars)          # for creating multidimensional spatio-temporal cubes
library(ggplot2)        # for visualization
library(tidyverse)      # for general data manipulation
library(viridis)        # visualization
```

### 1. Download and Crop Bioclimatic Data for Madagascar
Bioclimatic variables represent biologically meaningful climate summaries
- bio1 = annual mean temperature
- bio12 = annual precipitation
- bio15 = precipitation seasonality

``` r
# download present and future climate data
BIO_current <- worldclim_country("MG", var = "bio", res = 5, path = tempdir())
BIO_future_global <- cmip6_world(
  model = "CNRM-CM6-1", ssp = "585", time = "2061-2080",
  var = "bio", res = 5, version = "2.1", path = tempdir()
)

# crop by country
BIO_future <- crop(BIO_future_global, BIO_current)

# select 3 relevant variables
BIO_current_sub <- BIO_current[[c("wc2.1_30s_bio_1", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15")]]
BIO_future_sub  <- BIO_future[[c("bio01", "bio12", "bio15")]]

# rename for simplicity
names(BIO_current_sub) <- names(BIO_future_sub) <- c("bio1", "bio12", "bio15")
```
### 2. Create virtual species 
Let's create 3 virtual species: each of them has their own response functions to the variables. 

```
# Each species has a Gaussian response curve to each variable.
# These represent:
# - Species 1: Temperate generalist
# - Species 2: Warm-adapted species
# - Species 3: Cool-adapted, dry-specialist

curve_list <- list(
  formatFunctions(
    # species 1
    bio1 = c(fun = 'dnorm', mean = 15, sd = 2),
    bio12 = c(fun = 'dnorm', mean = 800, sd = 200),
    bio15 = c(fun = 'dnorm', mean = 60, sd = 10)
  ),
  formatFunctions(
    # species 2
    bio1 = c(fun = 'dnorm', mean = 22, sd = 3),
    bio12 = c(fun = 'dnorm', mean = 1200, sd = 300),
    bio15 = c(fun = 'dnorm', mean = 80, sd = 10)
  ),
  formatFunctions(
    # species 3
    bio1 = c(fun = 'dnorm', mean = 10, sd = 2),
    bio12 = c(fun = 'dnorm', mean = 400, sd = 100),
    bio15 = c(fun = 'dnorm', mean = 40, sd = 8)
  )
)
```
### 3. Generate Virtual Species Suitability (Current & Future)
Taking as input the bioclimatic variables we downloaded before and the response functions, `generateSpFromFun` creates the suitability map for each of them in the present and in the future. 
``` r
species_current <- list()
species_future  <- list()

for (i in 1:3) {
  species_current[[paste0("species_", i)]] <- generateSpFromFun(
    raster.stack = BIO_current_sub,
    parameters = curve_list[[i]],
    rescale = FALSE
  )
  species_future[[paste0("species_", i)]] <- generateSpFromFun(
    raster.stack = BIO_future_sub,
    parameters = curve_list[[i]],
    rescale = FALSE
  )
}
```
### 4. Create Regular Grid Over Madagascar 
We want to incorporate the 3 species into the same object.

To do so, we need to collapse the x and y spatial dimensions into one, creating a grid that represents individual cells within the study area. This grid allows us to reduce the dimensions and facilitates the transition to a vector data cube.

``` r
madagascar_shape <- gadm("Madagascar", level = 0, path = tempdir())

# Create a hexagonal grid with 0.9° spacing
mad_grid <- st_make_grid(
  madagascar_shape,
  cellsize = 0.7,
  what = "polygons",
  square = FALSE
) %>%
  st_as_sf() %>%
  mutate(cell_id = 1:n())
```
### 5. Aggregate Suitability Over Spatial Grid 

Species' suitability maps are grouped into two stars objects: the first is related to the present, the second to the future. Then, pixels in the suitability rasters are grouped based on their spatial intersection with a geometry (the grid). Each group is then reduced to a single value using an aggregation function, such as the mean or maximum. The result is a one-dimensional sequence of feature geometries defined in space.

``` r
suitability_aggregated <- list()

for (i in 1:3) {
  sp <- paste0("species_", i)
  
  # convert terra raster to stars object
  r_2020 <- st_as_stars(species_current[[sp]]$suitab.raster)
  r_2070 <- st_as_stars(species_future[[sp]]$suitab.raster)
  
  # aggregate mean suitability over each cell of the grid
  suitability_aggregated[[paste0(sp, "_2020")]] <- aggregate(r_2020, mad_grid, FUN = mean, na.rm = TRUE)
  suitability_aggregated[[paste0(sp, "_2070")]] <- aggregate(r_2070, mad_grid, FUN = mean, na.rm = TRUE)
}
```
### 6. Build Suitability Data Cube (stars)
Now that we have the suitability aggregated over polygons, we recombine into two data cubes (present and future). Finally, species dimension and time dimension are set. We have one object that contains everything we need to explore suitability changes in community. 
```  r
# combine species into stars objects for each time period
current_cube <- c(
  suitability_aggregated$species_1_2020,
  suitability_aggregated$species_2_2020,
  suitability_aggregated$species_3_2020
)

future_cube <- c(
  suitability_aggregated$species_1_2070,
  suitability_aggregated$species_2_2070,
  suitability_aggregated$species_3_2070
)

# Add "species" dimension
current_cube <- st_redimension(current_cube) %>%
  st_set_dimensions(2, values = paste0("species_", 1:3), names = "species")

future_cube <- st_redimension(future_cube) %>%
  st_set_dimensions(2, values = paste0("species_", 1:3), names = "species")

# Combine along the time dimension
suitability_cube <- c(current_cube, future_cube, along = list(time = c(2020, 2070)))

# Set final dimension names
suitability_cube <- st_set_dimensions(suitability_cube, 1, names = "cell")
names(suitability_cube) <- "suitability"
```
### 7. Explore and Visualize the Data Cube
 We have one object that contains everything we need to explore suitability changes in community. 
``` r
print(suitability_cube)

ggplot() +
  geom_stars(data = suitability_cube) +
  facet_grid(time ~ species) +
  scale_fill_viridis_c(name = "Suitability", na.value = "white") +
  theme_minimal() +
  labs(title = "Suitability of Virtual Species in Madagascar",
       subtitle = "Under Present (2020) and Future (2070) Climate Conditions")

```
