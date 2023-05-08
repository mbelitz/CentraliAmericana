# Set up pipeline ---------------------------------------------------------

## Step 1: Load in the necessary libraries ----
source('R/00_setup.R')
library(rnaturalearth)
library(rgeos)
library(dismo)
library(ENMeval)
library(stringr)
library(dplyr)
library(terra)
library(sf)

## Step 2: Load in the pipeline scripts ----
list.files(wd$fun, full.names = T) %>% 
  lapply(source)

## Step 3: Prepare the spatial data ----

# Load predictors
#cl <- list.files("ClimateOnly/", full.names = T)
#mod_vars <- rast(cl)
#mod_vars <- aggregate(mod_vars, 5) %>% 
#  terra::project(proj4_aea)
### 3C. load occurrence records
occs <- data.table::fread("data/BiologiaCentraliAmericanaGeoreferencedSpecies.csv") %>% 
  dplyr::rename(decimalLongitude = Longitude, decimalLatitude = Latitude) %>% 
  dplyr::filter(!is.na(decimalLongitude),
                !is.na(decimalLatitude)) %>% 
  dplyr::mutate(bw = Species)

### 3D. Read in basemap for visualizing
world <- ne_countries(scale = 110, returnclass = "sf")

### 3E. Choose projection and project data if necessary
study_proj <- proj4_aea #aea projection
world <- st_transform(world, crs = study_proj) 

occs_sf <- st_as_sf(occs, 
                    coords = c("decimalLongitude", "decimalLatitude"),
                    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 
occs_sf <- st_transform(occs_sf, crs = study_proj)
occs <- occs %>% 
  mutate(x = st_coordinates(occs_sf)[,1], y = st_coordinates(occs_sf)[,2])


############################ Part 3: Execute the pipeline #############################

rerunPipeline <- FALSE

#### Step 1: Set up species list
spp_list <- unique(occs$bw)

#### Step 2: Pass pipeline through the species_season list

sdm_pipeline <- function(x){
  
  print("=========================================================")
  print(paste("============== Working on", spp_list[x], "=============="))  
  print("=========================================================")
  binomial <- spp_list[x]
  species_df_raw <-   dplyr::filter(as.data.frame(occs), bw == binomial) %>% 
    filter(!is.na(x),
           !is.na(y))
  stopifnot(nrow(species_df_raw) >= 1)
  
  # Clip all potential predictor variables based on all cleaned occurrence
  # records from this species.
  
  aa_shp <- define_accessibleArea(species_df = species_df_raw, minBuff = 200e3,
                                  buff_prop = 0.80, projCRS = study_proj)
  
  mod_vars <- rast("climateStack_centraliAmericana.tif")
  # Clip environmental variable layers to the defined accessible area
  mymod_vars <- clip_variableLayers(rstack = mod_vars, accessibleArea = aa_shp)
  
  # Thin points based on the accessible area.
  ### 2A. Prepare the coordinates for rarefaction
  coordinates(species_df_raw) <- ~ x + y
  
  ### 2B. Perform the rarefaction
#  area_sqkm <- raster::area(aa_shp)*0.000001 
#  species_df <- thinPoints(
#    spp_df = species_df_raw, 
#    area_sqkm = area_sqkm, 
#    bio = mymod_vars[[2]], 
#    method = "complex",
#    simpleMult = 5
#  )
spp_df <- species_df_raw
  
  #### Step 3: Test and fine-tune the model
  
  ### 3A. Select top performing variables to reduce colinearity
  ## First run a test model
  print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
  max_model <- maxent(x = mymod_vars, p = coordinates(spp_df), progress = "text", silent = TRUE) 
  ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
  print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
  predictors <- select_sdmVariables(pred_vars = mymod_vars, maxent_mod = max_model, maxVIF = 5)
  
  ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
  print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
  
  eval1 <- ENMeval::ENMevaluate(
    occ = coordinates(spp_df), env = predictors,
    method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
    fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
    parallel = TRUE, numCores = 10, algorithm = 'maxent.jar',
    quiet = TRUE, updateProgress = TRUE
  )
  
  #### Step 4: Create the output for the best performing model
  
  ### 4A. Prepare the output path and coordinates
  bw <- spp_df$bw[1] # make a speciesName string without a space for better saving
  bw <- str_replace(bw, " ", "_")
  
  resultDir <-  wd$out
  if (!dir.exists(resultDir)) { dir.create(resultDir) }
  
  ### 4B. Output the model
  ## Output the model evaluation
  #  save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
  ## Output the best performing model, including both the actual model and the presence-absence model
  print(">>>>>>>>> Saving best model <<<<<<<<<")
  
 # sp_df2 <- as.data.frame(coordinates(species_df)) %>% 
  #  dplyr::rename(x = 1, y = 2)
  
  sp_df2 <- as.data.frame(coordinates(species_df_raw)) %>% 
    dplyr::rename(x = 1, y = 2)
  
  save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = wd$out,
                   spp = bw, occ_df = sp_df2)
  
  ### 4C. Visualize the model
  ## Load in the rasters
  r    <- raster(file.path(resultDir, paste0(bw,"_SDM.tif")))
  r_pa <- raster(file.path(resultDir, paste0(bw,"_SDM_PA.tif")))
  
  create_sdmFigure(
    spp = bw, r = r, r_pa = r_pa, occ_df = sp_df2, world = world, resultDir = resultDir, bw = bw
  )
  
  print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
  
  ## now let's classify rasters by a bunch of different thresholding values
  
  lpt_occs <- sp_df2 %>% 
    dplyr::select(x, y)
  lpt <- terra::extract(x = r, y = lpt_occs) 
  
  ## LPT 0  
  lptVal <- 0
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_000 <- classify(rast(r), pa_mat)
  
  ## LPT 0.01
  lptVal <- 0.01
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_001 <- classify(rast(r), pa_mat)
  
  ## LPT 0.025 
  lptVal <- 0.025
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_0025 <- classify(rast(r), pa_mat)
  
  ## LPT 0.05 
  lptVal <- 0.05
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_005 <- classify(rast(r), pa_mat)
  
  ## LPT 0.1  
  lptVal <- 0.1
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_01 <- classify(rast(r), pa_mat)
  
  ## LPT 0.15 
  lptVal <- 0.15
  lpt_10 <- quantile(lpt, probs = lptVal, na.rm = TRUE)
  pa <- c(0, lpt_10, 0, lpt_10, 1, 1)
  pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
  r_best_pa_015 <- classify(rast(r), pa_mat)
  
  
  p0 <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r_best_pa_000,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
              alpha = 0.8) +
    geom_point(sp_df2, mapping = aes(x = x, y = y),
               size = 0.5, alpha = 0.15, shape = 1) +
    coord_sf(xlim = c(min(sp_df2$x) - 1000000, max(sp_df2$x) + 1000000),
             ylim = c(min(sp_df2$y) - 1000000, max(sp_df2$y) + 1000000)) +
    labs(fill = "Presence") +
    scale_fill_viridis_d() +
    ggtitle("LPT 0%") +
    theme_classic()
  
  p001 <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r_best_pa_001,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
              alpha = 0.8) +
    geom_point(sp_df2, mapping = aes(x = x, y = y),
               size = 0.5, alpha = 0.15, shape = 1) +
    coord_sf(xlim = c(min(sp_df2$x) - 1000000, max(sp_df2$x) + 1000000),
             ylim = c(min(sp_df2$y) - 1000000, max(sp_df2$y) + 1000000)) +
    labs(fill = "Presence") +
    scale_fill_viridis_d() +
    ggtitle("LPT 1%") +
    theme_classic()
  
  p0025 <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r_best_pa_0025,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
              alpha = 0.8) +
    geom_point(sp_df2, mapping = aes(x = x, y = y),
               size = 0.5, alpha = 0.15, shape = 1) +
    coord_sf(xlim = c(min(sp_df2$x) - 1000000, max(sp_df2$x) + 1000000),
             ylim = c(min(sp_df2$y) - 1000000, max(sp_df2$y) + 1000000)) +
    labs(fill = "Presence") +
    scale_fill_viridis_d() +
    ggtitle("LPT 2.5%") +
    theme_classic()
  
  p005 <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r_best_pa_005,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
              alpha = 0.8) +
    geom_point(sp_df2, mapping = aes(x = x, y = y),
               size = 0.5, alpha = 0.15, shape = 1) +
    coord_sf(xlim = c(min(sp_df2$x) - 1000000, max(sp_df2$x) + 1000000),
             ylim = c(min(sp_df2$y) - 1000000, max(sp_df2$y) + 1000000)) +
    labs(fill = "Presence") +
    scale_fill_viridis_d() +
    ggtitle("LPT 5%") +
    theme_classic()
  
  p01 <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r_best_pa_01,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
              alpha = 0.8) +
    geom_point(sp_df2, mapping = aes(x = x, y = y),
               size = 0.5, alpha = 0.15, shape = 1) +
    coord_sf(xlim = c(min(sp_df2$x) - 1000000, max(sp_df2$x) + 1000000),
             ylim = c(min(sp_df2$y) - 1000000, max(sp_df2$y) + 1000000)) +
    labs(fill = "Presence") +
    scale_fill_viridis_d() +
    ggtitle("LPT 10%") +
    theme_classic()
  
  p015 <- ggplot() +
    geom_sf(world, mapping = aes()) +
    geom_tile(as.data.frame(r_best_pa_015,xy = T) %>% na.omit() %>% rename(ClogLog = 3),
              mapping = aes(x = x, y = y, fill = as.character(ClogLog)),
              alpha = 0.8) +
    geom_point(sp_df2, mapping = aes(x = x, y = y),
               size = 0.5, alpha = 0.15, shape = 1) +
    coord_sf(xlim = c(min(sp_df2$x) - 1000000, max(sp_df2$x) + 1000000),
             ylim = c(min(sp_df2$y) - 1000000, max(sp_df2$y) + 1000000)) +
    labs(fill = "Presence") +
    scale_fill_viridis_d() +
    ggtitle("LPT 15%") +
    theme_classic()
  
  # p <- ggplot() +
  #   geom_sf(world, mapping = aes()) +
  #   geom_point(sp_df2, mapping = aes(x = x, y = y),
  #              size = 0.5, alpha = 0.8, shape = 1) +
  #   coord_sf(xlim = c(min(sp_df2$x) - 1000000, max(sp_df2$x) + 1000000),
  #            ylim = c(min(sp_df2$y) - 1000000, max(sp_df2$y) + 1000000)) +
  #   labs(fill = "Presence") +
  #   scale_fill_viridis_d() +
  #   ggtitle("points only") +
  #   theme_classic()
  
  
  cp <- cowplot::plot_grid(p0, p001, p0025, p005, p01, p015, 
                           nrow = 3,
                           labels = c("A", "B", "C",
                                      "D", "E", "F"))
  
  ggsave(plot = cp, filename = paste0("thresholdingOut/",bw,"_thresh.png"),
         width = 8, height = 6)
  
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  unixtools::set.tempdir("tmp/")
  tmp_dir <- unixtools::set.tempdir("tmp/")
  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()
  
}

sdm_pipeline(39)


# looping through pipeline
# looping through pipeline
for(i in 39:304){
  tryCatch(sdm_pipeline(x = i),
           error = function(e) print(paste(spp_list[i], "Error in Code Skipping for now!")))
}  

