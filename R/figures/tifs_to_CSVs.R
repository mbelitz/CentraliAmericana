library(dplyr)
library(terra)
library(stringr)
library(data.table)
library(sf)

source("R/00_setup.R")

options(java.parameters = "-Xmx20000m")
raster::rasterOptions(tmpdir="tmp")
tmp_dir <- unixtools::set.tempdir("tmp/")
files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
file.remove(files)
gc()


# get list of species with SDMs PA
r_list <- list.files(path = "out", pattern = "*_PA.tif")

## Step 2: Load in the pipeline scripts ----
list.files(wd$fun, full.names = T) %>% 
  lapply(source)

## Step 3: Prepare the spatial data ----
# Load predictors
biorast <- rast("ClimateOnly/bio1.tif")
biorast <- aggregate(biorast, 5) %>% 
  terra::project(proj4_aea)

biorast_df <- terra::as.data.frame(biorast, xy = T)%>% 
  mutate(z = 1:nrow(.)) 

biorast_reg <- biorast_df %>% 
  dplyr::select(x,y,z) %>% 
  raster::rasterFromXYZ(.) 

fp <- file.path("out/")

#' function to genereate endemism tifs per species
sdm_to_csv <- function(x){
  
  bw <- str_replace(x, pattern = " ", replacement = "_")
  r <- rast(paste(fp,bw,"_SDM_PA.tif", sep = ""))
    
    rdf <- raster::as.data.frame(r, xy = T) %>% 
      na.omit() %>% 
      dplyr::rename(PA = 3)
    
    #regularlize
    spdf <- rdf
    sp::coordinates(spdf) <- ~ x + y
    
    e <- raster::extract(biorast_reg, spdf)
    
    rdf <- dplyr::mutate(rdf, z = e)
    
    lj <- left_join(rdf, biorast_df, by = "z") %>% 
      dplyr::filter(!is.na(z))
    
    lj <- lj %>% 
      dplyr::rename(x = x.y, y = y.y) %>% 
      dplyr::mutate(binomial = bw) %>% 
      dplyr::select(x, y, z, PA, binomial) %>% 
      dplyr::filter(PA > 0)
    
    fwrite(lj, 
           file = paste0("SDM_as_CSVs/", 
                         bw, ".csv"))
    
    print(x)
  
}

# looping through pipeline
for(i in seq_along(spp_list)){
  tryCatch(sdm_to_csv(x = spp_list[i]),
           error = function(e) print(paste(spp_list[i], "Error in Code Skipping for now!")))
}  



## read these all in and save as single CSV
ls <- list.files(path = "SDM_as_CSVs/", 
                 full.names = T)
library(tidyverse)
tbl_fread <- ls %>% 
  map_df(~read_csv(., col_types = cols(.default = "c"))) 
tbl_fread$x <- as.numeric(tbl_fread$x)
tbl_fread$y <- as.numeric(tbl_fread$y)
tbl_fread$PA <- as.numeric(tbl_fread$PA)

head(tbl_fread)

gb <- tbl_fread %>% 
  group_by(x,y) %>% 
  summarise(count = n())

ggplot(gb) +
  geom_tile(mapping = aes(x = x, y = y, fill = count, color = count))


gb2 <- tbl_fread %>% 
  group_by(x,y) %>% 
  summarise(richness = length(unique(binomial)))


library(sf)
proj4_aea <- "+proj=laea +lon_0=17.58 +lat_0=0 +datum=WGS84 +units=m +no_defs"

# read in data
world <- rnaturalearth::ne_countries(continent = "Africa", scale = 110, returnclass = 'sf') %>% 
  st_transform(crs = proj4_aea)

ggplot(gb2) +
  geom_sf(world, mapping = aes(), fill = NA) +
  geom_tile(mapping = aes(x = x, y = y, fill = richness, color = richness), alpha = 0.8) +
  scale_fill_viridis_c(option = "inferno") +
  scale_color_viridis_c(option = "inferno") + 
  coord_sf(xlim = c(-5000000,5000000),
            ylim = c(-4000000, 4500000)) +
  theme_void() 

ggsave(filename = "figs/richness.png")



