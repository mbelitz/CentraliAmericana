library(data.table)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
source("R/00_setup.R00_setup.R")

# read in data
df <- data.table::fread("data/BiologiaCentraliAmericanaGeoreferencedSpecies.csv")

# tidy afro
df <- df %>% 
  rename(scientificName = Species,
         decimalLongitude = Longitude,
         decimalLatitude = Latitude,
         Year = Year,
         Locality = Locality
  ) %>% 
  select(scientificName, 
         decimalLongitude, 
         decimalLatitude,
         Locality)


## combine into not cleaned occs data frame
nco <- afro %>% 
  filter(!is.na(decimalLongitude)) %>% 
  filter(!is.na(decimalLatitude))

## get mapping things
world <- ne_countries(scale = 10, returnclass = "sf")

# list of species
spp_list <- unique(nco$scientificName)


#map_coord_clean(dataframe = nco, spp = spp_list[3])

for(i in 1:length(spp_list)) {
  spp <- spp_list[[i]]
  print(paste("Working on", spp))
  tryCatch(
    map_coord_clean(dataframe = nco, spp),
    error = function(e) print(paste(spp[i], "Error, did not run, skipping"))
  )
}


# what species have we done already?
completed <- list.files(path = "AustroAsiaProject/bin/phil_coords/", pattern = ".png")
completed <- completed %>% 
  str_remove(".png") %>% 
  str_replace(pattern = "_", replacement = " ")

spp_list <- data.frame(spp = spp_list) %>% 
  filter(!spp_list %in% completed)

spp_list <-  spp_list$spp