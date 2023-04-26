library(data.table)
library(ggplot2)
library(rnaturalearth)
source("R/00_setup.R00_setup.R")

# read in data
afro <- data.table::fread("data/150_species_test_run_from_ODA.csv")

# tidy afro
afro <- afro %>% 
  rename(scientificName = Species,
         decimalLongitude = Longitude,
         decimalLatitude = Latitude,
         Year = Year,
         Locality = Locality
  ) %>% 
  select(scientificName, 
         decimalLongitude, 
         decimalLatitude,
         Year,
         Locality)


## combine into not cleaned occs data frame
nco <- afro %>% 
  filter(!is.na(decimalLongitude)) %>% 
  filter(!is.na(decimalLatitude))

## get mapping things
world <- ne_countries(scale = "medium", returnclass = "sf")

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