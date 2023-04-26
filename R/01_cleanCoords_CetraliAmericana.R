library(dplyr)
library(ggplot2)
library(CoordinateCleaner)
library(rnaturalearth)
library(sf)
library(stringr)
library(dbplyr)
library(RSQLite)
library(ggpubr)

# Function to clean coordinates from occurrence records.
# Input is a dataframe containing occurrence records for one species only (of name "spp").
map_coord_clean <- function(dataframe, spp){
  
  # remove NAs n decimalLongitude & Latitude and make sure decimalLong & Lats are distinct
  cs <- dataframe %>% 
    filter(scientificName == spp) %>% 
    filter(decimalLongitude >= -180 & decimalLongitude <= 180) %>% 
    filter(decimalLatitude >= -90 & decimalLatitude <= 90) %>% 
    distinct(decimalLongitude, decimalLatitude, .keep_all = T)
  
  cs2 <- clean_coordinates(cs, lon = "decimalLongitude", lat = "decimalLatitude",
                           species = "scientificName",
                           tests = c("capitals", "centroids", "equal",
                                     "gbif", "institutions", "outliers"),
                           verbose = T) %>% 
    dplyr::mutate(.summary = factor(.summary, levels = c(FALSE, TRUE)))
  
  mytab <- cs2 %>% 
    tidyr::pivot_longer(cols =c(".val", ".equ", ".cap", ".cen", ".otl", ".gbf", ".inst") ) %>%
    group_by(name, value) %>% 
    dplyr::summarise(n=n()) %>% 
    dplyr::filter(value == FALSE | name == ".gbf")
  
  bw <- str_replace(spp, " ", "_")
  write.csv(x = cs2, file = file.path(wd$phil_coords, paste0(bw, ".csv")))
  
  lonRange <- (max(cs$decimalLongitude)- min(cs$decimalLongitude))/10
  latRange <- (max(cs$decimalLatitude)- min(cs$decimalLatitude))/10
  
  a <- ggplot() +
    geom_sf(world, mapping = aes(), fill = "grey95", color = "grey80") +
    geom_point(cs, mapping = aes(x = decimalLongitude, y = decimalLatitude),
               color = "#18806d") + 
    coord_sf(xlim = c(min(cs$decimalLongitude) - lonRange, max(cs$decimalLongitude) + lonRange),
             ylim = c(min(cs$decimalLatitude) - latRange, max(cs$decimalLatitude) + latRange)) +
    ggtitle("Uncleaned") + 
    theme_bw() +
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
  
  b <- ggplot() +
    geom_sf(world, mapping = aes(), fill = "#f5f5f5", color = "grey80") +
    geom_point(
      dplyr::filter(cs2, .summary == "TRUE"), 
      mapping = aes(x = decimalLongitude, y = decimalLatitude),
      color = "black")  +
    geom_point(
      dplyr::filter(cs2, .summary == "FALSE"), 
      mapping = aes(x = decimalLongitude, y = decimalLatitude),
      shape = 21, fill = "yellow", color = "#786213")  +
    coord_sf(xlim = c(min(cs$decimalLongitude) - lonRange, max(cs$decimalLongitude) + lonRange),
             ylim = c(min(cs$decimalLatitude) - latRange, max(cs$decimalLatitude) + latRange)) +
    labs(color = "Not flagged") +
    ggtitle("Cleaned") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
  cp <- ggarrange(a, b, ggtexttable(mytab, rows = NULL),  nrow = 3, ncol = 1, heights = c(2,2,1))
  
  ggsave(plot = cp, 
         filename = file.path(wd$phil_coords, paste0(bw, ".png")),
         width = 4, height = 8)
}

quickPlot <- function(cs) {
  lonRange <- (max(cs$decimalLongitude)- min(cs$decimalLongitude))/10
  latRange <- (max(cs$decimalLatitude)- min(cs$decimalLatitude))/10
  
  ggplot() +
    geom_sf(world, mapping = aes(), fill = "grey95", color = "grey80") +
    geom_point(cs, mapping = aes(x = decimalLongitude, y = decimalLatitude),
               color = "#18806d") + 
    coord_sf(xlim = c(min(cs$decimalLongitude) - lonRange, max(cs$decimalLongitude) + lonRange),
             ylim = c(min(cs$decimalLatitude) - latRange, max(cs$decimalLatitude) + latRange)) +
    ggtitle("Uncleaned") + 
    theme_bw() +
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
}