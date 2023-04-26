library(magrittr)
library(dplyr)
library(terra)
library(ggplot2)

my_dir_path <- "/srv/duo/mbelitz/CentraliAmericana/"
wd <- list()
wd$R       <- file.path( my_dir_path, "R" )
wd$fun     <- file.path( my_dir_path, "functions" )
wd$bin     <- file.path( my_dir_path, "bin" )
wd$preds   <- file.path( my_dir_path, "ClimateOnly/" )
wd$data    <- file.path( my_dir_path, "data" )
wd$occs    <- file.path( wd$data, "BiologiaCentraliAmericanaGeoreferencedSpecies.csv")
wd$figs    <- file.path( my_dir_path, "figs" )
wd$out     <- file.path( my_dir_path, "out" )
invisible({ lapply(wd, function(i) if( dir.exists(i) != 1 ) dir.create(i) ) })

proj4_aea <- "+proj=laea +lon_0=-69.61 +lat_0=-10.36 +datum=WGS84 +units=m +no_defs"
proj4_wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"