## feems mix of samples 
{rm(list = ls())
library(dplyr)
library(sf)
library(raster)
library(sp)
library(fields)
library(nngeo)}

setwd("path/to/working/directory")

# Load all the shapefiles by locality
feems_files <- list.files(pattern = "\\.shp$", full.names = T)
feems_files_names <- list.files(pattern = "\\.shp$", full.names = F)
feems_files_names <- gsub(".bestCV.shp", "", feems_files_names)

list_feems <- vector('list', length(feems_files))
names(list_feems) <- feems_files_names
for (i in 1:length(list_feems)) {
  list_feems[[i]] <- st_read(feems_files[i])
}
# Add the inverse logarithm information
for (i in 1:length(list_feems)) {
  list_feems[[i]]$logW <- 10^(list_feems[[i]]$feems)
}
## Exctract all non gecko species
list_geckos <- list_feems[c(-7,-10,-11,-14)]
#3. Extract the geometry of a dataset with all the desired edges (in this case is the first one)
geometry <- st_geometry(list_geckos[[1]])
geometry <- st_as_sf(geometry) # convert it to a spatial object
# Craete a variable with the names of each species to append afterwards
names_df_geckos <- c('geometry', names(list_geckos))



##################################### PERCENTAGE ####################################

# standardise all values from each analysis to have each species edge value in a range of 0-1
# 1. Convert negative values to positive values by adding the lowest number of speciesX$feems
# 2. create a percentage dividing all craeted values by the max(dt$feems_sum)
## Loop it for all species
dt_pcnt <- data.frame (geometry = st_geometry(list_geckos[[1]])) ## We can leave the 1 because the dataset arno_gall has all the edges that will be used in the analysis
for (i in 1:length(list_geckos)) {
  dt = list_geckos[[i]]
  dt$feems_sum = dt$feems + (- min(dt$feems))
  dt$pcnt = dt$feems_sum / max(dt$feems_sum)
  dat <- dt %>% dplyr::select(pcnt)
  dt_pcnt <- left_join(dt_pcnt, dat, by = c("geometry" = "geometry"))
}
print(head(dt_pcnt))
names(dt_pcnt) <- names_df_geckos[1:length(names_df_geckos)]
print(head(dt_pcnt))
# create a new variable to store all average numbers
dt_pcnt$Average <- NA
dt_pcnt <- dt_pcnt[,2:length(dt_pcnt)]
# Average all values of each edge for only the species that are present in each edge
for (i in 1:NROW(dt_pcnt)) {
  row <- dt_pcnt[i,]
  row <- row[,!is.na(row[1,])]
  divided <- length(row)
  sum <- rowSums(row)
  dt_pcnt$Average[i] <- sum/divided
}
head(dt_pcnt)
#write.csv2(dt_pcnt, "percentage_data_geckos_highleands.csv")
# Append the resulting file to the previously stored geometry 
dt_pcnt_shp <- cbind(dt_pcnt, geometry)
# Write a shapefile
st_write(dt_pcnt_shp,"results/pcnt_averaged_feems_all_hjrs_all_gecck_log.shp")




################################## Logarithmic analysis ####################################


## Export all species with the inverse logarithm to check them with QGIS
for (i in 1:length(list_feems)) {
  Outfile <- paste("log_added/", feems_files_names[i], ".bestCV.shp", sep = "")
  st_write ( list_feems[[i]], Outfile)
}

# Get values from each analysis to have each species edge value in a range of non logarithmic values
## Loop it for all species
dt_log <- data.frame (geometry = st_geometry(list_geckos[[1]])) ## We can leave the 1 because the dataset arno_gall has all the edges that will be used in the analysis
for (i in 1:length(list_geckos)) {
  dt = list_geckos[[i]]
  dt$feems = dt$logW
  dat <- dt %>% dplyr::select(feems)
  dt_log <- left_join(dt_log, dat, by = c("geometry" = "geometry"))
}
print(head(dt_log))
names(dt_log) <- names_df_geckos[1:length(names_df_geckos)]
print(head(dt_log))
# create a new variable to store all average numbers
dt_log$Average <- NA
dt_log <- dt_log[,2:length(dt_log)]
# Average all values of each edge for only the species that are present in each edge
for (i in 1:NROW(dt_log)) {
  row <- dt_log[i,]
  row <- row[,!is.na(row[1,])]
  divided <- length(row)
  sum <- rowSums(row)
  dt_log$Average[i] <- sum/divided
}
head(dt_log)

## create_the log10 logarithm of each value
dt_log$feems_log <- log10(dt_log$Average)

#write.csv2(dt_log, "percentage_data_geckos_highleands.csv")
# Append the resulting file to the previously stored geometry 
dt_log_shp <- cbind(dt_log, geometry)
# Write a shapefile
st_write(dt_log_shp,"results/pcnt_averaged_feems_all_hjrs_all_gecck_log.shp")


################################## Averaged ####################################

#### Create an analysis with all the model values scaled to 0 but maintainging each standard deviation

# generate an empty data_frame 
dt_deff <- data.frame (geometry = st_geometry(list_geckos[[1]])) ## We can leave the 1 because the dataset arno_gall has all the edges that will be used in the analysis
## appen the geometry to all shapefiles
for (i in 1:length(list_geckos)) {
  dt <- dplyr::select(list_geckos[[i]], feems)
  dt_deff <- left_join(dt_deff, dt, by = c("geometry" = "geometry"))
}
#check that the result is in order
write.csv2(dt_deff,"dt_deff_test.csv")

names(dt_deff) <- names_df_geckos[1:length(names_df_geckos)]
head(dt_deff)
scaled.dat <- scale(dt_deff[,2:length(dt_deff)])
dtt <- as.data.frame(scaled.dat[,1:NCOL(scaled.dat)])
dtt$average <- NA
head(dtt)
for (i in 1:NROW(dtt)) {
  row <- dtt[i,]
  row <- row[,!is.na(row[1,])]
  divided <- length(row)
  sum <- rowSums(row)
  dtt$average[i] <- sum/divided
}
head(dtt)

dt_shp <- cbind(scaled.dat, dtt$average, geometry) 
st_write(dt_shp,"results/averaged_feems_all_hjrs_geckos.shp")


################################## Interpolate result ####################################


dt_pcnt_shp <- st_read("pcnt_averaged_feems_all_hjrs_all_gecck_log.shp") %>% st_transform(4326)
dt_avg_shp <- st_read("results/averaged_feems_all_hjrs_geckos.shp") %>% st_transform(4326)
plot(dt_avg_shp[(length(dt_avg_shp))-1])
r <- raster(dt_avg_shp, res = .01)
rr<- raster::rasterize(as(dt_avg_shp, "Spatial"), r, field = "dtt_vrg") 
plot(rr)
xyz <- rasterToPoints(rr)
tps <- Tps(xyz[,1:2], xyz[,3])

p <- raster(rr)
p <- interpolate(p, tps)
plot(rr)
plot(p)

sp <- rasterToPolygons(rr, dissolve = T) %>% st_as_sf() 
sp <- st_union(sp)
sp <- st_remove_holes(sp)

hjrs_30km <- st_read("../../BUFFERS_2_10_20_30km/hj_mnt_bf_30km.shp") 

dt_avg_shp_r <- mask(p, as(sp,"Spatial"))
dt_avg_shp_r_30 <- mask(p, as(hjrs_30km,"Spatial"))
library(RColorBrewer)
colors <- brewer.pal(8, 'YlOrRd')

plot(dt_avg_shp_r_30, col = heat.colors(15))
plot(dt_avg_shp[(length(dt_avg_shp))-1])

writeRaster(dt_avg_shp_r_30, "../../avg_raster_interpolated_all_geckos_0.1_TPS.tiff")