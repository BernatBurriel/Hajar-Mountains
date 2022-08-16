# building a Global grid

#install.packages('dggridR')
library(dggridR)
library(dplyr)
library(sf)
library(sp)


## Construct a grid that spaces around 10km between centroids
tst <- dgconstruct(topology = "TRIANGLE", res = 15, metric = T, aperture = 4, resround = "down")
dggetres(tst)
# From the starting Discrete Global Grid, get the resolution needed to have a spacing of around 10km
res <- dg_closest_res_to_spacing(tst,spacing=5,round='down',metric=T)
# Convert the DGG with the calculated resolution
dggs <- dgsetres(tst,res)


## Import coordinates from a species sample localities and for the outer margin of the area of study
coords = read.table("PATH/To/coords.coords", header = F)
outer = read.table("path/to/vertices.outer")

#Get the corresponding grid cells for each data_point and for the outer polygon (lat-long pair) converted into the produced dggs
dggcoords <- dgGEO_to_SEQNUM(dggs, in_lon_deg =  coords$V1, in_lat_deg =  coords$V2)$seqnum
dggouter <- dgGEO_to_SEQNUM(dggs, in_lon_deg =  outer$V1, in_lat_deg =  outer$V2)$seqnum

#Converting SEQNUM to GEO gives the center coordinates of the cells. This is optional because feems also works without moving the data points into the centroids of each grid cell. 
cellcenters_coords <- dgSEQNUM_to_GEO(dggs,dggcoords)
coords_def <- data.frame(Lon = cellcenters_coords$lon_deg, Lat = cellcenters_coords$lat_deg)
write.table (coords_def, 'DATA/temp/ptyo.coords', quote = F, row.names = F, col.names = F, sep = "\t")
cellcenters_outer <- dgSEQNUM_to_GEO(dggs,dggouter)
outer_def <- data.frame(Lon = cellcenters_outer$lon_deg, Lat = cellcenters_outer$lat_deg)
write.table (outer_def, 'DATA/temp/ptyo.outer', quote = F, row.names = F, col.names = F, sep = "\t")


## Save the final grid that will be the background for all the species analysis.

shapefile_to_grid <- dgshptogrid(dggs,       ## abstract DGGs to get the data from
                                 "DATA/temp/outline_for_grid_creation.shp",  ## Shapefile that outlines the margins of the resulting grid
                                 cellsize = 0.04545,  ## Distance, in degrees, between the sample points used to generate the grid. Small values yield long generation times while large values may omit cells
                                 savegrid = "DATA/temp/test_shapefile_to_grid_10km.shp") # produced grid.