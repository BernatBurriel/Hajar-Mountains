## Retrieve elevation data for the S1 table 

# 1. Set wd 

setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/bayestraits_elevation/")

library(dplyr);library(sf);library(raster);library(tidytree);library(treeio);library(phytools);library(geiger)


species_files <- list.files('../BSSVS/locations', pattern = "\\.txt$", full.names = T)
genera_names_F <- list.files('../BSSVS/locations', pattern = "\\.txt$", full.names = F)
genera_names <- gsub('.txt', '',genera_names_F)
genera_names <- gsub('[[:digit:]]+', '', genera_names)

elevation <- read.csv2('data/TableS1_ddrad_individuals.csv')
names(elevation)
elevation <- elevation %>% select(Name = SPECIMEN_CODE, elev = elevation, lineage = Lineage.name )
elevation$Name <- gsub('-','_', elevation$Name)
#elevation$Name <- gsub('AO054','_', elevation$Name)

genera_files <- vector('list', length(species_files))
names(genera_files) <- genera_names

for (i in 1:length(genera_files)) {
  genera_files[[i]] <- read.table(species_files[i], header = T)  
}
genera_files[[7]]$Name2 <- genera_files[[7]]$Name
genera_files[[7]]$Name <- genera_files[[7]]$New_name
genera_files[[7]]$Name <- gsub('_C' , '', genera_files[[7]]$Name)
genera_files[[7]]$Name <- gsub('_E' , '', genera_files[[7]]$Name)
genera_files[[7]]$Name <- gsub('_W' , '', genera_files[[7]]$Name)
genera_files[[7]]$Name <- gsub('_O' , '', genera_files[[7]]$Name)

elevaton_cont <- genera_files
for (i in 1:length(genera_files)) {
  genera_files[[i]]$elev <- elevation$elev[base::match(genera_files[[i]]$Name, elevation$Name)]
  genera_files[[i]]$lineage <- elevation$lineage[base::match(genera_files[[i]]$Name, elevation$Name)]
  elevaton_cont[[i]] <- genera_files[[i]] %>% dplyr::select(Name, elev, lineage)
}

for (i in 1:length(genera_files)) {
  genera_files[[i]]$low <- 0
  genera_files[[i]]$mid <- 0
  genera_files[[i]]$high <- 0
  for (j in 1:NROW(genera_files[[i]])) {
    if (as.numeric(genera_files[[i]]$elev[j]) < 300) {
      genera_files[[i]]$low[j] <- 1
    }
    if (as.numeric(genera_files[[i]]$elev[j]) >= 300 & as.numeric(genera_files[[i]]$elev[j]) < 1500) {
      genera_files[[i]]$mid[j] <- 1
    }
    if (as.numeric(genera_files[[i]]$elev[j]) >= 1500) {
      genera_files[[i]]$high[j] <- 1
    }
  }
}
tbl_data_simmap <- vector('list', length(genera_files))
names(tbl_data_simmap) <- genera_names_F
table_renaming <- vector('list', length(genera_files))
names(table_renaming) <- genera_names_F

for (i in 1:length(tbl_data_simmap)) {
  tbl_data_simmap[[i]] <- genera_files[[i]] %>% dplyr::select(Name, low,mid,high)
  table_renaming[[i]] <- genera_files[[i]] %>% dplyr::select(Name, New_name)
}
table_renaming$`07Prup.txt` <- genera_files$Prup %>% dplyr::select(Name2, Name)

for (i in 1:length(genera_files)) {
  Outfile_elev <- paste0('../Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/data/elevation/', names(genera_files[i]), "_elevation.txt")
  Outfile_simmmap <- paste0('../Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/data/species/', names(table_renaming[i]))
  Outfile_names <- paste0('../Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/data/names/', names(table_renaming[i]))
  write.table(elevaton_cont[[i]], Outfile_elev, quote = F, row.names = F)
  write.table(table_renaming[[i]], Outfile_names, quote = F, row.names = F)
  write.table(tbl_data_simmap[[i]], Outfile_simmmap, quote = F, row.names = F)
}


tree_file <- read.beast('../Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/trees/07Prup.tree')
tree <- treeio::rename_taxa(tree_file, data = table_renaming$`07Prup.txt`, key = Name2, value = Name)
write.beast(tree, '../Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/trees/07Prup.tree')

table_renaming <- read.table(file = '../Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/data/names/10Trachydactylus.txt', header = T)
tree_file <- read.beast('../Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/trees/10Trachidactylus.tree')
tree <- treeio::rename_taxa(tree_file, data = table_renaming, key = New_name, value = Name)
write.beast(tree, '../Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/trees/10Trachidactylus.tree')


