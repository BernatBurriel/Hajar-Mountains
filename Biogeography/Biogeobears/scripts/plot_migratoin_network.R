## create a dataframe with the dispersions from each block to another in different period times

rm(list = ls())
libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape","DescTools", "qgraph")
lapply(libs, require, character.only = TRUE)

setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/")
#setwd("C:/Users/User/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography")

# read the consensus states for all the species (except Pseudotrapelus and Pristurus gallagheri)
consensus_statesnoj <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/noJ/consensus_states.rds")
consensus_states <- consensus_statesnoj
# read the consensus_genera trees 
consensus_genera <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/genera_trees.rds")
genera = names(consensus_genera)

# 1. create a cumulative list of accumulation of lineages through time in each block
lin_t_time <- vector("list", length(genera))
names(lin_t_time) <- genera
for (i in 1:length(genera)) {
  states <- consensus_statesnoj[[i]]
  lin_t_time[[i]] <- tibble(genera = genera[i], node = c(1:length(states)), heigth_max = 0,heigth_min = 0, biogeo = states, biogeoP = 0, TIP_INNER = 0, parent = 0)
  }

for (i in 1:length(genera)) {
  for (j in 1:nrow(lin_t_time[[i]])){
    tipnodes <- 1:length(consensus_genera[[i]]@phylo$tip.label)    
    innernodes <- c(length(consensus_genera[[i]]@phylo$tip.label)+1:nrow(lin_t_time[[i]]))
    if (lin_t_time[[i]]$node[j] %in% tipnodes) {
      lin_t_time[[i]]$TIP_INNER[j] <- "TIP"
    } 
    if(lin_t_time[[i]]$node[j] %in% innernodes){
      lin_t_time[[i]]$TIP_INNER[j] <- "INNODE"
    }
    if (j %in% consensus_genera[[i]]@phylo$edge[, 2]) {
      lin_t_time[[i]]$parent[j] <- consensus_genera[[i]]@phylo$edge[consensus_genera[[i]]@phylo$edge[,2] == j, 1]   
    }
    lin_t_time[[i]]$heigth_max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[as_tibble(consensus_genera[[i]])$node == j][[1]][2]
    lin_t_time[[i]]$heigth_min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[as_tibble(consensus_genera[[i]])$node == j][[1]][1]
    }}  

for (i in 1:length(genera)) {
  for (j in 1:nrow(lin_t_time[[i]])){
    if(lin_t_time[[i]]$parent[j] == 0){ 
      lin_t_time[[i]]$parent[j] <- lin_t_time[[i]]$node[j]
      lin_t_time[[i]]$heightvec[[j]][2] <- lin_t_time[[i]]$heightvec[[j]][1]
    }
  }}
for (i in 1:length(genera)) {
  for (j in 1:nrow(lin_t_time[[i]])){
lin_t_time[[i]]$heightvec[j] <- list(c(lin_t_time[[i]]$heigth_max[lin_t_time[[i]][j,]$parent], lin_t_time[[i]][j,]$heigth_min))
lin_t_time[[i]]$biogeoP[j] <- lin_t_time[[i]]$biogeo[lin_t_time[[i]][j,]$parent]
}}



#### CALCULATE LINEAGE ACUMULATION IN EACH BLOCK ####
states <- c()
for (i in 1:length(consensus_states)) {
  states <- c(states, consensus_states[[i]])
}
states <- unique(states)

#### define central, east and west states ####
central <- c("C", "CI", "CW", "EW", "CE", "CEW")
east <- c("E",  "EW",  "CE",  "CEW" ,"EM")
west <- c( "W","IW","CW", "EW","CEW")

####define the time vectors to assess
time_vec <- list(c(80,40),c(40,39), c(30,29.9), c(20,19.9), c(10,9.9), c(5,4.9))

lineage_accumulation  <- data.frame(year = c("all", '80-40','40-39','30-29.9', '20-19.9', '10-9.9','5-4.9'), CENTRAL = 0, WEST = 0, EAST = 0)

for(i in 1:length(genera)){
  vec_matrix <- as.matrix(data.frame(A= sapply(lin_t_time[[i]]$heightvec, "[[", 1), B = sapply(lin_t_time[[i]]$heightvec, "[[", 2)), byrow = T, ncol = 2)         
  for (j in 1:length(lineage_accumulation$year)) {
    if (lineage_accumulation$year[j] == "all") {
      tbl_states <- as.data.frame(lin_t_time[[i]] %>% dplyr::filter(TIP_INNER=="TIP")%>% group_by(biogeo) %>% summarise(count = n()))
      CENTRAL <- sum(tbl_states$count[tbl_states$biogeo %in% central])
      EAST <- sum(tbl_states$count[tbl_states$biogeo %in% east])
      WEST <- sum(tbl_states$count[tbl_states$biogeo %in% west])
      year_name <- "all"
      lineage_accumulation[j,] <- data.frame(year = year_name, CENTRAL = sum(lineage_accumulation$CENTRAL[j], CENTRAL), WEST = sum(lineage_accumulation$WEST[j], WEST), EAST = sum(lineage_accumulation$EAST[j], EAST))
    }
    if(j > 1) {
      tbl_states <- as.data.frame(lin_t_time[[i]][Overlap(na.omit(vec_matrix),time_vec[[j-1]]) > 0, ]) %>% group_by(biogeo) %>% summarise(count = n())
      CENTRAL <- sum(tbl_states$count[tbl_states$biogeo %in% central])
      EAST <- sum(tbl_states$count[tbl_states$biogeo %in% east])
      WEST <- sum(tbl_states$count[tbl_states$biogeo %in% west])
      year_name <- paste(time_vec[[j-1]][1], "-", time_vec[[j-1]][2], sep = "")
      lineage_accumulation[j,] <- data.frame(year = year_name, CENTRAL = sum(lineage_accumulation$CENTRAL[j], CENTRAL), WEST = sum(lineage_accumulation$WEST[j], WEST), EAST = sum(lineage_accumulation$EAST[j], EAST))
    }}}

saveRDS(lineage_accumulation, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/lineage_accumulation_in_blocks.rds")

### Create the matrices of dispersion between blocks by each time period

# 1. create the variables of states to search
lineage_accumulation <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/lineage_accumulation_in_blocks.rds")
central = central
central_noE <- c("C", "CI", "CW")
central_noW <- c("C", "CI","CE")
east = east
east_noC <- c("E", "EM")
east_noW <- c("E",   "CE" ,"EM")
west = west 
west_noC <- c("W","IW")
west_noE  <- c("W","IW","CW")

# create a dataframe with the age range of the dispersion, the dispersion event and the time of the dispersion
dispersions <- vector('list', length = NROW(lineage_accumulation))
names(dispersions) <-  c("all", '80-40','40-30','30-20', '20-10', '10-5','5-0')
for (i in 1:length(dispersions)) {
  dispersions[[i]] <- data.frame(Dispersal_from = c('CENTRAL','EAST','WEST'), CENTRAL = 0, EAST = 0, WEST = 0) 
  }

time_vec2 <- time_vec <- list(c(80,40),c(40,30), c(30,20), c(20,10), c(10,5), c(5,0))
for (j in 1:length(dispersions)) {
if(j == 1){
    for (i in 1:length(genera)) {
    # Dispersion to central block
    EtoC <- nrow(lin_t_time[[i]][lin_t_time[[i]]$biogeo %in% central & lin_t_time[[i]]$biogeoP %in% east_noC,]) 
    WtoC <- nrow(lin_t_time[[i]][lin_t_time[[i]]$biogeo %in% central & lin_t_time[[i]]$biogeoP %in% west_noC,]) 
    dispersions[[j]][,2] <- c(lineage_accumulation$CENTRAL[j], sum(dispersions[[j]][2,2],EtoC), sum(dispersions[[j]][3,2],WtoC))
    # Dispersion to east block
    CtoE <- nrow(lin_t_time[[i]][lin_t_time[[i]]$biogeo %in% east & lin_t_time[[i]]$biogeoP %in% central_noE,]) 
    WtoE <-nrow(lin_t_time[[i]][lin_t_time[[i]]$biogeo %in% east & lin_t_time[[i]]$biogeoP %in% west_noE,]) 
    dispersions[[j]][,3] <- c( sum(dispersions[[j]][1,3],CtoE),lineage_accumulation$EAST[j], sum(dispersions[[j]][3,3],WtoE))
    # Dispersion to west block
    CtoW <- nrow(lin_t_time[[i]][lin_t_time[[i]]$biogeo %in% west & lin_t_time[[i]]$biogeoP %in% central_noW,]) 
    EtoW <- nrow(lin_t_time[[i]][lin_t_time[[i]]$biogeo %in% west & lin_t_time[[i]]$biogeoP %in% east_noW,])
    dispersions[[j]][,4] <- c(sum(dispersions[[j]][1,4],CtoW), sum(dispersions[[j]][2,4],WtoC), lineage_accumulation$WEST[j])
    }}
if (j > 1) {
for (i in 1:length(genera)) {
  vec_matrix <- as.matrix(data.frame(A= sapply(lin_t_time[[i]]$heightvec, "[[", 1), B = sapply(lin_t_time[[i]]$heightvec, "[[", 2)), byrow = T, ncol = 2)         
  dt <- lin_t_time[[i]]
  dt <- dt[Overlap(na.omit(vec_matrix),time_vec[[j-1]]) > 0,]
  # Dispersion to central block
  EtoC <- nrow(dt[dt$biogeo %in% central & dt$biogeoP %in% east_noC,]) 
  WtoC <- nrow(dt[dt$biogeo %in% central & dt$biogeoP %in% west_noC,]) 
  dispersions[[j]][,2] <- c(lineage_accumulation$CENTRAL[j], sum(dispersions[[j]][2,2],EtoC), sum(dispersions[[j]][3,2],WtoC))
  # Dispersion to east block
  CtoE <- nrow(dt[dt$biogeo %in% east & dt$biogeoP %in% central_noE,]) 
  WtoE <-nrow(dt[dt$biogeo %in% east & dt$biogeoP %in% west_noE,]) 
  dispersions[[j]][,3] <- c( sum(dispersions[[j]][1,3],CtoE),lineage_accumulation$EAST[j], sum(dispersions[[j]][3,3],WtoE))
  # Dispersion to west block
  CtoW <- nrow(dt[dt$biogeo %in% west & dt$biogeoP %in% central_noW,]) 
  EtoW <- nrow(dt[dt$biogeo %in% west & dt$biogeoP %in% east_noW,])
  dispersions[[j]][,4] <- c(sum(dispersions[[j]][1,4],CtoW), sum(dispersions[[j]][2,4],WtoC), lineage_accumulation$WEST[j])
}}}

class(dispersions[[1]])
input <- as.matrix(dispersions[[1]][,2:4])
input[1,1] <-input[2,2] <- input[3,3] <- 0
library(qgraph)
# plot the migration network
pops <- c("A1.1","A1.2","A2","A3","B1","B2","B3","B4","B5","O2.1","O2.2","O3")
df1 <- zip(pops,migRes[[3]])


colnames (migRes[[3]]) <- c("A1.1","A1.2","A2","A3","B1","B2","B3","B4","B5","O2.1","O2.2","O3")
nm <- migRes[[3]]
qgraph(input, 
       edge.labels = FALSE,
       layout="spring", theme="Hollywood",
       sampleSize = TRUE, edge.color = "navy",
       esize=3 , curveAll=TRUE,
       asize=2.7, curve=3,
       border.width=1, vsize=c(dispersions[[i]][1,2],dispersions[[i]][2,3],dispersions[[i]][3,4]),
       labels = colnames(input))

?qgraph

c(dispersions[[i]][1,1],dispersions[[i]][2,2],dispersions[[i]][3,3])


















