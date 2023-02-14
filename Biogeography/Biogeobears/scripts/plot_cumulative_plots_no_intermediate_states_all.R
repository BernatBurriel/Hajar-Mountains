
##### PLOT cumulative plot ###########

# Calculate how many colonization events have there been 
# Define  biogeographic scenarios 
# Establish the different scenarios. Each scenario will be a list of vectors of three elements:
# (state1, state2, state3), where state1 = state in parental node, state2 = state in descendentant 1,
# state3 = state in descendent 2.
libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape","DescTools")
lapply(libs, require, character.only = TRUE)

setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/")
#setwd("C:/Users/User/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography")

big_tree <- read.beast("Mountain_colonization_Simmaps/TREES/SQUAMATA_ALL_PF.tree")
big_tree@phylo$tip.label <- gsub("_", " ", big_tree@phylo$tip.label)
big_tree_phylo <- as.phylo(big_tree)

consensus_statesnoj <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/noJ/consensus_states.rds")


mountains_simmap <- readRDS("Mountain_colonization_Simmaps/simmaps_all_tree/objects/mountains_simmap_alltree.rds")
colonization <- list(c("Out", "Out", "Hajars"), c("Out", "Hajars", "Out"))
extirpation <- list (c("Hajars", "Out", "Hajars"), c("Hajars",  "Hajars", "Out"), c("Hajars", "Out", "Out"))


## biogeobears trees

# read the consensus_genera trees 
consensus_genera <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/genera_trees.rds")
genera = names(consensus_genera)

states <- c("C", "CE", "CI", "CW", "E", "CEW", "EM", "I", "IW", "M", "R", "W", "WC", "IW")
#### 
#C_col <-  list(c("CI", "C", "I"), c("CI", "I", "C"), c("CE", "C", "E"),
#              c("CE", "E", "C"), c("W", "W", "C"), c("W", "C", "W"),
#             c("E","C","E"), c("E","E","C"), c("CW", "C","W"), c("CW", "W","C"),
#            c("E","CE","C"), c("E","C","CE"))

#C_col <-  list(c("CI", "C"),  c("CE", "C"), c("W", "C"), c("E","C"), c("CW", "C"))
C_col <- list(c("W", "C"), c("E","C"), c("W","CW"), c("W","CE"),c("E","CE"),c("E","CW"), c("W","EW"), c("EM","EW"))
#C_div <- list(c("C", "C", "C"))
C_div <- list(c("C", "C", "C"), c("C", "CE", "C"), c("C", "C", "CE"),c("CE", "C", "CE"), c("CE", "CE", "C"), c("CE","CEW", "CEW"), c("CW", "C", "CW"),c("CW", "CW", "C"),
              c("CW", "C", "C"))

#E_col <- list(c("C","E","C"), c("C","C","E"), c("CW","E","CW"), c("CW","CW","E"), 
#             c("EW","CW","E"), c("EW", "E","EW"), c("EW", "E", "CW"), c("EW","CW","E"),
#            c("CE","C","E"), c("CE","E","C"), c("CE","E","CW"),c("CE","CW","E"),
#           c("CE","E","CW"), c("CE","CW","E"))
#E_col <- list(c("C","E"), c("CW","E"),c("EW","E"), c("CE","E")) 
E_col <- list(c("C","E"), c("CW","E"),c("CW","EW"), c("CW", "CE"), c("C","CE"), c("W","EW"), c("EM", "EC"), c("EM","CE"))
#E_div <- list(c("E", "E","E"))
E_div <- list(c("E", "E","E"),c("CE","CEW", "CEW"), c("E","CE","E"), c("E","E","CE"))

#W_col <- list(c("C", "W", "C"), c("C", "C", "W"), c("IW", "W", "I"), c("IW","I","W"), c("CW", "W", "C"), c("CW", "C", "W"),
#             c("CW", "CE", "W"), c("CW","W","CE"), c("CEW", "CE", "W"), c("CEW", "W", "CE"))
#W_col <- list(c("C", "W"), c("IW", "W"), c("CW", "W"), c("CEW", "W"))
W_col <- list(c("C", "W"), c("I", "IW"), c("C", "CW"), c("CE", "CEW"),c("CE", "CW"), c("EC","CW"), c("EM","EW"))
#W_div <- list(c("W","W","W"))
W_div <- list(c("W","W","W"),c("W","EW","W"),c("W","W","EW"), c("W","EW","W"))

vicariance_CE <- list(c("CE", "C","E"), c("CE","E","C"))
vicariance_CW <-list(c("CW", "C","W"), c("CW","W","C"), c("CW", "CE", "W"), c("CW", "W","CE"))
vicariance_EW <- list(c("EW", "W","E"), c("EW","E","W"), c("EW", "E", "CW"), c("EW", "CW", "E"), c("CEW", "W", "CE"), c("CEW", "CE", "W"))
vicariance_IW<- list(c("IW","I","W"), c("IW","W","I"))

# this are the ones to take into account when counting intermediate states
#vicariance_CE <- list(c("CW", "CE"), c("C","CE"), c("CEW", "CE"), c("E","CE"), c("EM", "CE")
#vicariance_CW <- list(c("W", "CW"), c("EW","CW"),c("CE", "CW"))
#vicariance_EW <- list(c("W","EW"), c("CW","EW"))
#vicariance_IW <- list(c("I", "IW"))

{nodes_hajarscol <- c()
  
  nodes_hajars_extirpation <- c()
  
  Ccol_nodes_list_cons <- vector("list", length(genera))
  names(Ccol_nodes_list_cons) <- genera
  Ccol_nodes_list_cons[1:length(genera)] <- 0
  
  Cdiv_nodes_list_cons <- vector("list", length(genera))
  names(Cdiv_nodes_list_cons) <- genera
  Cdiv_nodes_list_cons[1:length(genera)] <- 0
  
  Ecol_nodes_list_cons <- vector("list", length(genera))
  names(Ecol_nodes_list_cons) <- genera
  Ecol_nodes_list_cons[1:length(genera)] <- 0
  
  Ediv_nodes_list_cons <- vector("list", length(genera))
  names(Ediv_nodes_list_cons) <- genera
  Ediv_nodes_list_cons[1:length(genera)] <- 0
  
  Wcol_nodes_list_cons <- vector("list", length(genera))
  names(Wcol_nodes_list_cons) <- genera
  Wcol_nodes_list_cons[1:length(genera)] <- 0
  
  Wdiv_nodes_list_cons <- vector("list", length(genera))
  names(Wdiv_nodes_list_cons) <- genera
  Wdiv_nodes_list_cons[1:length(genera)] <- 0
  
  vicariance_CE_nodes_list_cons <- vector("list", length(genera))
  names(vicariance_CE_nodes_list_cons) <- genera
  vicariance_CE_nodes_list_cons[1:length(genera)] <- 0
  
  vicariance_CW_nodes_list_cons <- vector("list", length(genera))
  names(vicariance_CW_nodes_list_cons) <- genera
  vicariance_CW_nodes_list_cons[1:length(genera)] <- 0
  
  vicariance_EW_nodes_list_cons <- vector("list", length(genera))
  names(vicariance_EW_nodes_list_cons) <- genera
  vicariance_EW_nodes_list_cons[1:length(genera)] <- 0
  
  vicariance_IW_nodes_list_cons <- vector("list", length(genera))
  names(vicariance_IW_nodes_list_cons) <- genera
  vicariance_IW_nodes_list_cons[1:length(genera)] <- 0
  
}



#### Prepare all the data to wrap all colonization and diversification nodes in a single list for each biogeographic state
colonization_states <- list(C_col, E_col, W_col)
names_colonization <- c("Central colonization", "East colonization", "West colonization")
colonization_nodes <- list(Ccol_nodes_list_cons, Ecol_nodes_list_cons, Wcol_nodes_list_cons)
names(colonization_nodes) <- names_colonization

## All Colonization nodes (descendant nodes)
for (z in 1:length(colonization_nodes)) {
  for (i in 1:length(genera)){
    for (j in sort(as.phylo(consensus_genera[[i]])$edge[,2])){  
      parent <- consensus_statesnoj[[i]][as.phylo(consensus_genera[[i]])$edge[,1][as.phylo(consensus_genera[[i]])$edge[,2]==j]]
      son <- consensus_statesnoj[[i]][j]
      d <- c(parent, son)
      if (list(d) %in% colonization_states[[z]]){
        colonization_nodes[[z]][[i]] <- append(colonization_nodes[[z]][[i]], j)
      }}}}

## All Diversification nodes (parent nodes) 
diversification_states <- list(C_div, E_div, W_div,vicariance_CE,vicariance_CW, vicariance_EW, vicariance_IW)
names_diversification <- c("Central diversification", "East diversification", "West diversification", "Central-East vicariance", 
                           "Central-West vicariance", "East-West vicariance", "Iran-West vicariance")
diversification_nodes <- list(Cdiv_nodes_list_cons, Ediv_nodes_list_cons, Wdiv_nodes_list_cons,vicariance_CE_nodes_list_cons, 
                              vicariance_CW_nodes_list_cons, vicariance_EW_nodes_list_cons, vicariance_IW_nodes_list_cons)
names(diversification_nodes) <- names_diversification

for (z in 1:length(diversification_nodes)) {
  for (i in 1:length(genera)){
    for (j in min(as.phylo(consensus_genera[[i]])$edge[,1]):max(as.phylo(consensus_genera[[i]])$edge[,1])){  
      v <- c(consensus_statesnoj[[i]][j],
             consensus_statesnoj[[i]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][1],
             consensus_statesnoj[[i]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][2])
      if (list(v) %in% diversification_states[[z]]){
        diversification_nodes[[z]][[i]] <- append(diversification_nodes[[z]][[i]], j)
      }}}}

# Now we have 4 lists, 1. nodes_hajarscol (descendant node of HAJARS colonization event); 2.nodes_hajars_extirpation
# (descendant node extirpated from Hajars); 3. colonization_nodes (list containing colonization descendant nodes of each mountain block)
# 4. diversification_nodes (list containing diversification parent nodes of diversification in the blocks, or vicariance between them)
# Lists 1 and 2 contain data for all the species (BIG_tree)
# Lists 3 and 4 are lists of lists containing the each colonization (3) or diversification (4) event in each bloack and in each genus
# First, we're going to remove the zeros from the vectors with transition nodes.

for (z in 1:length(colonization_nodes)) {
  for (i in 1:length(genera)){
    if (length(colonization_nodes[[z]][[i]]) != 1){
      colonization_nodes[[z]][[i]] <- colonization_nodes[[z]][[i]][-1]
    }}}

for (z in 1:length(diversification_nodes)) {
  for (i in 1:length(genera)){
    if (length(diversification_nodes[[z]][[i]]) != 1){
      diversification_nodes[[z]][[i]] <- diversification_nodes[[z]][[i]][-1]
    }}}


#### Create a dataframe for each colonization event and afterwards join all of them together
colonization_nodes_df <- vector('list', length(colonization_nodes))
names(colonization_nodes_df) <- names(colonization_nodes)
for (i in 1:length(colonization_nodes_df)) {
  colonization_nodes_df[[i]] <- vector('list', length(genera))
  names(colonization_nodes_df[[i]]) <- genera
}
diversification_nodes_df <- vector('list', length(diversification_nodes))
names(diversification_nodes_df) <- names(diversification_nodes)
for (i in 1:length(diversification_nodes_df)) {
  diversification_nodes_df[[i]] <- vector('list', length(genera))
  names(diversification_nodes_df[[i]]) <- genera
}

for (z in 1:length(colonization_nodes_df)){
  for (i in 1:length(genera)){
    colonization_nodes_df[[z]][[i]] <- data.frame(node=colonization_nodes[[z]][[i]], genus=genera[i],
                                                  height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                                                  vicariance_EW = 0, vicariance_IW = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
    if (names(colonization_nodes_df[z]) == "Central colonization") {colonization_nodes_df[[z]][[i]]$Ccol = 1}
    if (names(colonization_nodes_df[z]) == "East colonization") {colonization_nodes_df[[z]][[i]]$Ecol = 1}
    if (names(colonization_nodes_df[z]) == "West colonization") {colonization_nodes_df[[z]][[i]]$Wcol = 1}
  }
}


for (z in 1:length(diversification_nodes_df)){
  for (i in 1:length(genera)){
    diversification_nodes_df[[z]][[i]] <- data.frame(node=diversification_nodes[[z]][[i]], genus=genera[i],
                                                     height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                                                     vicariance_EW = 0, vicariance_IW = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
    if (names(diversification_nodes_df[z]) == "Central diversification") {diversification_nodes_df[[z]][[i]]$Cdiv = 1}
    if (names(diversification_nodes_df[z]) == "East diversification") {diversification_nodes_df[[z]][[i]]$Ediv = 1}
    if (names(diversification_nodes_df[z]) == "West diversification") {diversification_nodes_df[[z]][[i]]$Wdiv = 1}
    if (names(diversification_nodes_df[z]) == "Central-East vicariance") {diversification_nodes_df[[z]][[i]]$vicariance_CE = 1}
    if (names(diversification_nodes_df[z]) == "Central-West vicariance") {diversification_nodes_df[[z]][[i]]$vicariance_CW = 1}
    if (names(diversification_nodes_df[z]) == "East-West vicariance") {diversification_nodes_df[[z]][[i]]$vicariance_EW = 1}
    if (names(diversification_nodes_df[z]) == "Iran-West vicariance") {diversification_nodes_df[[z]][[i]]$vicariance_IW = 1}
  }
}

#### JOIN ALL_TREE nodes and create also a dataframe. It will be attaach later on because we will have to keep processing the data in differentways. ##
hajar_nodes <- list(nodes_hajarscol, nodes_hajars_extirpation)

#### To continue generating the dataframe for this two datasets ####
event_list_cons_genera <- vector("list", length(genera))
names(event_list_cons_genera) <- genera

#### Generate a dataframe for all the biogeografic events separated by genus ####
for (i in 1:length(genera)){
  for (j in 1:length(colonization_nodes_df)) {
    event_list_cons_genera[[i]] <- rbind(event_list_cons_genera[[i]], colonization_nodes_df[[j]][[i]])
  }
  for (z in 1:length(diversification_nodes_df)) {
    event_list_cons_genera[[i]] <- rbind(event_list_cons_genera[[i]], diversification_nodes_df[[z]][[i]])
  }}
# Remove all rows with 0 in the node variable
for (i in 1:length(event_list_cons_genera)){ 
  event_list_cons_genera[[i]] <- event_list_cons_genera[[i]][event_list_cons_genera[[i]]$node != 0, ] 
}

#### Populate the dataframe with the min, max, and height columns ####

# Divergence times and 95% HPD (min and max) for each diversification node 
# and branch length for colonization and dispersal nodes.
for (i in 1:length(genera)){
  if (nrow(event_list_cons_genera[[i]]) > 0){
    for (j in 1:length(event_list_cons_genera[[i]]$node)){
      # diversification nodes
      if (event_list_cons_genera[[i]][j, "Cdiv"] == 1 | event_list_cons_genera[[i]][j, "Ediv"] == 1 | event_list_cons_genera[[i]][j, "Wdiv"] == 1
          | event_list_cons_genera[[i]][j, "vicariance_CE"] == 1 | event_list_cons_genera[[i]][j, "vicariance_CW"] == 1 | event_list_cons_genera[[i]][j, "vicariance_EW"] == 1
          | event_list_cons_genera[[i]][j, "vicariance_IW"] == 1){
        event_list_cons_genera[[i]]$min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[event_list_cons_genera[[i]]$node[j]]][1]
        event_list_cons_genera[[i]]$max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[event_list_cons_genera[[i]]$node[j]]][2]
        event_list_cons_genera[[i]]$height[j] <- as_tibble(consensus_genera[[i]])$height_median[[event_list_cons_genera[[i]]$node[j]]]
      }
      # Colonization or Dispersal nodes
      if (event_list_cons_genera[[i]][j, "Ccol"] == 1 | event_list_cons_genera[[i]][j, "Ecol"] == 1 | event_list_cons_genera[[i]][j, "Wcol"] == 1){
        dat <- as_tibble(consensus_genera[[i]])
        # Min is the age of the descendent node + its 95% HPD
        desc_age <- dat$height_0.95_HPD[dat$node == event_list_cons_genera[[i]]$node[j]][[1]][1]
        event_list_cons_genera[[i]]$min[j] <- desc_age
        # Max is the age of the parental node + its 95% HPD
        parent_age <- dat$height_0.95_HPD[dat$parent[dat$node == event_list_cons_genera[[i]]$node[j]]][[1]][2]
        event_list_cons_genera[[i]]$max[j] <- parent_age
        # let's set the height as the midpoint between the parental and the descendent nodes
        event_list_cons_genera[[i]]$height[j] <- mean(c(desc_age, parent_age))
      }}}}
         
# With this object we can represent the age and type of each biogeographic event for each clade.

# Before saving the object we will add the data of Pristurus gallagheri central diversification
pgall <- read.beast("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/species_only_two_tips/06Pgall.tree")
names_tbl <- read.table("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/species_only_two_tips/06Pgall_name_correspondance.txt", header = T)
pgall <- treeio::rename_taxa(pgall, data = names_tbl, key = Name, value = New_name)
tips_to_drop <- names_tbl[names_tbl$Out_In == "out",]
pgall <- treeio::drop.tip(pgall, as.character(tips_to_drop$New_name))
pgalldt <- as_tibble(pgall)

Pgall_df <- data.frame(node=3, genus="Pgall",
                       height=pgalldt$height_median[3], min=pgalldt$height_0.95_HPD[[3]][1], max=pgalldt$height_0.95_HPD[[3]][2], Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                       vicariance_EW = 0, vicariance_IW = 0, Cdiv = 1, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
rownames(Pgall_df) <- 1
event_list_cons_genera[[10]] <-  Pgall_df
names(event_list_cons_genera) <- c(genera, "Pgall")


### SAVE all events with 
saveRDS(event_list_cons_genera, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/event_list_cons_genera.rds")
event_list_cons_genera <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/event_list_cons_genera.rds")


# Get Hajar_mountain_colonization_nodes in the simmaps tree

#### MOUNTAIN COLONIZATION NODES ####
# get the first colonization nodes
mountains_simmap$ace 
nodes_tips <- data.frame(node = 1:length(big_tree_phylo$tip.label), name = big_tree_phylo$tip.label)

nodes_hajars <- c(rownames(mountains_simmap$ace[mountains_simmap$ace[,1] > 0.5,]), 
                  rownames(nodes_tips[nodes_tips[,2] %in% rownames(mountains_simmap$tips[mountains_simmap$tips[,1] > 0.5,]),]))
'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in 1:length(nodes_hajars)) {
  node <- as.numeric(nodes_hajars[i])
  if (big_tree_phylo$edge[big_tree_phylo$edge[,2] == node,1] %!in% nodes_hajars) {
    nodes_hajarscol <- c(nodes_hajarscol , node)  
  } 
}
length(nodes_hajarscol)

#### MOUNTAIN EXTIRPATION NODES #### 
# Get the descendant node of all extirpaton events  
nodes_hajarscolex <- as.numeric(nodes_hajars[as.numeric(nodes_hajars) > 285])
for (ex in 1:length(nodes_hajarscolex)) {
  node <- nodes_hajarscolex[ex]
  if(big_tree_phylo$edge[big_tree_phylo$edge[,1] == node,2][1] %!in% nodes_hajars){
    nodes_hajars_extirpation <- c(nodes_hajars_extirpation, big_tree_phylo$edge[big_tree_phylo$edge[,1] == node,2][1])
  }
  if(big_tree_phylo$edge[big_tree_phylo$edge[,1] == node,2][2] %!in% nodes_hajars){
    nodes_hajars_extirpation <- c(nodes_hajars_extirpation, node)
  }
}
first_col_nodes <- list(nodes_hajarscol, nodes_hajars_extirpation)
names(first_col_nodes) <- c("First colonization", "Extirpation")
for (i in 1:2) {
  first_col_nodes[[i]] <- data.frame(node=first_col_nodes[[i]], genus=names(first_col_nodes[i]),
                                     height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                                     vicariance_EW = 0, vicariance_IW = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
  if (names(first_col_nodes[i]) == "First colonization") {first_col_nodes[[i]]$colonization = 1}
  if (names(first_col_nodes[i]) == "Extirpation") {first_col_nodes[[i]]$extirpation = 1}
}

# Divergence times and branch length for colonization and extirpation nodes.
for (i in 1:length(first_col_nodes)) {
  if (nrow(event_list_cons_genera[[i]]) > 0){
    for (j in 1:length(first_col_nodes[[i]]$node)){
      # Colonization or Dispersal nodes
      if (first_col_nodes[[i]][j, "colonization"] == 1 | first_col_nodes[[i]][j, "extirpation"] == 1){
        dat <- as_tibble(big_tree)
        # Min is the age of the descendent node + its 95% HPD
        desc_age <- dat$height_0.95_HPD[dat$node == first_col_nodes[[i]]$node[j]][[1]][1]
        first_col_nodes[[i]]$min[j] <- desc_age
        # Max is the age of the parental node + its 95% HPD
        parent_age <- dat$height_0.95_HPD[dat$parent[dat$node == first_col_nodes[[i]]$node[j]]][[1]][2]
        first_col_nodes[[i]]$max[j] <- parent_age
        # let's set the height as the midpoint between the parental and the descendent nodes
        first_col_nodes[[i]]$height[j] <- mean(c(desc_age, parent_age))
      }}}}
dat[dat$node %in% first_col_nodes[[1]]$node,]
genera_order <- c("Omanosaura", "Asaccus", "Ptyodactylus", "Hemidactylus", "Prup", "Echis", "Pcele", "Pgall", "Pseudotrapelus", "Trachydactylus")
first_col_nodes_df <- first_col_nodes
first_col_nodes_df[[1]]$genus <- genera_order
first_col_nodes_df[[2]]$genus <- "Asaccus"

# Add each event to the each genera
event_list_cons_genera2 <- event_list_cons_genera
for (i in 1) {
  for (j in 1:nrow(first_col_nodes_df[[i]])) {
    for (z in 1:length(genera_order)) {
      name_z <- genera_order[z]
      if (first_col_nodes_df[[i]]$genus[j] == genera_order[z]) {
        first_col_nodes_df[[i]][j,]$min <- max(event_list_cons_genera2[[name_z]]$max)
        first_col_nodes_df[[i]][j,]$height <- median(c(first_col_nodes_df[[i]][j,]$min, first_col_nodes_df[[i]][j,]$max))
        event_list_cons_genera2[[name_z]] <- rbind(event_list_cons_genera2[[name_z]], first_col_nodes_df[[i]][j,]) }  
    }}}
event_list_cons_genera2[["Asaccus"]] <- rbind(event_list_cons_genera2[["Asaccus"]], first_col_nodes_df[[2]][1,])
event_list_cons_genera <- event_list_cons_genera2
### SAVE all events with 
saveRDS(event_list_cons_genera, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/event_list_cons_genera.rds")
event_list_cons_genera <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/event_list_cons_genera.rds")

########## Plot the results ########


# Now we will do the consensus by million years of each genus
# The deepest colonization event in the dataset is the one given by the Asaccus group and pristurus celerrimus. Its node is around 40-45 Mya ago. Just to be sure, we will reconstruct from the first knowledge of 
# mountain formation 80 Mya
myr <- 80
ma <- myr
genera2 <- names(event_list_cons_genera)

# We will create a list of vectors representing intervals of 1 Ma.
ma <- vector("list", ma)
for (i in 1:length(ma)){
  ma[[i]] <- c(i-1, i)
}

#### Create all the lists containing the data to plot 
{{event_myr_list_cons <- vector('list', length(event_list_cons_genera))
names(event_myr_list_cons) <- names(event_list_cons_genera)
vicariance_myr_list_cons <- vector("list", length(genera2))
names(vicariance_myr_list_cons) <- genera2
Ccol_myr_list_cons <- vector("list", length(genera2))
names(Ccol_myr_list_cons) <- genera2
Ecol_myr_list_cons <- vector("list", length(genera2))
names(Ecol_myr_list_cons) <- genera2
Wcol_myr_list_cons <- vector("list", length(genera2))
names(Wcol_myr_list_cons) <- genera2
Cdiv_myr_list_cons <- vector("list", length(genera2))
names(Cdiv_myr_list_cons) <- genera2
Ediv_myr_list_cons <- vector("list", length(genera2))
names(Ediv_myr_list_cons) <- genera2
Wdiv_myr_list_cons <- vector("list", length(genera2))
names(Wdiv_myr_list_cons) <- genera2
vicarianceCE_myr_list_cons <- vector("list", length(genera2))
names(vicarianceCE_myr_list_cons) <- genera2
vicarianceCW_myr_list_cons <- vector("list", length(genera2))
names(vicarianceCW_myr_list_cons) <- genera2
vicarianceEW_myr_list_cons <- vector("list", length(genera2))
names(vicarianceEW_myr_list_cons) <- genera2
vicarianceIW_myr_list_cons <- vector("list", length(genera2))
names(vicarianceIW_myr_list_cons) <- genera2
colonization_myr_list_cons <- vector("list", length(genera2))
names(colonization_myr_list_cons) <- genera2
extirpation_myr_list_cons <- vector("list", length(genera2))
names(extirpation_myr_list_cons) <- genera2
}

  
events_myr_names <- c("event_myr_list_cons","vicariance_myr_list_cons","Ccol_myr_list_cons", "Ecol_myr_list_cons", "Wcol_myr_list_cons",
                        "Cdiv_myr_list_cons","Ediv_myr_list_cons","Wdiv_myr_list_cons", 
                        "vicarianceCE_myr_list_cons", "vicarianceCW_myr_list_cons", "vicarianceEW_myr_list_cons",
                        "vicarianceIW_myr_list_cons", "colonization_myr_list_cons", "extirpation_myr_list_cons")
event_plots_genus <- list(event_myr_list_cons, vicariance_myr_list_cons, Ccol_myr_list_cons, Ecol_myr_list_cons, Wcol_myr_list_cons,
                            Cdiv_myr_list_cons,Ediv_myr_list_cons,Wdiv_myr_list_cons,
                            vicarianceCE_myr_list_cons, vicarianceCW_myr_list_cons, vicarianceEW_myr_list_cons,
                            vicarianceIW_myr_list_cons, colonization_myr_list_cons, extirpation_myr_list_cons)
names(event_plots_genus) <- events_myr_names
  
## GENERATE A VECTOR OF 80 FOR EACH EVENT AND EACH GENUS
for (j in 1:length(event_plots_genus)) {
  for (i in 1:length(genera2)){
    event_plots_genus[[j]][[i]] <- data.frame(Ma=c(1:80), N=0)
  }}
# Counting the maximum number of all events transition in each million year for each genus.
for (i in 1:length(genera2)){
    if (nrow(event_list_cons_genera[[i]]) > 0){
      for (j in 1:nrow(event_list_cons_genera[[i]])){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus[[1]][[i]][k,]$N <- event_plots_genus[[1]][[i]][k,]$N + 1
          }}}}}
# Counting the maximum number of all vicariances in each million year for each genus.
for (i in 1:length(genera2)){
    for (j in 1:nrow(event_list_cons_genera[[i]])){
      if (event_list_cons_genera[[i]]$vicariance_IW[j] > 0 | event_list_cons_genera[[i]]$vicariance_EW[j] > 0 |
          event_list_cons_genera[[i]]$vicariance_CW[j] > 0 | event_list_cons_genera[[i]]$vicariance_CE[j] > 0){
      vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
      for (k in 1:length(ma)){
        if (Overlap(na.omit(vi), ma[[k]]) != 0){
          event_plots_genus[[2]][[i]][k,]$N <- event_plots_genus[[2]][[i]][k,]$N + 1
        }}}}}
##### specific events PER MILLION YEAR FOR EACH GENUS (CONSENSUS) #####
# Counting the maximum number of events in each million year for each genus.
  for (i in 1:length(genera2)){
    for (j in 1:nrow(event_list_cons_genera[[i]])){
      if (event_list_cons_genera[[i]]$Ccol[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$Ccol_myr_list_cons[[i]][k,]$N <- event_plots_genus$Ccol_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$Ecol[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$Ecol_myr_list_cons[[i]][k,]$N <- event_plots_genus$Ecol_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$Wcol[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$Wcol_myr_list_cons[[i]][k,]$N <- event_plots_genus$Wcol_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$Cdiv[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$Cdiv_myr_list_cons[[i]][k,]$N <- event_plots_genus$Cdiv_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$Ediv[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$Ediv_myr_list_cons[[i]][k,]$N <- event_plots_genus$Ediv_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$Wdiv[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$Wdiv_myr_list_cons[[i]][k,]$N <- event_plots_genus$Wdiv_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$vicariance_CE[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$vicarianceCE_myr_list_cons[[i]][k,]$N <- event_plots_genus$vicarianceCE_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$vicariance_CW[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$vicarianceCW_myr_list_cons[[i]][k,]$N <- event_plots_genus$vicarianceCW_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$vicariance_EW[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$vicarianceEW_myr_list_cons[[i]][k,]$N <- event_plots_genus$vicarianceEW_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$vicariance_IW[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$vicarianceIW_myr_list_cons[[i]][k,]$N <- event_plots_genus$vicarianceIW_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$colonization[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$colonization_myr_list_cons[[i]][k,]$N <- event_plots_genus$colonization_myr_list_cons[[i]][k,]$N + 1
          }}}
      if (event_list_cons_genera[[i]]$extirpation[j] > 0){
        vi <- c(event_list_cons_genera[[i]]$min[j], event_list_cons_genera[[i]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_plots_genus$extirpation_myr_list_cons[[i]][k,]$N <- event_plots_genus$extirpation_myr_list_cons[[i]][k,]$N + 1
          }}}
    }}
  
}

# Sum all the transitions per Ma (of all the genera)
# Make a data frame with number of events in each 1 Mya period.
events_all_names <- c("event_all_cons","vicariance_all_cons","Ccol_all_cons", "Ecol_all_cons", "Wcol_all_cons",
                      "Cdiv_all_cons","Ediv_all_cons","Wdiv_all_cons", 
                      "vicarianceCE_all_cons", "vicarianceCW_all_cons", "vicarianceEW_all_cons",
                      "vicarianceIW_all_cons", "Colonization_all_cons", "Extirpation_all_cons")
event_plots_all <- vector("list", length(events_all_names))
names(event_plots_all) <- events_all_names

for (j in 1:length(event_plots_all)) {
  event_plots_all[[j]] <- data.frame(Ma=c(1:myr), N=0)
  for (i in 1:length(genera2)){
    event_plots_all[[j]]$N <- event_plots_all[[j]]$N + event_plots_genus[[j]][[i]]$N
  }}

saveRDS(event_plots_genus, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_genus.rds")
saveRDS(event_plots_all, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_all.rds")
event_plots_genus <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_genus.rds") 
event_plots_all <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_all.rds")


# All events
plot(1,type='n',xlim=c(myr,0),ylim=c(0,max(event_plots_all[[1]]$N)),xlab='Ma', ylab='N', main="All events (consensus)")
# smooth
lines(spline(event_plots_all[[1]]$Ma-1, event_plots_all[[1]]$N), type="l", col="red", lwd=5, pch=16, cex=1)
# not smooth

# Plot all smooth
plot(1,type='n',xlim=c(myr,0),ylim=c(0,max(event_plots_all[[1]]$N)),xlab='Ma', ylab='N', main="All events (consensus)")
lines(spline(event_plots_all[[1]]$Ma-1, event_plots_all[[1]]$N), type="l", col="grey", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[2]]$Ma-1, event_plots_all[[2]]$N), type="l", col="lightgreen", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[3]]$Ma-1, event_plots_all[[3]]$N), type="l", col="red", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[4]]$Ma-1, event_plots_all[[4]]$N), type="l", col="gold", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[5]]$Ma-1, event_plots_all[[5]]$N), type="l", col="blue", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[6]]$Ma-1, event_plots_all[[6]]$N), type="l", col="lightcoral", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[7]]$Ma-1, event_plots_all[[7]]$N), type="l", col="khaki1", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[8]]$Ma-1, event_plots_all[[8]]$N), type="l", col="lightblue", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[9]]$Ma-1, event_plots_all[[9]]$N), type="l", col="orange", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[10]]$Ma-1, event_plots_all[[10]]$N), type="l", col="purple", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[11]]$Ma-1, event_plots_all[[11]]$N), type="l", col="salmon", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[12]]$Ma-1, event_plots_all[[12]]$N), type="l", col="darkblue", lwd=5, pch=16, cex=1)
lines(spline(event_plots_all[[13]]$Ma-1, event_plots_all[[13]]$N), type="l", col="black", lwd=5, pch=16, cex=1)
legend_colors <- c("grey", "lightgreen", "red", "gold", "blue", "lightcoral", "khaki1", "lightblue", "orange", "purple", "salmon", "darkblue", "black")
legend_text <- events_all_names
legend_text <- gsub("_all_cons", "", legend_text)
legend_text <- gsub("event", "all events", legend_text)
legend("topleft", legend=legend_text, fill=legend_colors, bty="n")


# Stacked area plot
# https://www.r-graph-gallery.com/136-stacked-area-chart.html
Ccol_stacked <- cbind(event_plots_all$Ccol_all_cons, group="Ccol", group2="col")
Ecol_stacked <- cbind(event_plots_all$Ecol_all_cons, group="Ecol", group2="col")
Wcol_stacked <- cbind(event_plots_all$Wcol_all_cons, group="Wcol",group2="col")
Cdiv_stacked <- cbind(event_plots_all$Cdiv_all_cons, group="Cdiv", group2="div")
Ediv_stacked <- cbind(event_plots_all$Ediv_all_cons, group="Ediv", group2="div")
Wdiv_stacked <- cbind(event_plots_all$Wdiv_all_cons, group="Wdiv", group2="div")
vicCE_stacked <- cbind(event_plots_all$vicarianceCE_all_cons, group="vicCE", group2="vic")
vicCW_stacked <- cbind(event_plots_all$vicarianceCW_all_cons, group="vicCW", group2="vic")
vicEW_stacked <- cbind(event_plots_all$vicarianceEW_all_cons, group="vicEW", group2="vic")
vicIW_stacked <- cbind(event_plots_all$vicarianceIW_all_cons, group="vicIW", group2="vic")


data_col <- rbind(Ccol_stacked, Ecol_stacked, Wcol_stacked)
data_div <- rbind(Cdiv_stacked, Ediv_stacked, Wdiv_stacked)
data_vic <- rbind(vicCE_stacked, vicCW_stacked, vicEW_stacked, vicIW_stacked)
data <- rbind(data_col, data_div, data_vic)

ggplot(data_col, aes(x=Ma, y=N, fill=group)) +
  geom_area() +
  geom_vline(xintercept = 20, size=20, color = "blue", alpha = 0.4)


