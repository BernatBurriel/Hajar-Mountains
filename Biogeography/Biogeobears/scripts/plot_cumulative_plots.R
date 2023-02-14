
##### PLOT cumulative plot ###########

# Calculate how many colonization events have there been 
# Define  biogeographic scenarios 
# Establish the different scenarios. Each scenario will be a list of vectors of three elements:
# (state1, state2, state3), where state1 = state in parental node, state2 = state in descendentant 1,
# state3 = state in descendent 2.
libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape")
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
C_col <- list(c("W", "C"), c("E","C"), c("CW", "C"), c("EW", "CE"), c("EW", "CW"), c("W","CW"), c("W","CE"),c("E","CE"),c("E","CW"))
#C_div <- list(c("C", "C", "C"))
C_div <- list(c("C", "C", "C"), c("C", "CE", "C"), c("C", "C", "CE"),c("CE", "C", "CE"), c("CE", "CE", "C"), c("CE","CEW", "CEW") )

#E_col <- list(c("C","E","C"), c("C","C","E"), c("CW","E","CW"), c("CW","CW","E"), 
 #             c("EW","CW","E"), c("EW", "E","EW"), c("EW", "E", "CW"), c("EW","CW","E"),
  #            c("CE","C","E"), c("CE","E","C"), c("CE","E","CW"),c("CE","CW","E"),
   #           c("CE","E","CW"), c("CE","CW","E"))
#E_col <- list(c("C","E"), c("CW","E"),c("EW","E"), c("CE","E")) 
E_col <- list(c("C","E"), c("CW","E"),c("CW","EW"), c("CW", "CE"), c("C","CE"))
#E_div <- list(c("E", "E","E"))
E_div <- list(c("E", "E","E"),c("CE","CEW", "CEW"), c("E","CE","E"), c("E","E","CE"))

#W_col <- list(c("C", "W", "C"), c("C", "C", "W"), c("IW", "W", "I"), c("IW","I","W"), c("CW", "W", "C"), c("CW", "C", "W"),
 #             c("CW", "CE", "W"), c("CW","W","CE"), c("CEW", "CE", "W"), c("CEW", "W", "CE"))
#W_col <- list(c("C", "W"), c("IW", "W"), c("CW", "W"), c("CEW", "W"))
W_col <- list(c("C", "W"), c("I", "IW"), c("C", "CW"), c("CE", "CEW"),c("CE", "CW"))
#W_div <- list(c("W","W","W"))
W_div <- list(c("W","W","W"),c("W","EW","W"),c("W","W","EW"), c("W","EW","W"))

# this are the ones to take into account when counting intermediate states
#CE_col <- list(c("CW", "CE"), c("C","CE"), c("CEW", "CE"), c("E","CE"), c("EM", "CE")
#CW_col <- list(c("W", "CW"), c("EW","CW"),c("CE", "CW"))
#EW_col <- list(c("W","EW"), c("CW","EW"))
#IW_col <- list(c("I", "IW"))

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

CEcol_nodes_list_cons <- vector("list", length(genera))
names(CEcol_nodes_list_cons) <- genera
CEcol_nodes_list_cons[1:length(genera)] <- 0

CWcol_nodes_list_cons <- vector("list", length(genera))
names(CWcol_nodes_list_cons) <- genera
CWcol_nodes_list_cons[1:length(genera)] <- 0

EWcol_nodes_list_cons <- vector("list", length(genera))
names(EWcol_nodes_list_cons) <- genera
EWcol_nodes_list_cons[1:length(genera)] <- 0

IWcol_nodes_list_cons <- vector("list", length(genera))
names(IWcol_nodes_list_cons) <- genera
IWcol_nodes_list_cons[1:length(genera)] <- 0

}


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


#### Prepare all the data no wrap all colonization and diversification nodes in a single list for each biogeographic state
colonization_states <- list(C_col, E_col, W_col,CE_col,CW_col, EW_col, IW_col)
names_colonization <- c("Central colonization", "East colonization", "West colonization", "Central-East colonization", 
                        "Central-West colonization", "East-West colonization", "Iran-West colonization")
colonization_nodes <- list(Ccol_nodes_list_cons, Ecol_nodes_list_cons, Wcol_nodes_list_cons,CEcol_nodes_list_cons, CWcol_nodes_list_cons, EWcol_nodes_list_cons, IWcol_nodes_list_cons)
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
      }
    }
  }  
}

## All Diversification nodes (parent nodes)
diversification_states <- list(C_div, E_div, W_div)
names_diversification <- c("Central diversification", "East diversification", "West diversification")
diversification_nodes <- list(Cdiv_nodes_list_cons, Ediv_nodes_list_cons, Wdiv_nodes_list_cons)
names(diversification_nodes) <- names_diversification

for (z in 1:length(diversification_nodes)) {
  for (i in 1:length(genera)){
    for (j in min(as.phylo(consensus_genera[[i]])$edge[,1]):max(as.phylo(consensus_genera[[i]])$edge[,1])){  
      v <- c(consensus_statesnoj[[i]][j],
             consensus_statesnoj[[i]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][1],
             consensus_statesnoj[[i]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][2])
      if (list(v) %in% diversification_states[[z]]){
        diversification_nodes[[z]][[i]] <- append(diversification_nodes[[z]][[i]], j)
      }
    }
  }  
}

# Now we have 4 lists, 1. nodes_hajarscol (descendant node of HAJARS colonization event); 2.nodes_hajars_extirpation
# (descendant node extirpated from Hajars); 3. colonization_nodes (list containing colonization descendant nodes of each mountain block)
# 4. diversification_nodes (list containing diversification parent nodes of diversification in the blocks)
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
                                                  height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, CEcol=0, CWcol=0,
                                                  EWcol = 0, IWcol = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
    if (names(colonization_nodes_df[z]) == "Central colonization") {colonization_nodes_df[[z]][[i]]$Ccol = 1}
    if (names(colonization_nodes_df[z]) == "East colonization") {colonization_nodes_df[[z]][[i]]$Ecol = 1}
    if (names(colonization_nodes_df[z]) == "West colonization") {colonization_nodes_df[[z]][[i]]$Wcol = 1}
    if (names(colonization_nodes_df[z]) == "Central-East colonization") {colonization_nodes_df[[z]][[i]]$CEcol = 1}
    if (names(colonization_nodes_df[z]) == "Central-West colonization") {colonization_nodes_df[[z]][[i]]$CWcol = 1}
    if (names(colonization_nodes_df[z]) == "East-West colonization") {colonization_nodes_df[[z]][[i]]$EWcol = 1}
    if (names(colonization_nodes_df[z]) == "Iran-West colonization") {colonization_nodes_df[[z]][[i]]$IWcol = 1}
  }
}


for (z in 1:length(diversification_nodes_df)){
  for (i in 1:length(genera)){
    diversification_nodes_df[[z]][[i]] <- data.frame(node=diversification_nodes[[z]][[i]], genus=genera[i],
                                                     height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, CEcol=0, CWcol=0,
                                                     EWcol = 0, IWcol = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
    if (names(diversification_nodes_df[z]) == "Central diversification") {diversification_nodes_df[[z]][[i]]$Cdiv = 1}
    if (names(diversification_nodes_df[z]) == "East diversification") {diversification_nodes_df[[z]][[i]]$Ediv = 1}
    if (names(diversification_nodes_df[z]) == "West diversification") {diversification_nodes_df[[z]][[i]]$Wdiv = 1}
  }
}

#### JOIN ALL_TREE nodes and create also a dataframe. It will be attaach later on because we will have to keep processing the data in differentways. ##
 hajar_nodes <- list(nodes_hajarscol, nodes_hajars_extirpation)
 
 #### To continue generating the dataframe for this two datasets ####
 event_list_cons_genera <- vector("list", length(genera))
 names(event_list_cons_genera) <- genera
 
 #### Generate a dataframe for all the biogeografic events separated by genus ####
 for (i in 1:length(genera)){
     event_list_cons_genera[[i]] <- colonization_nodes_df[[1]][[i]]
     for (j in 2:length(colonization_nodes_df)) {
       event_list_cons_genera[[i]] <- rbind(event_list_cons_genera[[i]], colonization_nodes_df[[j]][[i]])
     }
     for (z in 1:length(diversification_nodes_df)) {
       event_list_cons_genera[[i]] <- rbind(event_list_cons_genera[[i]], diversification_nodes_df[[z]][[i]])
     }}
# Remove all rows with 0 in the node variable
 for (i in 1:length(event_list_cons_genera)){ 
 event_list_cons_genera[[i]] <- event_list_cons_genera[[i]][event_list_cons_genera[[i]]$node != 0, ] 
 }

#### populate the dataframe with the min, max, and height columns ####
 
 # Divergence times and 95% HPD (min and max) for each diversification node 
 # and branch length for colonization and dispersal nodes.
 for (i in 1:length(genera)){
   if (nrow(event_list_cons_genera[[i]]) > 0){
     for (j in 1:length(event_list_cons_genera[[i]]$node)){
       # diversification nodes
       if (event_list_cons_genera[[i]][j, "Cdiv"] == 1 | event_list_cons_genera[[i]][j, "Ediv"] == 1 | event_list_cons_genera[[i]][j, "Wdiv"] == 1){
         event_list_cons_genera[[i]]$min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[event_list_cons_genera[[i]]$node[j]]][1]
         event_list_cons_genera[[i]]$max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[event_list_cons_genera[[i]]$node[j]]][2]
         event_list_cons_genera[[i]]$height[j] <- as_tibble(consensus_genera[[i]])$height_median[[event_list_cons_genera[[i]]$node[j]]]
       }
       # Colonization or Dispersal nodes
       if (event_list_cons_genera[[i]][j, "Ccol"] == 1 | event_list_cons_genera[[i]][j, "Ecol"] == 1 | event_list_cons_genera[[i]][j, "Wcol"] == 1 
           | event_list_cons_genera[[i]][j, "CEcol"] == 1 | event_list_cons_genera[[i]][j, "CWcol"] == 1 | event_list_cons_genera[[i]][j, "EWcol"] == 1
           | event_list_cons_genera[[i]][j, "IWcol"] == 1){
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
 
 
 
 