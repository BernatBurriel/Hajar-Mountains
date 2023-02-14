##### PLOT cumulative plot in Simulations###########

# Calculate how many colonization events have there been 
# Define  biogeographic scenarios 
# Establish the different scenarios. Each scenario will be a list of vectors of three elements:
# (state1, state2, state3), where state1 = state in parental node, state2 = state in descendentant 1,
# state3 = state in descendent 2.

libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape","DescTools")
lapply(libs, require, character.only = TRUE)

setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/")
#setwd("C:/Users/User/Desktop/__MACOSX")
#setwd("C:/Users/User/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography")

sims_consensus_letters <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/sims_consensus_letters.rds")
consensus_genera <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/genera_trees.rds")
genera <- names(consensus_genera)
nsim = 1000

## Reading biogeographic states 
C_col <- list(c("W", "C"), c("E","C"), c("W","CW"), c("W","CE"),c("E","CE"),c("E","CW"), c("W","EW"), c("EM","EW"))
E_col <- list(c("C","E"), c("CW","E"),c("CW","EW"), c("CW", "CE"), c("C","CE"), c("W","EW"), c("EM", "EC"), c("EM","CE"))
W_col <- list(c("C", "W"), c("I", "IW"), c("C", "CW"), c("CE", "CEW"),c("CE", "CW"), c("EC","CW"), c("EM","EW"))

C_div <- list(c("C", "C", "C"), c("C", "CE", "C"), c("C", "C", "CE"),c("CE", "C", "CE"), c("CE", "CE", "C"), c("CE","CEW", "CEW"), c("CW", "C", "CW"),c("CW", "CW", "C"),
              c("CW", "C", "C"))
E_div <- list(c("E", "E","E"),c("CE","CEW", "CEW"), c("E","CE","E"), c("E","E","CE"))
W_div <- list(c("W","W","W"),c("W","EW","W"),c("W","W","EW"), c("W","EW","W"))

vicariance_CE <- list(c("CE", "C","E"), c("CE","E","C"))
vicariance_CW <-list(c("CW", "C","W"), c("CW","W","C"), c("CW", "CE", "W"), c("CW", "W","CE"))
vicariance_EW <- list(c("EW", "W","E"), c("EW","E","W"), c("EW", "E", "CW"), c("EW", "CW", "E"), c("CEW", "W", "CE"), c("CEW", "CE", "W"))
vicariance_IW<- list(c("IW","I","W"), c("IW","W","I"))

##### .-.-.-.-. COUNTING AND CATEGORIZING BIOGEOGRAPHIC EVENTS IN SIMULATIONS (CONSENSUS) #####

{Ccol_SIM_nodes_list_cons <- vector("list", length(genera))
names(Ccol_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){Ccol_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){Ccol_SIM_nodes_list_cons[[i]][[s]] <- 0}}

Cdiv_SIM_nodes_list_cons <- vector("list", length(genera))
names(Cdiv_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){Cdiv_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){Cdiv_SIM_nodes_list_cons[[i]][[s]] <- 0}}

Ecol_SIM_nodes_list_cons <- vector("list", length(genera))
names(Ecol_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){Ecol_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){Ecol_SIM_nodes_list_cons[[i]][[s]] <- 0}}

Ediv_SIM_nodes_list_cons <- vector("list", length(genera))
names(Ediv_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){Ediv_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){Ediv_SIM_nodes_list_cons[[i]][[s]] <- 0}}

Wcol_SIM_nodes_list_cons <- vector("list", length(genera))
names(Wcol_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){Wcol_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){Wcol_SIM_nodes_list_cons[[i]][[s]] <- 0}}

Wdiv_SIM_nodes_list_cons <- vector("list", length(genera))
names(Wdiv_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){Wdiv_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){Wdiv_SIM_nodes_list_cons[[i]][[s]] <- 0}}

vicariance_CE_SIM_nodes_list_cons <- vector("list", length(genera))
names(vicariance_CE_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){vicariance_CE_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){vicariance_CE_SIM_nodes_list_cons[[i]][[s]] <- 0 }}

vicariance_CW_SIM_nodes_list_cons <- vector("list", length(genera))
names(vicariance_CW_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){vicariance_CW_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){vicariance_CW_SIM_nodes_list_cons[[i]][[s]] <- 0 }}
  
vicariance_EW_SIM_nodes_list_cons <- vector("list", length(genera))
names(vicariance_EW_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){vicariance_EW_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
  for (s in 1:nsim){vicariance_EW_SIM_nodes_list_cons[[i]][[s]] <- 0 }}


vicariance_IW_SIM_nodes_list_cons <- vector("list", length(genera))
names(vicariance_IW_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){vicariance_IW_SIM_nodes_list_cons[[i]] <- vector("list", nsim) 
for (s in 1:nsim){ vicariance_IW_SIM_nodes_list_cons[[i]][[s]] <- 0 }}
}

#### Prepare all the data to wrap all colonization and diversification nodes in a single list for each biogeographic state
colonization_states <- list(C_col, E_col, W_col)
names_colonization <- c("Central colonization", "East colonization", "West colonization")
colonization_SIM_nodes <- list(Ccol_SIM_nodes_list_cons, Ecol_SIM_nodes_list_cons, Wcol_SIM_nodes_list_cons)
names(colonization_SIM_nodes) <- names_colonization

## All Colonization nodes (descendant nodes)
for (z in 1:length(colonization_SIM_nodes)) {
  for (i in 1:length(genera)){
    for (s in 1:nsim){
    for (j in sort(as.phylo(consensus_genera[[i]])$edge[,2])){  
      parent <- sims_consensus_letters[[i]][[s]][as.phylo(consensus_genera[[i]])$edge[,1][as.phylo(consensus_genera[[i]])$edge[,2]==j]]
      son <- sims_consensus_letters[[i]][[s]][j]
      d <- c(parent, son)
      if (list(d) %in% colonization_states[[z]]){
        colonization_SIM_nodes[[z]][[i]][[s]] <- append(colonization_SIM_nodes[[z]][[i]][[s]], j)
      }}}}}

## All Diversification and Vicariant nodes (parent nodes) 
diversification_states <- list(C_div, E_div, W_div,vicariance_CE,vicariance_CW, vicariance_EW, vicariance_IW)
names_diversification <- c("Central diversification", "East diversification", "West diversification", "Central-East vicariance", 
                           "Central-West vicariance", "East-West vicariance", "Iran-West vicariance")
diversification_SIM_nodes <- list(Cdiv_SIM_nodes_list_cons, Ediv_SIM_nodes_list_cons, Wdiv_SIM_nodes_list_cons,vicariance_CE_SIM_nodes_list_cons, 
                              vicariance_CW_SIM_nodes_list_cons, vicariance_EW_SIM_nodes_list_cons, vicariance_IW_SIM_nodes_list_cons)
names(diversification_SIM_nodes) <- names_diversification

for (z in 1:length(diversification_SIM_nodes)) {
  for (i in 1:length(genera)){
    for (s in 1:nsim){
    for (j in min(as.phylo(consensus_genera[[i]])$edge[,1]):max(as.phylo(consensus_genera[[i]])$edge[,1])){  
      v <- c(sims_consensus_letters[[i]][[s]][j],
             sims_consensus_letters[[i]][[s]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][1],
             sims_consensus_letters[[i]][[s]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][2])
      if (list(v) %in% diversification_states[[z]]){
        diversification_SIM_nodes[[z]][[i]][[s]] <- append(diversification_SIM_nodes[[z]][[i]][[s]], j)
      }}}}}
# Now we have 4 lists, 1. nodes_hajarscol (descendant node of HAJARS colonization event); 2.nodes_hajars_extirpation
# (descendant node extirpated from Hajars); 3. colonization_nodes (list containing colonization descendant nodes of each mountain block)
# 4. diversification_nodes (list containing diversification parent nodes of diversification in the blocks, or vicariance between them)
# Lists 1 and 2 contain data for all the species (BIG_tree)
# Lists 3 and 4 are lists of lists containing the each colonization (3) or diversification (4) event in each bloack and in each genus
# First, we're going to remove the zeros from the vectors with transition nodes.

for (z in 1:length(colonization_SIM_nodes)) {
  for (i in 1:length(genera)){
    for (s in 1:nsim){
    if (length(colonization_SIM_nodes[[z]][[i]][[s]]) != 1){
      colonization_SIM_nodes[[z]][[i]][[s]] <- colonization_SIM_nodes[[z]][[i]][[s]][-1]
    }}}}

for (z in 1:length(diversification_SIM_nodes)) {
  for (i in 1:length(genera)){
    for (s in 1:nsim){
    if (length(diversification_SIM_nodes[[z]][[i]][[s]]) != 1){
      diversification_SIM_nodes[[z]][[i]][[s]] <- diversification_SIM_nodes[[z]][[i]][[s]][-1]
    }}}}

#### Create a dataframe for each colonization event and afterwards join all of them together

colonization_SIM_nodes_df <- vector('list', length(colonization_SIM_nodes))
names(colonization_SIM_nodes_df) <- names(colonization_SIM_nodes)
for (i in 1:length(colonization_SIM_nodes_df)) {
  colonization_SIM_nodes_df[[i]] <- vector('list', length(genera))
  names(colonization_SIM_nodes_df[[i]]) <- genera
  for (j in 1:length(genera)){
    colonization_SIM_nodes_df[[i]][[j]] <- vector("list", nsim)
  }
}
diversification_SIM_nodes_df <- vector('list', length(diversification_SIM_nodes))
names(diversification_SIM_nodes_df) <- names(diversification_SIM_nodes)
for (i in 1:length(diversification_SIM_nodes_df)) {
  diversification_SIM_nodes_df[[i]] <- vector('list', length(genera))
  names(diversification_SIM_nodes_df[[i]]) <- genera
  for (j in 1:length(genera)) {
    diversification_SIM_nodes_df[[i]][[j]] <- vector("list", nsim)    
  }
}

for (z in 1:length(colonization_SIM_nodes_df)){
  for (i in 1:length(genera)){
    for (s in 1:nsim){
    colonization_SIM_nodes_df[[z]][[i]][[s]] <- data.frame(node=colonization_SIM_nodes[[z]][[i]][[s]], genus=genera[i],
                                                  height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                                                  vicariance_EW = 0, vicariance_IW = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
    if (names(colonization_SIM_nodes_df[z]) == "Central colonization") {colonization_SIM_nodes_df[[z]][[i]][[s]]$Ccol = 1}
    if (names(colonization_SIM_nodes_df[z]) == "East colonization") {colonization_SIM_nodes_df[[z]][[i]][[s]]$Ecol = 1}
    if (names(colonization_SIM_nodes_df[z]) == "West colonization") {colonization_SIM_nodes_df[[z]][[i]][[s]]$Wcol = 1}
  }}}


for (z in 1:length(diversification_SIM_nodes_df)){
  for (i in 1:length(genera)){
    for (s in 1:nsim){
    diversification_SIM_nodes_df[[z]][[i]][[s]] <- data.frame(node=diversification_SIM_nodes[[z]][[i]][[s]], genus=genera[i],
                                                     height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                                                     vicariance_EW = 0, vicariance_IW = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
    if (names(diversification_SIM_nodes_df[z]) == "Central diversification") {diversification_SIM_nodes_df[[z]][[i]][[s]]$Cdiv = 1}
    if (names(diversification_SIM_nodes_df[z]) == "East diversification") {diversification_SIM_nodes_df[[z]][[i]][[s]]$Ediv = 1}
    if (names(diversification_SIM_nodes_df[z]) == "West diversification") {diversification_SIM_nodes_df[[z]][[i]][[s]]$Wdiv = 1}
    if (names(diversification_SIM_nodes_df[z]) == "Central-East vicariance") {diversification_SIM_nodes_df[[z]][[i]][[s]]$vicariance_CE = 1}
    if (names(diversification_SIM_nodes_df[z]) == "Central-West vicariance") {diversification_SIM_nodes_df[[z]][[i]][[s]]$vicariance_CW = 1}
    if (names(diversification_SIM_nodes_df[z]) == "East-West vicariance") {diversification_SIM_nodes_df[[z]][[i]][[s]]$vicariance_EW = 1}
    if (names(diversification_SIM_nodes_df[z]) == "Iran-West vicariance") {diversification_SIM_nodes_df[[z]][[i]][[s]]$vicariance_IW = 1}
  }}}

# We want a list with one data frame per genus per simulation, including all types of events.
sim_event_consensus <- vector("list", length(genera))
names(sim_event_consensus) <- genera
for (i in 1:length(genera)){
  sim_event_consensus[[i]] <- vector("list", nsim)
}

#### Generate a dataframe for all the biogeografic events separated by genus and simulation####
for (i in 1:length(genera)){
  for (s in 1:nsim) {
    for (j in 1:length(colonization_SIM_nodes_df)) {
    sim_event_consensus[[i]][[s]] <- rbind(sim_event_consensus[[i]][[s]], colonization_SIM_nodes_df[[j]][[i]][[s]])
  }}
  for(k in 1:nsim){
    for (z in 1:length(diversification_SIM_nodes_df)) {
    sim_event_consensus[[i]][[k]] <- rbind(sim_event_consensus[[i]][[k]], diversification_SIM_nodes_df[[z]][[i]][[k]])
  }}}

sim_event_consensus$Asaccus[[453]]
sim_event_consensus$Trachydactylus[[231]]

# Remove all rows with 0 in the node variable
for (i in 1:length(sim_event_consensus)){ 
  for (s in 1:nsim) {
    sim_event_consensus[[i]][[s]] <- sim_event_consensus[[i]][[s]][sim_event_consensus[[i]][[s]]$node != 0, ] 
  }
}
sim_event_consensus$Asaccus[[453]]
sim_event_consensus$Trachydactylus[[231]]
##### INCLUDE INFO ON DIVERGENCE TIMES FOR SIMS #####
#### Populate the dataframe with the min, max, and height columns

# Divergence times and 95% HPD (min and max) for each diversification node 
# and branch length for colonization and dispersal nodes.
for (i in 1:length(genera)){
  for (s in 1:nsim) {
  if (nrow(sim_event_consensus[[i]][[s]]) > 0){
    for (j in 1:length(sim_event_consensus[[i]][[s]]$node)){
      # diversification nodes
      if (sim_event_consensus[[i]][[s]][j, "Cdiv"] == 1 | sim_event_consensus[[i]][[s]][j, "Ediv"] == 1 | sim_event_consensus[[i]][[s]][j, "Wdiv"] == 1
          | sim_event_consensus[[i]][[s]][j, "vicariance_CE"] == 1 | sim_event_consensus[[i]][[s]][j, "vicariance_CW"] == 1 | sim_event_consensus[[i]][[s]][j, "vicariance_EW"] == 1
          | sim_event_consensus[[i]][[s]][j, "vicariance_IW"] == 1){
        sim_event_consensus[[i]][[s]]$min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[sim_event_consensus[[i]][[s]]$node[j]]][1]
        sim_event_consensus[[i]][[s]]$max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[sim_event_consensus[[i]][[s]]$node[j]]][2]
        sim_event_consensus[[i]][[s]]$height[j] <- as_tibble(consensus_genera[[i]])$height_median[[sim_event_consensus[[i]][[s]]$node[j]]]
      }
      # Colonization or Dispersal nodes
      if (sim_event_consensus[[i]][[s]][j, "Ccol"] == 1 | sim_event_consensus[[i]][[s]][j, "Ecol"] == 1 | sim_event_consensus[[i]][[s]][j, "Wcol"] == 1){
        dat <- as_tibble(consensus_genera[[i]])
        # Min is the age of the descendent node + its 95% HPD
        desc_age <- dat$height_0.95_HPD[dat$node == sim_event_consensus[[i]][[s]]$node[j]][[1]][1]
        sim_event_consensus[[i]][[s]]$min[j] <- desc_age
        # Max is the age of the parental node + its 95% HPD
        parent_age <- dat$height_0.95_HPD[dat$parent[dat$node == sim_event_consensus[[i]][[s]]$node[j]]][[1]][2]
        sim_event_consensus[[i]][[s]]$max[j] <- parent_age
        # let's set the height as the midpoint between the parental and the descendent nodes
        sim_event_consensus[[i]][[s]]$height[j] <- mean(c(desc_age, parent_age))
      }}}}}


saveRDS(sim_event_consensus, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/sim_event_consensus.rds")
sim_event_consensus <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/sim_event_consensus.rds")


##### EVENTS PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
myr <- 80
ma <- vector("list", myr)
for(i in 1:myr){
  ma[[i]] <- c(i-1, i)
}



# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
sim_event_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  sim_event_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    sim_event_all_cons[[s]]$N <- sim_event_all_cons[[s]]$N + sim_event_myr_list_cons[[i]][[s]]$N
  }}

#### Create all the lists containing the data to plot  
{event_myr_SIM_list_cons <- vector('list', length(sim_event_consensus)) # All genera, all events, each simulation
names(event_myr_SIM_list_cons) <- names(sim_event_consensus)
vicariance_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, vicariance events, each simulation
names(vicariance_myr_SIM_list_cons) <- genera
Ccol_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, Ccol events, each simulation
names(Ccol_myr_SIM_list_cons) <- genera
Ecol_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, Ecol events, each simulation
names(Ecol_myr_SIM_list_cons) <- genera
Wcol_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, Wcol events, each simulation
names(Wcol_myr_SIM_list_cons) <- genera
Cdiv_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, Cdiv events, each simulation
names(Cdiv_myr_SIM_list_cons) <- genera
Ediv_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, Ediv events, each simulation
names(Ediv_myr_SIM_list_cons) <- genera
Wdiv_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, Wdiv events, each simulation
names(Wdiv_myr_SIM_list_cons) <- genera
vicarianceCE_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, VicCE events, each simulation
names(vicarianceCE_myr_SIM_list_cons) <- genera
vicarianceCW_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, VicCW events, each simulation
names(vicarianceCW_myr_SIM_list_cons) <- genera
vicarianceEW_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, VicEW events, each simulation
names(vicarianceEW_myr_SIM_list_cons) <- genera
vicarianceIW_myr_SIM_list_cons <- vector("list", length(genera)) #All genera, VicIW events, each simulation
names(vicarianceIW_myr_SIM_list_cons) <- genera
}
  
  
  events_myr_SIM_names <- c("event_myr_SIM_list_cons","vicariance_myr_SIM_list_cons","Ccol_myr_SIM_list_cons", "Ecol_myr_SIM_list_cons", "Wcol_myr_SIM_list_cons",
                        "Cdiv_myr_SIM_list_cons","Ediv_myr_SIM_list_cons","Wdiv_myr_SIM_list_cons", 
                        "vicarianceCE_myr_SIM_list_cons", "vicarianceCW_myr_SIM_list_cons", "vicarianceEW_myr_SIM_list_cons",
                        "vicarianceIW_myr_SIM_list_cons")
  event_SIM_plots_genus <- list(event_myr_SIM_list_cons, vicariance_myr_SIM_list_cons, Ccol_myr_SIM_list_cons, Ecol_myr_SIM_list_cons, Wcol_myr_SIM_list_cons,
                            Cdiv_myr_SIM_list_cons,Ediv_myr_SIM_list_cons,Wdiv_myr_SIM_list_cons,
                            vicarianceCE_myr_SIM_list_cons, vicarianceCW_myr_SIM_list_cons, vicarianceEW_myr_SIM_list_cons,
                            vicarianceIW_myr_SIM_list_cons)
  names(event_SIM_plots_genus) <- events_myr_SIM_names

  
## GENERATE A VECTOR OF 80 FOR EACH EVENT AND EACH GENUS AND EACH SIMULATION
for (j in 1:length(event_SIM_plots_genus)) {
    for (i in 1:length(genera)){
      event_SIM_plots_genus[[j]][[i]] <- vector("list", nsim)}
 }  

for (j in 1:length(event_SIM_plots_genus)) {
  for (i in 1:length(genera)){    
for (s in 1:nsim) {
    event_SIM_plots_genus[[j]][[i]][[s]] <- data.frame(Ma=c(1:80), N=0)
  }}}

event_SIM_plots_genus$Ccol_myr_SIM_list_cons$Asaccus[234]
  
  # Counting the maximum number of all events transition in each million year for each genus.
  for (i in 1:length(genera)){
    for (s in 1:nsim) {
          if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus[[1]][[i]][[s]][k,]$N <- event_SIM_plots_genus[[1]][[i]][[s]][k,]$N + 1
          }}}}}}
  
# Counting the maximum number of all vicariances in each million year for each genus.
  for (i in 1:length(genera)){
    for (s in 1:nsim) {
      if (nrow(sim_event_consensus[[i]][[s]]) > 0){
    for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
      if (sim_event_consensus[[i]][[s]]$vicariance_IW[j] > 0 | sim_event_consensus[[i]][[s]]$vicariance_EW[j] > 0 |
          sim_event_consensus[[i]][[s]]$vicariance_CW[j] > 0 | sim_event_consensus[[i]][[s]]$vicariance_CE[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus[[2]][[i]][[s]][k,]$N <- event_SIM_plots_genus[[2]][[i]][[s]][k,]$N + 1
          }}}}}}}
  
##### specific events PER MILLION YEAR FOR EACH GENUS (CONSENSUS) #####
  # Counting the maximum number of events in each million year for each genus.
  for (i in 1:length(genera)){
    for (s in 1:nsim) {
      if (nrow(sim_event_consensus[[i]][[s]]) > 0){
    for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
      if (sim_event_consensus[[i]][[s]]$Ccol[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$Ccol_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$Ccol_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$Ecol[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$Ecol_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$Ecol_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$Wcol[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$Wcol_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$Wcol_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$Cdiv[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$Cdiv_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$Cdiv_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$Ediv[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$Ediv_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$Ediv_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$Wdiv[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$Wdiv_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$Wdiv_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$vicariance_CE[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$vicarianceCE_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$vicarianceCE_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$vicariance_CW[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$vicarianceCW_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$vicarianceCW_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$vicariance_EW[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$vicarianceEW_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$vicarianceEW_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
      if (sim_event_consensus[[i]][[s]]$vicariance_IW[j] > 0){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            event_SIM_plots_genus$vicarianceIW_myr_SIM_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$vicarianceIW_myr_SIM_list_cons[[i]][[s]][k,]$N + 1
          }}}
#      if (sim_event_consensus[[i]][[s]]$colonization[j] > 0){
#        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
#        for (k in 1:length(ma)){
#          if (Overlap(na.omit(vi), ma[[k]]) != 0){
#            event_SIM_plots_genus$colonization_myr_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$colonization_myr_list_cons[[i]][[s]][k,]$N + 1
#          }}}
#      if (sim_event_consensus[[i]][[s]]$extirpation[j] > 0){
#        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
#        for (k in 1:length(ma)){
#          if (Overlap(na.omit(vi), ma[[k]]) != 0){
#            event_SIM_plots_genus$extirpation_myr_list_cons[[i]][[s]][k,]$N <- event_SIM_plots_genus$extirpation_myr_list_cons[[i]][[s]][k,]$N + 1
#          }}}
    }}}}

  
  # Sum all the transitions per Ma (of all the genera)
  # Make a data frame with number of events in each 1 Mya period.
  events_all_names <- c("event_all_cons","vicariance_all_cons","Ccol_all_cons", "Ecol_all_cons", "Wcol_all_cons",
                        "Cdiv_all_cons","Ediv_all_cons","Wdiv_all_cons", 
                        "vicarianceCE_all_cons", "vicarianceCW_all_cons", "vicarianceEW_all_cons",
                        "vicarianceIW_all_cons")#, "Colonization_all_cons", "Extirpation_all_cons")
  event_SIM_plots_all <- vector("list", length(events_all_names))
  names(event_SIM_plots_all) <- events_all_names
  
  for (j in 1:length(event_SIM_plots_all)) {
    event_SIM_plots_all[[j]] <- vector("list", nsim)
  for (s in 1:nsim){
    event_SIM_plots_all[[j]][[s]] <- data.frame(Ma=c(1:myr), N=0)
  }}
  
  for (j in 1:length(event_SIM_plots_all)) {
    for (s in 1:nsim) {
      for (i in 1:length(genera)){
    event_SIM_plots_all[[j]][[s]]$N <- event_SIM_plots_all[[j]][[s]]$N + event_SIM_plots_genus[[j]][[i]][[s]]$N
    }}}
  
  saveRDS(event_SIM_plots_genus, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_SIM_plots_genus.rds")
  saveRDS(event_SIM_plots_all, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_SIM_plots_all.rds")
  event_SIM_plots_genus <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_SIM_plots_genus.rds") 
  event_SIM_plots_all <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_SIM_plots_all.rds")
  event_plots_all <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_all.rds")
  event_plots_genus <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_genus.rds")
  event_plots_all_biogeo <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_all_biogeo.rds")
    ##### PLOT SIMULATIONS ####
  plot(1,type='n',xlim=c(myr,0),ylim=c(0,60),xlab='Ma', ylab='N', main="All simulated events (consensus)")
  # smooth
  for (i in 1:nsim){
    lines(spline(event_SIM_plots_all[[1]][[i]]$Ma-1, event_SIM_plots_all[[1]][[i]]$N),
          type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
  }
  lines(spline(event_plots_all[[1]]$Ma-1, event_plots_all[[1]]$N), type="l", col="black", lwd=5, pch=16, cex=1)

  
  # 95% CI ALL EVENTS
  ## quantiles and mean of each My.
  prob_qup <- 0.975
  prob_qlow <- 0.025
  
  # Let's do a list where each element is a vector with the
  # number of transitions per simulations for each Ma.
  # 60 vectors (60 Ma), and 1000 numbers in each vector.
  dlist <- vector("list", length(event_SIM_plots_all))
  names(dlist) <- names(event_SIM_plots_all)
  for (j in 1:length(dlist)) {
    dlist[[j]] <- vector("list", myr)  
    for (i in 1:myr){
      dlist[[j]][[i]] <- vector("numeric", nsim)
    }}
  
  
  dlist[[1]][[1]][2]
  
  for (j in 1:length(dlist)) {
   for (s in 1:nsim){
    for (i in 1:myr){
      dlist[[j]][[i]][s] <- event_SIM_plots_all[[j]][[s]]$N[i]
    }}}
  dlist[[1]][[1]]
  # Now we can create a dataframe with the quantiles and
  # the mean per Ma, which we will calculate with the different
  # elements of the list of distributions per Ma (dlist).
  Q_df <- vector("list", length(event_SIM_plots_all))
  names(Q_df) <- names(event_SIM_plots_all)
  for (j in 1:length(Q_df)) {
    Q_df[[j]] <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)  
  }
 
  for (j in 1:length(dlist)) {
  for (i in 1:myr){
    Q_df[[j]]$mean[i] <- mean(dlist[[j]][[i]])
    Q_df[[j]]$qlow[i] <- quantile(dlist[[j]][[i]], probs=c(prob_qlow, prob_qup))[1]
    Q_df[[j]]$qup[i] <- quantile(dlist[[j]][[i]], probs=c(prob_qlow, prob_qup))[2]
  }}
  
  # PLOT SMOOTH
  plot(1,type='n',xlim=c(40,0),ylim=c(0,30),xlab='Ma', ylab='N', main="Simulated transitions")
  for (i in 1:nsim){
    lines(spline(event_SIM_plots_all[[2]][[i]]$Ma-1, event_SIM_plots_all[[2]][[i]]$N),
          type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
  }
  #lines(as.data.frame(spline(Q_df$Ma-1, Q_df$mean)), type="l", col="purple", lwd=2, pch=16, cex=1)
  lines(spline(Q_df[[2]]$Ma-1, Q_df[[2]]$mean), type="l", col="purple", lwd=2, pch=16, cex=1)
  lines(spline(Q_df[[2]]$Ma-1, Q_df[[2]]$qlow), type="l", col="black", lwd=1, pch=16, cex=1)
  lines(spline(Q_df[[2]]$Ma-1, Q_df[[2]]$qup), type="l", col="black", lwd=1, pch=16, cex=1)

  # empiric
  lines(spline(event_plots_all_biogeo[[2]]$Ma-1, event_plots_all_biogeo[[2]]$N), type="l", col="red", lwd=5, pch=16, cex=1)
  
  # PLOT NO SMOOTH
  
  plot(1,type='n',xlim=c(myr,0),ylim=c(0,60),xlab='Ma', ylab='N', main="Simulated transitions")
  for (i in 1:nsim){
    lines(event_SIM_plots_all[[1]][[i]]$Ma-1, event_SIM_plots_all[[1]][[i]]$N,
          type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
  }
  #lines(as.data.frame(spline(Q_df$Ma-1, Q_df$mean)), type="l", col="purple", lwd=2, pch=16, cex=1)
  lines(Q_df$Ma-1, Q_df$mean, type="l", col="purple", lwd=2, pch=16, cex=1)
  lines(Q_df$Ma-1, Q_df$qlow, type="l", col="black", lwd=1, pch=16, cex=1)
  lines(Q_df$Ma-1, Q_df$qup, type="l", col="black", lwd=1, pch=16, cex=1)
  
  # empiric
  lines(event_plots_all_biogeo[[1]]$Ma-1, event_plots_all_biogeo[[1]]$N, type="l", col="red", lwd=5, pch=16, cex=1)
  