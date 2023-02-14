## Author Bernat Burriel
# Lineage accumulation plot

# load all necessary libraries
rm(list = ls())
libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape","DescTools")
lapply(libs, require, character.only = TRUE)

setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/")

#setwd("C:/Users/User/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography")

big_tree <- read.beast("Mountain_colonization_Simmaps/TREES/SQUAMATA_ALL_PF.tree")
big_tree@phylo$tip.label <- gsub("_", " ", big_tree@phylo$tip.label)
big_tree_phylo <- as.phylo(big_tree)

mountains_simmap <- readRDS("Mountain_colonization_Simmaps/simmaps_all_tree/objects/mountains_simmap_alltree.rds")
colonization <- list(c("Out", "Out", "Hajars"), c("Out", "Hajars", "Out"))
extirpation <- list (c("Hajars", "Out", "Hajars"), c("Hajars",  "Hajars", "Out"), c("Hajars", "Out", "Out"))
consensus_genera <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/genera_trees.rds")
genera = names(consensus_genera)

# Create a data frame for each genera to store each internal node
nodes_genera_df <- vector('list', length(genera))
names(nodes_genera_df) <- genera
for (i in 1:length(genera)){
  nodes_genera_df[[i]] <- data.frame(node=(length(consensus_genera[[i]]@phylo$tip.label)+1):(length(consensus_genera[[i]]@phylo$tip.label) + consensus_genera[[i]]@phylo$Nnode), 
  genus=genera[i],  height=0, min=0, max=0, diversification = 1, first_col = 0)
}

# Add node height information in each genera
for (i in 1:length(genera)){
  for (j in 1:length(nodes_genera_df[[i]]$node)){
    nodes_genera_df[[i]]$min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[nodes_genera_df[[i]]$node[j]]][1]
    nodes_genera_df[[i]]$max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[nodes_genera_df[[i]]$node[j]]][2]
    nodes_genera_df[[i]]$height[j] <- as_tibble(consensus_genera[[i]])$height_median[[nodes_genera_df[[i]]$node[j]]]
  }}

# delete non hajar nodes
nodes_genera_df[[1]] <- nodes_genera_df[[1]][nodes_genera_df[[1]]$node != 17, ]
nodes_genera_df[[1]] <- nodes_genera_df[[1]][nodes_genera_df[[1]]$node != 18, ]
nodes_genera_df[[8]] <- nodes_genera_df[[8]][nodes_genera_df[[8]]$node != 3, ]
 
# Now we will add the data on P.Gall
# Before saving the object we will add the data of Pristurus gallagheri central diversification
pgall <- read.beast("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/species_only_two_tips/06Pgall.tree")
names_tbl <- read.table("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/species_only_two_tips/06Pgall_name_correspondance.txt", header = T)
pgall <- treeio::rename_taxa(pgall, data = names_tbl, key = Name, value = New_name)
tips_to_drop <- names_tbl[names_tbl$Out_In == "out",]
pgall <- treeio::drop.tip(pgall, tips_to_drop$New_name)
pgalldt <- as_tibble(pgall)

Pgall_df <- data.frame(node=3, genus="Pgall",
                       height=pgalldt$height_median[3], min=pgalldt$height_0.95_HPD[[3]][1], max=pgalldt$height_0.95_HPD[[3]][2],
                       diversification= 1, first_col = 0)
rownames(Pgall_df) <- 1
nodes_genera_df[[9]] <-  Pgall_df
consensus_genera[[9]] <- pgall
names(consensus_genera) <- c(genera, "Pgall")
names(nodes_genera_df) <- c(genera, "Pgall")
genera2 = c(genera , "Pgall")

# Now we will add the data on Pseudotrapelus
# Before saving the object we will add the data of Pristurus gallagheri central diversification
Pseudotrapelus <- read.beast("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/species_only_two_tips/pseudotrapelus_annotatedR3.tree")
names_tbl <- read.table("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/old_files/08Pseudotrapelus.txt", header = T)
Pseudotrapelus <- treeio::rename_taxa(Pseudotrapelus, data = names_tbl, key = Name, value = New_name)
tips_to_drop <- names_tbl[names_tbl$Out_In == "out",]
Pseudotrapelus <- treeio::drop.tip(Pseudotrapelus, tips_to_drop$New_name)
Pseudotrapelusdt <- as_tibble(Pseudotrapelus)

Pseudotrapelus_df <- data.frame(node=3, genus="Pseudotrapelus",
                       height=pgalldt$height_median[3], min=pgalldt$height_0.95_HPD[[3]][1], max=pgalldt$height_0.95_HPD[[3]][2],
                       diversification= 1, first_col = 0)
rownames(Pseudotrapelus_df) <- 1
nodes_genera_df[[10]] <-  Pseudotrapelus_df
consensus_genera[[10]] <- Pseudotrapelus
names(consensus_genera) <- c(genera2, "Pseudotrapelus")
names(nodes_genera_df) <- c(genera2, "Pseudotrapelus")
genera2 = c(genera2 , "Pseudotrapelus")


#### MOUNTAIN COLONIZATION NODES ####
# get the first colonization nodes
nodes_hajarscol <- c()
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

genera_order <- c("Omanosaura", "Asaccus", "Ptyodactylus", "Hemidactylus", "Prup", "Echis", "Pcele", "Pgall", "Pseudotrapelus", "Trachydactylus")
first_col_nodes <- data.frame(node=nodes_hajarscol, genus=genera_order, height=0, min=0, max=0, diversification= 0, first_col = 1)


# Divergence times and branch length for colonization nodes.
      for (j in 1:length(first_col_nodes$node)){
      # Colonization or Dispersal nodes
        dat <- as_tibble(big_tree)
        # Min is the age of the descendent node + its 95% HPD
        first_col_nodes$min[j] <- dat$height_0.95_HPD[dat$node == first_col_nodes$node[j]][[1]][1]
        # Max is the age of the parental node + its 95% HPD
        first_col_nodes$max[j] <- dat$height_0.95_HPD[dat$parent[dat$node == first_col_nodes$node[j]]][[1]][2]
        # let's set the height as the midpoint between the parental and the descendent nodes
        first_col_nodes$height[j] <- mean(c(first_col_nodes$min[j], first_col_nodes$max[j]))
      }
dat[dat$node %in% first_col_nodes$node,]
first_col_nodes_df <- first_col_nodes

#for (i in 1:length(consensus_genera)) {
#  for (j in 1:length(first_col_nodes$node)) {
#    if (first_col_nodes$genus[j] == names(consensus_genera[i])) {
#    first_col_nodes_df$min[j] <- max(nodes_genera_df[[i]]$max) 
#  }}
#}

nodes_genera_df[[11]] <-  first_col_nodes_df
names(nodes_genera_df) <- c(genera2, "first_col")
genera2 <- c(genera2, "first_col")
# bind all genera in one file
events_div <- data.frame(  node = 0,        genus = 0,    height = 0,       min = 0,       max = 0, diversification = 0, first_col = 0)
for( i in 1:length(nodes_genera_df)){
  events_div <-  rbind(events_div, nodes_genera_df[[i]])
}
# remove 0 in first row
events_div <- events_div[-1,]
NROW(events_div)
saveRDS(events_div, "Mountain_colonization_Simmaps/objects/events_div.rds")

########## Plot the results ########

# Now we will do the consensus by million years of each genus
# The deepest colonization event in the dataset is the one given by the Asaccus group and pristurus celerrimus. Its node is around 40-45 Mya ago. Just to be sure, we will reconstruct from the first knowledge of 
# mountain formation 80 Mya
myr <- 80
ma <- myr
genera2

# We will create a list of vectors representing intervals of 1 Ma.
ma <- vector("list", ma)
for (i in 1:length(ma)){
  ma[[i]] <- c(i-1, i)
}

# Create a list containning the data to plot 

event_myr_list_cons <- vector('list', length(genera2))

## GENERATE A VECTOR OF 80 FOR EACH EVENT AND EACH GENUS
  for (i in 1:length(genera2)){
    event_myr_list_cons[[i]] <- data.frame(Ma=c(0:80), N=0)
  }
# Counting the maximum number of all events transition in each million year for each genus.
for (i in 1:length(genera2)){
  if (nrow(nodes_genera_df[[i]]) > 0){
    for (j in 1:nrow(nodes_genera_df[[i]])){
      vi <- c(nodes_genera_df[[i]]$min[j], nodes_genera_df[[i]]$max[j])
      for (k in 1:length(ma)){
        if (Overlap(na.omit(vi), ma[[k]]) != 0){
          event_myr_list_cons[[i]][k,]$N <- event_myr_list_cons[[i]][k,]$N + 1
        }}}}}

# Sum all the transitions per Ma (of all the genera)
# Make a data frame with number of events in each 1 Mya period.

event_plots_all <- data.frame(Ma=c(0:myr), N=0)
for (i in 1:length(genera2)){
  event_plots_all$N <- event_plots_all$N + event_myr_list_cons[[i]]$N
}

saveRDS(event_myr_list_cons, "Mountain_colonization_Simmaps/objects/genera_lineage_accumulation.rds")
saveRDS(event_plots_all, "Mountain_colonization_Simmaps/objects/genera_lineage_accumulation.rds")

#### CUMULATIVE PLOT ####
#sort(event_df$min)
event_df_sorted <- events_div[rev(order(events_div$max)),]
xx <- event_df_sorted[order(event_df_sorted$genus),]

##

time_vec <- seq(from=80, to=0, by=-1)

rmat <- matrix(NA, ncol = 2, nrow=length(time_vec))
colnames(rmat) <- c("Ma", "Ncum")
rmat[, "Ma"] <- time_vec

nrow(xx[xx$max >= time_vec[20], ])
for (ii in 1:length(time_vec)){
  
  tmp <- xx[xx$max >= time_vec[ii], ]
  rmat[ii, "Ncum"] <- nrow(tmp)
}

saveRDS(rmat, "Mountain_colonization_Simmaps/objects/rmat_first_col.rds")

libs <- c("tidyverse", "deeptime", "here", "cowplot")
lapply(libs, require, character.only = TRUE)
# Theme ----
# Set a customized theme for all the plots
theme_htc <- function(){
  theme_bw() +
    theme(panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          #        axis.text.x = element_blank(),
          #        axis.text.y = element_text(size = 5),
          #        axis.title.x = element_blank(),
          #        axis.title.y = element_blank(),
          plot.title = element_text(size = 13, vjust = 1, hjust = 0.5)
    )
}

# Geologic timescale ----
# Set the geologic period information 
data(periods)
data(epochs)

periods_htc <- periods
periods_htc$name[1] <- "Q"

epochs_htc <- epochs
epochs_htc$abbr[epochs_htc$abbr == "Plicn"] <- "Pli"
epochs_htc$abbr[epochs_htc$abbr == "Pls"] <- "Ple"
epochs_htc$abbr[epochs_htc$abbr == "Eo"] <- "Eocene"
epochs_htc$abbr[epochs_htc$abbr == "Mc"] <- "Miocene"
epochs_htc$abbr[epochs_htc$abbr == "Pal"] <- "Paleocene"
epochs_htc$abbr[epochs_htc$abbr == "Ol"] <- "Oligocene"

lwd_events <- 2


rmat_df <- as.data.frame(rmat)

# Create the lines to plot ----
# spline consensus
line_all <- data.frame(spline(rmat_df$Ma, rmat_df$Ncum)) 
colnames(line_all) <- c("Ma", "Ncum")
  
# create first colonization segments 
first_col <- data.frame(X0 = 0, X1 = 0) 

for (i in 1:nrow(xx[xx$first_col == 1,])){
  dt <- xx[xx$first_col == 1,][i,]
  first_col[i,]$X1 = dt$min
  first_col[i,]$X0 = dt$max
  }
  
first_col <- first_col[rev(order(first_col$X0)),]
first_col$Y <- seq(4, 40, 4)
first_col$event <- c(1:10)

plot_all <- ggplot() +
  

  # Add First colonization for each genera
  geom_segment(data = first_col, aes(x = X0, xend= X1, y=Y, yend = Y), lineend = "round",
               alpha = 0.4, color = "#92BA8C", size = 5 ) +
  # "#63A05A"
  # "#92BA8C"
  # "#B2D7AE"
  
  # Add colonization events
  geom_point(data = first_col, aes(x = (X0-0), y = (Y+1)), colour = "#D9EDD6", alpha = 0.9, size = 4) +
  geom_text(data = first_col, aes(x = (X0-0), y = (Y+1)), label = first_col$event, size = 2, fontface = "bold") +
  
  # Add orogeny events
  geom_rect(aes(xmin=70, xmax=50, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
  geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
  geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
  geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
  #geom_text(aes(x = 62, y = max(rmat_df$Ncum - 12), angle = 90), label = "Exhumation of high pressure rocks", size = 3.5) +
 # geom_text(aes(x = 61, y = max(rmat_df$Ncum - 12), angle = 90), label = "Proto-Jebel Akhdar and proto-Saih", size = 3, fontface = "bold") +
#  geom_text(aes(x = 59.5, y = max(rmat_df$Ncum - 12), angle = 90), label = "Hatat culminations", size = 3, fontface = "bold") +
#  geom_text(aes(x = 36, y = max(rmat_df$Ncum - 12), angle = 90), label = "Hajar Mountain's main", size = 3, fontface = "bold") +
#  geom_text(aes(x = 34.5, y = max(rmat_df$Ncum - 12), angle = 90), label = "uplift activity", size = 3, fontface = "bold")   +
#  geom_text(aes(x = 18.3, y = max(rmat_df$Ncum - 12), angle = 90), label = "Northern Hajar's", size = 3, fontface = "bold")   +
#  geom_text(aes(x = 16.8, y = max(rmat_df$Ncum - 12), angle = 90), label = "secondary uplift", size = 3, fontface = "bold")   +
#  geom_text(aes(x = 7.5, y = max(rmat_df$Ncum - 12), angle = 90), label = "Jebel Akhdar", size = 3, fontface = "bold")   +
#  geom_text(aes(x = 6, y = max(rmat_df$Ncum - 12), angle = 90), label = "secondary uplift", size = 3, fontface = "bold")   +
  # Add climate aridification
  geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
 # geom_text(aes(x = 1.5, y = max(rmat_df$Ncum - 20), angle = 90), label = "Arabia aridification", size = 3, fontface = "bold") +
  
  
  #observed events
    geom_line(data = line_all, aes(x = Ma, y = Ncum), color = '#EE6A50',
            size = lwd_events) + 
  
  xlim(80,0) +
  labs(x = "Time before present (Ma)", y = "Number of lineages") +
  ggtitle("Lineage accumulation through time") +
  # Insert geologic scale
  coord_geo(xlim = c(80, 0), ylim = c(0,max(rmat_df$Ncum + 5)), pos = as.list(rep("bottom", 2)),
            dat = list(epochs_htc, periods_htc),
            height = list(unit(1, "lines"), unit(1, "line")),
            rot = list(0, 0), size = list(2, 3), abbrv = list(TRUE, FALSE), 
            center_end_labels = T, 
            skip = c('Holocene'), 
            lab = TRUE) +
  # Set the theme 
  theme_htc()
ggsave("Mountain_colonization_Simmaps/simmaps_all_tree/plots/plot_lineage_accumulation.pdf", width = 7, height = 5)
}
