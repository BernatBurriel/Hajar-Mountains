##### THE BEGINNING #####
# Simmaps analysis to learn how many times the Hajar Mountains have been colonized in each genera
setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/Mountain_colonization_Simmaps/")
rm(list = ls())

libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape")
lapply(libs, require, character.only = TRUE)

#### Import the tree in beast format for each genus
big_tree <- read.beast("TREES/SQUAMATA_ALL_PF.tree")
big_tree_phylo <- as.phylo(big_tree)
big_tree_tibble <- as_tibble(big_tree)
# import genera names 
species.files <- list.files(path = "GENUS", full.names = T)
names.genera <- gsub(".txt", "", species.files)
genera <- gsub("GENUS/", "", names.genera)

nodelist.consensus <- vector(mode = "numeric", length(names.genera))
names(nodelist.consensus) <- names.genera  

for (i in 1:length(genera)) {
  species <- read.table(species.files[i], header = T)
  nodelist.consensus[i] <- MRCA(big_tree_tibble, species$Species)$node
  #species$Species %in% big_tree@phylo$tip.label # In case the for loop does not work, run this to check all names are te same
      }
consensus_genera <- vector("list", length(genera))
names(consensus_genera) <- genera
for (i in 1:length(nodelist.consensus)){
  consensus_genera[[i]] <- tree_subset(big_tree, nodelist.consensus[i], levels_back=0)
}
mountain_data <- vector("list", length(genera))
names(mountain_data) <- genera
for (i in 1:length(mountain_data)) {
  dt <- read.table(species.files[i], header = T, row.names = 1)
  mountain_data[[i]] <- dt[,1]
  names(mountain_data[[i]]) <- rownames(dt)
}

saveRDS(mountain_data, "objects/mountain_data.rds")
saveRDS(consensus_genera, "objects/consensus_genera.rds")
mountain_data <- readRDS("objects/mountain_data.rds")
consensus_genera <- readRDS("objects/consensus_genera.rds")

##### ANCESTRAL RECONSTRUCTION OF PRESENCE IN THE MOUNTAINS #####

registerDoParallel(cores=7)
?fitDiscrete
?multi2di

{fit_ER_dnamountains <- vector("list", length(genera))
fit_SYM_dnamountains <- vector("list", length(genera))
fit_ARD_dnamountains <- vector("list", length(genera))
names(fit_ER_dnamountains) <- genera
names(fit_SYM_dnamountains) <- genera
names(fit_ARD_dnamountains) <- genera
}
for (i in 1:length(genera)) {
  td <- geiger::treedata(consensus_genera[[i]]@phylo, mountain_data[[i]], sort = TRUE)
  td_data <- td$data[,1]
  # Fit ER
  fit_ER_dnamountains[[i]] <- fitDiscrete(consensus_genera[[i]]@phylo, td_data, model = 'ER', ncores = 7)
  # Fit SYM
  fit_SYM_dnamountains[[i]] <- fitDiscrete(consensus_genera[[i]]@phylo, td_data, model = 'SYM', ncores = 7)
  # Fit ARD
  fit_ARD_dnamountains[[i]] <- fitDiscrete(consensus_genera[[i]]@phylo, td_data, model = 'ARD', ncores = 7)
}

saveRDS(fit_ER_dnamountains, "objects/fitted_models/fit_ER_dnamountains.rds")
saveRDS(fit_SYM_dnamountains, "objects/fitted_models/fit_SYM_dnamountains.rds")
saveRDS(fit_ARD_dnamountains, "objects/fitted_models/fit_ARD_dnamountains.rds")
fit_ER_dnamountains <- readRDS("objects/fitted_models/fit_ER_dnamountains.rds")
fit_SYM_dnamountains <- readRDS("objects/fitted_models/fit_SYM_dnamountains.rds")
fit_ARD_dnamountains <- readRDS("objects/fitted_models/fit_ARD_dnamountains.rds")

fitted_models <- data.frame(ER = NA, SYM = NA, ARD = NA)
for (i in 1:length(genera)) {
  fitted_models[i,"ER"] <- fit_ER_dnamountains[[i]]$opt$aicc
  fitted_models[i,"SYM"] <- fit_SYM_dnamountains[[i]]$opt$aicc
  fitted_models[i,"ARD"] <- fit_ARD_dnamountains[[i]]$opt$aicc
}

rownames(fitted_models) <- genera
best_models <- names(fitted_models)[apply(fitted_models, MARGIN = 1, FUN = which.min)]

#### Make simmap
?make.simmap
mountains_simmap1000 <- vector("list", length(genera))
names(mountains_simmap1000) <- genera
mountains_simmap <- vector("list", length(genera))
names(mountains_simmap) <- genera
for (i in 1:length(genera)) {
  mountains_simmap1000[[i]] <- make.simmap(consensus_genera[[i]]@phylo, mountain_data[[i]], 
                                           model = best_models[i], nsim = 1000)
  mountains_simmap[[i]] <- describe.simmap(mountains_simmap1000[[i]])
}  

saveRDS(mountains_simmap1000, "objects/mountains_simmap1000.rds")
saveRDS(mountains_simmap, "objects/mountains_simmap.rds")
mountains_simmap1000 <- readRDS("objects/mountains_simmap1000.rds")
mountains_simmap <- readRDS("objects/mountains_simmap.rds")

# Explore the data
colnames(mountains_simmap$Asaccus$count)
summary(mountains_simmap$Asaccus$count)
mountains_simmap$Asaccus$ace

# Set the colours
levels(as.factor(mountain_data[[1]]))
mountain_colors <- c("#a31429", "#c2c2cc")
names(mountain_colors) <- levels(as.factor(mountain_data[1]))
plot(1:2, cex = 5 , pch = 16, col = mountain_colors)

#### PLOT SIMMAP ####

for (i in 1:length(genera)) {
  outfile <- paste("plots/", genera[i], "_reconstruction.pdf",  sep ="")
  v_breaks <- seq(from = 0, to = max(node.depth.edgelength(consensus_genera[[i]]@phylo)[1]), by = 5)
  pdf(outfile, height=8.3, width=5.8)
  plotTree(consensus_genera[[i]]@phylo, fsize = .6, lwd = 1, type = "phylogram", 
           color = "darkgrey", ftype = "i", offset = .6, mar = c(3, 1, 0.5, 1))
  tiplabels(pie = mountains_simmap[[i]]$tips, piecol = mountain_colors, cex = 0.56, lty = par(lty="blank"))
  nodelabels(pie = mountains_simmap[[i]]$ace, piecol = mountain_colors, cex = 0.5, lty = par(lty="blank"))
#  add.simmap.legend(colors=mountain_colors, fsize=0.5, prompt = F, x=1, y=-2)
  abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
  axisPhylo(side = 1, lwd = .8, cex.axis = .8)
  dev.off()
}


######## SIMMAPS RECONSTRUCTION FOR THE WHOLE TREE ########
#### Import the tree in beast format for each genus
rm(list = ls())
{setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/Mountain_colonization_Simmaps/simmaps_all_tree/")
big_tree <- read.beast("../TREES/SQUAMATA_ALL_PF.tree")
big_tree_phylo <- as.phylo(big_tree)
big_tree_phylo$tip.label <- gsub("_", " ", big_tree_phylo$tip.label)
big_tree@phylo$tip.label <- gsub("_", " ", big_tree@phylo$tip.label)
big_tree_tibble <- as_tibble(big_tree)
# import species data
sps <- read.csv2("data/species_tree.csv")
sps$Species_names_pasted <- gsub("_"," ", sps$Species_names_pasted)
rownames(sps) <- sps$Species_names_pasted
name.check(big_tree_phylo, sps)
sps$Hajars <- "Out"
for (i in 1:NROW(sps)) {
if (sps$H[i] == 1) {
  sps$Hajars[i] = "Hajars"
  } else {sps$Hajars[i] = "Out"}  
  
}
}
saveRDS(big_tree, "objects/big_tree.rds")
saveRDS(sps, "objects/sps.rds")

##### ANCESTRAL RECONSTRUCTION #####
big_tree <- readRDS("objects/big_tree.rds")
sps <- readRDS("objects/sps.rds")
mountain_data <- sps$Hajars
names(mountain_data) <- rownames(sps)
registerDoParallel(cores=8)

head(sps)
nrow(sps) # number of species = 284

# Fit models of discrete character evolution to the general tree
td <- geiger::treedata(big_tree@phylo, mountain_data, sort = TRUE)
td_data <- td$data[,1]
#fit ER
fit_ER_mountains <- fitDiscrete(big_tree@phylo, td_data, model = 'ER', ncores = 8)
saveRDS(fit_ER_mountains, 'objects/fit_ER_alltree.rds')
# fit SYM
fit_SYM_mountains <- fitDiscrete(big_tree@phylo, td_data, model = 'SYM', ncores = 24)
saveRDS(fit_SYM_mountains, 'objects/fit_SYM_alltree.rds')
# AICc ARD
fit_ARD_mountains <- fitDiscrete(big_tree@phylo, td_data, model = 'ARD', ncores = 24)
saveRDS(fit_ARD_mountains, 'objects/fit_ARD_alltree.rds')

fit_ER_mountains <- readRDS('objects/fit_ER_alltree.rds')
fit_SYM_mountains <- readRDS('objects/fit_SYM_alltree.rds')
fit_ARD_mountains <- readRDS('objects/fit_ARD_alltree.rds')
aicw(c(fit_ER_mountains$opt$aicc, fit_SYM_mountains$opt$aicc, 
       fit_ARD_mountains$opt$aicc))

# best fit model for mountain occupancy is ARD (the lowest one)
# Make simmap
mountains_simmap1000 <- make.simmap(big_tree@phylo, mountain_data, model = "ARD", nsim = 1000)
mountains_simmap <- describe.simmap(mountains_simmap1000)
saveRDS(mountains_simmap1000, "objects/mountains_simmap1000_alltree.rds")
saveRDS(mountains_simmap, "objects/mountains_simmap_alltree.rds")

mountains_simmap <- readRDS("objects/mountains_simmap_alltree.rds")
colnames(mountains_simmap$count)

# Set the colours
levels(as.factor(mountain_data))
mountain_colors <- c("#a31429", "#c2c2cc")
names(mountain_colors) <- levels(as.factor(mountain_data))
plot(1:2, cex = 5 , pch = 16, col = mountain_colors)

# PLOT SIMMAP ----
# create a checklist with all the species names for the piechart to be plotted
species.files <- list.files(path = "../GENUS", full.names = T)
names.genera <- gsub(".txt", "", species.files)
genera <- gsub("../GENUS/", "", names.genera)
species_checklist <- data.frame(Species = NA)
for (i in 1:length(genera)) {
  species <- read.table(species.files[i], header = T)
  species$Species <- gsub("_", " ", species$Species)
  species <- species%>% dplyr::select(Species)
  species_checklist <- rbind(species_checklist, species)
}
species_checklist <- species_checklist[-1,]
species_checklist <- species_checklist[species_checklist != "Podarcis lilfordi"]
species_checklist <- species_checklist[species_checklist != "Podarcis pityusensis"]

nodes_simmap <- as.data.frame(mountains_simmap$ace)
nodes_cal <- read.table('data/calibration_nodes.txt', header = T)
nodes_cal$num <- 1:NROW(nodes_cal)
nodes_cal$node <- as.numeric(nodes_cal$node)
node_mountains <- as.numeric(rownames(nodes_simmap[nodes_simmap$Out<.80,]))
nodes_simmap_matrix <- as.matrix(nodes_simmap)
tips_simmap <- as.data.frame(mountains_simmap$tips)
tips_simmap_matrix <- as.matrix(tips_simmap)
tip_mountains <- as.numeric(which(tips_simmap$Out!=1))
tips_mountains_genus <- list()
for (i in 1:NROW(species_checklist)) {
  tips_mountains_genus <- (rbind(tips_mountains_genus, as.numeric(which(rownames(tips_simmap) == species_checklist[i]))))
}  
tips_mountains_genus <- as.numeric(tips_mountains_genus)

big_tree@phylo$tip.label
node_mountain_genus <- as.numeric(c(302:309,338:342,392:399,425:456, 460:464, 470:521, 536:566))

#tsta <-  ape::rotate(phy=big_tree@phylo, node = 286)


# Plot all big
pdf("plots/reconstruction_big.pdf", height=10, width=15)
plotTree(big_tree@phylo, fsize=0.1, lwd=0.1, type="fan", color="gray")
tiplabels(pie=mountains_simmap$tips, piecol=mountain_colors, cex=0.1, lty=par(lty="blank"))
nodelabels(pie=mountains_simmap$ace, piecol=mountain_colors, cex=0.1, lty=par(lty="blank") )
axisPhylo()
dev.off()

#plot node numbers
pdf("plots/reconstruction_big_nodes.pdf", height=10, width=15)
plotTree(big_tree@phylo, fsize=0.4, lwd=1, type="fan", color="gray", offset = 5)
tiplabels(pie=mountains_simmap$tips, piecol=mountain_colors, cex=0.15, lty=par(lty="blank"))
nodelabels(pie=mountains_simmap$ace, piecol=mountain_colors, cex=0.12, lty=par(lty="blank") )
nodelabels(text=rownames(nodes_simmap), cex=0.4, bg = "lightgreen", frame = "circle")
#axisPhylo()
dev.off()

# Plot only mountain big (small piecharts)
pdf("plots/mountain_reconstruction_big_ARD.pdf", height=10, width=10)
plotTree(big_tree@phylo, fsize=0.3, ftype='i', lwd=1, type="fan", color="gray", offset=5)
tiplabels(tip=tips_mountains_genus, pie=tips_simmap_matrix[tips_mountains_genus,], piecol=mountain_colors, cex=0.15, lty=par(lty="blank"))
nodelabels(node=node_mountain_genus, pie=nodes_simmap_matrix[as.character(node_mountain_genus),], piecol=mountain_colors, cex=0.1, lty=par(lty="blank"))
#add.simmap.legend(colors=mountain_colors, fsize=1.8, prompt = F, x=30, y=104)
axisPhylo()
dev.off()

# Plot only mountain big (small piecharts) final
pdf("plots/mountain_reconstruction_big_ARD_noaxis.pdf", height=10, width=10)
plotTree(ladderize(big_tree@phylo, right = F), fsize=0.3, ftype='i', lwd=1, type="fan", color="gray", offset=5, part = 0.98)
tiplabels(tip=tips_mountains_genus, pie=tips_simmap_matrix[tips_mountains_genus,], piecol=mountain_colors, cex=0.15, lty=par(lty="blank"))
nodelabels(node=node_mountain_genus, pie=nodes_simmap_matrix[as.character(node_mountain_genus),], piecol=mountain_colors, cex=0.12, lty=par(lty="blank"))
nodelabels(node=nodes_cal$node, text=nodes_cal$num, bg="violet", cex=0.14, frame = "circle")
add.simmap.legend(colors=mountain_colors, fsize=1.2, prompt = F, x=60, y=-70, vertical = T)
#axisPhylo()
line<-axis(1,pos=-8,at=c(0:233),cex.axis=1.5,labels=FALSE, col = "grey", lwd.ticks = 0)
obj <- axis(1, pos = -8, at=c(seq(0,200,by=50),233), lwd = 0, lwd.ticks = 1.5, labels = F, col = "darkgrey")
text(obj,rep(-19,length(obj)),rev(obj),cex=0.4)
text(mean(obj),-35,"time (mya)",cex=0.5)
dev.off()

# Plot only mountain big (small piecharts) final
{pdf("plots/mountain_reconstruction_big_ARD_noaxis_notips.pdf", height=12, width=12)
par(mar=c(20,20,20,20))
plotTree(ladderize(big_tree@phylo, right = F), fsize=0.00001, ftype='i', lwd=1, type="fan", color="gray", offset=5, part = 0.95, mar = c(6,6,6,6))
nodes_groups <- as.numeric(c(290,299,327,364,374,416))
groups_names <- c('SCINCOIDEA', 'LACERTOIDEA', 'IGUANIA','ANGUIMORPHA', 'SERPENTES','GEKKOTA')
nodes_genera <- as.numeric(c(470, 442, 426,309,464,337,536,392))
genera_names <- c('Hemidactylus arid clade', 'Ptyodactylus', 'Asaccus', 'Omanosaura', 'Trachydactylus', 'Pseudotrapelus', 'Pristurus', 'Echis')


nulo<-mapply(arc.cladelabels,text=genera_names,node=nodes_genera,col='black',
             mark.node=F,lwd=3.5,cex=1.1,lab.offset=1.1, ftype = 'i',
             ln.offset=1.05)

#nulo<-mapply(arc.cladelabels,text=groups_names,node=nodes_groups,col='black',
 #            mark.node=F,lwd=2.5,cex=0.8,lab.offset=1.28,
#             ln.offset=1.14)



tiplabels(tip=tips_mountains_genus, pie=tips_simmap_matrix[tips_mountains_genus,], piecol=mountain_colors, cex=0.3, lty=par(lty="blank"))
nodelabels(node=node_mountain_genus, pie=nodes_simmap_matrix[as.character(node_mountain_genus),], piecol=mountain_colors, cex=0.3, lty=par(lty="blank"))
nodelabels(node=nodes_cal$node, text=c("","","","","","","","","","","","",""),bg="violet", cex=0.3, frame = "circle")
nodelabels(node=nodes_cal$node, text=nodes_cal$num, bg="violet", cex=0.25, frame = "circle")
add.simmap.legend(colors=mountain_colors, fsize=1.2, prompt = F, x=60, y=-70, vertical = T)
#axisPhylo()
line<-axis(1,pos=-8,at=c(0:233),cex.axis=1.5,labels=FALSE, col = "grey", lwd.ticks = 0)
obj <- axis(1, pos = -8, at=c(seq(0,200,by=50),233), lwd = 0, lwd.ticks = 1.5, labels = F, col = "darkgrey")
text(obj,rep(-19,length(obj)),rev(obj),cex=0.4)
text(mean(obj),-35,"time (mya)",cex=0.5)

dev.off()
}










for(i in 1:length(nodes_groups)){
  arc.cladelabels(tree=big_tree@phylo,groups_names[i],nodes_groups[i],
                  orientation="curved",mark.node=F, offset = 3, col = 'black', MoreArgs=list(mark.node=FALSE,lwd=8,lend=2,cex=2.2,lab.offset=1.12))}
for(i in 1:length(nodes_genera)){
  arc.cladelabels(tree=big_tree@phylo,genera_names[i],nodes_genera[i],
                  orientation="curved",mark.node=FALSE)}




col=cols[5:6],
MoreArgs=list(mark.node=FALSE,lwd=8,lend=2,cex=2.2,lab.offset=1.12
              








