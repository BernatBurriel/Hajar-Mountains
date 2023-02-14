## plot simmap eelvation blocks and violin plots

setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/Mountain_colonization_Simmaps/Geologic_regions/")
rm(list = ls())

libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape", 'ggplot2')
lapply(libs, require, character.only = TRUE)


## import data 

#elevation
consensus_genera <- readRDS('lowland_midland_highland/objects/consensus_genera.rds')
genera <- names(consensus_genera)
elev_data <- readRDS('lowland_midland_highland/objects/elev_data.rds')
elev_simmap <- readRDS('lowland_midland_highland/objects/elev_simmap.rds')

# continuous elevation
elev_cont_files <- list.files(path = 'lowland_midland_highland/data/elevation', pattern = '\\.txt$', full.names = T)
elev_cont <- elev_cont300 <- elev_cont1500 <- elev_cont1501 <- vector('list', length(genera))
names(elev_cont) <- names(elev_cont300) <- names(elev_cont1500) <- names(elev_cont1501) <- genera
for (i in 1:length(genera)) {
  temp <- read.table(elev_cont_files[i], header = T)
  elev_data[[i]]$cont_elev <- temp$elev[match(elev_data[[i]]$Name, temp$Name)]
  elev_data[[i]]$lineage <- temp$lineage[match(elev_data[[i]]$Name, temp$Name)]
  elev_cont300[[i]] <- elev_data[[i]][elev_data[[i]]$elev == 'low',]
  elev_cont1500[[i]] <- elev_data[[i]][elev_data[[i]]$elev == 'mid',]
  elev_cont1501[[i]] <- elev_data[[i]][elev_data[[i]]$elev == 'high',]
}
#{elev_data[[1]] %>% ggplot(aes(x = lineage, y = cont_elev)) + 
  geom_violin(trim = F) + 
  geom_boxplot(width = 0.1) + 
  #geom_dotplot(binaxis = 'y', show.legend = F, stackdir = 'center', binwidth = 25) +
  theme(axis.title.x = element_blank())#}

#mountain block 
block_data <- readRDS('mountain_blocks/objects/block_data.rds')
block_simmap <- readRDS('mountain_blocks/objects/block_simmap.rds')

### prepare data for plotting trees
## set the colors 
elevation_colors <- c( high="#694F5D", low="#41D3BD", mid="#E07A5F")


block_colors <- c(C="#a83545", E="#f9c87a", I="#bbbbbb",
                  M="#bbbbbb", O = "#bbbbbb", W="#437ba8")
plot(1:3, cex = 5 , pch = 16, col = elevation_colors)
plot(1:6, cex = 5 , pch = 16, col = block_colors)

nodes_simmap <- vector('list', length(genera))
for (i in 1:length(genera)) {
  nodes_simmap[[i]] <- as.data.frame(elev_simmap[[i]]$ace)  
}



nodes_1 <- nodes_3 <- nodes_2 <- vector('list', length(genera))
names(nodes_1) <- names(nodes_3) <- names(nodes_2) <- genera

nodes_1[[1]]<- c(127:131, 146:148, 199,224, 225, 190);nodes_2[[1]]<- c(132,190,134, 222, 208, 191, 179,149)
nodes_1[[2]] <- elev_simmap[[2]]$ace;  nodes_2[[1]]<- elev_simmap[[2]]$tips
nodes_1[[3]] <- c(54,55,77); nodes_2[[3]] <- c(87,94,91,85,68,88,104,95:98,89,78,79,84,76,63,64,56,59,105,90,65,57)
nodes_1[[4]] <- c(20,21,27);nodes_2[[4]] <- c(22:26,28,34)
nodes_1[[5]] <- c(33); nodes_2[[5]] <- c(34,35,50,62)
nodes_1[[6]] <- c(25, 47,26,27); nodes_2[[6]] <- c(NA) 
nodes_1[[7]] <- c(174:176,334); nodes_2[[7]]<- c(335,336,343,340,285,177:181,211,226,240,241,242,262,271)
nodes_1[[8]] <- c(23,24,42) ; nodes_2[[8]] <- c(25,34:36,40,41,26,30)
nodes_1[[9]] <- c(73,74);nodes_2[[9]] <- c(75,108,133); #plot(rotateNodes(consensus_genera[[9]]@phylo,c(108:112,79)))
nodes_1[[10]] <- c(44,45,63) ; nodes_2[[10]] <- c(64:78,45:62)

for (i in 1:length(genera)) {
  mask <-as.numeric(consensus_genera[[i]]@data$node) %in% c(nodes_1[[i]], nodes_2[[i]], 1:length(consensus_genera[[i]]@phylo$tip.label))

  nodes_3[[i]] <- as.numeric(consensus_genera[[i]]@data$node)[!mask]
  
}        




for (i in 1:length(genera)) {
  outfile <- paste0('plots/', genera[i], '_blocks.pdf')
  v_breaks <- seq(from = 0, to = max(node.depth.edgelength(consensus_genera[[i]]@phylo)[1]), by = 5)
  colors_plot_elev <- elevation_colors[names(elevation_colors) %in% unique(colnames(elev_simmap[[i]]$tips))]
  pdf(outfile, height=8.3, width=5.8)
  par(mfrow = c(1,2))
  plotTree(consensus_genera[[i]]@phylo, fsize = .00000001, lwd = 1, type = "phylogram",
           color = "darkgrey", ftype = "i", offset = .6, mar = c(3, 1, 1, 0.25))
  tiplabels(pie = elev_simmap[[i]]$tips, piecol = colors_plot_elev, cex = 0.6, lty = par(lty="blank"))
  level_1_nodes <- elev_simmap[[i]]$ace[rownames(elev_simmap[[i]]$ace) %in% nodes_1[[i]],]
#  nodelabels(pie = elev_simmap[[i]]$ace, 
#             piecol = colors_plot_elev, cex = 0.9, lty = par(lty="blank"))
  nodelabels(node = as.numeric(rownames(level_1_nodes)), pie = elev_simmap[[i]]$ace[rownames(elev_simmap[[i]]$ace) %in% nodes_1[[i]],], 
             piecol = colors_plot_elev, cex = 0.9, lty = par(lty="blank"))
  nodelabels(node = as.numeric(rownames(elev_simmap[[i]]$ace[rownames(elev_simmap[[i]]$ace) %in% nodes_2[[i]],])), 
             pie = elev_simmap[[i]]$ace[rownames(elev_simmap[[i]]$ace) %in% nodes_2[[i]],], piecol = colors_plot_elev, 
             cex = 0.7, lty = par(lty="blank"))
  nodelabels(node = as.numeric(rownames(elev_simmap[[i]]$ace[rownames(elev_simmap[[i]]$ace) %in% nodes_3[[i]],])), 
             pie = elev_simmap[[i]]$ace[rownames(elev_simmap[[i]]$ace) %in% nodes_3[[i]],], 
             piecol = colors_plot_elev, cex = 0.7, lty = par(lty="blank"))
  
  #  add.simmap.legend(colors=elevation_colors, fsize=0.5, prompt = F, x=1, y=-2)
  #abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
  axisPhylo(side = 1, lwd = .8, cex.axis = .8)
  
  colors_plot_block <- block_colors[names(block_colors) %in% unique(colnames(block_simmap[[i]]$tips))]
  par(lty=1)
  plotTree(consensus_genera[[i]]@phylo, fsize = .000001, lwd = 1, type = "phylogram",
           color = "darkgrey", direction = 'leftwards',ftype = "i", offset = 0, mar = c(3, 0.25, 1, 1))
  tiplabels(pie = block_simmap[[i]]$tips, piecol = colors_plot_block, cex = 0.6, lty = par(lty="blank"))
#  nodelabels(pie = block_simmap[[i]]$ace, 
#             piecol = colors_plot_block, cex = 0.9, lty = par(lty="blank"))
  nodelabels(node = as.numeric(rownames(level_1_nodes)), pie = block_simmap[[i]]$ace[rownames(block_simmap[[i]]$ace) %in% nodes_1[[i]],], 
             piecol = colors_plot_block, cex = 0.9, lty = par(lty="blank"))
  nodelabels(node = as.numeric(rownames(block_simmap[[i]]$ace[rownames(block_simmap[[i]]$ace) %in% nodes_2[[i]],])), 
             pie = block_simmap[[i]]$ace[rownames(block_simmap[[i]]$ace) %in% nodes_2[[i]],], piecol = colors_plot_block, 
             cex = 0.7, lty = par(lty="blank"))
  nodelabels(node = as.numeric(rownames(block_simmap[[i]]$ace[rownames(block_simmap[[i]]$ace) %in% nodes_3[[i]],])), 
             pie = block_simmap[[i]]$ace[rownames(block_simmap[[i]]$ace) %in% nodes_3[[i]],], 
             piecol = colors_plot_block, cex = 0.7, lty = par(lty="blank"))
  #nodelabels(pie = block_simmap[[i]]$ace, piecol = colors_plot_block, cex = 0.4, lty = par(lty="blank"))
  #mask <- rownames(block_simmap[[i]]$ace) %in% rownames(block_simmap[[i]]$tips)
  #nodelabels(text = rownames(block_simmap[[i]]$ace[!mask]), cex = 0.4, bg = "lightgreen")
  #  add.simmap.legend(colors=elevation_colors, fsize=0.5, prompt = F, x=1, y=-2)
  #abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
  axisPhylo(side = 1, lwd = .8, cex.axis = .8)
  par(lty=1)
  dev.off()
}
  