
# This is to plot the maximum likelihood biogeographic states on the trees 
# as yielded by BioGeoBEARS analyses.

setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/")

# Load packages ----
libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "MCMCtreeR", "viridis", 'eply')
lapply(libs, require, character.only = TRUE)

# import tables with biogeographic state
location_files <- list.files("location/locations_corrected", full.names=TRUE)
#location_files[5]

# Import list of clade trees ----
consensus_genera <- readRDS("objects/genera_trees.rds")
genera <- names(consensus_genera)

tipranges_list <- vector("list", length(genera))
names(tipranges_list) <- genera

for(i in 1:length(location_files)){
  tipranges_list[[i]] = getranges_from_LagrangePHYLIP(lgdata_fn=location_files[i])
}


# Import list of tip and node states for each clade ----
consensus_probs <- readRDS("objects/noJ/consensus_probs.rds")
consensus_states <- readRDS("objects/noJ/consensus_states.rds")
consensus_states$Omanosaura[[5]] <- "CEW"
saveRDS(consensus_states , "objects/noJ/consensus_states.rds")
consensus_states <- readRDS("objects/noJ/consensus_states.rds")

# Define clade names 
clade_names <- names(consensus_probs)

# Set the colors of each biogeographic state ----
state_colors <- c(C="#b22125",
                  CE="#ef8354",
                  CEW="#CCF5AC",
                  CI="#b2838a",
                  CW="#33F1FF",
                  E="#fec250",
                  EC="#ef8354",
                  ECW = "#CCF5AC",
                  EM="#f9d9a5",
                  EW = "#CCF5AC",
                  I="#bbbbbb",
                  IW="#748fa5",
                  M="#bbbbbb",
                  R="#206ab4" ,
                  W="#206ab4",
                  WC="#33F1FF"
                  )
names_colors <- c("C", "CE", "CEW", "CI", "CW", "E", "EC", "ECW", "EM","EW", "I", "IW", "M", "R", "W", "WC")
levels(state_colors) <- levels(as.factor(names_colors))
names(state_colors) <- levels(as.factor(names_colors))


plot(1:16, 1:16, col=state_colors, cex=5, pch=16)
text(1:16, 1:16, labels = names(state_colors), cex = 1, pos = 4, col = "black")



#pdf("plots/Biogeobears_consensus_plots.pdf", height=8.3, width=5.8)
for (i in 1:length(consensus_probs)) {
  # prepare data for the plot
  probs = consensus_probs[[i]]
  probs <- probs[,2:length(probs)]
  nodes <- c(1:NROW(probs))
  tipnodes <- nodes[1:length(consensus_genera[[i]]@phylo$tip.label)]
  innernodes <- nodes[(length(consensus_genera[[i]]@phylo$tip.label)+1):NROW(consensus_probs[[i]])]
  
  
  # Create a new matrix with the proper assignation to each color
  df_nodes <- as.data.frame(matrix(nrow = NROW(consensus_probs[[i]]), ncol = 16))
  names(df_nodes) <- names_colors
  
  df_nodes$C = if (is.null(consensus_probs[[i]]$C)){NA} else {consensus_probs[[i]]$C}; df_nodes$CE =  if (is.null(consensus_probs[[i]]$CE)){NA} else {consensus_probs[[i]]$CE}; 
  df_nodes$CI =  if (is.null(consensus_probs[[i]]$CI)){NA} else {consensus_probs[[i]]$CI}; df_nodes$CW = if (is.null(consensus_probs[[i]]$CW)){NA} else {consensus_probs[[i]]$CW}; 
  df_nodes$E =  if (is.null(consensus_probs[[i]]$E)){NA} else {consensus_probs[[i]]$E}; df_nodes$EC = if (is.null(consensus_probs[[i]]$EC)){NA} else {consensus_probs[[i]]$EC};
  df_nodes$ECW =if (is.null(consensus_probs[[i]]$CEW)){NA} else {consensus_probs[[i]]$CEW}; df_nodes$EM = if (is.null(consensus_probs[[i]]$EM)){NA} else {consensus_probs[[i]]$EM}; 
  df_nodes$I = if (is.null(consensus_probs[[i]]$I)){NA} else {consensus_probs[[i]]$I}; df_nodes$IW = if (is.null(consensus_probs[[i]]$IW)){NA} else {consensus_probs[[i]]$IW}; 
  df_nodes$M = if (is.null(consensus_probs[[i]]$M)){NA} else {consensus_probs[[i]]$M}; df_nodes$R = if (is.null(consensus_probs[[i]]$R)){NA} else {consensus_probs[[i]]$R}; 
  df_nodes$W = if (is.null(consensus_probs[[i]]$W)){NA} else {consensus_probs[[i]]$W}; df_nodes$WC = if (is.null(consensus_probs[[i]]$WC)){NA} else {consensus_probs[[i]]$WC};
  df_nodes$WM = if (is.null(consensus_probs[[i]]$WM)){NA} else {consensus_probs[[i]]$WM}; df_nodes$EW = if (is.null(consensus_probs[[i]]$EW)){NA} else {consensus_probs[[i]]$EW}
  df_nodes[is.na(df_nodes)] <- 0
  
  # PREPARE COLORS 
  states <- consensus_states[[i]]
  states_df <- data.frame(node = c(1:length(states)), biogeo = states)
  states_df$biogeo <- as.factor(states_df$biogeo)
  colors_df <- state_colors[levels(state_colors) %in% levels(states_df$biogeo)]
  
  # Prepare line breaks for the time scale
  v_breaks <- seq(from = 0, to = max(node.depth.edgelength(consensus_genera[[i]]@phylo)[1]), by = 5)
  # Prepare the outfile name 
  outfile <- paste("plots_noJ/Biogeobears_", genera[i], "_reconstruction.pdf",  sep ="")
  
  pdf(outfile, height=8.3, width=5.8)
  plotTree(consensus_genera[[i]]@phylo, fsize = .6, lwd = 1, type = "phylogram", 
           color = "darkgrey", ftype = "i", offset = .6, mar = c(3, 1, 0.5, 1))
  abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
  tiplabels(node = tipnodes, pie = df_nodes[1:length(tipnodes),], piecol = state_colors, cex = 0.9, lty = par(lty="blank"))
  nodelabels(node = innernodes, pie = df_nodes[(length(tipnodes)+1):max(innernodes),], piecol = state_colors, cex = 0.9, lty = par(lty="blank"))
  #add.simmap.legend(colors=state_colors, fsize=0.5, leg = names(state_colors), prompt = F, x=1, y=-2)
  axisPhylo(side = 1, lwd = .8, cex.axis = .8 )
  dev.off() 
  outfile2 <- paste("plots_noJ/node_numbers/Biogeobears_", genera[i], "_NODES_reconstruction.pdf",  sep ="")
  
  pdf(outfile2, height=8.3, width=5.8)
  plotTree(consensus_genera[[i]]@phylo, fsize = .6, lwd = 1, type = "phylogram", 
           color = "darkgrey", ftype = "i", offset = .6, mar = c(3, 1, 0.5, 1))
  abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
  tiplabels(text = tipnodes, bg = colors_df[states_df[tipnodes, 2]], cex = 0.5, frame = "circle")
  nodelabels(text = innernodes, bg = colors_df[states_df[innernodes, 2]], cex = 0.5, frame = "circle")
  #add.simmap.legend(colors=state_colors, fsize=0.5, leg = names(state_colors), prompt = F, x=1, y=-2)
  axisPhylo(side = 1, lwd = .8, cex.axis = .8 )
  dev.off() 
  
  outfile3 <- paste("plots_noJ/states/Biogeobears_", genera[i], "_STATES_reconstruction.pdf",  sep ="")
  pdf(outfile3, height=8.3, width=5.8)
  plotTree(consensus_genera[[i]]@phylo, fsize = .6, lwd = 1, type = "phylogram", 
           color = "darkgrey", ftype = "i", offset = .6, mar = c(3, 1, 0.5, 1))
  abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
  tiplabels(text = states_df[tipnodes,2],  bg = colors_df[states_df[tipnodes, 2]], cex = 0.5, frame = "circle")
  nodelabels(text = states_df[innernodes,2], bg = colors_df[states_df[innernodes, 2]], cex = 0.5, frame = "circle")
  #add.simmap.legend(colors=state_colors, fsize=0.5, leg = names(state_colors), prompt = F, x=1, y=-2)
  axisPhylo(side = 1, lwd = .8, cex.axis = .8 )
  dev.off() 
}
  #dev.off()

  














