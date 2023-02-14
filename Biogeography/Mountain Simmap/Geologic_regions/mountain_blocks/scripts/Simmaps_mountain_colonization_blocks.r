##### THE BEGINNING #####
# Simmaps analysis to learn how many times the Hajar Mountains have been colonized in each genera
setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/Mountain_colonization_Simmaps/Geologic_regions/mountain_blocks/")
rm(list = ls())

libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape")
lapply(libs, require, character.only = TRUE)

#### Import the tree in beast format for each genus
tree_files <- list.files(path = "trees/", pattern = '\\.tree$',full.names = T)
tree_names <- list.files(path = "trees/", pattern = '\\.tree$', full.names = F)
tree_names <- gsub('_all', '', tree_names)
tree_names <- gsub('.tree','', tree_names)
genera <- gsub('[[:digit:]]+', '', tree_names)

consensus_genera <-  vector("list", length(genera))
names(consensus_genera) <- genera
for (i in 1:length(genera)){
  consensus_genera[[i]] <- read.beast(tree_files[i])
}


# import associated table names 
species.files <- list.files(path = "data/names/", full.names = T)
block_data <- vector("list", length(genera))
for (i in 1:length(block_data)) {
  block_data[[i]] <- read.table(species.files[i], header = T)
  names_row <- block_data[[i]]$Name
  block_data[[i]] <- block_data[[i]] %>% dplyr::select(Name,block=name)
  rownames(block_data[[i]]) <- names_row
  }

block_data[[10]] <- block_data[[10]][block_data[[10]]$Name != "CN10658",]
block_data[[10]] <- block_data[[10]][block_data[[10]]$Name != "CN10695",]
block_data[[10]] <- block_data[[10]][block_data[[10]]$Name != "TMHC_2013_10_404",]
block_data[[10]] <- block_data[[10]][block_data[[10]]$Name != "TMHC_2013_10_405",]

check_names <- vector('list', length(genera))
for (i in 1:length(genera)) {
 check_names[[i]] <- name.check(consensus_genera[[i]]@phylo,block_data[[i]])  
}



saveRDS(block_data, "objects/block_data.rds")
saveRDS(consensus_genera, "objects/consensus_genera.rds")
block_data <- readRDS("objects/block_data.rds")
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
  td <- geiger::treedata(consensus_genera[[i]]@phylo, block_data[[i]], sort = TRUE)
  td_data <- td$data[,2]
  # Fit ER
  fit_ER_dnamountains[[i]] <- fitDiscrete(consensus_genera[[i]]@phylo, td_data, model = 'ER', ncores = 30)
  # Fit SYM
  fit_SYM_dnamountains[[i]] <- fitDiscrete(consensus_genera[[i]]@phylo, td_data, model = 'SYM', ncores = 30)
  # Fit ARD
  fit_ARD_dnamountains[[i]] <- fitDiscrete(consensus_genera[[i]]@phylo, td_data, model = 'ARD', ncores = 30)
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
best_models <- data.frame(genera = genera, best_model = names(fitted_models)[apply(fitted_models, MARGIN = 1, FUN = which.min)])
write.table(best_models, 'best_models.txt', row.names = F, quote = F, col.names = T)
best_models = names(fitted_models)[apply(fitted_models, MARGIN = 1, FUN = which.min)]
#### Make simmap
?make.simmap
block_simmap1000 <- vector("list", length(genera))
names(block_simmap1000) <- genera
block_simmap <- vector("list", length(genera))
names(block_simmap) <- genera

block_data_def <- block_data
for (i in 1:length(genera)) {
  block_data_def[[i]] <- block_data[[i]]$block 
  names(block_data_def[[i]]) <- rownames(block_data[[i]])
}

for (i in 1:length(genera)) {
  block_simmap1000[[i]] <- make.simmap(consensus_genera[[i]]@phylo, block_data_def[[i]], 
                                           model = best_models[i], nsim = 1000)
  block_simmap[[i]] <- describe.simmap(block_simmap1000[[i]])
}  

saveRDS(block_simmap1000, "objects/block_simmap1000.rds")
saveRDS(block_simmap, "objects/block_simmap.rds")
block_simmap1000 <- readRDS("objects/block_simmap1000.rds")
block_simmap <- readRDS("objects/block_simmap.rds")

# Explore the data
colnames(block_simmap$Asaccus$count)
summary(block_simmap$Asaccus$count)
block_simmap$Asaccus$ace

# Set the colours
block_data <- readRDS("objects/block_data.rds")
block_data_def <- block_data
for (i in 1:length(block_data)){
  block_data_def[[i]] <- as.data.frame(block_data_def[[i]][,2])
  rownames(block_data_def[[i]]) <- rownames(block_data[[i]])
}
write_rds(block_data_def, 'objects/block_data_def.rds')
consensus_genera <- readRDS("objects/consensus_genera.rds")


block_colors <- c(C="#b21020",
                  E="#fec14d",
                  I="#bbbbbb",
                  M="#bbbbbb",
                  O = "#bbbbbb",
                  W="#1868b5")


plot(1:6, cex = 5 , pch = 16, col = block_colors)

#### PLOT SIMMAP ####

for (i in 1:length(genera)) {
    outfile <- paste("plots/", genera[i], "_reconstruction.pdf",  sep ="")
    v_breaks <- seq(from = 0, to = max(node.depth.edgelength(consensus_genera[[i]]@phylo)[1]), by = 5)
    colors_plot <- block_colors[names(block_colors) %in% unique(colnames(block_simmap[[i]]$tips))]
    pdf(outfile, height=8.3, width=5.8)
    plotTree(consensus_genera[[i]]@phylo, fsize = .3, lwd = 1, type = "phylogram",
             color = "darkgrey", ftype = "i", offset = .6, mar = c(3, 1, 0.5, 1))
    tiplabels(pie = block_simmap[[i]]$tips, piecol = colors_plot, cex = 0.3, lty = par(lty="blank"))
    nodelabels(pie = block_simmap[[i]]$ace, piecol = colors_plot, cex = 0.3, lty = par(lty="blank"))
  #  add.simmap.legend(colors=elevation_colors, fsize=0.5, prompt = F, x=1, y=-2)
    abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
    axisPhylo(side = 1, lwd = .8, cex.axis = .8)
    dev.off()
  }

