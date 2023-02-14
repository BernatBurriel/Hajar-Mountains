##### THE BEGINNING #####
# Simmaps analysis to learn how many times the Hajar Mountains have been colonized in each genera
setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/continuous_trait/")
rm(list = ls())

libs <- c("treeio", "phytools", "geiger", "tidytree", "dplyr", "tidyverse", "doParallel",
          "RColorBrewer", "ggtree", "scico", "ape")
lapply(libs, require, character.only = TRUE)

#### Import the tree in beast format for each genus
tree_names <- list.files(path = "../trees/", pattern = '\\.tree$', full.names = F)
tree_files <- list.files(path = "../trees/", pattern = '\\.tree$',full.names = T)
tree_names <- gsub('_all', '', tree_names)
tree_names <- gsub('.tree','', tree_names)
genera <- gsub('[[:digit:]]+', '', tree_names)

consensus_genera <-  vector("list", length(genera))
names(consensus_genera) <- genera
for (i in 1:length(genera)){
  consensus_genera[[i]] <- read.beast(tree_files[i])
}


# import associated table names 
species.files <- list.files(path = "data/elevation/", full.names = T)
elev_data <- vector("list", length(genera))
for (i in 1:length(elev_data)) {
  elev_data[[i]] <- read.table(species.files[i], header = T)
  elev_data[[i]] <- elev_data[[i]] %>% dplyr::select(Name, elev)
  rownames(elev_data[[i]]) <- elev_data[[i]]$Name
  }


check_names <- vector('list', length(genera))
for (i in 1:length(genera)) {
 check_names[[i]] <- name.check(consensus_genera[[i]]@phylo,elev_data[[i]])  
}

elev_data[[10]] <- elev_data[[10]][elev_data[[10]]$Name != "CN10658",]
elev_data[[10]] <- elev_data[[10]][elev_data[[10]]$Name != "CN10695",]
elev_data[[10]] <- elev_data[[10]][elev_data[[10]]$Name != "TMHC_2013_10_404",]
elev_data[[10]] <- elev_data[[10]][elev_data[[10]]$Name != "TMHC_2013_10_405",]


saveRDS(elev_data, "objects/elev_data.rds")
saveRDS(consensus_genera, "objects/consensus_genera.rds")
elev_data <- readRDS("objects/elev_data.rds")
names(elev_data) <- genera

species<- elev_data$Prup[,2]
tree <- consensus_genera$Prup
names(species) <- rownames(elev_data$Prup)

fit <- fastAnc(consensus_genera$Prup@phylo, species, vars = T, CI = T)
obj<-contMap(consensus_genera$Prup@phylo, species,plot=FALSE)

plot(obj,legend=0.7*max(nodeHeights(consensus_genera$Prup@phylo)),
     sig=2,fsize=c(0.7,0.9))

ggtree(obj$tree, aes(color=maps)) +
  scale_color_continuous(low='darkgreen', high='red') +
  theme(legend.position="right")


dev.off()


## load `tree_anole` and `df_svl` from 'TDbook'
tree <- consensus_genera$Prup@phylo
fit <- fastAnc(tree, species, vars = T, CI = T)

elev <- as.matrix(as.data.frame(species))[,1]

td <- data.frame(node = nodeid(tree, names(elev)),
                 trait = elev)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree, d, by = 'node')

ggtree(tree, aes(color=trait), layout = 'circular', 
             ladderize = FALSE, continuous = 'colour', size=2) +
  scale_color_gradient(colours=c('#3dd3bd', '#e1785a', '#694e5d')) +
  geom_tiplab(hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(.05, .85))

ggtree(tree, aes(color=trait), layout = 'rectangular', ladderize = T, conitnuous = 'colour', size = 1) +
  scale_color_gradientn(colours=c('#3dd3bd', '#e1785a', '#694e5d')) 
   scale_colour_gradient2(low='#e1785a', mid = '#3dd3bd', high ='#694e5d') +
  geom_tiplab(hjust = -.1) + 
  theme(legend.position = c(.05, .85))




## load `tree_anole` and `df_svl` from 'TDbook'
svl <- as.matrix(df_svl)[,1]
fit <- phytools::fastAnc(tree_anole, svl, vars=TRUE, CI=TRUE)

td <- data.frame(node = nodeid(tree_anole, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree_anole, d, by = 'node')

ggtree(tree, aes(color=trait), layout = 'rectangular', 
             ladderize = FALSE, continuous = 'colour', size=2) +
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(.05, .85)) 








phenogram(consensus_genera$Prup@phylo,species,fsize=0.6,spread.costs=c(1,0))


## simulate a tree & some data
tree<-pbtree(n=26,scale=1,tip.label=LETTERS[26:1])
## simulate with ancestral states
x<-fastBM(tree,internal=TRUE)
## ancestral states
a<-x[as.character(1:tree$Nnode+Ntip(tree))]
## tip data
x<-x[tree$tip.label]
fit<-fastAnc(tree,x,CI=TRUE)
fit
plot(a,fit$ace,xlab="true states",ylab="estimated states")
lines(range(c(x,a)),range(c(x,a)),lty="dashed",col="red") ## 1:1 line
## custom function that conducts a simulation, estimates ancestral
## states, & returns the fraction on 95% CI
foo<-function(){
  tree<-pbtree(n=100)
  x<-fastBM(tree,internal=TRUE)
  fit<-fastAnc(tree,x[1:length(tree$tip.label)],CI=TRUE)
  mean(((x[1:tree$Nnode+length(tree$tip.label)]>=fit$CI95[,1]) +
          (x[1:tree$Nnode+length(tree$tip.label)]<=fit$CI95[,2]))==2)
}
## conduct 100 simulations
pp<-replicate(100,foo())
mean(pp)

dev.off()

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
  td <- geiger::treedata(consensus_genera[[i]]@phylo, elev_data[[i]], sort = TRUE)
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
best_models <- names(fitted_models)[apply(fitted_models, MARGIN = 1, FUN = which.min)]

#### Make simmap
?make.simmap
elev_simmap1000 <- vector("list", length(genera))
names(elev_simmap1000) <- genera
elev_simmap <- vector("list", length(genera))
names(elev_simmap) <- genera


for (i in 1:length(genera)) {
  elev_data_def[[i]] <- elev_data[[i]]$elev 
  names(elev_data_def[[i]]) <- rownames(elev_data[[i]])
}

for (i in 1:length(genera)) {
  elev_simmap1000[[i]] <- make.simmap(consensus_genera[[i]]@phylo, elev_data_def[[i]], 
                                           model = best_models[i], nsim = 1000)
  elev_simmap[[i]] <- describe.simmap(elev_simmap1000[[i]])
}  

saveRDS(elev_simmap1000, "objects/elev_simmap1000.rds")
saveRDS(elev_simmap, "objects/elev_simmap.rds")
elev_simmap1000 <- readRDS("objects/elev_simmap1000.rds")
elev_simmap <- readRDS("objects/elev_simmap.rds")

# Explore the data
colnames(elev_simmap$Asaccus$count)
summary(elev_simmap$Asaccus$count)
elev_simmap$Asaccus$ace

# Set the colours
elev_data <- readRDS("objects/elev_data.rds")
elev_data_def <- elev_data
for (i in 1:length(elev_data)) {
  elev_data_def[[i]] <- elev_data[[i]]$elev 
  names(elev_data_def[[i]]) <- rownames(elev_data[[i]])
}
consensus_genera <- readRDS("objects/consensus_genera.rds")

levels(as.factor(elev_data_def[[1]]))
elevation_colors <- c( "#694F5D", "#41D3BD", "#E07A5F")
names(elevation_colors) <- levels(as.factor(elev_data_def[[1]]))
plot(1:3, cex = 5 , pch = 16, col = elevation_colors)

#### PLOT SIMMAP ####
plotBranchbyTrait(consensus_genera[[i]]@phylo, x  = consensus_genera[[i]]@data$location, mode = 'nodes')

for (i in 1:length(genera)) {
    outfile <- paste("plots/", genera[i], "_reconstruction.pdf",  sep ="")
    v_breaks <- seq(from = 0, to = max(node.depth.edgelength(consensus_genera[[i]]@phylo)[1]), by = 5)
    colors_plot <- elevation_colors[names(elevation_colors) %in% unique(colnames(elev_simmap[[i]]$tips))]
    pdf(outfile, height=8.3, width=5.8)
    plotTree(consensus_genera[[i]]@phylo, fsize = .3, lwd = 1, type = "phylogram",
             color = "darkgrey", ftype = "i", offset = .6, mar = c(3, 1, 0.5, 1))
    tiplabels(pie = elev_simmap[[i]]$tips, piecol = colors_plot, cex = 0.3, lty = par(lty="blank"))
    nodelabels(pie = elev_simmap[[i]]$ace, piecol = colors_plot, cex = 0.3, lty = par(lty="blank"))
    
  #  add.simmap.legend(colors=elevation_colors, fsize=0.5, prompt = F, x=1, y=-2)
    abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
    axisPhylo(side = 1, lwd = .8, cex.axis = .8)
    dev.off()
  }

#### PLOT SIMMAP with BSSVS ####

## Read simmap-files 
elev_data <- readRDS("objects/elev_data.rds")
elev_data_def <- elev_data
for (i in 1:length(elev_data)) {
  elev_data_def[[i]] <- elev_data[[i]]$elev 
  names(elev_data_def[[i]]) <- rownames(elev_data[[i]])
}
consensus_genera <- readRDS("objects/consensus_genera.rds")
elev_simmap <- readRDS("objects/elev_simmap.rds")

tree_files_bssvs <- list.files(path = "/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/BSSVS/consensus_trees", pattern = "\\.tree$", full.names = TRUE)
tree_posterior_files <- list.files(path = "/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/BSSVS/tree_files/trees_renamed", pattern = "\\.tree$", full.names = TRUE) # This are the trees with the associated posterior probability


consensus_genera_bssvs <- vector("list", length(genera))
posterior_trees <- vector("list", length(genera))
names(consensus_genera_bssvs) <- names(posterior_trees) <- genera

for (i in 1:length(consensus_genera)) {
  consensus_genera_bssvs[[i]] <- read.beast(file = tree_files_bssvs[i])
  posterior_trees[[i]] <- read.beast(file = tree_posterior_files[i])
}
write_rds(consensus_genera_bssvs, "objects/consensus_genera_bssvs.rds")
write_rds(posterior_trees, "objects/posterior_trees.rds")
consensus_genera_bssvs <- readRDS("objects/consensus_genera_bssvs.rds")
posterior_trees <- readRDS("objects/posterior_trees.rds")

#Get a table with the posterior probabilities and pass it to the new tree
bssvs_location <- vector("list", length(genera))
names(bssvs_location) <- genera


for (i in 1:length(bssvs_location)) {
  bssvs_location[[i]] <- as.data.frame(data.frame(pp = as_tibble(consensus_genera_bssvs[[i]])$location, node = as_tibble(consensus_genera_bssvs[[i]])$node))

  tst_cons_gen <- match(consensus_genera_bssvs[[i]]@data$node,consensus_genera[[i]]@data$node)
  
  consensus_genera[[i]]@data$location <- NA
  consensus_genera[[i]]@data$location.prob <- NA
  consensus_genera[[i]]@data$location.set <- NA
  consensus_genera[[i]]@data$location.set.prob <- NA
  consensus_genera[[i]]@data$location[tst_cons_gen[!is.na(tst_cons_gen)]] <- consensus_genera_bssvs[[i]]@data$location  
  consensus_genera[[i]]@data$location[tst_cons_gen[is.na(tst_cons_gen)]] <- "I"
  consensus_genera[[i]]@data$location[is.na(consensus_genera[[i]]@data$location)] <- "I"
  consensus_genera[[i]]@data$location.prob[tst_cons_gen[!is.na(tst_cons_gen)]] <- consensus_genera_bssvs[[i]]@data$location.prob
  consensus_genera[[i]]@data$location.prob[tst_cons_gen[is.na(tst_cons_gen)]] <- 1
  consensus_genera[[i]]@data$location.prob[is.na(consensus_genera[[i]]@data$location.prob)] <- 1
  consensus_genera[[i]]@data$location.set[tst_cons_gen[!is.na(tst_cons_gen)]] <- consensus_genera_bssvs[[i]]@data$location.set
  consensus_genera[[i]]@data$location.set[tst_cons_gen[is.na(tst_cons_gen)]] <- "I"
  consensus_genera[[i]]@data$location.set[is.na(consensus_genera[[i]]@data$location.set)] <- "I"
  consensus_genera[[i]]@data$location.set.prob[tst_cons_gen[!is.na(tst_cons_gen)]] <- consensus_genera_bssvs[[i]]@data$location.set.prob
  consensus_genera[[i]]@data$location.set.prob[tst_cons_gen[is.na(tst_cons_gen)]] <- '1'
  consensus_genera[[i]]@data$location.set.prob[is.na(consensus_genera[[i]]@data$location.set.prob)] <- '1'
  }  

bssvs_probs <- vector("list", length(genera))
names(bssvs_probs) <- genera

for (i in 1:length(bssvs_probs)) {
  bssvs_probs[[i]] <- data.frame(node = as.numeric(as_tibble(consensus_genera[[i]])$node), number_of_states=0, C=0,E=0,W=0,I=0,M=0)
  for (j in 1:NROW(bssvs_probs[[i]])) {
    loc.set <- as_tibble(consensus_genera[[i]])$location.set[j][[1]]
    bssvs_probs[[i]][j, "number_of_states"] <-length(loc.set)
    if (length(loc.set) > 1) {
      df <- data.frame(node = as.numeric(as_tibble(consensus_genera[[i]])$node[j]), number_of_states=length(loc.set), C=0,E=0,W=0,I=0,M=0)
      for (k in 1:length(loc.set)) {
        df[,names(df) %in% loc.set[k]] <- as.numeric(as_tibble(consensus_genera[[i]])$location.set.prob[j][[1]][k])
      }
      bssvs_probs[[i]][j,] <- df
    }
    if (length(loc.set) == 1){
      bssvs_probs[[i]][j,names(bssvs_probs[[i]]) %in% loc.set[1]] <- as.numeric(as_tibble(consensus_genera[[i]])$location.set.prob[j][[1]])
    }
  }}

bssvs_probs[[6]]$I = 1-bssvs_probs[[6]]$C
# Get only the percents below 90% of assignation to a group
bssvs_probs_reduced <- vector("list", length(genera))
names(bssvs_probs_reduced) <- genera
for (i in 1:length(bssvs_probs_reduced)) {
  bssvs_probs_reduced[[i]] <- bssvs_probs[[i]][bssvs_probs[[i]]$number_of_states > 1, ]
  bssvs_probs_reduced[[i]] <- bssvs_probs_reduced[[i]][between(bssvs_probs_reduced[[i]]$C, 0, 0.9),]
  bssvs_probs_reduced[[i]] <- bssvs_probs_reduced[[i]][between(bssvs_probs_reduced[[i]]$E, 0, .9),]
  bssvs_probs_reduced[[i]] <- bssvs_probs_reduced[[i]][between(bssvs_probs_reduced[[i]]$W, 0, .9),]
  bssvs_probs_reduced[[i]] <- bssvs_probs_reduced[[i]][between(bssvs_probs_reduced[[i]]$M, 0, .9),]
  bssvs_probs_reduced[[i]] <- bssvs_probs_reduced[[i]][between(bssvs_probs_reduced[[i]]$I, 0, .9),]
}


# Get the probs for tiplabels and nodelabels separately
bssvs_probs_tips <- bssvs_probs_nodes <- vector("list", length(genera))
names(bssvs_probs_nodes) <- names(bssvs_probs_tips) <- genera


as_tibble(consensus_genera[[i]])$label


# Set the colors of each biogeographic state ----
levels(as.factor(elev_data_def[[1]]))
elevation_colors <- c( "#694F5D", "#41D3BD", "#E07A5F")
names(elevation_colors) <- levels(as.factor(elev_data_def[[1]]))

colors_plot <- elevation_colors[names(elevation_colors) %in% unique(colnames(elev_simmap[[i]]$tips))]




# Set the colors of each biogeographic state ----
state_colors <- c(C="#b21020",
                  E="#fec14d",
                  I="#bbbbbb",
                  M="#bbbbbb",
                  W="#1868b5")
names_colors <- c("C", "E", "W","I", "M")
levels(state_colors) <- levels(as.factor(names_colors))
names(state_colors) <- levels(as.factor(names_colors)) 
colors_plot <- state_colors[names(state_colors) %in% unique(consensus_genera[[i]]@data$location)]
probs_plot <- bssvs_probs_reduced[[i]][,c('node','number_of_states',names(colors_plot))]

cols<-setNames(c("#b21020","#fec14d","bbbbbb",
                 "bbbbbb","1868b5"),c('C','E','I','M','W'))
for(i in 1:length(genera)){
  consensus_genera[[i]]@data$location_cols = NA
  for (j in 1:length(state_colors)) {
    consensus_genera[[i]]@data$location_cols[consensus_genera[[i]]@data$location == names(state_colors[j])]  <- state_colors[j]
  }
}

#### PLOT Asaccus  ####
#1 set the variables 
colors_plot1 <- state_colors[names(state_colors) %in% unique(consensus_genera[[1]]@data$location)]
probs_plot <- bssvs_probs_reduced[[i]][,c('node','number_of_states',names(colors_plot1))]



v_breaks <- seq(from = 0, to = max(node.depth.edgelength(consensus_genera[[i]]@phylo)[1]), by = 5)
colors_plot <- elevation_colors[names(elevation_colors) %in% unique(colnames(elev_simmap[[i]]$tips))]
outfile = 'test.pdf'
pdf(outfile, height=8.3, width=5.8)
par(mfrow = c(1,2))
plotTree(consensus_genera[[i]]@phylo, fsize = .3, lwd = 1, type = "phylogram", 
         color = 'darkgrey', ftype = "i", offset = .6, mar = c(3, 1, 0.5, 1))
tiplabels(pie = elev_simmap[[i]]$tips, piecol = colors_plot, cex = 0.3, lty = par(lty="blank"))
nodelabels(pie = elev_simmap[[i]]$ace, piecol = colors_plot, cex = 0.3, lty = par(lty="blank"))

#  add.simmap.legend(colors=elevation_colors, fsize=0.5, prompt = F, x=1, y=-2)
abline(v = v_breaks[-length(v_breaks)] +(max(node.depth.edgelength(consensus_genera[[i]]@phylo) - max(v_breaks))), lty = 2, lwd = 0.7, col = "lightgrey")
axisPhylo(side = 1, lwd = .8, cex.axis = .8)


plotTree(consensus_genera[[i]]@phylo, fsize = 0.0000001, lwd = 1, type = "phylogram", 
          color = 'darkgrey', ftype = "i", offset = .6, direction = 'leftwards', mar = c(3, 1, 0.5, 1))
tiplabels(pie = bssvs_probs[[i]][bssvs_probs[[i]]$node %in% (1:length(consensus_genera[[i]]@phylo$tip.label)), 3:length(probs_plot)], 
          piecol = colors_plot1, cex = 0.3, lty = par(lty="blank"))
nodelabels(pie = bssvs_probs[[i]][bssvs_probs[[i]]$node %in% (length(consensus_genera[[i]]@phylo$tip.label)+1):NROW(as_tibble(consensus_genera[[i]])),
                                  3:length(probs_plot)], piecol = colors_plot1, cex = 0.3, lty = par(lty="blank"))
dev.off()


