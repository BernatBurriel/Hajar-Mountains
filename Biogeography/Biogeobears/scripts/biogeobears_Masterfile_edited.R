##### THE BEGINNING #####
# Useful page to understand BioGeoBEARS output: http://phylo.wikidot.com/example-biogeobears-scripts#toc13

setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/")
# Direct the output to a file as well as printing it to the screen (split=T)
#sink(file="screen_output.txt", append = FALSE, split=TRUE)
###### PACKAGES #####
#library(devtools)
#devtools::install_github(repo="YuLab-SMU/tidytree")
#devtools::install_github(repo="YuLab-SMU/treeio")
#devtools::install_github(repo="YuLab-SMU/ggtree")
#devtools::install_github(repo="wilkelab/cowplot")
#devtools::install_github('nmatzke/BioGeoBears')
rm(list = ls())
{library(dplyr)
library(DescTools)
#library(ips)
#library(cowplot) #package 'cowplot' requires R >= 3.5.0
library(ggplot2)
library(ggtree)
library(tidytree)
library(treeio)
library(phytools)
library(geiger)
#library(treeman)
library(optimx)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(devtools)
library(BioGeoBEARS)}

#source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)# slight speedup hopefully

##### Import the tree (beast format) #####
tree_files <- list.files(path = "tree_files", pattern = "\\.tree$", full.names = TRUE)
names_genera0 <- list.files(path = "tree_files", pattern = "\\.tree$", full.names = FALSE)
names_genera <- gsub(".tree", "", names_genera0)
names_genera <- gsub('[[:digit:]]+', '', names_genera)
genera <- names_genera
# gsub("01"|"02"|"03", "", names_genera)  LOOK FOR EXTRACTING NUMERIC VALUES WITH GSUB

consensus_genera <- vector("list", length(names_genera))
names(consensus_genera) <- names_genera

for (i in 1:length(consensus_genera)) {
  consensus_genera[[i]] <- read.beast(file = tree_files[i])
}

# modify tip labels 
name_correspondance_files <- list.files(path = "name_correspondance", pattern = "\\.txt$", full.names = TRUE)
for (i in 1:length(genera)) {
  names_tbl <- read.table(name_correspondance_files[i], header = T)
  tree <- consensus_genera[[i]]
  tree <- treeio::rename_taxa(tree, data = names_tbl, key = Name, value = New_name)
  tips_to_drop <- names_tbl[names_tbl$Out_In == "out",]
    tree_sub <- treeio::drop.tip(tree, as.character(tips_to_drop$New_name))
    consensus_genera[[i]] <- tree_sub
    }
  


saveRDS(consensus_genera, "objects/genera_trees.rds")


consensus_genera <- readRDS("objects/genera_trees.rds")

##### READING LOCATION FILES AND STORING THEM IN A LIST (EACH ELEMENT IS CLASS "tipranges" from BioGeoBEARS) #####
location_files <- list.files("location/locations_corrected", full.names=TRUE)
#location_files[5]

tipranges_list <- vector("list", length(genera))
names(tipranges_list) <- genera

for(i in 1:length(location_files)){
  tipranges_list[[i]] = getranges_from_LagrangePHYLIP(lgdata_fn=location_files[i])
}
#class(tipranges_list$Acanthodactylus)
#tipranges_list
# Now we have a list, tipranges_list, with the locations for each genus.
# BioGeoBEARS will call each element of that list to reconstruct the biogeography of each genus.

# Set the maximum number of areas any species may occupy; this cannot be larger
# than the number of areas you set up, but it can be smaller.
max_range_size = 2

# Since BioGeoBEARS doesn't work with actual R objects, but with file paths, we'll use the location_files list.
# Furthermore, each time in the loop, the correspondent tree will be exported and then the file path will be used.

##### BIOGEOGRAPHIC RECONSTRUCTION (BioGeoBEARS)... #####
## ...on genera from consensus tree.
{DECconsensus_output_biogeo <- vector("list", length(genera))
names(DECconsensus_output_biogeo) <- genera
DIVAconsensus_output_biogeo <- vector("list", length(genera))
names(DIVAconsensus_output_biogeo) <- genera
BAYAREAconsensus_output_biogeo <- vector("list", length(genera))
names(BAYAREAconsensus_output_biogeo) <- genera

DECjconsensus_output_biogeo <- vector("list", length(genera))
names(DECjconsensus_output_biogeo) <- genera
DIVAjconsensus_output_biogeo <- vector("list", length(genera))
names(DIVAjconsensus_output_biogeo) <- genera
BAYAREAjconsensus_output_biogeo <- vector("list", length(genera))
names(BAYAREAjconsensus_output_biogeo) <- genera

BESTMODELconsensus <- vector("list", length(genera))
names(BESTMODELconsensus) <- genera
Qmat_consensus <- vector("list", length(genera))
names(Qmat_consensus) <- genera
COOmat_consensus <- vector("list", length(genera))
names(COOmat_consensus) <- genera
rootstate_consensus <- vector("list", length(genera))
names(rootstate_consensus) <- genera
consensus_genera_ape <- vector("list", length(genera))
names(consensus_genera_ape) <- genera

resDEC_list <- vector("list", length(genera))
names(resDEC_list) <- genera
resDIVA_list <- vector("list", length(genera))
names(resDIVA_list) <- genera
resBAYAREA_list <- vector("list", length(genera))
names(resBAYAREA_list) <- genera

resDECj_list <- vector("list", length(genera))
names(resDECj_list) <- genera
resDIVAj_list <- vector("list", length(genera))
names(resDIVAj_list) <- genera
resBAYAREAj_list <- vector("list", length(genera))
names(resBAYAREAj_list) <- genera


restable_list <- vector("list", length(genera))
names(restable_list) <- genera

table_bestmodels <- data.frame(genera=genera, bestmodel=NA)
rownames(table_bestmodels) <- genera
}

wd <- getwd()
trfn <- np(paste(wd,"/biogeo_tmp/","tree_tmp.nex",sep=""))
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
cores <- 6 # number of cores for BioGeoBEARS analysis

# Modified functions from BioGeoBEARS to allow a tree in NEXUS format
{source("functions/bears_optim_run_nexus.R")
source("functions/check_BioGeoBEARS_run_nexus.R")
source("functions/check_trfn_nexus.R")
source("functions/get_Qmat_COOmat_from_BioGeoBEARS_run_object_nexus.R")
source("functions/get_Qmat_COOmat_from_res_nexus.R")
source("functions/readfiles_BioGeoBEARS_run_nexus.R")
source("functions/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("functions/modsel.R")
}
#################################################
##### BIOGEOBEARS LOOP FOR CONSENSUS GENERA #####
#################################################
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    write.beast(consensus_genera[[i]], "biogeo_tmp/tree_tmp.nex")
    tr <- read.nexus("biogeo_tmp/tree_tmp.nex")
    consensus_genera_ape[[i]] <- tr
    #trfn <- np(paste(wd,"/biogeo_tmp/","tree_tmp.tre",sep=""))
    print(genera[i])
    geogfn <- location_files[i]
    tipranges <- tipranges_list[[i]]
    # DEC MODEL
    print(paste(genera[i], "DEC"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.0000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc. (see Massana et al.)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = cores
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    runslow = TRUE
    resfn = "biogeo_tmp/output_DEC.RData"
    if (runslow)
    {
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res
      
      save(res, file=resfn)
      resDEC = res
    } else {
      # Loads to "res"
      load(resfn)
      resDEC = res
    }
    
    #resDEC$condlikes_of_each_state
    # store the output table
    resDEC_list[[i]] <- resDEC
    DECconsensus_output_biogeo[[i]] <- resDEC$ML_marginal_prob_each_state_at_branch_top_AT_node
    #resDEC$outputs@params_table["d","est"] #rate of dispersal (range expansion)
    #resDEC$outputs@params_table["e", "est"] # rate of "extinction" (range contraction or extirpation)
    
    # Plot
    #plot(tr, cex=.5)
    #nodelabels(frame="circle", cex=0.4)
    #tiplabels(frame="circle", cex=0.4)
    analysis_titletxt = paste(genera[i], "DEC")
    results_object = resDEC
    resDEC$inputs
    # States
    pdf(paste("DEC_plots/", genera[i],"_dec.pdf", sep=""))
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    # DEC+J model
    print(paste(genera[i], "DEC+J"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O???Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDEC$outputs@params_table["d","est"]
    estart = resDEC$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 6
    
    # Add j as a free parameter
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    resfn = "biogeo_tmp/output_DECJ.RData"
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDECj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDECj = res
    }
    
    resDECj_list[[i]] <- resDECj
    DECjconsensus_output_biogeo[[i]] <- resDECj$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    #######################################################
    # Plot ancestral states - DECJ
    #######################################################
    analysis_titletxt = paste(genera[i], "DEC+J")
    
    # Setup
    results_object = resDECj
    #scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
    
    pdf(paste("DECJ_plots/", genera[i],"_decj.pdf", sep=""))
    # States
    res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    
    dev.off()  # Turn off PDF
    
    # DIVA MODEL
    print(paste(genera[i], "DIVA"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc. (see Massana et al.)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = cores
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = "biogeo_tmp/output_DIVA.RData"
    if (runslow)
    {
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res
      
      save(res, file=resfn)
      resDIVALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKE = res
    }
    
    # store the output table
    resDIVA_list[[i]] <- resDIVALIKE
    DIVAconsensus_output_biogeo[[i]] <- resDIVALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    # Plot
    analysis_titletxt = paste(genera[i], "DIVALIKE")
    results_object = resDIVALIKE
    # States
    pdf(paste("DIVA_plots/", genera[i],"_diva.pdf", sep=""), width = 10, height = 15)
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    # DIVALIKE+J MODEL
    print(paste(genera[i], "DIVA+J"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O???Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDIVALIKE$outputs@params_table["d","est"]
    estart = resDIVALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 6
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # Add jump dispersal/founder-event speciation
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
    
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    resfn = "biogeo_tmp/output_DIVAJ.RData"
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDIVALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKEj = res
    }
    
    # store the output table
    resDIVAj_list[[i]] <- resDIVALIKEj
    DIVAjconsensus_output_biogeo[[i]] <- resDIVALIKEj$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    # Plot
    analysis_titletxt = paste(genera[i], "DIVALIKE+J")
    results_object = resDIVALIKEj
    # States
    pdf(paste("DIVAJ_plots/", genera[i],"_divaj.pdf", sep=""), width = 10, height = 15)
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    
    # BAYAREA MODEL
    print(paste(genera[i], "BAYAREA"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc. (see Massana et al.)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = cores
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = "biogeo_tmp/output_BAYAREA.RData"
    if (runslow)
    {
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res
      
      save(res, file=resfn)
      resBAYAREALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKE = res
    }
    # store the output table
    resBAYAREA_list[[i]] <- resBAYAREALIKE
    BAYAREAconsensus_output_biogeo[[i]] <- resBAYAREALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    # Plot
    analysis_titletxt = paste(genera[i], "BAYAREALIKE")
    results_object = resBAYAREALIKE
    pdf(paste("BAYAREA_plots/", genera[i],"_bayarea.pdf", sep=""))
    # States
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    # BAYAREA+J model
    print(paste(genera[i], "BAYAREA+J"))
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O???Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = "GenSA"
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resBAYAREALIKE$outputs@params_table["d","est"]
    estart = resBAYAREALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 6
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
    # machines. I can't replicate this on my Mac machines, but it is almost certainly
    # just some precision under-run issue, when optim/optimx tries some parameter value 
    # just below zero.  The "min" and "max" options on each parameter are supposed to
    # prevent this, but apparently optim/optimx sometimes go slightly beyond 
    # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
    # slightly for each parameter:
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-12
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-12
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 6
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 1e-12
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    
    check_BioGeoBEARS_run_nexus(BioGeoBEARS_run_object)
    
    resfn = "biogeo_tmp/output_BAYAREAJ.RData"
    runslow = TRUE
    if (runslow)
    {
      res = bears_optim_run_nexus(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resBAYAREALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKEj = res
    }
    # store the output table
    resBAYAREAj_list[[i]] <- resBAYAREALIKEj
    BAYAREAjconsensus_output_biogeo[[i]] <- resBAYAREALIKEj$ML_marginal_prob_each_state_at_branch_top_AT_node
    
    
    # Plot
    analysis_titletxt = paste(genera[i], "BAYAREALIKE+J")
    results_object = resBAYAREALIKEj
    # States
    pdf(paste("BAYAREAJ_plots/", genera[i],"_bayareaj.pdf", sep=""), width = 10, height = 15)
    res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    # Pie chart
    plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
    dev.off()
    
    
    ##### ... Model selection ####
    # Once we have run all the models, we can choose the best fit model according to the AICc
    
    LnL_DEC = get_LnL_from_BioGeoBEARS_results_object(resDEC)
    numparams_DEC = 2
    cols_DEC = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_DECJ = get_LnL_from_BioGeoBEARS_results_object(resDECj)
    numparams_DECJ = 3
    cols_DECJ = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_DIVA = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
    numparams_DIVA = 2
    cols_DIVA = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_DIVAJ = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
    numparams_DIVAJ = 3
    cols_DIVAJ = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_BAYAREA = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
    numparams_BAYAREA = 2
    cols_BAYAREA = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    LnL_BAYAREAJ = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
    numparams_BAYAREAJ = 3
    cols_BAYAREAJ = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
    
    
    # Table with LnL
    restable = rbind(cols_DEC, cols_DECJ, cols_DIVA, cols_DIVAJ, cols_BAYAREA, cols_BAYAREAJ)
    row.names(restable) = c("DEC", "DECJ", "DIVALIKE", "DIVALIKEJ", "BAYAREALIKE", "BAYAREALIKEJ")
    restable
    # Add AIC and AIC weight
    AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
    restable = cbind(restable, AICtable)
    #restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
    #restable_AIC_rellike
    
    if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 4){
      # Add AICc and AICc weight
      samplesize = length(tr$tip.label)
      AICctable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
      restable = cbind(restable, AICctable)
      #restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
      #restable_AICc_rellike
      
      restable_list[[i]] <- restable
      
      # Model selection (min AICc): 1 is DEC, 2 is DEC+J, 3 is DIVA, 4 is DIVA+J,
      # 5 is BAYAREA, 6 is BAYAREA+J
      # and store the Q matrix of the best model for posterior simulations
      if (restable$AICc[1] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- DECconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DEC"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDEC)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDEC)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDEC$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[2] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- DECjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DECj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDECj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDECj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDECj$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[3] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- DIVAconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DIVA"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKE)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKE)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDIVALIKE$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[4] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- DIVAjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DIVAj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKEj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKEj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDIVALIKEj$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[5] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- BAYAREAconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "BAYAREA"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKE)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKE)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resBAYAREALIKE$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AICc[6] == min(restable$AICc)){
        BESTMODELconsensus[[i]] <- BAYAREAjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "BAYAREAj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKEj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKEj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resBAYAREALIKEj$relative_probs_of_each_state_at_bottom_of_root_branch)
      }
      
      # For groups with less than 4 species, AICc calculation returns an error. 
      # So we use AIC instead.
    } else {
      restable_list[[i]] <- restable
      
      # Model selection (min AIC): 1 is DEC, 2 is DEC+J, 3 is DIVA, 4 is DIVA+J,
      # 5 is BAYAREA, 6 is BAYAREA+J
      # and store the Q matrix of the best model for posterior simulations
      if (restable$AIC[1] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- DECconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DEC"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDEC)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDEC)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDEC$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[2] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- DECjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DECj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDECj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDECj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDECj$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[3] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- DIVAconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DIVA"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKE)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKE)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDIVALIKE$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[4] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- DIVAjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "DIVAj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKEj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resDIVALIKEj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resDIVALIKEj$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[5] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- BAYAREAconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "BAYAREA"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKE)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKE)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resBAYAREALIKE$relative_probs_of_each_state_at_bottom_of_root_branch)
        
      } else if (restable$AIC[6] == min(restable$AIC)){
        BESTMODELconsensus[[i]] <- BAYAREAjconsensus_output_biogeo[[i]]
        table_bestmodels[genera[i], ]$bestmodel <- "BAYAREAj"
        Qmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKEj)$Qmat
        COOmat_consensus[[i]] <- get_Qmat_COOmat_from_res_nexus(resBAYAREALIKEj)$COO_weights_columnar
        rootstate_consensus[[i]] <- which.max(resBAYAREALIKEj$relative_probs_of_each_state_at_bottom_of_root_branch)
      }
      
      
      
    }
    
  } else {
    print(paste(genera[i], " has less than 3 species! =(", sep=""))
  }
  
}




##### END OF THE CONSENSUS BIOGEOBEARS LOOP ####

#DECconsensus_output_biogeo
#DIVAconsensus_output_biogeo
#BAYAREAconsensus_output_biogeo

saveRDS(resDEC_list, "objects/resDEC_list.rds")
saveRDS(resDECj_list, "objects/resDECj_list.rds")
saveRDS(resDIVA_list, "objects/resDIVA_list.rds")
saveRDS(resDIVAj_list, "objects/resDIVAj_list.rds")
saveRDS(resBAYAREA_list, "objects/resBAYAREA_list.rds")
saveRDS(resBAYAREAj_list, "objects/resBAYAREAj_list.rds")

resDEC_list <- readRDS("objects/resDEC_list.rds")
resDIVA_list <- readRDS("objects/resDIVA_list.rds")
resBAYAREA_list <- readRDS("objects/resBAYAREA_list.rds")

resDECj_list <- readRDS("objects/resDECj_list.rds")
resDIVAj_list <- readRDS("objects/resDIVAj_list.rds")
resBAYAREAj_list <- readRDS("objects/resBAYAREAj_list.rds")

consensus_genera <- readRDS("objects/genera_trees.rds")
source("functions/modsel.R")
source("functions/get_Qmat_COOmat_from_res_nexus.R")
source("functions/get_Qmat_COOmat_from_BioGeoBEARS_run_object_nexus.R")
source("functions/check_trfn_nexus.R")
source("functions/readfiles_BioGeoBEARS_run_nexus.R")

source("functions/bears_optim_run_nexus.R")
source("functions/check_BioGeoBEARS_run_nexus.R")
source("functions/BioGeoBEARS_extract_Qmat_COOmat_v1.R")



# MODEL SELECTION ----
res_bgb_list_all <- list(resDEC_list, resDECj_list, resDIVA_list, resDIVAj_list,
                         resBAYAREA_list, resBAYAREAj_list)
names(res_bgb_list_all) <- c('DEC', 'DECJ', 'DIVA', 'DIVAJ', 'BAYAREA', 'BAYAREAJ')


# Run the model selection function
a <- modsel(results_list = res_bgb_list_all, treedata_list = consensus_genera)
length(res_bgb_list_all$DEC)
names(res_bgb_list_all$DEC)
names(consensus_genera)
length(consensus_genera)

# Save output:
BESTMODELconsensus <- a$BESTMODELconsensus
restable_list <- a$restable_list
table_bestmodels <- a$table_bestmodels
Qmat_consensus <- a$Qmat_consensus
COOmat_consensus <- a$COOmat_consensus
rootstate_consensus <- a$rootstate_consensus


saveRDS(BESTMODELconsensus, "objects/BESTMODELconsensus.rds")
saveRDS(restable_list, "objects/restable_list.rds")
saveRDS(table_bestmodels, "objects/table_bestmodels.rds")
saveRDS(Qmat_consensus, "objects/Qmat_consensus.rds")
saveRDS(COOmat_consensus, "objects/COOmat_consensus.rds")
saveRDS(rootstate_consensus, "objects/rootstate_consensus.rds")


table_bestmodels <- readRDS("objects/table_bestmodels.rds")
write.table(table_bestmodels, "objects/table_bestmodels.csv", sep=";", 
            row.names = FALSE, quote = FALSE)
BESTMODELconsensus <- readRDS("objects/BESTMODELconsensus.rds")

# MODEL SELECTION without J models----
res_bgb_list_all <- list(resDEC_list, resDIVA_list, 
                         resBAYAREA_list)
names(res_bgb_list_all) <- c('DEC', 'DIVA', 'BAYAREA')


# Run the model selection function
a <- modsel(results_list = res_bgb_list_all, treedata_list = consensus_genera)
length(res_bgb_list_all$DEC)
names(res_bgb_list_all$DEC)
names(consensus_genera)
length(consensus_genera)

# Save output:
BESTMODELconsensusnoJ <- a$BESTMODELconsensus
restable_listnoJ <- a$restable_list
table_bestmodelsnoJ <- a$table_bestmodels
Qmat_consensusnoJ <- a$Qmat_consensus
COOmat_consensusnoJ <- a$COOmat_consensus
rootstate_consensusnoJ <- a$rootstate_consensus


saveRDS(BESTMODELconsensusnoJ , "objects/noJ/BESTMODELconsensus.rds")
saveRDS(restable_listnoJ , "objects/noJ/restable_list.rds")
saveRDS(table_bestmodelsnoJ , "objects/noJ/table_bestmodels.rds")
saveRDS(Qmat_consensusnoJ , "objects/noJ/Qmat_consensus.rds")
saveRDS(COOmat_consensusnoJ , "objects/noJ/COOmat_consensus.rds")
saveRDS(rootstate_consensusnoJ , "objects/noJ/rootstate_consensus.rds")


table_bestmodelsnoJ <- readRDS("objects/noJ/table_bestmodels.rds")
write.table(table_bestmodelsnoJ, "objects/noJ/table_bestmodels.csv", sep=";", 
            row.names = FALSE, quote = FALSE)
BESTMODELconsensusnoJ <- readRDS("objects/noJ/BESTMODELconsensus.rds")


# In BESTMODELconsensus there are matrices with the probabilities of each biogeographic state
# in each tip and node.

##### TRANSFORM PROBABILITIES MATRICES INTO VECTORS OF STATES #####
for(i in 1:length(location_files)){
  tipranges_list[[i]] = getranges_from_LagrangePHYLIP(lgdata_fn=location_files[i])
}

# In consensus genera
consensus_states <- vector("list", length(genera))
consensus_statesnoj <- vector("list", length(genera))
names(consensus_states) <- genera
names(consensus_statesnoj) <- genera
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    areas = getareas_from_tipranges_object(tipranges_list[[i]])
    statenames = areas_list_to_states_list_new(areas, maxareas = max_range_size,
                                               include_null_range = TRUE,
                                               split_ABC = FALSE)
    statenames = unlist(statenames)
    relprobs_matrix_noj <- BESTMODELconsensusnoJ[[i]] # NoJ
    relprobs_matrix <- BESTMODELconsensus[[i]] #With J
    MLprobs = get_ML_probs(relprobs_matrix)
    MLprobs_noj = get_ML_probs(relprobs_matrix_noj)
    MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames,
                                           returnwhat = "states", if_ties = "takefirst")
    MLstates_noj = get_ML_states_from_relprobs(relprobs_matrix_noj, statenames,
                                           returnwhat = "states", if_ties = "takefirst")
    consensus_states[[i]] <- MLstates
    consensus_statesnoj[[i]] <- MLstates_noj
  }
}

#In consensus genera 
consensus_probs <- vector("list", length(genera))
consensus_probsnoj <- vector("list", length(genera))
names(consensus_probs) <- genera
names(consensus_probsnoj) <- genera
for (i in 1:length(genera)) {
  areas = getareas_from_tipranges_object(tipranges_list[[i]])
  statenames = areas_list_to_states_list_new(areas, maxareas = max_range_size,
                                             include_null_range = TRUE,
                                             split_ABC = FALSE)
  statenames = unlist(statenames)
  relprobs_matrixnoj <- BESTMODELconsensusnoJ[[i]] %>% round(digits = 2) %>% as.data.frame() # NO J
  relprobs_matrix <- BESTMODELconsensus[[i]] %>% round(digits = 2) %>% as.data.frame() # With J
  names(relprobs_matrix) <- statenames
  names(relprobs_matrixnoj) <- statenames
  consensus_probs[[i]] <- relprobs_matrix
  consensus_probsnoj[[i]] <- relprobs_matrixnoj
}

# consensus_states is a list with vectors of states for each genera (consensus)

saveRDS(consensus_statesnoj, "objects/noJ/consensus_states.rds")
saveRDS(consensus_probsnoj, "objects/noJ/consensus_probs.rds")
saveRDS(consensus_states, "objects/consensus_states.rds")
saveRDS(consensus_probs, "objects/consensus_probs.rds")


consensus_statesnoj <- readRDS("objects/noJ/consensus_states.rds")
consensus_probsnoj <- readRDS("objects/noJ/consensus_probs.rds")
consensus_states <- readRDS("objects/consensus_states.rds")
consensus_probs <- readRDS("objects/consensus_probs.rds")

consensus_all <- c()
for (i in 1:length(consensus_statesnoj)) {
  consensus_all <- c(consensus_all, consensus_statesnoj[[i]])
  }
consensus_all <- unique(consensus_all)


##### DEFINE BIOGEOGRAPHIC SCENARIOS #####
## To do so we can move to the script plot_cumulative_plots_no_intermediate_states.R where all biogeobears states will be taken into account but also all the 
# Simmaps first colonization events and the Pristurus gallagheri colonization and speciation which can't be caluclated through Biogeobears. 



#######################
##### SIMULATIONS #####
#######################
file.sources = list.files("functions/", pattern="*.R", full.names=T)
#file.sources <- file.sources[-6]
source("functions/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
  for (i in 1:length(file.sources)){
  source(file.sources[10])
}

#saveRDS(BESTMODELconsensus, "objects/BESTMODELconsensus.rds")
#saveRDS(restable_list, "objects/restable_list.rds")
#saveRDS(table_bestmodels, "objects/table_bestmodels.rds")
#saveRDS(Qmat_consensus, "objects/Qmat_consensus.rds")
#saveRDS(COOmat_consensus, "objects/COOmat_consensus.rds")
#saveRDS(rootstate_consensus, "objects/rootstate_consensus.rds")

Qmat_consensus <- readRDS("objects/noJ/Qmat_consensus.rds")
COOmat_consensus <- readRDS("objects/noJ/COOmat_consensus.rds")
rootstate_consensus <- readRDS("objects/noJ/rootstate_consensus.rds")


##### Simulations with consensus tree #####
sims_consensus <- vector("list", length(genera))
names(sims_consensus) <- genera
nsim <- 1000
for (i in 1:length(genera)){
  if (Ntip(consensus_genera[[i]]) > 2){
    for (s in 1:nsim){
      sims_consensus[[i]][[s]] <- simulate_biogeog_history(phy=as.phylo(consensus_genera[[i]]), Qmat=Qmat_consensus[[i]],
                                                           COO_probs_columnar=COOmat_consensus[[i]],
                                                           index_Qmat_0based_of_starting_state=rootstate_consensus[[i]]-1)
    }
  }
}

# Add 1 to all the states because the function doesn't count the null state (state 1).
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (!is.null(sims_consensus[[i]])){
      sims_consensus[[i]][[s]] <- sims_consensus[[i]][[s]] + 1
    }
  }
}

# Names of the states for each genus
consensus_states_names <- vector("list", length(genera))
names(consensus_states_names) <- genera
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    areas = getareas_from_tipranges_object(tipranges_list[[i]])
    consensus_states_names[[i]] = areas_list_to_states_list_new(areas, maxareas = max_range_size,
                                                                include_null_range = TRUE,
                                                                split_ABC = FALSE)
    consensus_states_names[[i]] = unlist(consensus_states_names[[i]])
    
  }
}

# Now we can transform the numeric states of sims_consensus into the names (F, R, etc)
sims_consensus_letters <- vector("list", length(genera))
names(sims_consensus_letters) <- genera
for (i in 1:length(genera)){
  sims_consensus_letters[[i]] <- vector("list", nsim)
}

for (i in 1:length(genera)){
  for (s in 1:nsim){
    sims_consensus_letters[[i]][[s]] <- vector(mode="character")
  }
}

for (i in 1:length(genera)){
  if (Ntip(consensus_genera[[i]]) > 2){
    for (s in 1:nsim){
      for (j in 1:length(sims_consensus[[i]][[s]])){
        sims_consensus_letters[[i]][[s]][j] <- consensus_states_names[[i]][sims_consensus[[i]][[s]][j]]
      }
    }
  }
}

saveRDS(sims_consensus_letters, "objects/sims_consensus_letters.rds")
sims_consensus_letters$Asaccus[25]


##### .-.-.-.-. COUNTING AND CATEGORIZING BIOGEOGRAPHIC EVENTS IN SIMULATIONS (CONSENSUS) #####
#consensus_genera$Acanthodactylus
#consensus_states[[1]][60]

Af2Ar_SIM_nodes_list_cons <- vector("list", length(genera))
names(Af2Ar_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  Af2Ar_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    Af2Ar_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

Ar2Af_SIM_nodes_list_cons <- vector("list", length(genera))
names(Ar2Af_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  Ar2Af_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    Ar2Af_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

vicariance_SIM_nodes_list_cons <- vector("list", length(genera))
names(vicariance_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  vicariance_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    vicariance_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

extirpfromAf_SIM_nodes_list_cons <- vector("list", length(genera))
names(extirpfromAf_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  extirpfromAf_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    extirpfromAf_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

extirpfromAr_SIM_nodes_list_cons <- vector("list", length(genera))
names(extirpfromAr_SIM_nodes_list_cons) <- genera
for (i in 1:length(genera)){
  extirpfromAr_SIM_nodes_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    extirpfromAr_SIM_nodes_list_cons[[i]][[s]] <- 0
  }
}

# Get vicariance nodes (the parental node)
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    for (s in 1:nsim){
      for (j in min(as.phylo(consensus_genera[[i]])$edge[,1]):max(as.phylo(consensus_genera[[i]])$edge[,1])){
        v <- c(sims_consensus_letters[[i]][[s]][j],
               sims_consensus_letters[[i]][[s]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][1],
               sims_consensus_letters[[i]][[j]][as.phylo(consensus_genera[[i]])$edge[,2][as.phylo(consensus_genera[[i]])$edge[,1]==j]][2])
        if (list(v) %in% vicariance){
          vicariance_SIM_nodes_list_cons[[i]][[s]] <- append(vicariance_SIM_nodes_list_cons[[i]][[s]], j)
        }}}}}

# Get dispersal and extirpation nodes (the descendent node)
'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in 1:length(genera)){
  if (length(as.phylo(consensus_genera[[i]])$tip.label) >= 3){
    for (s in 1:nsim){
      for (j in sort(as.phylo(consensus_genera[[i]])$edge[,2])){
        parent <- sims_consensus_letters[[i]][[s]][as.phylo(consensus_genera[[i]])$edge[,1][as.phylo(consensus_genera[[i]])$edge[,2]==j]]
        son <- sims_consensus_letters[[i]][[s]][j]
        d <- c(parent, son)
        if (list(d) %in% dispersal_Af2Ar){
          Af2Ar_SIM_nodes_list_cons[[i]][[s]] <- append(Af2Ar_SIM_nodes_list_cons[[i]][[s]], j)
          
        }
        if (list(d) %in% dispersal_Ar2Af){
          Ar2Af_SIM_nodes_list_cons[[i]][[s]] <- append(Ar2Af_SIM_nodes_list_cons[[i]][[s]], j)
        }
        
        # extirpation (if the parent node is not a vicariant node)
        if (as.phylo(consensus_genera[[i]])$edge[,1][as.phylo(consensus_genera[[i]])$edge[,2]==j] %!in% vicariance_SIM_nodes_list_cons[[i]][[s]]){
          if (list(d) %in% extirpation_from_Ar){
            extirpfromAr_SIM_nodes_list_cons[[i]][[s]] <- append(extirpfromAr_SIM_nodes_list_cons[[i]][[s]], j)
          }
          if (list(d) %in% extirpation_from_Af){
            extirpfromAf_SIM_nodes_list_cons[[i]][[s]] <- append(extirpfromAf_SIM_nodes_list_cons[[i]][[s]], j)
          }
        }  
      }
    }
  }
}





Af2Ar_SIM_nodes_list_cons[[10]][[s]]

# Remove zeros.

for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (length(Af2Ar_SIM_nodes_list_cons[[i]][[s]]) != 1){
      Af2Ar_SIM_nodes_list_cons[[i]][[s]] <- Af2Ar_SIM_nodes_list_cons[[i]][[s]][-1]
    }
    if (length(Ar2Af_SIM_nodes_list_cons[[i]][[s]]) != 1){
      Ar2Af_SIM_nodes_list_cons[[i]][[s]] <- Ar2Af_SIM_nodes_list_cons[[i]][[s]][-1]
    }
    if (length(vicariance_SIM_nodes_list_cons[[i]][[s]]) != 1){
      vicariance_SIM_nodes_list_cons[[i]][[s]] <- vicariance_SIM_nodes_list_cons[[i]][[s]][-1]
    }
    if (length(extirpfromAf_SIM_nodes_list_cons[[i]][[s]]) != 1){
      extirpfromAf_SIM_nodes_list_cons[[i]][[s]] <- extirpfromAf_SIM_nodes_list_cons[[i]][[s]][-1]
    }
    if (length(extirpfromAr_SIM_nodes_list_cons[[i]][[s]]) != 1){
      extirpfromAr_SIM_nodes_list_cons[[i]][[s]] <- extirpfromAr_SIM_nodes_list_cons[[i]][[s]][-1]
    }
  }
}

# Make data frames
Af2Ar_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(Af2Ar_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  Af2Ar_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}

Ar2Af_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(Ar2Af_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  Ar2Af_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}

vicariance_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(vicariance_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  vicariance_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}

extirpfromAf_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(extirpfromAf_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  extirpfromAf_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}

extirpfromAr_SIM_nodes_dfs_cons <- vector("list", length(genera))
names(extirpfromAr_SIM_nodes_dfs_cons) <- genera
for (i in 1:length(genera)){
  extirpfromAr_SIM_nodes_dfs_cons[[i]] <- vector("list", nsim)
}


for (i in 1:length(genera)){
  for (s in 1:nsim){
    Af2Ar_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=Af2Ar_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                     height=0, min=0, max=0, Af2Ar=1, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=0)
    Ar2Af_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=Ar2Af_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                     height=0, min=0, max=0, Af2Ar=0, Ar2Af=1, Vic=0, ExtAf=0, ExtAr=0)
    vicariance_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=vicariance_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                          height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=1, ExtAf=0, ExtAr=0)
    extirpfromAf_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=extirpfromAf_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                            height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=1, ExtAr=0)
    extirpfromAr_SIM_nodes_dfs_cons[[i]][[s]] <- data.frame(node=extirpfromAr_SIM_nodes_list_cons[[i]][[s]], genus=genera[i],
                                                            height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=1)
  }
}

# We want a list with one data frame per genus per simulation, including all types of events.
sim_event_consensus <- vector("list", length(genera))
names(sim_event_consensus) <- genera
for (i in 1:length(genera)){
  sim_event_consensus[[i]] <- vector("list", nsim)
}

for (i in 1:length(genera)){
  for (s in 1:nsim){
    sim_event_consensus[[i]][[s]] <- rbind(Af2Ar_SIM_nodes_dfs_cons[[i]][[s]], Ar2Af_SIM_nodes_dfs_cons[[i]][[s]],
                                           vicariance_SIM_nodes_dfs_cons[[i]][[s]], extirpfromAf_SIM_nodes_dfs_cons[[i]][[s]],
                                           extirpfromAr_SIM_nodes_dfs_cons[[i]][[s]])
    
  }
}
sim_event_consensus$Mesalina[[729]]

# Now we already have one data frame per genus per sim with the transition nodes. We can remove the zeros:
for (i in 1:length(genera)){
  for (s in 1:nsim){
    sim_event_consensus[[i]][[s]] <- sim_event_consensus[[i]][[s]][sim_event_consensus[[i]][[s]]$node!=0,]
  }
}

# Now we need to include information about the divergence times: the height of the nodes with transitions.
##### INCLUDE INFO ON DIVERGENCE TIMES FOR SIMS #####

# Divergence times and 95% HPD (min and max) for each node
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:length(sim_event_consensus[[i]][[s]]$node)){
        
        # Vicariance nodes:
        if (sim_event_consensus[[i]][[s]][j, "Vic"] == 1){
          sim_event_consensus[[i]][[s]]$min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[sim_event_consensus[[i]][[s]]$node[j]]][1]
          sim_event_consensus[[i]][[s]]$max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[sim_event_consensus[[i]][[s]]$node[j]]][2]
          sim_event_consensus[[i]][[s]]$height[j] <- as_tibble(consensus_genera[[i]])$height_median[[sim_event_consensus[[i]][[s]]$node[j]]]
        }
        
        # Dispersal and extirpation nodes:
        if (sim_event_consensus[[i]][[s]][j, "Af2Ar"] == 1 | sim_event_consensus[[i]][[s]][j, "Ar2Af"] == 1 | sim_event_consensus[[i]][[s]][j, "ExtAf"] == 1 | sim_event_consensus[[i]][[s]][j, "ExtAr"] == 1){
          dat <- as_tibble(consensus_genera[[i]])
          # Min is the age of the descendent node
          desc_age <- dat$height[dat$node == sim_event_consensus[[i]][[s]]$node[j]]
          sim_event_consensus[[i]][[s]]$min[j] <- desc_age
          
          # Max is the age of the parental node
          parent_age <- dat$height[dat$parent[dat$node == sim_event_consensus[[i]][[s]]$node[j]]]
          sim_event_consensus[[i]][[s]]$max[j] <- parent_age
          
          # let's set the height as the midpoint between the parental and the descendent nodes
          sim_event_consensus[[i]][[s]]$height[j] <- mean(c(desc_age, parent_age))
        }
        
      }
    }
  }
}





#for (i in 1:length(genera)){
#  for (s in 1:nsim){
#    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
#      for (j in 1:length(sim_event_consensus[[i]][[s]]$node)){
#        sim_event_consensus[[i]][[s]]$min[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[sim_event_consensus[[i]][[s]]$node[j]]][1]
#        sim_event_consensus[[i]][[s]]$max[j] <- as_tibble(consensus_genera[[i]])$height_0.95_HPD[[sim_event_consensus[[i]][[s]]$node[j]]][2]
#        sim_event_consensus[[i]][[s]]$height[j] <- as_tibble(consensus_genera[[i]])$height_median[[sim_event_consensus[[i]][[s]]$node[j]]]
#      }
#    }
#  }
#}

saveRDS(sim_event_consensus, "objects/sim_event_cons.rds")

##### EVENTS PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
myr <- 60
sim_event_myr_list_cons <- vector("list", length(genera))
names(sim_event_myr_list_cons) <- genera
for (i in 1:length(genera)){
  sim_event_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    sim_event_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}
sim_event_myr_list_cons$Acanthodactylus[[1]][1,]$N

# We also can create a list of vectors representing intervals of 1 Ma.
ma <- vector("list", myr)
for(i in 1:myr){
  ma[[i]] <- c(i-1, i)
}

# Counting the maximum number of transition in each million year for each genus and each simulation.
sim_event_myr_list_cons0 <- sim_event_myr_list_cons
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
        for (k in 1:length(ma)){
          if (Overlap(na.omit(vi), ma[[k]]) != 0){
            sim_event_myr_list_cons[[i]][[s]][k,]$N <- sim_event_myr_list_cons[[i]][[s]][k,]$N + 1
          }}}}}}


# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
sim_event_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  sim_event_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    sim_event_all_cons[[s]]$N <- sim_event_all_cons[[s]]$N + sim_event_myr_list_cons[[i]][[s]]$N
  }
}

##### AFRICA TO ARABIA PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
Af2Ar_sim_myr_list_cons <- vector("list", length(genera))
names(Af2Ar_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  Af2Ar_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    Af2Ar_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}
Af2Ar_sim_myr_list_cons$Acanthodactylus[[1]][1,]$N

# Counting the maximum number of Af -> Ar dispersals in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$Af2Ar[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              Af2Ar_sim_myr_list_cons[[i]][[s]][k,]$N <- Af2Ar_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
Af2Ar_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  Af2Ar_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    Af2Ar_sim_all_cons[[s]]$N <- Af2Ar_sim_all_cons[[s]]$N + Af2Ar_sim_myr_list_cons[[i]][[s]]$N
  }
}

##### ARABIA TO AFRICA PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
Ar2Af_sim_myr_list_cons <- vector("list", length(genera))
names(Ar2Af_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  Ar2Af_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    Ar2Af_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}
Ar2Af_sim_myr_list_cons$Acanthodactylus[[1]][1,]$N

# Counting the maximum number of Ar -> Af dispersals in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$Ar2Af[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              Ar2Af_sim_myr_list_cons[[i]][[s]][k,]$N <- Ar2Af_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
Ar2Af_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  Ar2Af_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    Ar2Af_sim_all_cons[[s]]$N <- Ar2Af_sim_all_cons[[s]]$N + Ar2Af_sim_myr_list_cons[[i]][[s]]$N
  }
}

##### VICARIANCE PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
vicariance_sim_myr_list_cons <- vector("list", length(genera))
names(vicariance_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  vicariance_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    vicariance_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}
vicariance_sim_myr_list_cons$Acanthodactylus[[1]][1,]$N

# Counting the maximum number of vicariance in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$Vic[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              vicariance_sim_myr_list_cons[[i]][[s]][k,]$N <- vicariance_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
vicariance_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  vicariance_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    vicariance_sim_all_cons[[s]]$N <- vicariance_sim_all_cons[[s]]$N + vicariance_sim_myr_list_cons[[i]][[s]]$N
  }
}

##### EXTIRPATION FROM AFRICA PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
extirp_from_Af_sim_myr_list_cons <- vector("list", length(genera))
names(extirp_from_Af_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  extirp_from_Af_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    extirp_from_Af_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}
extirp_from_Af_sim_myr_list_cons$Acanthodactylus[[1]][1,]$N

# Counting the maximum number of extirpation from Africa in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$ExtAf[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              extirp_from_Af_sim_myr_list_cons[[i]][[s]][k,]$N <- extirp_from_Af_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
extirp_from_Af_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  extirp_from_Af_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    extirp_from_Af_sim_all_cons[[s]]$N <- extirp_from_Af_sim_all_cons[[s]]$N + extirp_from_Af_sim_myr_list_cons[[i]][[s]]$N
  }
}

##### EXTIRPATION FROM ARABIA PER MILLION YEAR FOR EACH GENUS AND SIMULATION (CONSENSUS) #####
extirp_from_Ar_sim_myr_list_cons <- vector("list", length(genera))
names(extirp_from_Ar_sim_myr_list_cons) <- genera
for (i in 1:length(genera)){
  extirp_from_Ar_sim_myr_list_cons[[i]] <- vector("list", nsim)
}
for (i in 1:length(genera)){
  for (s in 1:nsim){
    extirp_from_Ar_sim_myr_list_cons[[i]][[s]] <- data.frame(Ma=c(1:60), N=0)
  }
}
extirp_from_Ar_sim_myr_list_cons$Acanthodactylus[[1]][1,]$N

# Counting the maximum number of Extirpation from Arabia in each million year for each genus.
for (i in 1:length(genera)){
  for (s in 1:nsim){
    if (nrow(sim_event_consensus[[i]][[s]]) > 0){
      for (j in 1:nrow(sim_event_consensus[[i]][[s]])){
        if (sim_event_consensus[[i]][[s]]$ExtAr[j] > 0){
          vi <- c(sim_event_consensus[[i]][[s]]$min[j], sim_event_consensus[[i]][[s]]$max[j])
          for (k in 1:length(ma)){
            if (Overlap(na.omit(vi), ma[[k]]) != 0){
              extirp_from_Ar_sim_myr_list_cons[[i]][[s]][k,]$N <- extirp_from_Ar_sim_myr_list_cons[[i]][[s]][k,]$N + 1
            }
          }
        }
      }
    }
  }
}

# Sum all the transitions per Ma (of all the genera) for each simulation
# Make a data frame with number of events in each 1 Mya period.
extirp_from_Ar_sim_all_cons <- vector("list", nsim)
for (s in 1:nsim){
  extirp_from_Ar_sim_all_cons[[s]] <- data.frame(Ma=c(1:myr), N=0)
}

for (s in 1:nsim){
  for (i in 1:length(genera)){
    extirp_from_Ar_sim_all_cons[[s]]$N <- extirp_from_Ar_sim_all_cons[[s]]$N + extirp_from_Ar_sim_myr_list_cons[[i]][[s]]$N
  }
}


##### PLOT SIMULATIONS CONSENSUS #####
# All events
plot(1,type='n',xlim=c(myr,0),ylim=c(0,60),xlab='Ma', ylab='N', main="All simulated events (consensus)")
# smooth
for (i in 1:nsim){
  lines(spline(sim_event_all_cons[[i]]$Ma-1, sim_event_all_cons[[i]]$N),
        type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
}

# not smooth
for (i in 1:nsim){
  lines(sim_event_all_cons[[i]]$Ma-1, sim_event_all_cons[[i]]$N, type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
}

#for (i in 1:length(genera)){
#  lines(event_myr_list_cons[[i]]$Ma-1, event_myr_list_cons[[i]]$N, type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
#}


# Africa -> Arabia dispersal
plot(1,type='n',xlim=c(myr,0),ylim=c(0,40),xlab='Ma', ylab='N', main="Africa to Arabia simulated dispersal (consensus)")
# smooth
for (i in 1:nsim){
  lines(spline(Af2Ar_sim_all_cons[[i]]$Ma-1, Af2Ar_sim_all_cons[[i]]$N), type="l", col="green", lwd=0.25, pch=16, cex=0.7)
}
lines(smooth.spline(Af2Ar_all_cons$Ma-1, Af2Ar_all_cons$N, spar=0.2), type="l", col="red", lwd=5, pch=16, cex=1)

plot(1,type='n',xlim=c(myr,0),ylim=c(0,40),xlab='Ma', ylab='N', main="Africa to Arabia simulated dispersal (consensus)")
# not smooth
for (i in 1:nsim){
  lines(Af2Ar_sim_all_cons[[i]]$Ma-1, Af2Ar_sim_all_cons[[i]]$N, type="l", col="green", lwd=0.25, pch=16, cex=0.7)
}
lines(Af2Ar_all_cons$Ma-1, Af2Ar_all_cons$N, type="l", col="red", lwd=5, pch=16, cex=1)




# Arabia -> Africa dispersal
plot(1,type='n',xlim=c(myr,0),ylim=c(0,35),xlab='Ma', ylab='N', main="Arabia to Africa simulated dispersal (consensus)")
# smooth
for (i in 1:nsim){
  lines(spline(Ar2Af_sim_all_cons[[i]]$Ma-1, Ar2Af_sim_all_cons[[i]]$N), type="l", col="purple", lwd=0.25, pch=16, cex=0.7)
}

# not smooth
for (i in 1:nsim){
  lines(Ar2Af_sim_all_cons[[i]]$Ma-1, Ar2Af_sim_all_cons[[i]]$N, type="l", col="purple", lwd=0.25, pch=16, cex=0.7)
}
lines(Ar2Af_all_cons$Ma-1, Ar2Af_all_cons$N, type="l", col="red", lwd=5, pch=16, cex=1)


# Vicariance
plot(1,type='n',xlim=c(myr,0),ylim=c(0,10),xlab='Ma', ylab='N', main="Simulated vicariance (consensus)")
# smooth
for (i in 1:nsim){
  lines(spline(vicariance_sim_all_cons[[i]]$Ma-1, vicariance_sim_all_cons[[i]]$N), type="l", col="gold", lwd=0.25, pch=16, cex=0.7)
}
# not smooth
for (i in 1:nsim){
  lines(vicariance_sim_all_cons[[i]]$Ma-1, vicariance_sim_all_cons[[i]]$N, type="l", col="gold", lwd=0.25, pch=16, cex=0.7)
}
lines(vicariance_all_cons$Ma-1, vicariance_all_cons$N, type="l", col="red", lwd=5, pch=16, cex=1)


# Extirpation from Africa
plot(1,type='n',xlim=c(myr,0),ylim=c(0,10),xlab='Ma', ylab='N', main="Simulated extirpation from Africa (consensus)")
# smooth
for (i in 1:nsim){
  lines(spline(extirp_from_Af_sim_all_cons[[i]]$Ma-1, extirp_from_Af_sim_all_cons[[i]]$N), type="l", col="blue", lwd=0.25, pch=16, cex=0.7)
}
# not smooth
for (i in 1:nsim){
  lines(extirp_from_Af_sim_all_cons[[i]]$Ma-1, extirp_from_Af_sim_all_cons[[i]]$N, type="l", col="blue", lwd=0.25, pch=16, cex=0.7)
}
lines(extirp_from_Af_all_cons$Ma-1, extirp_from_Af_all_cons$N, type="l", col="red", lwd=5, pch=16, cex=1)



# Extirpation from Arabia
plot(1,type='n',xlim=c(myr,0),ylim=c(0,10),xlab='Ma', ylab='N', main="Simulated extirpation from Arabia (consensus)")
# smooth
for (i in 1:nsim){
  lines(spline(extirp_from_Ar_sim_all_cons[[i]]$Ma-1, extirp_from_Ar_sim_all_cons[[i]]$N), type="l", col="orange", lwd=0.25, pch=16, cex=0.7)
}
# not smooth
for (i in 1:nsim){
  lines(extirp_from_Ar_sim_all_cons[[i]]$Ma-1, extirp_from_Ar_sim_all_cons[[i]]$N, type="l", col="orange", lwd=0.25, pch=16, cex=0.7)
}
lines(extirp_from_Ar_all_cons$Ma-1, extirp_from_Ar_all_cons$N, type="l", col="red", lwd=5, pch=16, cex=1)

# Plot all
legend_colors <- c("red", "green", "purple", "gold", "blue", "orange")
legend_text <- c("All events", "Af -> Ar", "Ar -> Af", "Vicariance", "Extirp. from Af", "Extirp. from Ar")
legend("topleft", legend=legend_text, fill=legend_colors, bty="n")
?legend

# plot all less smooth
legend("topleft", legend=legend_text, fill=legend_colors, bty="n")

# plot all no smooth
legend("topleft", legend=legend_text, fill=legend_colors, bty="n")


# 95% CI ALL EVENTS
## quantiles and mean of each My.
prob_qup <- 0.975
prob_qlow <- 0.025

# Let's do a list where each element is a vector with the
# number of transitions per simulations for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist <- vector("list", myr)
for (i in 1:myr){
  dlist[[i]] <- vector("numeric", nsim)
}
dlist[[1]][2]

for (i in 1:nsim){
  for (j in 1:myr){
    dlist[[j]][i] <- sim_event_all_cons[[i]]$N[j]
  }
}

# Now we can create a dataframe with the quantiles and
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df$mean[i] <- mean(dlist[[i]])
  Q_df$qlow[i] <- quantile(dlist[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df$qup[i] <- quantile(dlist[[i]], probs=c(prob_qlow, prob_qup))[2]
}

# PLOT SMOOTH
plot(1,type='n',xlim=c(myr,0),ylim=c(0,60),xlab='Ma', ylab='N', main="Simulated transitions")
for (i in 1:nsim){
  lines(spline(sim_event_all_cons[[i]]$Ma-1, sim_event_all_cons[[i]]$N),
        type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
}

# mean and quantiles
# as.data.frame for ggplot2 (geom_line())
#lines(as.data.frame(spline(Q_df$Ma-1, Q_df$mean)), type="l", col="purple", lwd=2, pch=16, cex=1)
lines(spline(Q_df$Ma-1, Q_df$mean), type="l", col="purple", lwd=2, pch=16, cex=1)
lines(spline(Q_df$Ma-1, Q_df$qlow), type="l", col="black", lwd=1, pch=16, cex=1)
lines(spline(Q_df$Ma-1, Q_df$qup), type="l", col="black", lwd=1, pch=16, cex=1)

# empiric
lines(spline(event_all_cons$Ma-1, event_all_cons$N), type="l", col="red", lwd=5, pch=16, cex=1)

#for (i in 1:nsim){
#  lines(smooth.spline(sim_event_all_cons[[i]]$Ma-1, sim_event_all_cons[[i]]$N, spar=0.1, cv=F),
#        type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
#}
#lines(smooth.spline(Q_df$Ma-1, Q_df$mean), type="l", col="purple", lwd=2, pch=16, cex=1)
#lines(smooth.spline(Q_df$Ma-1, Q_df$qlow), type="l", col="black", lwd=1, pch=16, cex=1)
#lines(smooth.spline(Q_df$Ma-1, Q_df$qup), type="l", col="black", lwd=1, pch=16, cex=1)

# PLOT NO SMOOTH
plot(1,type='n',xlim=c(myr,0),ylim=c(0,60),xlab='Ma', ylab='N', main="Simulated transitions consensus")
for (i in 1:nsim){
  lines(sim_event_all_cons[[i]]$Ma-1, sim_event_all_cons[[i]]$N, type="l", col="grey", lwd=0.25, pch=16, cex=0.7)
}
lines(event_all_cons$Ma-1, event_all_cons$N, type="l", col="red", lwd=5, pch=16, cex=1)

lines(Q_df$Ma-1, Q_df$mean, type="l", col="purple", lwd=2, pch=16, cex=1)
lines(Q_df$Ma-1, Q_df$qlow, type="l", col="black", lwd=1, pch=16, cex=1)
lines(Q_df$Ma-1, Q_df$qup, type="l", col="black", lwd=1, pch=16, cex=1)

lines(event_all_cons$Ma-1, event_all_cons$N, type="l", col="red", lwd=5, pch=16, cex=1)

# Plot to save
plot(1,type='n',xlim=c(myr,0),ylim=c(0,60),xlab='Ma', ylab='N', main="Simulated transitions")


##### Save event objects #####
saveRDS(sim_event_all_cons, "objects/sim_event_all_cons.rds")
saveRDS(event_all_cons, "objects/event_all_cons.rds")
saveRDS(Af2Ar_sim_all_cons, "objects/Af2Ar_sim_all_cons.rds")
saveRDS(Af2Ar_all_cons, "objects/Af2Ar_all_cons.rds")
saveRDS(Ar2Af_sim_all_cons, "objects/Ar2Af_sim_all_cons.rds")
saveRDS(Ar2Af_all_cons, "objects/Ar2Af_all_cons.rds")
saveRDS(vicariance_sim_all_cons, "objects/vicariance_sim_all_cons.rds")
saveRDS(vicariance_all_cons, "objects/vicariance_all_cons.rds")
saveRDS(extirp_from_Af_sim_all_cons, "objects/extirp_from_Af_sim_all_cons.rds")
saveRDS(extirp_from_Af_all_cons, "objects/extirp_from_Af_all_cons.rds")
saveRDS(extirp_from_Ar_sim_all_cons, "objects/extirp_from_Ar_sim_all_cons.rds")
saveRDS(extirp_from_Ar_all_cons, "objects/extirp_from_Ar_all_cons.rds")

############ THE END #################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
