


rm(list = ls())
libs <- c("tidyverse", "deeptime", "here", "cowplot")
library(RCurl)
library(gridExtra)
library(plyr)
lapply(libs, require, character.only = TRUE)
# https://github.com/willgearty/deeptime


setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography")




# Import event objects ----


sim_event_all_cons <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_SIM_plots_all_no_diversification.rds")
event_all_cons <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_all_biogeo_no_diversification.rds")



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
epochs_htc$abbr[epochs_htc$abbr == "Mc"] <- "Miocene"
epochs_htc$abbr[epochs_htc$abbr == "Ol"] <- "Oligocene"
epochs_htc$abbr[epochs_htc$abbr == "Palcn"] <- "P"

# Quantiles information ----
prob_qup <- 0.975
prob_qlow <- 0.025
nsim <- 1000
myr <- 40


# Parameters (color, size, transparency) ----
# Colors
{color_sim <- "gray93"
col_mean <- 'black'
col_q <- 'black'
color_Ccol <- "#b22125"
color_Ecol <- "#fec250"
color_Wcol <- "#206ab4"
color_Cdiv <- "salmon"
color_Ediv <- "gold"
color_Wdiv <- "lightblue"
color_vicCE <- "#ef8354"
color_vicCW <- "#33F1FF"
color_vicWE <- "#CCF5AC"
color_vicIW <- "#bbbbbb"
color_all_events <- "#EE6A50"
color_all_vicariances <- "#52ddb9"
}

# Line width (size)
lwd_q <- 0.5
lwd_mean <- 0.3
lwd_events <- 2
lwd_sim <- 0.3

# Transparency (alpha)
alpha_sim <- 0.9

######### :::::::::::::::::::::::::::::::::::::::::::::::::::::: #############
# ALL BIOGEOGRAPHIC EVENTS ----

# Create the lines to plot ----
# spline consensus
line_all <- vector("list", length(event_all_cons)) 
names(line_all) <- names(event_all_cons)
for(i in 1:length(line_all)){   
   line_all[[i]] <- data.frame(spline(event_all_cons[[i]]$Ma-1, event_all_cons[[i]]$N))
colnames(line_all[[i]]) <- c("time", "events")
}
# spline simulations
lines_sim_all <- vector("list", length(sim_event_all_cons))
names(lines_sim_all) <- names(sim_event_all_cons)
for(j in 1:length(lines_sim_all)){
   lines_sim_all[[j]] <- vector("list", nsim)
   for (i in 1:nsim){
      lines_sim_all[[j]][[i]] <- data.frame(spline(sim_event_all_cons[[j]][[i]]$Ma-1,
                                              sim_event_all_cons[[j]][[i]]$N))
      colnames(lines_sim_all[[j]][[i]]) <- c("time", "events")
      
   }}

## spline maximum and minimimum 
lines_sim_max_min <- vector("list", length(sim_event_all_cons))
names(lines_sim_max_min) <- names(sim_event_all_cons)
for(j in 1:length(lines_sim_max_min)){
  lines_sim_max_min[[j]] <- vector("list", 2)
  lines_sim_max_min[[j]][[1]] <- lines_sim_all[[j]][[1]]
  lines_sim_max_min[[j]][[2]] <- lines_sim_all[[j]][[1]]
  names(lines_sim_max_min[[j]]) <- c ("max", "min")
}
 
for(j in 1:length(lines_sim_max_min)){
   for (z in 1:nrow(lines_sim_all[[j]][[2]])) {
      val = data.frame(time = NA, events = NA) 
        for (i in 1:length(lines_sim_all[[j]])){
        val = rbind(val,lines_sim_all[[j]][[i]][z,]) 
        }
      lines_sim_max_min[[j]][[1]][z,] <- c(val$time[2], max(val$events, na.rm = T))
      lines_sim_max_min[[j]][[2]][z,] <- c(val$time[2], min(val$events,  na.rm = T))
  }}

saveRDS(lines_sim_max_min, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/lines_sim_min_max.rds")
lines_sim_max_min <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/lines_sim_min_max.rds")
lines_sim_max_min[[1]][[2]][1,]
lines_sim_max_min$event_all_cons

# Quantiles and mean ----
# 95% CI ALL EVENTS
# Let's do a list where each element is a vector with the
# number of transitions per simulations for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist_all <- vector("list", length(lines_sim_all))
names(dlist_all) <- names(lines_sim_all)
for (i in 1:length(dlist_all)) {
   dlist_all[[i]] <- vector("list",  myr)
   for (j in 1:myr){
      dlist_all[[i]][[j]] <- vector("numeric", nsim)
   }}

for(k in 1:length(dlist_all)){
for (i in 1:nsim){
  for (j in 1:myr){
    dlist_all[[k]][[j]][i] <- sim_event_all_cons[[k]][[i]]$N[j]
  }}}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df_all <- vector("list", length(lines_sim_all))
names(Q_df_all) <- names(lines_sim_all)
for(j in 1:length(Q_df_all)){
   Q_df_all[[j]] <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)   
   for (i in 1:myr){
      Q_df_all[[j]]$mean[i] <- mean(dlist_all[[j]][[i]])
      Q_df_all[[j]]$qlow[i] <- quantile(dlist_all[[j]][[i]], probs=c(prob_qlow, prob_qup))[1]
      Q_df_all[[j]]$qup[i] <- quantile(dlist_all[[j]][[i]], probs=c(prob_qlow, prob_qup))[2]
   }}


mean_line_all <- vector("list", length(lines_sim_all))
qup_line_all <- vector("list", length(lines_sim_all))
qlow_line_all <- vector("list", length(lines_sim_all))
names(mean_line_all) <- names(qup_line_all) <- names(qlow_line_all) <- names(lines_sim_all)
for (i in 1:length(mean_line_all)) {
   mean_line_all[[i]]   <- data.frame(spline(Q_df_all[[i]]$Ma-1, Q_df_all[[i]]$mean))   
   qup_line_all[[i]] <- data.frame(spline(Q_df_all[[i]]$Ma-1, Q_df_all[[i]]$qup))
   qlow_line_all[[i]] <- data.frame(spline(Q_df_all[[i]]$Ma-1, Q_df_all[[i]]$qlow))
   colnames(mean_line_all[[i]]) <- colnames(qup_line_all[[i]]) <- colnames(qlow_line_all[[i]]) <- 
      c("time", "events")
}


# Plot Everything ----
for (i in 1:length(lines_sim_all)) {
  
  
  
   # Plot all events ----
   if (names(lines_sim_all[i]) == "event_all_cons") {
   plot_all <- ggplot() +
  
  # simulations
     geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
               mapping = aes(x=time, y=events, group=sim), 
               color=color_sim, alpha=alpha_sim, size=lwd_sim) +
     # Mean and quantiles
      geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
      geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
      geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
 #     geom_line(climate_data_line, mapping=aes(x=Time, y=Ts), color='red', size=1) +
      # Add orogeny events
      geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
      geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
      geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
      #geom_text(aes(x = 36, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
      #geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
      #geom_text(aes(x = 18.3, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
      #geom_text(aes(x = 16.8, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
      #geom_text(aes(x = 7.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
      #geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
      # Add climate aridification
      geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
      geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N + 7), angle = 90), label = "Arabia aridification", size = 3.5) +
      
      
 #observed events
  geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_all_events,
             size = lwd_events) + 
  
  xlim(30,0) +
  labs(x = "Time before present (Mya)", y = "Number of events") +
  ggtitle("All biogeographic events") +
    # Insert geologic scale
  coord_geo(xlim = c(30, 0), ylim = c(0,38), pos = as.list(rep("bottom", 2)),
            dat = list(epochs_htc, periods_htc),
            height = list(unit(1, "lines"), unit(1, "line")),
            rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
            center_end_labels = T, 
            skip = c('Holocene'), 
            lab = TRUE) +
    # Set the theme 
  theme_htc()
  ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_all.pdf")
}

   # ALL Vicariances  ----
if (names(lines_sim_all[i]) == "vicariance_all_cons") {
  
   plot_All_vic <- ggplot() +
  
   # simulations
   geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
             mapping = aes(x=time, y=events, group=sim), 
             color=color_sim, alpha=alpha_sim, size=lwd_sim) +
   
      # Mean and quantiles
      geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
      geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
      geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
      # Add orogeny events
      geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
      geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
      geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
      geom_text(aes(x = 36, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
      geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
      geom_text(aes(x = 18.3, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
      geom_text(aes(x = 16.8, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
      geom_text(aes(x = 7.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
      geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
       # Add climate aridification
      geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
      geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
      
      #observed events
   geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_all_vicariances,
             size = lwd_events) +
      xlim(30,0) +
   labs(x = "Time before present (Mya)", y = "Number of events") +
   ggtitle("All Vicariance Events") +
   # Insert geologic scale
   coord_geo(xlim = c(30, 0), ylim = c(0,20), pos = as.list(rep("bottom", 2)),
             dat = list(epochs_htc, periods_htc),
             height = list(unit(1, "lines"), unit(1, "line")),
             rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
             skip = c('Holocene'), 
             center_end_labels = T, 
                          lab = TRUE) +
   # Set the theme 
   theme_htc() 
   # Save it
   ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_All_vic.pdf")
   }

   # Central Colonization  ----
      if (names(lines_sim_all[i]) == "Ccol_all_cons") {
      plot_Ccol <- ggplot() +
          # simulations
          geom_line(bind_rows(lines_sim_all[[i]], .id="sim"), 
                    mapping = aes(x=time, y=events, group=sim), 
                    color=color_sim, alpha=alpha_sim, size=lwd_sim) +
         # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
         #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_Ccol,
                    size = lwd_events) +
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("Central Colonization Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'),
                    center_end_labels = T, 
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_Ccol.pdf")
      }
   # East Colonization  ----
   if (names(lines_sim_all[i]) == "Ecol_all_cons") {
      plot_Ecol <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
         
           #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_Ecol,
                    size = lwd_events) +
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("East Colonization Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T, 
                    lab = TRUE) +
          # Set the theme 
          theme_htc()
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_Ecol.pdf")
      }
   # West Colonization  ----
   if (names(lines_sim_all[i]) == "Wcol_all_cons") {
      plot_Wcol <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
         
         #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_Wcol,
                    size = lwd_events) +
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("West Colonization Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T,
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_Wcol.pdf")
      }
   # Central Diversification  ----
   if (names(lines_sim_all[i]) == "Cdiv_all_cons") {
      plot_Cdiv <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
         #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_Cdiv,
                    size = lwd_events) +
          
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("Central Diversification Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T,
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_Cdiv.pdf")
      }
   # East Diversification  ----
   if (names(lines_sim_all[i]) == "Ediv_all_cons") {
      plot_Ediv <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
          #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_Ediv,
                    size = lwd_events) +
         
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("East Diversification Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T,
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_Ediv.pdf")
      }
   # West Diversification  ----
   if (names(lines_sim_all[i]) == "Wdiv_all_cons") {
      plot_Wdiv <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
            #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_Wdiv,
                    size = lwd_events) +
         
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("West Diversification Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T,
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_Wdiv.pdf")
      }
   # CE Vicariance  ----
   if (names(lines_sim_all[i]) == "vicarianceCE_all_cons") {
      plot_vic_CE <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
         
            #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_vicCE,
                    size = lwd_events) +
          
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("Central-East Vicariant Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T,
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_CEvic.pdf")
      }
   # CW Vicariance  ----
   if (names(lines_sim_all[i]) == "vicarianceCW_all_cons") {
      plot_vic_CW <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
          #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_vicCW,
                    size = lwd_events) +
         
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("Central-West Vicariant Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T,
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_CWvic.pdf")
      }
   # EW Vicariance  ----
   if (names(lines_sim_all[i]) == "vicarianceEW_all_cons") {
      plot_vic_EW <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        
         # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
         
          #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_vicWE,
                    size = lwd_events) +
          
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("East-West Vicariant Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T,
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_EWvic.pdf")
      }
   
   # IW Vicariance  ----
   if (names(lines_sim_all[i]) == "vicarianceIW_all_cons") {
      plot_vic_IW <- ggplot() +
          # simulations
        geom_line(bind_rows(lines_sim_max_min[[i]], .id="sim"), 
                  mapping = aes(x=time, y=events, group=sim), 
                  color=color_sim, alpha=alpha_sim, size=lwd_sim) +
        # Mean and quantiles
         geom_line(mean_line_all[[i]], mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
         geom_line(qup_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         geom_line(qlow_line_all[[i]], mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
         # Add orogeny events
         geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
         geom_vline(xintercept=c(5), linetype="dashed", colour = "lightskyblue3", lwd = 0.8) +
         geom_text(aes(x = 35.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Hajar Mountain's main", size = 3.5) +
         geom_text(aes(x = 34.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "uplift event", size = 3.5)   +
         geom_text(aes(x = 18, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Western Hajar's", size = 3.5)   +
         geom_text(aes(x = 17, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         geom_text(aes(x = 7, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Jabal Akhdar", size = 3.5)   +
         geom_text(aes(x = 6, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "secondary uplift", size = 3.5)   +
         # Add climate aridification
         geom_rect(aes(xmin=3, xmax=0.02, ymin=-Inf, ymax=Inf), fill="#e5ea88", alpha=0.3, inherit.aes = FALSE) +
         geom_text(aes(x = 1.5, y = max(event_all_cons[[i]]$N +2), angle = 90), label = "Arabia aridification", size = 3.5) +
         
         
         #observed events
          geom_line(data = line_all[[i]], aes(x = time, y = events), color = color_vicIW,
                    size = lwd_events) +
          
          xlim(30,0) +
          labs(x = "Time before present (Mya)", y = "Number of events") +
          ggtitle("Iran-West Vicariant Events") +
          # Insert geologic scale
          coord_geo(xlim = c(30, 0), ylim = c(0,(max(event_all_cons[[i]]$N + 5))), pos = as.list(rep("bottom", 2)),
                    dat = list(epochs_htc, periods_htc),
                    height = list(unit(1, "lines"), unit(1, "line")),
                    rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
                    skip = c('Holocene'), center_end_labels = T,
                    lab = TRUE) +
          # Set the theme 
          theme_htc() 
          # Save it
          ggsave(width = 7, height = 5, units = 'in',"Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots/plot_IWvic.pdf")
      }
   
}
