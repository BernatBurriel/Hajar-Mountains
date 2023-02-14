

rm(list = ls())
libs <- c("tidyverse", "deeptime", "here", "cowplot")
lapply(libs, require, character.only = TRUE)

events <- read.csv2("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/Mountain_colonization_Simmaps/Geologic_regions/lowland_midland_highland/elevation_events/elevation_events.csv")

myr <- 30
ma <- myr

upwards_elevation_plots <- data.frame(M = c(1:myr), N = 0)
downward_migration <- data.frame(M = c(1:myr), N = 0)

# We will create a list of vectors representing intervals of 1 Ma.
ma <- vector("list", ma)
for (i in 1:length(ma)){
  ma[[i]] <- c(i-1, i)
}
for (i in 1:NROW(events)){
    if (events$upwards[i] > 0){
      vi <- as.numeric(c(events$min[i], events$max[i]))
      for (k in 1:length(ma)){
        if (Overlap(na.omit(vi), ma[[k]]) != 0){
          upwards_elevation_plots[k,]$N <- upwards_elevation_plots[k,]$N + 1
        }}}
  if (events$downwards[i] > 0){
    vi <- as.numeric(c(events$min[i], events$max[i]))
    for (k in 1:length(ma)){
      if (Overlap(na.omit(vi), ma[[k]]) != 0){
        downward_migration[k,]$N  <- downward_migration[k,]$N + 1
        }}}
}


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

# Prepare the lines

  line_up <- data.frame(spline(upwards_elevation_plots$M-1, upwards_elevation_plots$N))
  line_dw <- data.frame(spline(downward_migration$M-1, downward_migration$N))
  colnames(line_up) <-  colnames(line_dw)<- c("time", "events")

saveRDS(upwards_elevation_plots, "/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/paleoclimate/data/upwards_elevation_plots.rds")

p <- ggplot() +
  
  #geom_rect(aes(xmin=40, xmax=30, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.15, inherit.aes = FALSE) +
  geom_rect(aes(xmin=20, xmax=15, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.40, inherit.aes = FALSE) +
geom_rect(aes(xmin=6, xmax=4, ymin=-Inf, ymax=Inf), fill="lightskyblue3", alpha=0.40, inherit.aes = FALSE) +
  geom_rect(aes(xmin=5.9, xmax=3.3, ymin=-Inf, ymax=Inf), fill="#e3ccb1", alpha=0.60, inherit.aes = FALSE) +
  geom_rect(aes(xmin=6.20, xmax=6.30, ymin=-Inf, ymax=Inf), fill="#e3ccb1", alpha=0.6, inherit.aes = FALSE) +
  geom_rect(aes(xmin=7.40, xmax=7.90, ymin=-Inf, ymax=Inf), fill="#e3ccb1", alpha=0.60, inherit.aes = FALSE) +
  geom_rect(aes(xmin=8.70, xmax=8.8, ymin=-Inf, ymax=Inf), fill="#e3ccb1", alpha=0.60, inherit.aes = FALSE) +
  #observed events
  
  geom_line(data = line_up, aes(x = time, y = events), color = '#e81785',
            size = 2) + 
  geom_line(data = line_dw, aes(x = time, y = events), color = '#35b44a',
            size = 2) + 
  xlim(20,0) +
  labs(x = "Time before present (Mya)", y = "Number of events") +
  ggtitle("Elevation shifts") +
  # Insert geologic scale
  coord_geo(xlim = c(20, 0), ylim = c(0,12), pos = as.list(rep("bottom", 2)),
            dat = list(epochs_htc, periods_htc),
            height = list(unit(1, "lines"), unit(1, "line")),
            rot = list(0, 0), size = list(3, 4), abbrv = list(TRUE, FALSE), 
            center_end_labels = T, 
            skip = c('Holocene'), 
            lab = TRUE) +
  
  theme_htc()
ggsave(width = 7, height = 5, units = 'in',"/Volumes/DROPBOX/Dropbox/Hajar_Mountains/figures/suplementary/Elevation.pdf")
