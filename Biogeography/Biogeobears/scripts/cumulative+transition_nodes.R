
# PLOT CUMULATIVE NUMBER OF EVENTS (EMPIRICAL AND SIMULATED)
#setwd("C:/Users/User/Desktop/__MACOSX/")
setwd("/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/")
rm(list=ls())
sim_event_all_cons <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_SIM_plots_all_no_diversification.rds")
event_all_cons <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_all_biogeo_no_diversification.rds")

nsim <- 1000
myr <- 80

climate_data_line <- read.table('paleoclimate/data/global_temperature_50kyr.txt', header = T)
climate_data_point <- read.csv2('paleoclimate/paleoclimatic_temperature_and_sea_level.csv')


# Import empirical
event_list_cons <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/event_list_cons_genera_biogeo_no_diversification.rds")
event_df <- data.frame(node=0, genus=0,
                       height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                       vicariance_EW = 0, vicariance_IW = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)

for (i in 1:length(event_list_cons)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      event_df <- rbind(event_df, event_list_cons[[i]][j,])
    }
  }
}

event_df <- event_df[-1,]

# take out diversification information 
event_df$Wdiv = event_df$Cdiv = event_df$Ediv = 0
nrow(event_df)
event_df <- event_df[rowSums(event_df[,6:length(event_df)]) > 0,]
nrow(event_df)
head(event_df)

#sort(event_df$min)
event_df_sorted <- event_df[rev(order(event_df$max)),]
xx <- event_df_sorted

# Import simulations
sim_event_cons <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/sim_event_consensus_no_diversification.rds")
sim_event_df <- vector("list", length(nsim))
for (i in 1:nsim){
  sim_event_df[[i]] <- data.frame(node=0, genus=0,
                                  height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                                  vicariance_EW = 0, vicariance_IW = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)
}

head(sim_event_df)

for (i in 1:length(sim_event_cons)){
  for (j in 1:nsim){
    if (nrow(sim_event_cons[[i]][[j]]) > 0){
      for (r in 1:nrow(sim_event_cons[[i]][[j]])){
        sim_event_df[[j]] <- rbind(sim_event_df[[j]], sim_event_cons[[i]][[j]][r,])
      }
    }
  }
}

for (i in 1:length(sim_event_df)){
  sim_event_df[[i]] <- sim_event_df[[i]][-1,]
  sim_event_df[[i]]$Cdiv = sim_event_df[[i]]$Ediv = sim_event_df[[i]]$Wdiv = 0
  sim_event_df[[i]] = sim_event_df[[i]][rowSums(sim_event_df[[i]][,6:length(sim_event_df[[i]])]) > 0,]
}

# sort each dataframe
event_df_sorted <- event_df[rev(order(event_df$max)),]

sim_event_df_sorted <- vector("list", nsim)
for (i in 1:nsim){
  sim_event_df_sorted[[i]] <- sim_event_df[[i]][rev(order(sim_event_df[[i]]$max)),]
}

head(sim_event_df_sorted)
sim_xx <- sim_event_df_sorted


##

time_vec <- seq(from=50, to=0, by=-1)

rmat <- matrix(NA, ncol = 2, nrow=length(time_vec))
colnames(rmat) <- c("Ma", "Ncum")
rmat[, "Ma"] <- time_vec

sim_rmat <- vector("list", nsim)
for (i in 1:nsim){
  sim_rmat[[i]] <- matrix(NA, ncol=2, nrow=length(time_vec))
  colnames(sim_rmat[[i]]) <- c("Ma", "Ncum")
  sim_rmat[[i]][, "Ma"] <- time_vec
}


nrow(xx[xx$max >= time_vec[20], ])
for (ii in 1:length(time_vec)){
  
  tmp <- xx[xx$max >= time_vec[ii], ]
  rmat[ii, "Ncum"] <- nrow(tmp)
}

nrow(sim_xx[[1]][sim_xx[[1]]$max >= time_vec[20], ])
for (s in 1:nsim){
  for (ii in 1:length(time_vec)){
    tmp <- sim_xx[[s]][sim_xx[[s]]$max >= time_vec[ii], ]
    sim_rmat[[s]][ii, "Ncum"] <- nrow(tmp)
  }
}

head(sim_rmat)

saveRDS(sim_rmat, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/sim_rmat_biogeo_no_diversification.rds")
saveRDS(rmat, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/rmat_col_no_diversification.rds")

## PLOT CUMULATIVE NUMBER OF EVENTS (OBSERVED AND SIM) ----
sim_rmat <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/sim_rmat_biogeo_no_diversification.rds")
rmat <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/rmat_col_no_diversification.rds")

##### 95% CI CALCULATION #####
# quantiles and mean of each My.
prob_qup <- 0.975
prob_qlow <- 0.025
nsim <- 1000
myr <- 50

##### 95% CI ALL EVENTS #####
# Let's do a list where each element is a vector with the
# number of transitions per simulation for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist <- vector("list", myr+1)
for (i in 1:length(dlist)){
  dlist[[i]] <- vector("numeric", nsim)
}
#dlist[[1]][2]

for (i in 1:nsim){
  for (j in 1:length(dlist)){
    dlist[[j]][i] <- as.data.frame(sim_rmat[[i]])$Ncum[j]
  }
}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df <- data.frame(Ma=c(myr:0), qlow=0, qup=0, mean=0)
for (i in 1:length(dlist)){
  Q_df$mean[i] <- mean(dlist[[i]])
  Q_df$qlow[i] <- quantile(dlist[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df$qup[i] <- quantile(dlist[[i]], probs=c(prob_qlow, prob_qup))[2]
}

mean_line <- spline(Q_df$Ma, Q_df$mean)
qup_line <- spline(Q_df$Ma, Q_df$qup)
qlow_line <- spline(Q_df$Ma, Q_df$qlow)


# Set colors for observed and simulated lines ----
color_all <- "#EE6A50"
color_sim <- "gray93"


pdf("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots_noJ/cumulative_events_OK.pdf", paper="a4", height=20, width=10)
plot(1,type='n',xlim=c((myr-10),0),ylim=c(0,350),xlab='Ma', ylab='Cumulative N', main="Cumulative biogeographic events")
for (s in 1:nsim){
  lines(spline(sim_rmat[[s]][, "Ma"]-1, sim_rmat[[s]][, "Ncum"]), type="l", col=color_sim, lwd=0.25, pch=16, cex=0.7)
}
lines(spline(rmat[, "Ma"]-1, rmat[, "Ncum"]), type="l", col=color_all, lwd=3)
lines(mean_line, type="l", col="black", lwd=1, pch=16, cex=1)
lines(qlow_line, type="l", col="black", lwd=2, pch=16, cex=1)
lines(qup_line, type="l", col="black", lwd=2, pch=16, cex=1)

axis(side=4)
axis(side = 1, at = seq(5,50, 10), labels = F, lwd.ticks = 0.5)
dev.off()


# Transition nodes ----
event_list_cons <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/event_list_cons_genera_biogeo_no_diversification.rds")
event_df <- data.frame(node=0, genus=0,
                       height=0, min=0, max=0, Ccol=0, Ecol=0, Wcol=0, vicariance_CE=0, vicariance_CW=0,
                       vicariance_EW = 0, vicariance_IW = 0, Cdiv = 0, Ediv = 0, Wdiv = 0,colonization = 0, extirpation = 0)

for (i in 1:length(event_list_cons)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      event_df <- rbind(event_df, event_list_cons[[i]][j,])
    }
  }
}

event_df <- event_df[-1,]
head(event_df)

myr <- 50

# Set color
color_all <- "#EE6A50"

# Plot all events in same color ----
pdf("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots_noJ/event_nodes.pdf", paper="a4r")
plot(1,type='n',xlim=c(myr,0),ylim=c(0,nrow(event_df)+5),xlab='Ma', ylab= "genus", main="Transition nodes", yaxt="n")
axis(2, at=c(nrow(event_df):1), labels=event_df$genus, col.axis="black", las=1, cex.axis=0.5)
j <- nrow(event_df)
for (i in 1:nrow(event_df)){
  segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_all, lwd=3)
  j <- j-1
}
dev.off()


# Define colors for different types of event ----
{color_Ccol <- "#b22125"
color_Ecol <- "#fec250"
color_Wcol <- "#206ab4"
color_Cdiv <- "salmon"
color_Ediv <- "gold"
color_Wdiv <- "lightblue"
color_vicCE <- "#ef8354"
color_vicCW <- "#33F1FF"
color_vicWE <- "#CCF5AC"
color_vicIW <- "#bbbbbb"
color_first_col <- "forestgreen"
}

#event_cols <- c(color_Ccol, color_Ecol, color_Wcol, color_Cdiv, color_Ediv, color_Wdiv, color_vicCE,
 #               color_vicCW, color_vicWE, color_vicIW, color_first_col)

event_cols <- c(color_Ccol, color_Ecol, color_Wcol, color_vicCE,
                color_vicCW, color_vicWE, color_vicIW)

#names(event_cols) <- c("Ccol", "Ecol", "Wcol", "Cdiv", "Ediv", "Wdiv", "vicCE", "vicCW", "vicEW", "vic_IW", "first_col")
names(event_cols) <- c("Ccol", "Ecol", "Wcol", "vicCE", "vicCW", "vicEW", "vic_IW")

plot(c(1:11), c(1:11), col=event_cols, pch=16, cex=2)
?segments




# Plot with different colors for different types of event ----
pdf("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots_noJ/event_nodes_colors_def.pdf", paper="a4", height=20, width=10)
plot(1,type='n',xlim=c(myr-10,0),ylim=c(0,nrow(event_df)+5),xlab='Ma', ylab= "genus", main="Transition nodes", yaxt="n")
axis(2, at=c(nrow(event_df):1), labels=event_df$genus, col.axis="black", las=2, cex.axis=0.5)
j <- nrow(event_df[rowSums(event_df[,6:length(event_df)] ) ])
for (i in 1:nrow(event_df)){
  if (event_df[i,]$Ccol==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Ccol, lwd=3)
  }
  if (event_df[i,]$Ecol==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Ecol, lwd=3)
  }
  if (event_df[i,]$Wcol==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Wcol, lwd=3)
  }
  if (event_df[i,]$Cdiv==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Cdiv, lwd=3)
  }
  if (event_df[i,]$Ediv==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Ediv, lwd=3)
  }
  if (event_df[i,]$Wdiv==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Wdiv, lwd=3)
  }
  if (event_df[i,]$vicariance_CE==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_vicCE, lwd=3)
  }
  if (event_df[i,]$vicariance_CW==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_vicCW, lwd=3)
  }
  if (event_df[i,]$vicariance_EW==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_vicWE, lwd=3)
  }
  if (event_df[i,]$vicariance_IW==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_vicIW, lwd=3)
  }
  if (event_df[i,]$colonization==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_first_col, lwd=3)
  }
  j <- j-1
}

segments(x0=68, x1=65, y0=seq(61, 71, 1), y1=seq(61, 71, 1), col=event_cols, lwd=5)
text(x=65, y=seq(61, 71, 1), labels = names(event_cols), cex=0.5, pos = 4)


dev.off()

saveRDS(event_df, "Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_df_no_diversification.rds")
event_df <- readRDS("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_df_no_diversification.rds")
colSums(event_df[,6:16])
sum(colSums(event_df[,6:16]))


#Merge all clades + species for plotting 
event_df_trimmed <- event_df
event_df_trimmed$Cdiv <- event_df_trimmed$Ediv <- event_df_trimmed$Wdiv <- 0
event_df_trimmed <- event_df_trimmed[rowSums(event_df_trimmed[,6:16]) > 0,]
event_df_trimmed$clade = NA
for (i in 1:length(unique(event_df_trimmed$genus))){
  gen = unique(event_df_trimmed$genus)[i]
  if (gen == "Asaccus" | gen == "Prup" | gen == "Hemidactylus" | gen ==  'Pcele'| gen == 'Ptyodactylus' | gen ==  'Trachydactylus') {
    event_df_trimmed$clade[event_df_trimmed$genus == gen] = 'Gekkota'
  }
  if (gen == "Echis"){
    event_df_trimmed$clade[event_df_trimmed$genus == gen] = 'Serpentes'
  }
  if (gen == "Omanosaura"){
    event_df_trimmed$clade[event_df_trimmed$genus == gen] = 'Lacertoidea'
  }}
event_df_trimmed = arrange(event_df_trimmed, clade)
rownames(event_df_trimmed) = c(1:nrow(event_df_trimmed))

# Plot cumulative in the background and the nodes on top (base R) ----

{pdf("Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/plots_noJ/cumulative_nodes.pdf", paper="a4", height=20, width=10)

# Plot cumulative number ----
par(mar=c(5,5,2,5))
plot(1,type='s',xlim=c((myr+1)-20,0),ylim=c(0,nrow(event_df_trimmed)+6),xlab='Mya', yaxt="n", ylab = "",
      yaxt="n", 
     frame = F)
axis(side=4, pos = -2, at = seq(35,0,-5), las = 1)
text(x = par('usr')[3]-5, y = par('usr')[4]/2, "Cumulative number of biogeographic events", srt = -90, xpd = NA)
mtext("Biogeographic events through time", side = 3, padj = 2, xpd = NA, cex = 1.2, font = 2)


# plot axis 2
position_species = c()
position_clades = c()
ablines <- c()
y0seg <- c()
y1seg <- c()
event_df_trimmed$genus = gsub('Pcele', 'Pristurus celerrimus', event_df_trimmed$genus)
event_df_trimmed$genus = gsub('Prup', 'Pristurus rupestris', event_df_trimmed$genus)

for (i in 1:length(unique(event_df_trimmed$genus))) {
  position_species = c(position_species, -1+ mean(as.numeric(rownames(event_df_trimmed[event_df_trimmed$genus == unique(event_df_trimmed$genus)[i],]))))
  ablines = c(ablines, nrow(event_df_trimmed) +1.5 - min(as.numeric(rownames(event_df_trimmed[event_df_trimmed$genus == unique(event_df_trimmed$genus)[i],]))))
}
for (i in 1:length(unique(event_df_trimmed$clade))) {
position_clades = c(position_clades, -1 + mean(as.numeric(rownames(event_df_trimmed[event_df_trimmed$clade == unique(event_df_trimmed$clade)[i],]))))
y0seg <- c(y0seg, nrow(event_df_trimmed)  +1.25 - min(as.numeric(rownames(event_df_trimmed[event_df_trimmed$clade == unique(event_df_trimmed$clade)[i],]))))
y1seg <- c(y1seg, nrow(event_df_trimmed) +.75 - max(as.numeric(rownames(event_df_trimmed[event_df_trimmed$clade == unique(event_df_trimmed$clade)[i],]))))
}
text(x = 32, y=nrow(event_df_trimmed)-(position_species), labels=unique(event_df_trimmed$genus),font = 3, cex=0.8, pos = 4)
text(y=nrow(event_df_trimmed)-(position_clades), x = (par("usr")[1] + .75), labels=unique(event_df_trimmed$clade),font = 1, cex=0.8, pos = 2, srt = 45, xpd = NA)
abline(h = ablines, col = 'grey80', lty = 1, lwd = 1)

segments(x0=par('usr')[1] + 0.5, x1=par('usr')[1] + 0.5, y0=y0seg, y1=y1seg, col='grey75', lwd=5, xpd = NA)



for (s in 1:nsim){
 lines(spline(sim_rmat[[s]][, "Ma"], sim_rmat[[s]][, "Ncum"]), type="l", col=color_sim, lwd=0.25, pch=16, cex=0.7)
}
rect(xleft= 20, xright = 15 , ytop =31.5, ybottom =  -3, col= rgb(141, 182, 205,alpha = 60, maxColorValue = 255), border = NA)
rect(xleft= 3, xright = 0.1 , ytop =31.5, ybottom =  -3, col= rgb(229, 224, 136, alpha = 60, maxColorValue = 255), border = NA)
lines(x = c(5,5),y= c(31.5,0), lwd=2, lty = 2, col = rgb(141, 182, 205,alpha = 60, maxColorValue = 255))

lines(spline(rmat[, "Ma"], rmat[, "Ncum"]), type="l", col=color_all, lwd=2)
lines(mean_line, type="l", col="black", lwd=1, pch=16, cex=1)
lines(qlow_line, type="l", col="black", lwd=0.5, pch=16, cex=1)
lines(qup_line, type="l", col="black", lwd=0.5, pch=16, cex=1)


# Plot event nodes with different colors



j <- nrow(event_df_trimmed)
for (i in 1:nrow(event_df_trimmed)){
  if (event_df_trimmed[i,]$Ccol==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_Ccol, lwd=3)
  }
  if (event_df_trimmed[i,]$Ecol==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_Ecol, lwd=3)
  }
  if (event_df_trimmed[i,]$Wcol==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_Wcol, lwd=3)
  }
  if (event_df_trimmed[i,]$Cdiv==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_Cdiv, lwd=3)
  }
  if (event_df_trimmed[i,]$Ediv==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_Ediv, lwd=3)
  }
  if (event_df_trimmed[i,]$Wdiv==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_Wdiv, lwd=3)
  }
  if (event_df_trimmed[i,]$vicariance_CE==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_vicCE, lwd=3)
  }
  if (event_df_trimmed[i,]$vicariance_CW==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_vicCW, lwd=3)
  }
  if (event_df_trimmed[i,]$vicariance_EW==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_vicWE, lwd=3)
  }
  if (event_df_trimmed[i,]$vicariance_IW==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_vicIW, lwd=3)
  }
  if (event_df_trimmed[i,]$colonization==1){
    segments(x0=event_df_trimmed[i,]$min, y0=j, x1=event_df_trimmed[i,]$max, y1=j, col=color_first_col, lwd=3)
  }
  j <- j-1
}


segments(x0=31, x1=32, y0=seq(35,33, -1), y1=seq(35,33, -1), col=event_cols[1:3], lwd=5)
text(x=31, y=seq(35,33, -1), labels = c('Central colonization', 'East colonization', 'West colonization'), cex=0.8, pos = 4)

segments(x0=20, x1=21, y0=seq(35,33, -1), y1=seq(35,33, -1), col=event_cols[4:6], lwd=5)
text(x=20, y=seq(35,33, -1), labels = c('Central-East vicariance', 'Central-West vicariance', 'East-West vicariance'), cex=0.8, pos = 4)

segments(x0=8, x1=9, y0=seq(35,35, -1), y1=seq(35,35, -1), col=event_cols[7], lwd=5)
text(x=8, y=seq(35,35, -1), labels = c('Iran-West vicariance'), cex=0.8, pos = 4)

#abline(h = 0, col = "black", lty = 2)

dev.off()
}





