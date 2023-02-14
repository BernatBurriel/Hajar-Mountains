# PLOT CUMULATIVE NUMBER OF EVENTS (EMPIRICAL AND SIMULATED)
setwd("~/Dropbox/AFRICA_ARABIA_H/BioGeoBEARS_AfAr/NEW2")


sim_event_all_cons <- readRDS("objects/sim_event_all_cons.rds")
event_all_cons <- readRDS("objects/event_all_cons.rds")

nsim <- 1000
myr <- 60

# Import empirical
event_list_cons <- readRDS("objects/event_list_cons.rds")
event_df <- data.frame(node=0, genus=NA, height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=0)

for (i in 1:length(event_list_cons)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      event_df <- rbind(event_df, event_list_cons[[i]][j,])
    }
  }
}

event_df <- event_df[-1,]
head(event_df)

#sort(event_df$min)
event_df_sorted <- event_df[rev(order(event_df$max)),]
xx <- event_df_sorted

# Import simulations
sim_event_cons <- readRDS("objects/sim_event_cons.rds")
sim_event_df <- vector("list", length(nsim))
for (i in 1:nsim){
  sim_event_df[[i]] <- data.frame(node=0, genus=NA, height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=0)
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

time_vec <- seq(from=60, to=1, by=-1)

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

saveRDS(sim_rmat, "objects/sim_rmat.rds")
saveRDS(rmat, "objects/rmat.rds")

## PLOT CUMULATIVE NUMBER OF EVENTS (OBSERVED AND SIM) ----
sim_rmat <- readRDS("objects/sim_rmat.rds")
rmat <- readRDS("objects/rmat.rds")

# Set colors for observed and simulated lines ----
color_all <- "#EE6A50"
color_sim <- "gray93"


pdf("plots/cumulative_events_OK.pdf", paper="a4", height=20, width=10)
plot(1,type='n',xlim=c(myr,0),ylim=c(0,80),xlab='Ma', ylab='Cumulative N', main="Cumulative biogeographic events")
for (s in 1:nsim){
  lines(spline(sim_rmat[[s]][, "Ma"]-1, sim_rmat[[s]][, "Ncum"]), type="l", col=color_sim, lwd=0.25, pch=16, cex=0.7)
}
lines(spline(rmat[, "Ma"]-1, rmat[, "Ncum"]), type="l", col=color_all, lwd=3)
axis(side=4)
dev.off()





