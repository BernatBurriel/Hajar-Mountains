##### Plot transitions nodes in consensus #####

setwd("~/Dropbox/AFRICA_ARABIA_H/BioGeoBEARS_AfAr/NEW2")

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

myr <- 60

# Set color
color_all <- "#EE6A50"

# Plot all events in same color ----
pdf("plots/event_nodes.pdf", paper="a4r")
plot(1,type='n',xlim=c(myr,0),ylim=c(0,nrow(event_df)+5),xlab='Ma', ylab= "genus", main="Transition nodes", yaxt="n")
axis(2, at=c(nrow(event_df):1), labels=event_df$genus, col.axis="black", las=1, cex.axis=0.5)
j <- nrow(event_df)
for (i in 1:nrow(event_df)){
  segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_all, lwd=3)
  j <- j-1
}
dev.off()


# Define colors for different types of event ----
color_Af2Ar <- "#A8CB66"
color_Ar2Af <- "#2E8B57"
color_Vic <- "#FFC125"
color_ExtAf <- "#27408B"
color_ExtAr <- "#5CACEE"
event_cols <- c(color_Af2Ar, color_Ar2Af, color_Vic, color_ExtAf, color_ExtAr)
names(event_cols) <- c("Af2Ar", "Ar2Af", "Vic", "ExtAf", "ExtAr")

plot(c(1:5), c(1:5), col=event_cols, pch=16, cex=5)
?segments

# Plot with different colors for different types of event ----
pdf("plots/event_nodes_colors.pdf", paper="a4", height=20, width=10)
plot(1,type='n',xlim=c(myr,0),ylim=c(0,nrow(event_df)+5),xlab='Ma', ylab= "genus", main="Transition nodes", yaxt="n")
axis(2, at=c(nrow(event_df):1), labels=event_df$genus, col.axis="black", las=1, cex.axis=0.5)
j <- nrow(event_df)
for (i in 1:nrow(event_df)){
  if (event_df[i,]$Af2Ar==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Af2Ar, lwd=3)
  }
  if (event_df[i,]$Ar2Af==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Ar2Af, lwd=3)
  }
  if (event_df[i,]$Vic==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Vic, lwd=3)
  }
  if (event_df[i,]$ExtAf==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_ExtAf, lwd=3)
  }
  if (event_df[i,]$ExtAr==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_ExtAr, lwd=3)
  }
  j <- j-1
}

segments(x0=60, x1=55, y0=c(70,75,80,85,90)-10, y1=c(70, 75, 80, 85, 90)-10, col=event_cols, lwd=4)
text(x=55, y=c(70, 75, 80, 85, 90)-10, labels = names(event_cols), cex=0.7, pos = 4)

dev.off()

saveRDS(event_df, "objects/event_df.rds")

colSums(event_df[,6:10])



