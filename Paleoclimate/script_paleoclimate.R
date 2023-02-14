setwd('/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/paleoclimate/')
rm(list = ls())
library(dplyr) ; library(sf); library(phytools); library(tidyr)

# 1. get the dataset with the delta 18_O data Extracted fro Hansen et al. (2013). 
climate_data <- read.csv2('paleoclimatic_temperature_and_sea_level.csv')

# Create a new dataset summarising the data with 500kyr windows for the final line 
kyr500_df <- data.frame(Ts = NA, Time= seq(0,65.5228, 0.5))
seq_kyr500 <- seq(0,65.5228, 0.5)


for (i in seq_kyr500) {
  row_name <- rownames(kyr500_df[kyr500_df$Time == i,])
  kyr500_df[as.numeric(row_name),1] <- mean(as.numeric(climate_data$Ts[as.numeric(climate_data$Time.1) > i & as.numeric(climate_data$Time.1) < i+0.3 ]))
}

pdf("results/temperature_through_time_65_3000.pdf", paper="a4r", height=20, width=40)
plot(1,type='s',xlim=c(65,0),ylim=c(min(as.numeric(climate_data$Ts))-5 , max(as.numeric(climate_data$Ts))+10),xlab='Mya', ylab = "Surface mean temperature",
     frame = F, yaxt = 'n')
axis(side=2, at = seq(34,4,-2.), las = 1)
#points(rev(climate_data$Ts) ~ rev(climate_data$Time.1), pch = 16, cex = 0.4, col = 'grey')
lines(rev(spline(kyr1500_df$Ts ~ kyr1500_df$Time)), type="l", col='red', lwd=1)
dev.off()

write.table(kyr50_df, 'data/global_temperature_50kyr.txt', row.names = F, col.names = T, quote = F)
