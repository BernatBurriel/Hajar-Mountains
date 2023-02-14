### Testing correlatinos with Biogeography, Vicariance and climate

# set working directory
setwd('/Volumes/DROPBOX/Dropbox/BERNAT/PHD/01.HAJAR_MOUNTAINS/COLONIZATION/03.Biogeography/')

library(car);library(dplyr);library(nlme);library(RRPP)

# import the data
kyr50_df <- read.table('paleoclimate/data/global_temperature_50kyr.txt', header = T)  # paleoclimate temperatures every 50 Kyr
kyr1000_df <- read.table('paleoclimate/data/global_temperature_1000kyr.txt', header = T)  # paleoclimate temperatures every 1000 Kyr
T21myr <- kyr1000_df[kyr1000_df$Time <= 21,]
T24myr <- kyr1000_df[kyr1000_df$Time <= 24,]


## Import biogeographic data
biogeo <- readRDS('Mountain_colonization_area_Biogeobears/3_mountain_blocks_2state/objects/cumulative_plots/event_plots_all_biogeo_no_diversification.rds')
all_biogeo <- biogeo$event_all_cons[1:length(T24myr$Time),] # all biogeographic events


## Build a data.frame with number of events and time
biogeo_vic <- biogeo$vicariance_all_cons[1:length(T21myr$Time),] #all vicariant events
names(biogeo_vic) <- c("Time", "N")
dt_vic <- left_join(T21myr, biogeo_vic)
dt_vic <- dt_vic[,c("N","Ts","Time")]
dt_vic = dt_vic[-1,]


#Here we account for temporal autocorrelation structure in our linear model. Two approaches
 # 2 approaches: REML-based gls (nlme) & RRPP

library(nlme); library(RRPP)

# analyses

#1: REML-based gls
# Fit AR(1) model
fit.gls = gls(N ~ Ts, data = dt_vic,
                            correlation = corAR1(form = ~Time))
res <- summary(fit.gls)
res
res$coefficients

#2: RRPP
phi = 0.8871927  #obtained from REML

## correlation structure
ar1 <- corAR1(form = ~dt_vic$Time, value = phi)
## initialize this constructor against our data
AR1 <- Initialize(ar1, data = data.frame(dt_vic$Time))
## generate a correlation matrix
V <- corMatrix(AR1) 

rownames(V) = colnames(V) = rownames(dt_vic)
dt_vic$V=V

# Permutation estimation
fit.rrpp <- lm.rrpp(N ~ Ts data = dt_vic, Cov = V)

anova(fit.rrpp)

##coefficients (identical btwn methods)
fit.rrpp$LM$gls.coefficients
res$coefficients

plot(dt_vic$Ts, dt_vic$N)
abline(fit.rrpp)

