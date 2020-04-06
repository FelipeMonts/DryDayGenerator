###############################################################################################################
#                             Tell the program where the package libraries are stored                        
###############################################################################################################


#  Tell the program where the package libraries are  #####################

.libPaths("C:/Felipe/SotwareANDCoding/R_Library/library")  ;


###############################################################################################################
#                             Setting up working directory  Loading Packages and Setting up working directory                        
###############################################################################################################


#      set the working directory

# readClipboard()
setwd("C:\\Felipe\\Students Projects\\MitchHunter\\DryDayGenerator\\Fw__Next_steps_on_AWC_drought_analysis") ;   # 



###############################################################################################################
#                            Install the packages that are needed                       
###############################################################################################################


###############################################################################################################
#                           load the libraries that are neded   
###############################################################################################################

###############################################################################################################################################################
# 
#  Using Sys.time, proc.time() and system.time() to check how much the code take.
#  
#   Tracking memory can be done with gc(), object.size(), memory.profile(), memory.size ()
# 
##############################################################################################################################################################
M1<-memory.size() ;
T1<-Sys.time() ;
L1<-41;

Time.Memory.use<-data.frame(Memory=M1, Time=T1) ;


ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages = c(# "agricolae", 
             "nlme", 
             # "car", 
             # "faraway", 
             "ggplot2", 
             "lattice", 
             "reshape2", 
             "multcomp", 
             # "vegan", 
             "lme4", 
             # "AICcmodavg", 
             # "pgirmess", 
             # "boot", 
             # "arm", 
             # "coefplot", 
             # "gridExtra", 
             # "devtools", 
             # "nls2", 
             # "piecewiseSEM", 
             "cowplot", 
             # "Rmisc", 
             "emmeans",        
             "multcompView", 
             "dplyr", 
             # "MASS", 
             "MuMIn",
             "randomForest")          

ipak(packages)  # Determine if packages are installed; install if needed; load all
T2<-Sys.time() ;
M2<-memory.size() ;
L2<-84 ;




# setwd("C:/Users/mhunter/Google Drive/UMN/Personal/AWC Analysis/") # Set the working directory; not needed when working in a Project

# rm(list=ls())  # remove all objects from workspace
# rm(rf.TM)  # remove specific objects
# rm(rf.EA)  # remove specific objects
# rm(rf.SOLAR)  # remove specific objects

# gc()  # clear unused memory

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","gray30")



#### Load data and calculate new variables ####

orig<-read.csv("WoosterOriginal.csv", header = T)
orig$YEAR = factor(orig$YEAR)

orig$TM = (orig$TX + orig$TN) / 2                                 # mean temperature
orig$DT = orig$TX - orig$TN                                       # difference between max and min temp
orig$ES.TX = .611 * exp((17.502 * orig$TX) / (orig$TX + 240.97))  # saturation vapor pressure at TX
orig$ES.TN = .611 * exp((17.502 * orig$TN) / (orig$TN + 240.97))  # saturation vapor pressure at TN
orig$EA.TX = orig$ES.TX * orig$RHN / 100                          # ambient vapor pressure at TX
orig$EA.TN = orig$ES.TN * orig$RHX / 100                          # ambient vapor pressure at TN
orig$DE = orig$EA.TX - orig$EA.TN                                 # difference between max and min vapor pressure
orig$EA = (orig$EA.TX + orig$EA.TN) / 2                           # mean ambient vapor pressure



#### ~~ Calculate PSOLAR - Potential Solar Radiation ####

lat = 40.772289       # direct input latitude
lat = lat * pi / 180  # convert degrees to radians

# Define function for calculating potential solar radiation give DOY and latitude
pot.solar <- function(DOY, lat){  # both numeric, lat in radians
  d = DOY
  sindec = 0.39785 * sin(4.869 + 0.0172 * d + 0.03345 * sin(6.224 + 0.0172 * d))
  D2 = 1 + 0.0334 * cos(0.01721 * d - 0.0552)
  cosdec = sqrt(1 - sindec * sindec)
  SinLatSinDec = sin(lat) * sindec    # lat must be defined outside of function
  CosLatCosDec = cos(lat) * cosdec
  coshs = -SinLatSinDec / CosLatCosDec

  if(coshs < -1){coshs = -0.9999}
  
  sinhs = sqrt(1 - coshs * coshs)
  hs = pi / 2 - atan(coshs / sinhs)
  potentialSolarGS = 117.5 * D2 * (hs * SinLatSinDec + CosLatCosDec * sinhs) / pi
  
  if(coshs > 1){potentialSolarGS = 0}
  potentialSolarGS
}


## Apply pot.solar function to calculate PSOLAR
orig = orig %>% 
  rowwise() %>% 
  mutate(PSOLAR = pot.solar(DOY, lat))  # create new column of PSOLAR dat using pot.solar function
orig = as.data.frame(orig)

plot(orig$PSOLAR[1:365])  # plot daily PSOLAR values

T3<-Sys.time() ;
M3<-memory.size() ;
L3<-154 ;


#### Fit random forest models ####

# Model for predicting mean temperature

rf.TM = randomForest(formula = TM ~ EA + PP + PSOLAR + DOY + WIND, 
                     data = orig, ntree = 300, mtry = 2, 
                     importance = T, proximity = T) 

orig$TM.pred = predict(rf.TM, orig)        # predicted values from original file
orig$TM.resid = orig$TM - orig$TM.pred     # residuals


partialPlot(rf.TM, orig, EA)
partialPlot(rf.TM, orig, PP)

T4<-Sys.time() ;
M4<-memory.size() ;
L4<-174;


# Model for predicting mean vapor pressure

rf.EA = randomForest(formula = EA ~ TM + PP + PSOLAR + DOY + WIND, 
                     data = orig, ntree = 300, mtry = 2, 
                     importance = T, proximity = T) 

orig$EA.pred = predict(rf.EA, orig)        # predicted values from original file
orig$EA.resid = orig$EA - orig$EA.pred     # residuals


partialPlot(rf.EA, orig, TM)
partialPlot(rf.EA, orig, PP)


T5<-Sys.time() ;
M5<-memory.size() ;
L5<-193;

# Model for predicting solar radiation

rf.SOLAR = randomForest(formula = SOLAR ~ TX + TN + PP + RHN + RHX + PSOLAR + DOY + WIND, 
                        data = orig, ntree = 300, mtry = 3, 
                        importance = T, proximity = T) 

orig$SOLAR.pred = predict(rf.SOLAR, orig)           # predicted values from original file
orig$SOLAR.resid = orig$SOLAR - orig$SOLAR.pred     # residuals


plot(orig$PP, orig$SOLAR.resid)  # check dependence of solar residuals on PP amount
# The model seems to under-predict low solar radiation events on dry days



T6<-Sys.time() ;
M6<-memory.size() ;
L6<-210;

#### Create zero-precip scenario ####
{
dry = orig
dry$PP = 0

dry$TM.pred = predict(rf.TM, dry)                    # predict TM for PP = 0
dry$TM = dry$TM.pred + orig$TM.resid                 # add residuals back to TM

dry$TX = dry$TM + orig$DT / 2                        # calc TX using original DT
dry$TN = dry$TM - orig$DT / 2                        # calc TN using original DT


dry$EA.pred = predict(rf.EA, dry)                    # predict EA for PP = 0
dry$EA = dry$EA.pred + orig$EA.resid                 # add residuals back to EA

dry$EA.TX = dry$EA + orig$DE / 2                     # calc EA.TX using original DE
dry$EA.TN = dry$EA - orig$DE / 2                     # calc EA.TN using original DE

dry$ES.TX = .611 * exp((17.502 * dry$TX) / (dry$TX + 240.97))  # saturation vapor pressure at TX
dry$ES.TN = .611 * exp((17.502 * dry$TN) / (dry$TN + 240.97))  # saturation vapor pressure at TX

dry$RHN = dry$EA.TX / dry$ES.TX * 100                # calc RHX
dry$RHX = dry$EA.TN / dry$ES.TN * 100                # calc RHN


dry$SOLAR.pred = predict(rf.SOLAR, dry)              # predict SOLAR for PP = 0
dry$SOLAR = dry$SOLAR.pred + orig$SOLAR.resid        # add residuals back to SOLAR
}


T7<-Sys.time() ;
M7<-memory.size() ;
L7<-244;




#### Functions for creating drier weather files ####

orig.test = orig[1:365,]  # To show functionality in examples below; for full implementation, use all values
dry.test = dry[1:365,]   # To show functionality in examples below; for full implementation, use all values

# sum(orig.test$PP)
# ggplot(orig.test, aes(x = PP)) + geom_density()  # plot distribution of PP values


## drier.threshold - Function to remove all PP events below a THRESHOLD of cm
drier.threshold = function(threshold, orig.func, dry.func){  
  new = orig.func[,1:18]                  # set up new dataframe
  to.dry = orig.func$PP < threshold       # identify days with PP > threshold
  new[to.dry,] = dry.func[to.dry, 1:18]   # swap in dry days
  new
}
# temp1 = drier.threshold(.5, orig.test, dry.test)  # remove events lower than .5 cm
# sum(temp1$PP) / sum(orig.test$PP)


## drier.quantile - Function to remove all PP events below a QUANTILE of cm
drier.quantile = function(quantile, orig.func, dry.func){  
  new = orig.func[,1:18]                                    # set up new dataframe
  to.dry = orig.func$PP < quantile(orig.func$PP, quantile)  # identify days with PP > quantile
  new[to.dry,] = dry.func[to.dry, 1:18]                     # swap in dry days
  new
}
# temp2 = drier.quantile(.5, orig.test, dry.test)   # remove values below the 50th percentile (.21 cm); percentile includes zeros
# sum(temp2$PP) / sum(orig.test$PP)


## drier.pct.low - Function to remove SMALL events until a specified PERCENT of total precip is removed
drier.pct.low = function(pct, orig.func, dry.func){                         # pct in decimal form
  new = orig.func[,1:18]                                                    # set up new dataframe
  tot.PP = sum(orig.func$PP)                                                # calc total PP
  PP.df = data.frame(index = 1:length(orig.func$PP), PP = orig.func$PP)     # set up DF with index and PP values
  PP.df = PP.df[order(PP.df$PP),]                                           # order DF by PP values
  PP.df$cum.PP = cumsum(PP.df$PP)                                           # calc cumulative sum of PP values
  cutoff = which(PP.df$cum.PP > tot.PP * pct)[1]                            # calc cutoff value for percentage of total PP
  to.dry = PP.df$index[1:cutoff-1]                                          # identify days to swap in based on index
  new[to.dry,] = dry.func[to.dry, 1:18]                                     # swap in dry days
  new
}
# temp3 = drier.pct.low(.3, orig.test, dry.test)   # remove 30% of precip on low end
# sum(temp3$PP) / sum(orig.test$PP)


## drier.pct.hi - Function to remove LARGE events until a specified PERCENT of total precip is removed
drier.pct.hi = function(pct, orig.func, dry.func){                    # pct in decimal form
  new = orig.func[,1:18]                                              # set up new dataframe
  tot.PP = sum(orig.func$PP)                                          # calc total PP
  PP.df = data.frame(index = 1:length(orig.func$PP), PP = orig.func$PP)       # set up DF with index and PP values
  PP.df = PP.df[order(PP.df$PP),]                                     # order DF by PP values
  PP.df$cum.PP = cumsum(PP.df$PP)                                     # calc cumulative sum of PP values
  cutoff = which(PP.df$cum.PP > tot.PP * (1-pct))[1]                  # calc cutoff value for percentage of total PP
  to.dry = PP.df$index[cutoff:nrow(PP.df)]                            # identify days to swap in based on index
  new[to.dry,] = dry.func[to.dry, 1:18]                               # swap in dry days
  new
}
# temp4 = drier.pct.hi(.3, orig.test, dry.test)   # remove 30% of precip on high end
# sum(temp4$PP) / sum(orig.test$PP)


## drier.random - Function to remove events RANDOMLY until a specified percent of total precip is removed
drier.random = function(pct, orig.func, dry.func){  # pct in decimal form
  new = orig.func[,1:18]                                                  # set up new dataframe
  tot.PP = sum(orig.func$PP)                                              # calc total PP
  PP.df = data.frame(index = 1:length(orig.func$PP), PP = orig.func$PP)   # set up DF with index and PP values
  set.seed(123)
  PP.df = PP.df[sample(nrow(PP.df)),]                                     # order DF randomly
  PP.df$cum.PP = cumsum(PP.df$PP)                                         # calc cumulative sum of PP values
  cutoff = which(PP.df$cum.PP > tot.PP * pct)[1]                          # calc cutoff value for percentage of total PP
  to.dry = PP.df$index[1:cutoff-1]                                        # identify days to swap in based on index
  new[to.dry,] = dry.func[to.dry, 1:18]                                   # swap in dry days
  new
}
# temp5 = drier.random(.3, orig.test, dry.test)   # remove 30% of precip from random events
# sum(temp5$PP) / sum(orig.test$PP)

T8<-Sys.time() ;
M8<-memory.size() ;
L8<-330;

#### Create drier weather files ####

# Output dataframes represent:
# - 30% or 50% reduction OVER ALL YEARS by removing SMALLER precip events
# - 30% or 50% reduction PER YEAR by removing SMALLER precip events
# - 30% or 50% reduction OVER ALL YEARS by removing LARGER precip events 
# - 30% or 50% reduction PER YEAR by removing LARGER precip events 
# - 30% or 50% reduction OVER ALL YEARS by removing RANDOM precip events
# - 30% or 50% reduction PER YEAR by removing RANDOM precip events


# 30% or 50% reduction OVER ALL YEARS by removing SMALLER precip events
dry.30.low = drier.pct.low(.3, orig, dry)
write.csv(dry.30.low, "dry.30.low", row.names = F)

dry.50.low = drier.pct.low(.5, orig, dry)
write.csv(dry.50.low, "dry.50.low", row.names = F)

# 30% or 50% reduction PER YEAR by removing SMALLER precip events
dry.30.low.annual = data.frame(matrix(nrow = 0, ncol = ncol(orig)))
colnames(dry.30.low.annual) = colnames(orig)
year.list = levels(orig$YEAR)
for(i in 1:length(year.list)){ 
  temp = drier.pct.low(.3, 
                       filter(orig, orig$YEAR == year.list[i]), 
                       filter(dry, dry$YEAR == year.list[i]))
  dry.30.low.annual = rbind(dry.30.low.annual, temp)
}
write.csv(dry.30.low.annual, "dry.30.low.annual", row.names = F)

dry.50.low.annual = data.frame(matrix(nrow = 0, ncol = ncol(orig)))
colnames(dry.50.low.annual) = colnames(orig)
year.list = levels(orig$YEAR)
for(i in 1:length(year.list)){ 
  temp = drier.pct.low(.5, 
                       filter(orig, orig$YEAR == year.list[i]), 
                       filter(dry, dry$YEAR == year.list[i]))
  dry.50.low.annual = rbind(dry.50.low.annual, temp)
}
write.csv(dry.50.low.annual, "dry.50.low.annual", row.names = F)


# - 30% or 50% reduction OVER ALL YEARS by removing LARGER precip events 
dry.30.hi = drier.pct.hi(.3, orig, dry)
write.csv(dry.30.hi, "dry.30.hi", row.names = F)

dry.50.hi = drier.pct.hi(.5, orig, dry)
write.csv(dry.50.hi, "dry.50.hi", row.names = F)


# - 30% or 50% reduction PER YEAR by removing LARGER precip events 
dry.30.hi.annual = data.frame(matrix(nrow = 0, ncol = ncol(orig)))
colnames(dry.30.hi.annual) = colnames(orig)
year.list = levels(orig$YEAR)
for(i in 1:length(year.list)){ 
  temp = drier.pct.hi(.3, 
                       filter(orig, orig$YEAR == year.list[i]), 
                       filter(dry, dry$YEAR == year.list[i]))
  dry.30.hi.annual = rbind(dry.30.hi.annual, temp)
}
write.csv(dry.30.hi.annual, "dry.30.hi.annual", row.names = F)

dry.50.hi.annual = data.frame(matrix(nrow = 0, ncol = ncol(orig)))
colnames(dry.50.hi.annual) = colnames(orig)
year.list = levels(orig$YEAR)
for(i in 1:length(year.list)){ 
  temp = drier.pct.hi(.5, 
                       filter(orig, orig$YEAR == year.list[i]), 
                       filter(dry, dry$YEAR == year.list[i]))
  dry.50.hi.annual = rbind(dry.50.hi.annual, temp)
}
write.csv(dry.50.hi.annual, "dry.50.hi.annual", row.names = F)


# - 30% or 50% reduction OVER ALL YEARS by removing RANDOM precip events
dry.30random = drier.random(.3, orig, dry)
write.csv(dry.30random, "dry.30random", row.names = F)

dry.50random = drier.random(.5, orig, dry)
write.csv(dry.50random, "dry.50random", row.names = F)


# - 30% or 50% reduction PER YEAR by removing RANDOM precip events
dry.30random.annual = data.frame(matrix(nrow = 0, ncol = ncol(orig)))
colnames(dry.30random.annual) = colnames(orig)
year.list = levels(orig$YEAR)
for(i in 1:length(year.list)){ 
  temp = drier.random(.3, 
                      filter(orig, orig$YEAR == year.list[i]), 
                      filter(dry, dry$YEAR == year.list[i]))
  dry.30random.annual = rbind(dry.30random.annual, temp)
}
write.csv(dry.30random.annual, "dry.30random.annual", row.names = F)

dry.50random.annual = data.frame(matrix(nrow = 0, ncol = ncol(orig)))
colnames(dry.50random.annual) = colnames(orig)
year.list = levels(orig$YEAR)
for(i in 1:length(year.list)){ 
  temp = drier.random(.5, 
                      filter(orig, orig$YEAR == year.list[i]), 
                      filter(dry, dry$YEAR == year.list[i]))
  dry.50random.annual = rbind(dry.50random.annual, temp)
}
write.csv(dry.50random.annual, "dry.50random.annual", row.names = F)

T9<-Sys.time() ;
M9<-memory.size() ;
L9<-439;


#### Create Monthly Versions of Dry Weather Files ####

dry.30.low


dry.30.low.mly = cbind(
  aggregate(dry.30.low["TN"], by = list(dry.30.low$yr, dry.30.low$month), FUN = mean, na.rm = T),
  aggregate(dry.30.low["TN"], by = list(dry.30.low$yr, dry.30.low$month), FUN = var, na.rm = T)[,3],
  aggregate(dry.30.low["TX"], by = list(dry.30.low$yr, dry.30.low$month), FUN = mean, na.rm = T)[,3],
  aggregate(dry.30.low["tTX"], by = list(dry.30.low$yr, dry.30.low$month), FUN = var, na.rm = T)[,3],
  aggregate(dry.30.low["ref.et"], by = list(dry.30.low$yr, dry.30.low$month), FUN = sum, na.rm = T)[,3],
  aggregate(dry.30.low["ref.et"], by = list(dry.30.low$yr, dry.30.low$month), FUN = var, na.rm = T)[,3],
  aggregate(dry.30.low["pp"], by = list(dry.30.low$yr, dry.30.low$month), FUN = sum, na.rm = T)[,3],
  aggregate(dry.30.low["pp"], by = list(dry.30.low$yr, dry.30.low$month), FUN = var, na.rm = T)[,3])
names(w.hist.mly) = c("yr", "month", "tn.mn", "tn.var", "tx.mn", "tx.var", "ref.et.sum", "ref.et.var", "pp.sum", "pp.var")
w.hist.mly = w.hist.mly[order(w.hist.mly$yr),]


Time.Memory.use<-data.frame(Memory_MB=c(M1,M2,M3,M4,M5,M6,M7,M8,M9), Time_s=c(T1,T2,T3,T4,T5,T6,T7,T8,T9), CodeLine=c(L1,L2,L3,L4,L5,L6,L7,L8,L9));

Time.Memory.use$Memory_consumed_MB<-c(0,diff(Time.Memory.use$Memory)) ;

Time.Memory.use$Time_used_s<-c(0,diff(Time.Memory.use$Time)) ;


write.csv(Time.Memory.use, file="Time.Memory.use.csv") ;


#### Calc SPEI for Dry Weather Files ####

w.hist.spei.weather = w.hist.mly$pp.sum - w.hist.mly$ref.et.sum
w.hist.spei.weather = ts(w.hist.spei.weather, end=c(2005,12), frequency=12)  # transform into time series
w.hist.spei.weather = as.matrix(cbind(w.hist.spei.weather, w.hist.spei.weather))  # create matrix consisting of two time series so that the result of fitted() will be a simple time series of results, rather than a matrix (bad functionality)

w.hist.spei.3 = spei(data = w.hist.spei.weather, scale = 3)  # create spei object with 3-month period
w.hist.spei.3.vals = as.data.frame(fitted(w.hist.spei.3))[,1]
w.hist.sumspei = w.hist.spei.3.vals[seq(8, 668, 12)]   # August 3-month SPEI values



w.spei.1 = spei(data = w.spei.weather, ref.start = c(1950, 1), ref.end = c(2005,12), scale = 1)  # create spei object with 1-month period







w.hist.spei = data.frame(cbind(w.hist.junspei, w.hist.julspei, w.hist.augspei, w.hist.jagspei, w.hist.sumspei))
names(w.hist.spei) = c("junspei", "julspei", "augspei", "jagspei", "sumspei")
write.csv(w.hist.spei, "w.hist.spei.csv", row.names = F)

T10<-Sys.time() ;
M10<-memory.size() ;



#### Compare original and dry.30pct.low values ####

# #### Precipitation
# plot(orig$PP[1:365], type = 'l', lwd = 1, main = "PP", xlab = "DOY in 1980", ylab = "cm")
# lines(dry.30pct.low$PP[1:365], col = "red", lwd = 1)
# 
# plot(dry.30pct.low$PP[1:365] - orig$PP[1:365], main = "Precipitation Difference (dry minus original)", xlab = "DOY in 1980", ylab = "cm")
# 
# 
# 
# #### Mean Temperature
# plot(orig$TM[1:365], type = 'l', lwd = 1, main = "Mean Temperature", xlab = "DOY in 1980", ylab = "C")
# lines(dry.30pct.low$TM[1:365], col = "red", lwd = 1)
# 
# plot(dry.30pct.low$TM[1:365] - orig$TM[1:365], main = "Mean Temperature Difference (dry minus original)", xlab = "DOY in 1980", ylab = "C")
# 
# 
# 
# #### Solar Radiation
# plot(orig$SOLAR[1:365], type = 'l', lwd = 1, main = "Solar Radiation", xlab = "DOY in 1980", ylab = "MJ/m2/d")
# lines(dry.30pct.low$SOLAR[1:365], col = "red", lwd = 1)
# 
# plot(dry.30pct.low$SOLAR[1:365] - orig$SOLAR[1:365], main = "Solar Radiation Difference (dry minus original)", xlab = "DOY in 1980", ylab = "MJ/m2/d")
# 
# # zoomed-in version
# plot(orig$SOLAR[150:200], type = 'l', lwd = 1, main = "Solar Radiation", xlab = "DOY in 1980", ylab = "MJ/m2/d")
# lines(dry.30pct.low$SOLAR[150:200], col = "red", lwd = 1)
# 
# plot(dry.30pct.low$SOLAR[150:200] - orig$SOLAR[150:200], main = "Solar Radiation Difference (dry minus original)", xlab = "DOY in 1980", ylab = "MJ/m2/d")
# 
# 
# 
# #### EA
# plot(orig$EA[1:365], type = 'l', lwd = 1, main = "EA", xlab = "DOY in 1980", ylab = "%")
# lines(dry.30pct.low$EA[1:365], col = "red", lwd = 1)
# 
# plot(dry.30pct.low$EA[1:365] - orig$EA[1:365], main = "EA Difference (dry minus original)", xlab = "DOY in 1980", ylab = "%")
# 
# 
# 
# #### RHX
# plot(orig$RHX[1:365], type = 'l', lwd = 1, main = "RHX", xlab = "DOY in 1980", ylab = "%")
# lines(dry.30pct.low$RHX[1:365], col = "red", lwd = 1)
# 
# plot(dry.30pct.low$RHX[1:365] - orig$RHX[1:365], main = "RHX Difference (dry minus original)", xlab = "DOY in 1980", ylab = "%")
# 
# 
# 
# #### RHN
# plot(orig$RHN[1:365], type = 'l', lwd = 1, main = "RHN", xlab = "DOY in 1980", ylab = "%")
# lines(dry.30pct.low$RHN[1:365], col = "red", lwd = 1)
# 
# plot(dry.30pct.low$RHN[1:365] - orig$RHN[1:365], main = "RHN Difference (dry minus original)", xlab = "DOY in 1980", ylab = "%")







