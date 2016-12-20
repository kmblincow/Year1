#Data analysis for Harvest Preference Project---KB 
#last edited 12/20/2016

setwd("C:/Users/Kayla/SkyDrive/PhDzNuts/Fall2016/Sugihara")

vesseldata<- read.table("harvestvd.txt",header=T)
combodata<-read.table("harvestd.txt",header=T)

#convert date values to R recognized date values
vesseldata$Date<-as.Date(vesseldata$Date, "%m/%d/%Y")
combodata$Date<-as.Date(combodata$Date, "%m/%d/%Y")

vesseldata<-vesseldata[order(vesseldata$Date, decreasing=FALSE),] 
combodata<-combodata[order(combodata$Date, decreasing=FALSE),] 
#orders data based on date

newdata<-data.frame(
  Date=seq(as.Date("2014-01-01"), as.Date("2016-12-31"),"days"),
  RFHarvest=NA,
  KBHarvest=NA)

#cross reference the dtaes and fill in data from original data set
idx<-match(combodata$Date,newdata$Date)
newdata[idx,c(2,3)]<-combodata[,c(2,3)]
newdata<-newdata[-which(newdata$Date==as.Date("2016-02-29")),]
na_idx <- is.na(newdata$KBHarvest)
newdata$KBHarvest[na_idx]<-NA

plot1<-plot(newdata$Date, newdata$RFHarvest, type="l",col="red", ylab="Fish Harvest", xlab="Date", xlim=c(as.Date("2014-03-01"),as.Date("2016-11-01")))
lines(newdata$Date,newdata$KBHarvest,col="blue")
legend("topleft", legend = c("Kelp Bass Harvest", "Rockfish Harvest"), 
       col = c("blue", "red"), lwd = 1)

library(rEDM)


#### univariate analysis ####

# use simplex to find optimal Embedding Dimesion
ts <- newdata$RFHarvest
lib <- c(1, length(ts))
pred <- c(1, length(ts))
simplex_output <- simplex(ts, lib, pred)

# plot simplex output
par(mfrow = c(1,1)) # reset to single panel plot
plot(simplex_output$E, simplex_output$rho, type = "l",ylab="rho",xlab="E")

# use smap to identify nonlinearity
smap_output <- s_map(ts, lib, pred, E = c(2,9))

# plot s-map output
par(mfrow = c(1,2))
plot(smap_output$theta[smap_output$E == 2], smap_output$rho[smap_output$E == 2], type = "l",xlab = "Nonlinearity (theta)", 
     ylab = "Forecast Skill (rho)")
plot(smap_output$theta[smap_output$E == 9], smap_output$rho[smap_output$E == 9], type = "l",xlab = "Nonlinearity (theta)", 
     ylab = "Forecast Skill (rho)")

plot(smap_output$theta, smap_output$rho, type = "l", xlab = "Nonlinearity (theta)", 
     ylab = "Forecast Skill (rho)")
#### use CCM to identify causality between Rockfish harvest and Kelp Bass Harvest ####

KB_xmap_RF <- ccm(newdata, lib, pred, E = 9, lib_sizes = seq(50, 500, by = 50), lib_column = "KBHarvest", target_column = "RFHarvest", num_samples = 100)
RF_xmap_KB <- ccm(newdata, lib, pred, E = 9, lib_sizes = seq(50, 500, by = 50), lib_column = "RFHarvest", target_column = "KBHarvest", num_samples = 100)

#get single rho value for each cross map
KB_xmap_RF1 <- ccm(newdata, lib, pred, E = 9, lib_sizes = length(ts), lib_column = "KBHarvest", target_column = "RFHarvest",  random_libs = FALSE)
RF_xmap_KB1 <- ccm(newdata, lib, pred, E = 9, lib_sizes = length(ts), lib_column = "RFHarvest", target_column = "KBHarvest",  random_libs = FALSE)



# get average cross map skill for each lib_size
xmap_means <- list(ccm_means(KB_xmap_RF), 
                   ccm_means(RF_xmap_KB))


# plot ccm between KB and RF
par(mfrow = c(1,1)) # reset to single panel plot
plot_limits <- range(xmap_means[[1]]$rho, xmap_means[[2]]$rho)
plot(xmap_means[[1]]$lib_size, xmap_means[[1]]$rho, type = "l", col = "red", ylim = plot_limits, ylab="rho", xlab="Library Size")
lines(xmap_means[[2]]$lib_size, xmap_means[[2]]$rho, col = "blue")
legend("right", legend = c("KB xmap RF", "RF xmap KB"), 
       col = c("red", "blue"), lwd = 1)


####Correlation Analysis (hoping for negative correlation)####
correlation<-cor(combodata$KBHarvest,combodata$RFHarvest)
correlation


#### generate surrogate temperature time series to test for significance ####
#create new data frame that is a list of all possible dates in the three year period


newdata$KBHarvest[na_idx]<-0
num_surr <- 500
surr_KB <- make_surrogate_data(newdata$KBHarvest, method = "seasonal", num_surr = num_surr, T_period = 365)
surr_KB[na_idx,] <- NA

newdata$KBHarvest[na_idx]<-NA

# plot surrogate time series
par(mfrow = c(4,1), mar = c(1,4,1,1))
plot(newdata$Date, newdata$KBHarvest, type = "l", xaxt = "n", col = "blue", ylab="Actual KB Harvest")
plot(newdata$Date, surr_KB[,1], type = "l", xaxt = "n",ylab="Surrogate 1")
plot(newdata$Date, surr_KB[,2], type = "l", xaxt = "n", ylab="Surrogate 2")
plot(newdata$Date, surr_KB[,3], type = "l", ylab= "Surrogate 3")



# compute ccm for actual KB harvest
RF_xmap_KB <- ccm(newdata, lib, pred, E = 6, lib_sizes = NROW(newdata), lib_column = "RFHarvest", target_column = "KBHarvest", random_libs = FALSE)

RF_xmap_surr_KB <- data.frame()
for(i in 1:num_surr)
{
  KB_block <- newdata
  KB_block$KBHarvest <- surr_KB[,i]
  RF_xmap_surr_KB <- rbind(RF_xmap_surr_KB, 
                                 ccm(KB_block, lib, pred, E = 6, lib_sizes = NROW(newdata), lib_column = "RFHarvest", target_column = "KBHarvest", random_libs = FALSE))
}


# compute a p-value
p_value <- (sum(RF_xmap_surr_KB$rho > RF_xmap_KB$rho) + 1) / (num_surr + 1)
pval<-signif(p_value,3)

# plot null distribution for RF xmap KB
par(mfrow = c(1,1), mar = c(4,4,1,1))
boxplot(RF_xmap_surr_KB$rho, ylim=c(0.2,0.7), ylab="rho")
points(1, RF_xmap_KB$rho, col = "blue", pch = 8, cex = 2)
legend("topleft", legend = c("Null (RF xmap surr KB)", "Actual (RF xmap KB)"), pch = c(NA, 8), col = c("black", "blue"), lwd = c(1, NA))
text(x=1.15, y=0.6, labels=italic("p = 0.002"))



#try with weekly grouping of data... LATER
weekdata<-data.frame(
  Date=seq(as.Date("2014-03-02"), as.Date("2016-10-29"),"weeks"),
  RFHarvest=sum(),
  KBHarvest=NA)
#cross reference the dtaes and fill in data from original data set
idx<-match(combodata$Date,weekdata$Date)
weekdata[idx,c(2,3)]<-combodata[,c(2,3)]
weekdata<-weekdata[-which(weekdata$Date==as.Date("2016-02-29")),]
na_idx <- is.na(weekdata$KBHarvest)
weekdata$KBHarvest[na_idx]<-0


##3Bs: bass, barracuda, bonito