#################################################################
#################################################################
#Title: Data Processing for Bayesian Hierarchical Model of 
# 		Atlantic Tropical Cyclones
#Author: Lindsey Dietz (diet0146@umn.edu)
#Objective: Process data to produce final data sets for analysis
#Created: 2/24/16
#Last updated: 3/17/2021 (by Sakshi Arya)
#################################################################
#################################################################
#Clear out any R junk
#rm(list=ls())

#Loads all necessary libraries
library(ggplot2)
library(grid)
library(reshape2)
library(plyr)
library(tseries)
library(stringr)
library(ncdf4)

file.original.21<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Data Sets/Original Data Sets 2021'
file.derived.21<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Data Sets/Derived Data Sets 2021'
#################################################################
#Useful Functions
#################################################################

#Function to melt data into correct format
melting_data<-function(data,name){
	data <-melt(data,id='YR')
	names(data)<-c('YR','MON',name)
	return(data)
}

#Correlation function
panel.cor <- function(x, y, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- (cor(x, y))
    txt <- format(round(r,2), digits = 2)[1]
    txt <- paste0(prefix, txt)
    text(0.5, 0.5, txt, cex = 2)
}
#################################################################
#################################################################
#Data Set 1: ICAT Damage Data
#http://www.icatdamageestimator.com/viewdata
#Note that this is not using modification for NOAA top 30 storms
#################################################################
#################################################################
Damages.21 <- read.csv(file.path(file.original.21,'stormDataICATDamages.csv'))

#Subsetting to include 1960-2017 (only goes up to 2017)
Damages.21$Year <- as.integer(substring(Damages.21$LANDFALL.DATE, 8,11))
Damages.21$Month <- substring(Damages.21$LANDFALL.DATE, 1,3)
Damages.21$Month <- as.factor(Damages.21$Month)
month.names <- substring(month.name, 1,3)
Damages.21$Month <- match(Damages.21$Month, month.names)

Damages_data.21<-Damages.21[Damages.21$Year>1959,]


#Adding a variable to help with merging to Hurdat data
Damages_data.21$first.letter<-toupper(substr(Damages_data.21$STORM.NAME,1,3))

#Modifying those below hurrican force to have category 0, i.e. TDs and TS
Damages_data.21$CATEGORY <- as.factor(Damages_data.21$CATEGORY)
levels(Damages_data.21$CATEGORY)[6]<-'0'

#Changing it so the categories are ordered
Damages_data.21$CATEGORY <-relevel(Damages_data.21$CATEGORY,'0')

#Modifying the names so they conform to R style
names(Damages_data.21)[4:5]<-c('CURRENT.DAMAGE.2021','BASE.DAMAGE')
Damages_data.21$STORM.NAME<-as.character(Damages_data.21$STORM.NAME)


# Removing commas from damage values
Damages_data.21$CURRENT.DAMAGE.2021<- as.numeric(gsub(",","",Damages_data.21$CURRENT.DAMAGE.2021))
Damages_data.21$BASE.DAMAGE<- as.numeric(gsub(",","",Damages_data.21$BASE.DAMAGE))



# Aggregating damage by storm; 
# some storms have separate damages for 
# each state a storm might make landfall in
Damages_Aggregate_Storm.21 <-ddply(Damages_data.21, ~ Year+Month +STORM.NAME+ first.letter ,summarize, BASE.DAMAGE =sum(BASE.DAMAGE), CURRENT.DAMAGE.2021 =sum(CURRENT.DAMAGE.2021), CATEGORY=mean(as.numeric(CATEGORY)-1), WINDS.MPH=mean(WINDS.MPH.))

#Manually Modifying the names of some storms based on manual inspection of Hurdat naming
Damages_Aggregate_Storm.21[Damages_Aggregate_Storm.21$STORM.NAME== 'Not Named 1960',]$STORM.NAME<-'Unnamed 1960'
Damages_Aggregate_Storm.21[Damages_Aggregate_Storm.21$STORM.NAME== 'Subtrop 1 1974',]$STORM.NAME<-'Unnamed 1974'
Damages_Aggregate_Storm.21[Damages_Aggregate_Storm.21$STORM.NAME== 'Subtrop 1 1982',]$STORM.NAME<-'Unnamed 1982'

Damages_Aggregate_Storm.21$first.letter<-toupper(substr(Damages_Aggregate_Storm.21$STORM.NAME,1,3))

#################################################################
#################################################################
#Data Set 2: Hurdat North Atlantic Data
#http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html
#Updated through Nov 2019
#################################################################
#################################################################
#Note Max wind speed is in knots
hurdat_NA<-read.table(file.path(file.derived.21,'summarized_atl_data_1960_2019.txt'),header=T)

#Missing values are described by -99 so replacing them with -Inf
hurdat_NA[hurdat_NA$maxWS==-99,]$maxWS<--Inf

#Excuding 2020 because the damages data doesn't go that far
hurdat_NA <-hurdat_NA[hurdat_NA$Year<2020,]

hurdat_NA$first.letter<-toupper(substr(str_trim(hurdat_NA$name),1,3))
hurdat_NA$maxWS <-round(hurdat_NA$maxWS)
hurdat_NA$minCP <-round(hurdat_NA$minCP)
hurdat_NA$avgLat<-floor(hurdat_NA$avgLat)
hurdat_NA$avgLon<-floor(hurdat_NA$avgLon)

#################################################################
#################################################################
#Data Set 3: Merged Hurdat and ICAT data
#1960-2013
#################################################################
#################################################################

#Damages merged with all hurricanes
final_data<-merge(Damages_Aggregate_Storm.21, hurdat_NA,by=c('Year','first.letter'),all = TRUE)

#Manual fix for the subtropical storms in 1974, 1982, 1987 
final_data[c(216:223,225:228),c('Month','STORM.NAME', 'BASE.DAMAGE', 'CURRENT.DAMAGE.2015', 'Deflation_2014_ICATWunderground1', 'Deflation_2014_ICATWunderground2','CATEGORY')]<-NA
final_data[c(377,379),c('Month','STORM.NAME', 'BASE.DAMAGE', 'CURRENT.DAMAGE.2015', 'Deflation_2014_ICATWunderground1', 'Deflation_2014_ICATWunderground2','CATEGORY')]<-NA
final_data[c(436:437,439:443),c('Month','STORM.NAME', 'BASE.DAMAGE', 'CURRENT.DAMAGE.2015', 'Deflation_2014_ICATWunderground1', 'Deflation_2014_ICATWunderground2','CATEGORY')]<-NA

#Setting up Saffir-Simpson categories in the data
final_data$max_Cat<-with(final_data,cut(maxWS,breaks=c(-Inf,64,83,96,113,137,250),right=F,include.lowest = T))
levels(final_data$max_Cat)<-c(0:5)

#Setting up Low (TS-2) and High (3-5) intensity categories in both data sets
final_data$Landfall_Cat_ICAT <-with(final_data,cut(CATEGORY,breaks=c(0,3,5),right=F,include.lowest = T))
levels(final_data$Landfall_Cat_ICAT)<-c('TS-2','3-5')

final_data$Cat_HURDAT <-with(final_data,cut(as.numeric(final_data$max_Cat)-1,breaks=c(0,3,5),right=F,include.lowest = T))
levels(final_data$Cat_HURDAT)<-c('TS-2','3-5')

#Setting a damage indicator by storm
final_data$DamageInd<-'Yes'
final_data[is.na(final_data$STORM.NAME),]$DamageInd<-'No'
final_data$DamageInd<-factor(final_data$DamageInd)
final_data$DamageInd <-relevel(final_data$DamageInd,'Yes')

#################################################################
#################################################################
#Data Set 4: Teleconnection Indices- Includes NAO
#https://psl.noaa.gov/data/correlation/nao.data
# Dowloaded: 03/15/2021
#################################################################
#################################################################
#Read in Teleconnections data set
indices.21 <- read.csv(file.path(file.original.21,'nao.data.txt'),header=F, skip= 1, na.strings=c('-99.90',NA), sep = "")
dim(indices.21)
names(indices.21) <- c("Year", c(1:12))
ind.21 <- indices.21[-c(1:2,75:77),]
ind.21.long <- melt(ind.21, id = "Year", na.rm = TRUE)
names(ind.21.long) <- c("YR","MON","NAO")
ind.21.long <- ind.21.long[order(ind.21.long$YR),]
ind.21.long$MON <- as.integer(ind.21.long$MON)
summary(ind.21.long)
ind.21.long$YR <- as.integer(ind.21.long$YR)
ind.21.long$NAO <- as.numeric(ind.21.long$NAO)

#Plot of time series with red lines for the portion used in analysis
##pdf("Monthly_NAO_Index_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(ind.21.long$YR+ind.21.long$MON/12,ind.21.long$NAO,type='l',xlab='Year',ylab='Monthly NAO Index',pch=15,cex = 0.5, cex.axis=1.5, cex.lab = 1.5,xlim=c(1950,2021))
abline(v=1960,lwd=2,col='red')
abline(v=2020,lwd=2,col='red')
###dev.off()

#################################################################
#################################################################
#Data Set 5: Southern Oscillation Index (SOI)
#ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/soi
#################################################################
#################################################################
#Read in SOI data set
soi.21 <-read.csv(file.path(file.original.21,'soi'),header=F,na.strings='-999.9', sep = "")
soi.21 <- soi.21[1:74,]
soi.21 <- soi.21[-c(1:3),]
names(soi.21)<-c('YR',1:12)
soi.21[71, 3]  <- 2.5


#Formatting data to merge
soi.21.long <-melting_data(soi.21,'SOI')
soi.21.long<-soi.21.long[order(soi.21.long$YR,as.numeric(soi.21.long$MON)),]
soi.21.long$YR <- as.integer(soi.21.long$YR)
soi.21.long$MON <- as.integer(soi.21.long$MON)
soi.21.long$SOI <- as.numeric(soi.21.long$SOI)

#Plot of time series with red lines for the portion used in analysis
##pdf("Monthly_SOI_Index_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(soi.21.long$YR+ as.numeric(soi.21.long$MON)/12, soi.21.long$SOI,xlab='Year',ylab='Monthly SOI Index',type='l',pch=15,cex = 0.5, cex.axis=1.5, cex.lab = 1.5,xlim=c(1950,2021))
abline(v=1960,lwd=2,col='red')
abline(v=2020,lwd=2,col='red')
###dev.off()

#Merging with previous data set
cov1<-merge(soi.21.long, ind.21.long[c('YR','MON','NAO')],by=c('YR','MON'))

#################################################################
#################################################################
#Data Set 6: AMO instead of Atlantic SST (unsmoothed)
#http://www.esrl.noaa.gov/psd/data/timeseries/AMO/
#################################################################
#################################################################
#Read in AMO data set
amo.21 <-read.csv(file.path(file.original.21,'AMO_2021.txt'),header=F,skip=1,na.strings='-99.990', sep = "")
amo.21 <- amo.21[1:166,]

#Formatting data to merge
names(amo.21) <- c('YR', 1:12)
amo.21 <-melting_data(amo.21,'AMO')
amo.21 <-na.omit(amo.21[order(amo.21$YR,as.numeric(amo.21$MON)),])
amo.21$YR <- as.integer(amo.21$YR)
amo.21$MON <- as.integer(amo.21$MON)
amo.21$AMO <- as.numeric(amo.21$AMO)

#Plot of time series with red lines for the portion used in analysis
##pdf("Monthly_AMO_Index_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(amo.21$YR+ as.numeric(amo.21$MON)/12, amo.21$AMO,type='l',xlab='Year',ylab='Monthly AMO Index',pch=15,cex = 0.5, cex.axis=1.5, cex.lab = 1.5,xlim=c(1850,2021), xaxs="i")
abline(v=1960,lwd=2,col='red')
abline(v=2020,lwd=2,col='red')
##dev.off()

#Merging with previous data set
cov2a<-merge(cov1, amo.21,by=c('YR','MON'))

#################################################################
#################################################################
#Data Set 7: Nino Indices 
#http://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii
#################################################################
#################################################################
#Read in Nino Indices data set
ninos.21 <-read.table(file.path(file.original.21,'ersst5.nino.mth.91-20.ascii.txt'),header=T)

#Merging with previous data set (only keeping Nino 3.4 Anomaly)
cov2b<-merge(cov2a, ninos.21[,c('YR','MON','ANOM.3')],by=c('YR','MON'))

#################################################################
#################################################################
#Data Set 8: Atlantic SST
#https://psl.noaa.gov/cgi-bin/db_search/DBSearch.pl?Dataset=NOAA+Extended+Reconstructed+SST+V5&Variable=Sea+Surface+Temperature&group=0&submit=Search
#All data not just atlantic
#################################################################
#################################################################
#Read in Atlantic SST data set- data is monthly from Jan 1960- Sep 2015
#which gives 669 months
SST.21 <- nc_open(file.path(file.original.21, 'sst.mnmean.nc'))
SST_array.21 <- ncvar_get(SST.21)
#time <- ncvar_get(AtlSST.21,  'time')
lat.21 <- ncvar_get(SST.21,  'lat')
lon.21 <- ncvar_get(SST.21,  'lon')

## Lat and Lon for the NA basin (10–25N, 80–20W area)
lat.21.Atl <- which(lat.21 > 9 & lat.21 < 26)
lon.21.Atl <- which(lon.21 > 279 & lon.21 < 341)

monthly_mean_ATL_SST.21 <-vector()

# The original data is from 1854 to 2021 so we get rid of the data from 1854 to 1960:
number.of.months.21 <- seq(1272,dim(SST_array.21)[3],by = 1)

#Derives the monthly mean over the 10–25N, 80–20W area
for(i in 1:length(number.of.months.21)){
  monthly_mean_ATL_SST.21[i]<-mean(SST_array.21[lon.21.Atl,lat.21.Atl,number.of.months.21[i]],na.rm=T)	
}

#Putting the monthly mean values into a data frame similar to the others
sst.21 <-data.frame(expand.grid(c(1:12),c(1960:2021))[1:length(number.of.months.21),], monthly_mean_ATL_SST.21)
names(sst.21)<-c('MON','YR','Atl_SST')

#Merging with previous data set 
cov2<-merge(cov2b, sst.21,by=c('YR','MON'))
cov2 <-cov2[order(cov2$YR,as.numeric(cov2$MON)),]

#Plot of time series with red lines for the portion used in analysis
##pdf("Monthly_Atl_SST_new.pdf", width = 12, height = 6)
par(mar=c(3,3,0,0)+0.1, mgp =c(2,1,0))#sets margins of plotting area
plot(cov2$YR+ as.numeric(cov2$MON)/12, cov2$Atl_SST,type='l',xlab='Year',ylab='Monthly Atlantic SST Index',pch=15,cex=0.5,xlim=c(1960,2021))
abline(v=1960,lwd=2,col='red')
abline(v=2020,lwd=2,col='red')
#dev.off()

# Plot for Nino 3.4 anomaly:
##pdf("Monthly_Anom34_new.pdf", width = 12, height = 6)
par(mar=c(3,3,0,0)+0.1, mgp =c(2,1,0))#sets margins of plotting area
plot(cov2$YR+ as.numeric(cov2$MON)/12, cov2$ANOM.3,type='l',xlab='Year',ylab='Monthly Nino 3.4 Anomaly',pch=15,cex=0.5,xlim=c(1960,2021))
abline(v=1960,lwd=2,col='red')
abline(v=2020,lwd=2,col='red')
#dev.off()

#################################################################
#################################################################
#Data Set 9: Sunspots
#http://www.sidc.be/silso/datafiles
#################################################################
#################################################################
#Read in Sunspot number datafile
sunspots.21 <-read.csv(file.path(file.original.21,'SN_m_tot_V2.0.csv'),sep=';',header=F)
names(sunspots.21)[1:7]<-c('YR','MON','YRMON','Sunspot_Avg','Sunspot_SD','Num_Obs_Mon','Revision')

#Plot of time series with red lines for the portion used in analysis
#pdf("Monthly_Sunspots_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
with(sunspots.21 ,plot(YRMON, Sunspot_Avg,type='l',xlab='Year',ylab='Monthly Average Sunspots',cex = 0.5, cex.axis=1.5, cex.lab = 1.5))
abline(v=1960,lwd=2,col='red')
abline(v=2020,lwd=2,col='red')
#dev.off()
# Creating July to June averages for each year
# we set all MON= 5 for the derivation in the next part
# %/% keeps the integer part (year in this case)
# sunspots[sunspots$YRMON> 1749.455,] => starts from July 1749
julyjuneavgs.21<-data.frame(YR=1750:2021, MON=5, Sunspots=with(sunspots.21[sunspots.21$YRMON> 1749.455,], 
                                                            sapply(
                                                              split(Sunspot_Avg, (seq_along(Sunspot_Avg) - 1) %/% 12), 
                                                              mean)))

#Merging with previous data set 
sunspots2.21<-merge(sunspots.21, julyjuneavgs.21,by=c('YR','MON'))

#Keeping only those from 1960 forward 
sunspots_1960.21 <-sunspots2.21[sunspots2.21$YR>1959,]

#Merging with previous data set 
cov2<-merge(cov2, sunspots_1960.21,by=c('YR'))
cov2 <-cov2[order(cov2 $YR,as.numeric(cov2 $MON.x)),]
names(cov2)[2]<-'MON'
#################################################################
#################################################################
#Final Derived Data sets
#################################################################
#################################################################

#Keeping only May/June averages of covariates
cov3b<-cov2[cov2$MON %in% c(5:6) & cov2$YR %in% c(1960:2019),]
cov.val<-cov2[cov2$MON %in% c(5:6) & cov2$YR %in% c(2020:2021),]
cov4<-ddply(cov3b,~ YR,summarize, NAO = mean(NAO), SOI =sum(SOI), AMO =mean(AMO), ANOM.3.4=mean(ANOM.3),Atl_SST =mean(Atl_SST), Sunspots =mean(Sunspots))
cov.val2<-ddply(cov.val,~ YR,summarize, NAO = mean(NAO), SOI =sum(SOI), AMO =mean(AMO), ANOM.3.4=mean(ANOM.3), Atl_SST =mean(Atl_SST), Sunspots =mean(Sunspots))

#Time series of each index using only May/June Avg.
#pdf("NAO_May_June_Avg_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(cov4$YR, cov4$NAO,type='l',pch=15,cex = 0.5, cex.axis=1.5, cex.lab = 1.5,xlab='Year',#main='May/June Average NAO',
ylab='NAO Index May/June Avg.')
#dev.off()

#pdf("SOI_May_June_Avg_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(cov4$YR, cov4$SOI,type='l',pch=15,cex = 0.5, cex.axis=1.5, cex.lab = 1.5,xlab='Year',#main='May/June Average SOI',
ylab='SOI Index May/June Avg.')
#dev.off()

#pdf("AMO_May_June_Avg_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(cov4$YR, cov4$AMO,type='l',pch=15,cex = 0.5, cex.axis=1.5, cex.lab = 1.5,xlab='Year',#main='May/June Average AMO',
ylab='AMO Index May/June Avg.')
#dev.off()

#pdf("ANOM34_May_June_Avg_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(cov4$YR, cov4$ANOM.3.4,type='l',pch=15,cex=0.5,xlab='Year',#main='May/June Average ANOM 3.4',
ylab='Nino 3.4 Anomaly May/June Avg.')
#dev.off()

#pdf("Atl_SST_May_June_Avg_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(cov4$YR, cov4$Atl_SST,type='l',pch=15,cex=0.5,xlab='Year',#main='May/June Average Atlantic SST',
ylab='Atlantic SST May/June Avg.')
#dev.off()

#pdf("Sunspot_Sept_Avg_new.pdf", width = 12, height = 6)
par(mar=c(5,5,0,0)+0.1, mgp =c(2,0.5,0))#sets margins of plotting area
plot(cov4$YR, cov4$Sunspots,type='l',pch=15,cex = 0.5, cex.axis=1.5, cex.lab = 1.5,xlab='Year',ylab='Sun Prev. July to June Avg.')
#dev.off()

#Correlations among all of the possible covariates
#pdf("Correlation_AMO_SOI_NAO_new.pdf")
pairs(cov4[,c(2,3,5,4,6,7)],pch=16,cex=.7,lower.panel = panel.smooth, upper.panel = panel.cor)
#dev.off()

#Correlations among the selected covariates only
pairs(cov4[,c(2,3,4,7)],pch=16,cex=.7,lower.panel = panel.smooth, upper.panel = panel.cor)

#Calculate annual damages by category (TS-2,3-5)
damage=ddply(final_data,~Year+ Cat_HURDAT,summarize, damage=sum(BASE.DAMAGE,na.rm=T))

#Calculate annual number of storms by category (TS-2,3-5)
storms= cbind(Year=c(1960:2019,1960:2019),rbind(data.frame(Freq= with(final_data,table(Year, Cat_HURDAT))[,1], Cat_HURDAT ='TS-2'), data.frame(Freq=with(final_data,table(Year, Cat_HURDAT))[,2], Cat_HURDAT ='3-5')))

#Calculate annual total number of storms
annual_storms<-ddply(storms,~ Year,summarize,total=sum(Freq))

#Add annual total number of storms to the annual covariates
cov5<-data.frame(annual_storms, cov4)

#Calculate annual number of landfalling storms by category (TS-2,3-5)
landfalling= cbind(Year=c(1960:2019,1960:2019),rbind(data.frame( Landfall= with(final_data,table(Year, DamageInd =='Yes',Cat_HURDAT))[,,1][,2], Cat_HURDAT ='TS-2'),data.frame(Landfall=with(final_data,table(Year, DamageInd =='Yes',Cat_HURDAT))[,,2][,2], Cat_HURDAT ='3-5')))

#Merge categorized landfalls, storms, and damage into one data set
cat_data<-rbind(merge(landfalling, merge(storms, damage, by=c('Year', 'Cat_HURDAT')),by=c('Year', 'Cat_HURDAT')),data.frame(Year=c(1968, 1972, 1986, 1994, 2019), Cat_HURDAT='3-5', Landfall=0, Freq=0,damage=0))

#Output Data Sets for analysis
write.table(final_data,
	file= file.path(file.derived.21,'Categorized_Storm_1960_2019.csv'),row.names=F,sep=',')
write.table(cat_data,
            file=file.path(file.derived.21,'Categorized_Annual_1960_2019.csv'),row.names=F,sep=',')
write.table(cov5,
	file=file.path(file.derived.21,'Annual_Covariates_1960_2019.csv'),row.names=F,sep=',')





