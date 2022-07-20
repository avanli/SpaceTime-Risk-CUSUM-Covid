# Tract level COVID Case count surveillance using CUSUM
# Author:   Arda Vanli
# Version: 2/20/2022
# see the helper scripts: plot_correlationSVIandCovid.R and plot_CDCThemes.R 
#### libraries
library(tigris)
library(rgeos)
library(sp)
library(sf)
library(spdep)
library(spatialreg)
library(MMWRweek)
library(censusapi)
library(plotrix)
library(scales)
library(lattice)
#### paths, flags, sourced scripts
flagTuneMethods=FALSE  
flagCreatePNG=FALSE
#paths
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path1='data/'
path2=''

source(paste(path2,"adjustForCorrelSAR.R",sep=""))
source(paste(path2,"getLikelihoodSAR.R",sep=""))
source(paste(path2,"get_spaceTimeCUSUMNormalStandard.R",sep=""))
source(paste(path2,"get_ARL_SCUSUMCorrel_NeighborCovarGiven.R",sep=""))
source(paste(path2,"simulateSAR_NBCovarGiven.R",sep=""))

#### shape file for FL tracts ####
options(tigris_use_cache = TRUE) #ZCTAs can take several minutes to download. To cache the data use this.
fl_tracts=as_Spatial(tracts("Florida",year=2018,cb=TRUE))
#### Florida zip code populations and tracts ####
# Retrieve tract-level data for any variable within a specific state/county 
# census api: https://github.com/hrecht/censusapi#api-key-setup
# new key obtained from http://api.census.gov/data/key_signup.html
Sys.setenv(CENSUS_KEY="0e95aa62e853f18e1d8099dbbc0ca054cdf14f81")
readRenviron("~/.Renviron")
Sys.getenv("CENSUS_KEY")
# uncomment the following identify variables in sf1. Here we use TOTAL POPULATION IN OCCUPIED HOUSING UNITS 
# mytab=listCensusMetadata(name = "dec/sf1", vintage = 2010,type='variables')
fl2010 <- getCensus(name = "dec/sf1", vintage = 2010,
                    vars = c("H010001"),
                    region = "tract:*", regionin = "state:12")
# read "zcta to tract relationship" file 
#zctaToTractFile=read.csv(paste(path2,"zcta_tract_rel_10.csv",sep =''))
zctaToTractFile=read.csv("https://www2.census.gov/geo/docs/maps-data/data/rel/zcta_tract_rel_10.txt")
zctaToTractFileFL=zctaToTractFile[which(zctaToTractFile$STATE=='12'),]
(noFLTracts=dim(fl_tracts)[1])
(dim(zctaToTractFileFL))

# how many tracts does each zip code have?
noOfTractsZipCodes=as.matrix(table(factor(zctaToTractFileFL$ZCTA5))) # rownames(noOfTractsZipCodes)[1]
(noFLZipCodes=dim(noOfTractsZipCodes)[1])
FLZipCodes=unique(zctaToTractFileFL$ZCTA5)  
# to match the GEOID field in zctaToTractFileFL, create geoID in fl2010
fl2010$geoID=NA
for (i in 1: dim(fl2010)[1]){
  fl2010$geoID[i]=paste(fl2010$state[i],fl2010$county[i],fl2010$tract[i],sep="")
}
matchTractToPop=match(fl_tracts$GEOID,fl2010$geoID) 
fl_tracts$pop=fl2010$H010001[matchTractToPop]
fl_tracts$validateGEOID=fl2010$geoID[matchTractToPop]

matchTractToZipCodes=match(fl_tracts$GEOID,as.character(zctaToTractFileFL$GEOID)) 
fl_tracts$ZipCode=zctaToTractFileFL$ZCTA5[matchTractToZipCodes]

#### FL zip codes weekly case counts ####
covid_case=read.csv(paste(path1,"historical_testingbyzip.csv",sep =''))
covid_case$Date=as.Date(covid_case$date,"%m/%d/%y")
covid_case$Week=format(covid_case$Date,"%W")
uniqueWeeks=unique(covid_case$Week)
uniqueZipCodesC=unique(covid_case$zip)
indCovidCaseofFLzipcodes=match(FLZipCodes,uniqueZipCodesC)
uniqueWeeksMonday=data.frame('week'=NA)
for (i in 1:length(uniqueWeeks)){
  uniqueWeeksMonday[i]=MMWRweek2Date(MMWRyear = 2020,MMWRweek = as.integer(uniqueWeeks[i]),MMWRday = 2)
}

covid_weeklyCase=matrix(NA,noFLZipCodes,length(uniqueWeeks))
for (i in 1:noFLZipCodes){ # in the order as FLZipCodes
  for (j in 1:length(uniqueWeeks)){
    covid_weeklyCase[i,j]=
      sum(as.numeric(covid_case$number_of_positives[covid_case$zip==uniqueZipCodesC[indCovidCaseofFLzipcodes[i]]
                                                    &covid_case$Week==uniqueWeeks[j]]))    
  }
}
# replace NA's by 0. NA means no case was observed that week.
covid_weeklyCase[is.na(covid_weeklyCase)]=0
rownames(covid_weeklyCase)=FLZipCodes
noWeeks=length(uniqueWeeks)
weeknames=1:noWeeks
fl_tracts=cbind(fl_tracts, setNames(lapply(weeknames, function(x) x=NA), weeknames) )


#### FL tracts weekly case counts ####
for (i in 1:length(fl_tracts)){
  i1=which(fl_tracts$ZipCode[i]==FLZipCodes)
  if(!is.na(fl_tracts$ZipCode[i])){
    # covid weekly counts is apportioned proportional to population
    for (j in 1:noWeeks){
      fl_tracts[i,12+j]=covid_weeklyCase[i1,j]/noOfTractsZipCodes[i1]  # to be correct, multiply with fl_tracts$pop[i]/totpop. But assume same counts in all tracts
    }
  }
}

#### infection rates ####
# non zero pop census tracts
# indNonZeroTract=which(fl_tracts$pop!=0)
# non zero pop and panhandle tracts
indNonZeroTract=which(fl_tracts$pop!=0&coordinates(fl_tracts)[,1]<(-83.347))
panhandleTracts=fl_tracts[indNonZeroTract,]
# population of 12-th tract is 1, creating outlier. replace it with average pop.
fl_tracts$pop[indNonZeroTract][12]=mean(fl_tracts$pop[indNonZeroTract])
# number of tracts
ns2=length(indNonZeroTract)
ns=noFLTracts
# total counts
total_weekly=colSums(fl_tracts@data[indNonZeroTract,13:(12+noWeeks)])
# baseline infection rate per person-first 4 weeks constitute the baseline window
lambda0_hat_time=NA
n_bl=4
for (i in 1:noWeeks){
  # weekly incidence rate per 100000 for the entire study area over time
  lambda0_hat_time[i]=total_weekly[i]/sum(fl_tracts$pop[indNonZeroTract])
}
# average incidence rate for the baseline period
lambda0_hat=mean(lambda0_hat_time[1:n_bl])
# define transformed variable- on nonzero population ZCs - should be normally distributed
y_weekly=log((fl_tracts@data[indNonZeroTract,13:(12+noWeeks)]+1)/fl_tracts$pop[indNonZeroTract])
matplot(t(y_weekly),lty=1,col='grey',type="l",main="transformed weekly tract counts")
mu0_hat=apply(y_weekly[,1:n_bl],1,mean)
sigSq0_hat=apply(y_weekly[,1:n_bl],1,var)
# assume constant variance for all zip codes, replicate average variance
sig0_hat=rep(sqrt(mean(sigSq0_hat[sigSq0_hat!=0])),ns2)
# plot infection rate
if(flagCreatePNG){png("infectionRate.png", width = 7, height =5, units = 'in', res = 600)}
par(mfrow=c(1,1),mar=c(3,3,0.5,0)+0.1,mgp=c(2,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
plot(lambda0_hat_time*100000,type='b',xlab='week',ylab='weekly Incidence per 100000',pch=1,
     cex.lab=1.2,cex.main=1,cex.axis=1.2,cex=1,xaxt="none")
axis(1,at=1:length(uniqueWeeksMonday),labels=format(uniqueWeeksMonday,"%m/%d"))
lines(c(1,n_bl),c(lambda0_hat*100000,lambda0_hat*100000),col='red',lty=2,las=2,lwd=2)
if(flagCreatePNG){dev.off()}

# ODmatrix, coordinates
ODmatrix=spDists(coordinates(fl_tracts)) 
ODmatrix2=ODmatrix[indNonZeroTract,indNonZeroTract]
X0=coordinates(fl_tracts[indNonZeroTract,])
# plot panhandle tracts
if(flagCreatePNG){png("panhandleTracts.png", width = 7, height =5, units = 'in', res = 600)}
par(mfrow=c(1,1),pty="s",mar=c(3,3,0.5,0)+0.1,mgp=c(1.5,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
plot(fl_tracts[,1],col='white',border='grey',main="",axes = FALSE)
plot(fl_tracts[coordinates(fl_tracts)[,1]<(-83.347),],border='grey',col='navajowhite2',add=TRUE)
if(flagCreatePNG){dev.off()}
# show centroids
if(flagCreatePNG){png("centroids.png", width = 7, height =5, units = 'in', res = 600)}
par(mfrow=c(1,1),mar=c(3,3,0.5,0)+0.1,mgp=c(2,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
plot(fl_tracts[coordinates(fl_tracts)[,1]<(-83.347),],border='grey',col='white')
points(coordinates(fl_tracts),cex=0.5)
if(flagCreatePNG){dev.off()}
# show scan circles
if(flagCreatePNG){png("scanCircles.png", width = 7, height =5, units = 'in', res = 600)}
par(mfrow=c(1,1),mar=c(3,3,0.5,0)+0.1,mgp=c(2,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
plot(fl_tracts[coordinates(fl_tracts)[,1]<(-83.347),],border='grey',col='white')
#identify(coordinates(fl_tracts)[,1],coordinates(fl_tracts)[,2],fl_tracts$GEOID)
rad_cluster=c(0,.25,.50)# to contain half of the counties
# from plotrix package
draw.circle(coordinates(fl_tracts[fl_tracts$GEOID=='12013010200',])[1],
            coordinates(fl_tracts[fl_tracts$GEOID=='12013010200',])[2],
            rad_cluster,
            border="blue",lty=3,lwd=2)
points(coordinates(fl_tracts[fl_tracts$GEOID=='12013010200',]),cex=0.5)
if(flagCreatePNG){dev.off()}
#### find ARL0 - non risk or correlation adjusted ####
# radius values
step_size=0.25
r=seq(0,2*step_size,by=step_size) 
nr=length(r)

n_MC=10
# baseline mean, standard deviation and risk adjustment coefficients
b0=mu0_hat # intercept equal to subregion means
b1=0 # non risk adjusted
sd_err=sig0_hat[1]
delt=0 # in control ARL
delt_c=sd_err*3  # control chart parameter
ind_inject=1  # immaterial for ARL0
flag_plot=FALSE
flag_moran=FALSE
time_inj=1
rho=0
l0=0 # mean of risk adjusted data

nb4rt2 <- poly2nb(fl_tracts[indNonZeroTract,])
xMatrix=matrix(0,ns2,1)  # should be 0 for non risk adjusted

if (flagTuneMethods){
  h_s=13.8  #h=13.8 gives ARL(seARL)=8.100000
  flagRiskCorrAdjust=TRUE  # with xMatrix=0,rho=0 and b1=0, this is equivalent to non risk adjusted
  set.seed(98765)
  res=get_ARL_SCUSUMCorrel_NeighborCovarGiven(r, X0, b0,b1,rho,sd_err,l0,delt_c, delt, h_s, 
                                              time_inj, n_MC,ind_inject,flag_plot,flag_moran,
                                              nb4rt2,xMatrix,flagRiskCorrAdjust)
  print(c("limit, delta,ARL,seRL"))
  print(c(h_s,delt,res$ARL,res$seRL))
}
#### CUSUM surveillance of infection data-non risk or correlation adjusted ####
S_prev=matrix(0,nr,ns2)
tau=matrix(0,nr,ns2)
S_max=0
Smax_vec=0
wk_begin=3
wk_end=11
wk=wk_begin
# cusum parameter: 3 sigma shift
h_s=13.8
while (S_max<h_s&wk<wk_end){
  # when sigSq0_hat=0 the numerator is also 0, so make yAdj=0 when that happens
  yAdj=(y_weekly[,wk]-mu0_hat)/sqrt(sigSq0_hat); yAdj[is.na(yAdj)]=0; yAdj[is.infinite(yAdj)]=0

  res_S=get_spaceTimeCUSUMNormalStandard(yAdj,X0,S_prev,delt_c/2,r)
  
  S_prev=res_S$S_new
  tau[which(S_prev==0)]=wk 
  # update cusum statistic
  S_max=res_S$S_max_t
  Smax_vec=c(Smax_vec,S_max)
  print(wk)
  wk=wk+1
}
# change point estimator
ind_max=which(res_S$S_new==S_max)
tau_est=tau[ind_max[1]]
print(paste('change point=',tau_est,sep=""))
# plot CUSUM
if(flagCreatePNG){png("cusumUnAdj.png", width = 7, height =5, units = 'in', res = 600)}
par(mfrow=c(1,1),mar=c(3,3,0.5,0)+0.1,mgp=c(2,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
plot(seq(wk_begin,wk),Smax_vec,type='b',xlab='Week',ylab='CUSUM',xaxt='none')
axis(1,at=seq(wk_begin,wk),labels=format(uniqueWeeksMonday[wk_begin:wk],"%m/%d"))
abline(h=h_s,lty=2)
if(flagCreatePNG){dev.off()}

# plot identified cluster - zero pop ZCs removed 
clust_center=as.character(fl_tracts$GEOID[indNonZeroTract[res_S$region_max]])
print(paste("cluster center is at tract",clust_center,sep=""))
rad_cluster=r[res_S$radius_max]
# second and third largest clusters
reg_max2=res_S$RadiusRegionTop3[2,2]; reg_max3=res_S$RadiusRegionTop3[2,3]
rad_max2=r[res_S$RadiusRegionTop3[1,2]]; rad_max3=r[res_S$RadiusRegionTop3[1,3]]
# find tracts included in the top cluster- zero pop ZCs not removed
ind_InCluster=which(ODmatrix[indNonZeroTract[res_S$region_max],]<=rad_cluster)
# find tracts included in the second and third largest clusters
ind_InCluster2=which(ODmatrix[indNonZeroTract[reg_max2],]<=rad_max2)
ind_InCluster3=which(ODmatrix[indNonZeroTract[reg_max3],]<=rad_max3)
colorTracts=rep(NA,ns)
colorTracts[indNonZeroTract]='navajowhite2'
colorTracts[ind_InCluster2]="red"       # uncomment to show second and third largest clusters
colorTracts[ind_InCluster3]="blue"
colorTracts[ind_InCluster]="seagreen4"
ylim1=c(29.753158,30.532798)
xlim1=c(-88.463172,-82.272418)
if(flagCreatePNG){png("ClusterUnAdj.png", width = 7, height =5, units = 'in', res = 600)}
plot(panhandleTracts,border='grey',main="") #main= paste('Clusters identified in wk ',wk))
plot(fl_tracts,border=NA,add=TRUE,col=alpha(colorTracts,0.5),#xlim=xlim1,ylim=ylim1,
     axes = FALSE)
text(coordinates(fl_tracts)[indNonZeroTract[res_S$region_max],1]+0.1,coordinates(fl_tracts)[indNonZeroTract[res_S$region_max],2]+.05,"1")
text(coordinates(fl_tracts)[indNonZeroTract[reg_max2],1]+0.1,coordinates(fl_tracts)[indNonZeroTract[reg_max2],2]+.05,"2")
text(coordinates(fl_tracts)[indNonZeroTract[reg_max3],1]+0.1,coordinates(fl_tracts)[indNonZeroTract[reg_max3],2]-0.1,"3")
if(flagCreatePNG){dev.off()}
# plot counts, normalized counts and transformed data for all tracts and cluster
par(mfrow=c(1,1),mar=c(3,3,0.5,0)+0.1,mgp=c(2,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
matplot(t(exp(y_weekly))*100000,type='l',
        col='gray',ylab='Incidence per 100000, zipcode',xlab='week',lty=1,log='y',main="")
#### fit SAR model between infection counts of week 3 and the SVI Theme 1, 2, 3 or 4 as a single covariate####
# read CDC SVI at tract level for FL 
CDC_SVI=read.csv(paste(path1,"cdc-sovi-tract 2018-csv-FL.csv",sep =''))

# reorder dataframe to match the shapefile order of the tracts
indToReorderCDC_SVI=match(fl_tracts$GEOID,as.character(CDC_SVI$FIPS))
CDC_SVIr=CDC_SVI[indToReorderCDC_SVI,]
#replace with na'a rows with -999 value in any of the vulnerability variable
CDC_SVIr$SPL_THEME1[which(CDC_SVIr$SPL_THEME1==-999)]=NA
CDC_SVIr$SPL_THEME2[which(CDC_SVIr$SPL_THEME2==-999)]=NA
CDC_SVIr$SPL_THEME3[which(CDC_SVIr$SPL_THEME3==-999)]=NA
CDC_SVIr$SPL_THEME4[which(CDC_SVIr$SPL_THEME4==-999)]=NA

# find average of k nearest neighbors for each theme
CDC_SVIr$RPL_THEME1Smooth=NA
kNeighbors=10
for (i in 1:noFLTracts){
  # indicators
  M1=cbind(ODmatrix[i,],0.25*CDC_SVIr$SPL_THEME1,0.25*CDC_SVIr$SPL_THEME2,0.5*CDC_SVIr$SPL_THEME3,0.2*CDC_SVIr$SPL_THEME4)
  sortedThemesByDist=M1[order(M1[,1],decreasing=FALSE),]
  
  CDC_SVIr$RPL_THEME1Smooth[i]=mean(sortedThemesByDist[1:kNeighbors,2],na.rm=TRUE)
  CDC_SVIr$RPL_THEME2Smooth[i]=mean(sortedThemesByDist[1:kNeighbors,3],na.rm=TRUE)
  CDC_SVIr$RPL_THEME3Smooth[i]=mean(sortedThemesByDist[1:kNeighbors,4],na.rm=TRUE)
  CDC_SVIr$RPL_THEME4Smooth[i]=mean(sortedThemesByDist[1:kNeighbors,5],na.rm=TRUE)
}
# smoothing - use for SAR modeling
tractSVI=cbind(CDC_SVIr$RPL_THEME1Smooth,CDC_SVIr$RPL_THEME2Smooth,
               CDC_SVIr$RPL_THEME3Smooth,CDC_SVIr$RPL_THEME4Smooth)
# no smoothing - use for plotting
#tractSVI=cbind(CDC_SVIr$SPL_THEME1/4,CDC_SVIr$SPL_THEME2/4,
#               CDC_SVIr$SPL_THEME3/2,CDC_SVIr$SPL_THEME4/5)
# manual calculation of indicators for panhandle
len1=length(indNonZeroTract)
x1_manual=0.25*(rank(CDC_SVIr$EP_POV[indNonZeroTract])/len1+rank(CDC_SVIr$EP_UNEMP[indNonZeroTract])/len1+
  1-rank(CDC_SVIr$E_PCI[indNonZeroTract])/len1+rank(CDC_SVIr$E_NOHSDP[indNonZeroTract])/len1)
x2_manual=0.25*(rank(CDC_SVIr$EP_AGE65[indNonZeroTract])/len1+rank(CDC_SVIr$EP_AGE17[indNonZeroTract])/len1+
                  rank(CDC_SVIr$EP_DISABL[indNonZeroTract])/len1+rank(CDC_SVIr$EP_SNGPNT[indNonZeroTract])/len1)
x3_manual=0.5*(rank(CDC_SVIr$EP_MINRTY[indNonZeroTract])/len1+rank(CDC_SVIr$EP_LIMENG[indNonZeroTract])/len1)
x4_manual=0.2*(rank(CDC_SVIr$EP_MUNIT[indNonZeroTract])/len1+rank(CDC_SVIr$EP_MOBILE[indNonZeroTract])/len1+
                  rank(CDC_SVIr$EP_CROWD[indNonZeroTract])/len1+rank(CDC_SVIr$EP_NOVEH[indNonZeroTract])/len1+
                 rank(CDC_SVIr$EP_GROUPQ[indNonZeroTract])/len1 )
tractSVI_Pan=cbind(x1_manual,x2_manual,x3_manual,x4_manual)
# weekWindow=seq(3,5)  # weeks 3,4,5
weekWindow=3           # weeks 3. estimate SAR from 1 week of data. replications don't work with SAR
themeNo=2
x1=tractSVI[indNonZeroTract,themeNo] #to use the published indicator values - ranked for entire state
#x1=tractSVI_Pan[,themeNo]             #to use manually calculated indicator values - ranked for only panhandle
xMatrix=x1 
yMat=as.matrix(y_weekly[,weekWindow])
nn=length(weekWindow)
# one zip code has a NA SVI value
x=rep(x1,nn)
# coordinates
X02=X0
# make it a one column matrix
y=cbind(c(yMat))
plot(y~x,xlab='Vulnerability Indicator',ylab='Transformed Covid counts in Week 3 ')
mod1=lm(y~x) # linear model
abline(mod1,col='red')
summary(mod1)
# Plot Covid count wk 3 vs SVI themes
# use the published indicator values
dataTable=rbind(data.frame(count=y,VulnIndic=tractSVI[indNonZeroTract,1],theme=1),# to use manual indicators, use VulnIndic=tractSVI_Pan[,1]
               data.frame(count=y,VulnIndic=tractSVI[indNonZeroTract,2],theme=2),
               data.frame(count=y,VulnIndic=tractSVI[indNonZeroTract,3],theme=3),
               data.frame(count=y,VulnIndic=tractSVI[indNonZeroTract,4],theme=4))
dataTable$theme=factor(dataTable$theme)
if(flagCreatePNG){png("CovidCount_Vuln3.png", width = 7, height =5, units = 'in', res = 600)}
par(mfrow=c(1,1),pty="s",mar=c(3,3,0.5,0)+0.1,mgp=c(2,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
xyplot(count ~ VulnIndic|theme, dataTable,
       xlab = 'Vulnerability Indicator',
       ylab = "Transformed Covid Counts",
       type=c('p','r'),
       between = list(x = 0.1),
       #auto.key = list(columns = 5,title=expression(rho)),cex=0.5,
       strip = strip.custom(strip.names = TRUE, var.name = "Theme"),
       index.cond=list(c(3,4,1,2)))

if(flagCreatePNG){dev.off()}


# fit spatial error model between average covid count from weeks 3, 4, 5 and RPL_THEME2 or RPL_THEME3
# (see SEM model in Ch 16, https://keen-swartz-3146c4.netlify.app/)
nb4rt <- poly2nb(fl_tracts)
W11=nb2mat(nb4rt,zero.policy = TRUE)
W12=W11[indNonZeroTract,indNonZeroTract]

# # SAR model - use manual fit
# Wn=W12; Wnew=NA;i=1
# while (i<nn){Wnew=bdiag(Wn,W12); Wn=Wnew; i=i+1}
# sol <- optimize(getLikelihoodSAR, interval=c(-1, 1),y=y,W=as.matrix(Wn),X=x)
# rho_hat=sol$minimum; B=rho_hat*Wn
# m=length(y); A=diag(m)-B
# X=matrix(c(rep(1,m),x),m,2)
# beta_hat=solve(t(A%*%X)%*%(A%*%X),t(A%*%X)%*%(A%*%y)) 
# yhat=X%*%beta_hat; e_vec=y-yhat
# sig2_hat=(1/m)*t(e_vec)%*%(t(A)%*%A)%*%e_vec

# SAR model - use errorsarlm
nb4rt2 <- poly2nb(fl_tracts[indNonZeroTract,])
weightMat=nb2listw(nb4rt2, zero.policy=TRUE)
y_avg=rowMeans(yMat)                      # replications don't work with SAR
sar_mg=errorsarlm(y_avg~1+x1, listw=weightMat,na.action =na.omit,zero.policy=TRUE)  #risk & correlation adjusted
#sar_mg=errorsarlm(y_avg~1, listw=weightMat,na.action =na.omit,zero.policy=TRUE)    #non risk adjusted; only correlation adjusted
summary(sar_mg)
yhat=sar_mg$fitted.values
par(mfrow=c(1,1),mar=c(3,3,0.5,0)+0.1,mgp=c(2,1,0),cex.lab=1.2,cex.main=1,cex.axis=1.2,cex=1)
plot(y_avg~yhat,ylab='observed',xlab='predicted')
abline(0,1,col='red')

e_vec=y_avg-yhat
# moran test on adjusted data
(moran.test(as.numeric(e_vec), listw=weightMat,zero.policy=TRUE)$p.value)
# moran test on raw data
(moran.test(y_avg, listw=weightMat,zero.policy=TRUE)$p.value)
# R-sq of the model- also try cor(y_avg,yhat)^2
(Rsq=1-sar_mg$SSE/(var(y_avg)*(ns2-1)))
#### find ARL0 - risk adjusted ####

# # use errorsarlm fit
b0=as.numeric(sar_mg$coefficients[1])
b1=as.numeric(sar_mg$coefficients[2])
rho=as.numeric(sar_mg$lambda)
sd_err=sqrt(sar_mg$s2)

n_MC=10
delt=0        # in control ARL
ind_inject=1  # immaterial for ARL0
flag_plot=FALSE
flag_moran=FALSE
time_inj=1
step_size=0.25
r=seq(0.05,2*step_size+0.05,by=step_size)
nr=length(r)
# cusum parameter: 3 sigma shift(sustained shift CUSUM)
delt_c=sd_err*3 
if (flagTuneMethods){
  h_s=4.85
  # alarm limits for 3 sigma shift, same for all Themes 
  # h=4.85 , ARL0=13.700000  (3.296631) for T1, T3 and T4, ARL0=20.900000  (4.552045) for T2 
  # with manual indicator calculation
  # h=4.85 , ARL0=8.600000,20.900000,  for T1,T2 
  set.seed(98765)
  flagRiskCorrAdjust=TRUE
  res=get_ARL_SCUSUMCorrel_NeighborCovarGiven(r, X02, b0,b1,rho,sd_err,l0,delt_c, delt, h_s, 
                                              time_inj, n_MC,ind_inject,flag_plot,flag_moran,nb4rt2,xMatrix,flagRiskCorrAdjust)
  
  print(c("limit, delta,ARL,seRL"))
  print(c(h_s,delt,res$ARL,res$seRL))
}
#### CUSUM surveillance of infection data- risk or correlation adjusted ####
S_prev=matrix(0,nr,ns2)
tau=matrix(0,nr,ns2)
pvalMoranAdj=NA
pvalMoran=NA
flag_RiskAdjust=TRUE; xL=NA;xU=NA;
nnOL=1
S_max=0
Smax_vec=0
wk_begin=3
wk_end=11
wk=wk_begin
r=seq(0.05,2*step_size+0.05,by=step_size)
h_s=4.85 
# storage matrix
yAdjMat=matrix(NA,ns2,1)
while (S_max<h_s&wk<wk_end){
  # risk and correlation adjusted cusum
  y=y_weekly[,wk]
  yAdj=adjustForCorrelSAR(ns2,rho,nnOL,b0,b1,sd_err,y,xMatrix,W12,xL,xU,flag_plot,flag_RiskAdjust)$yAdj
  #yAdj=adjustForCorrelSAR(ns2,rho,nnOL,b0,0,sd_err,y,rep(0,len1),W12,xL,xU,flag_plot,flag_RiskAdjust)$yAdj # non risk adjusted
  yAdjMat=cbind(yAdjMat,yAdj)
  res_S=get_spaceTimeCUSUMNormalStandard(yAdj,X02,S_prev,delt_c/2,r) 
  
  S_prev=res_S$S_new
  tau[which(S_prev==0)]=wk 
  # update cusum statistic
  S_max=res_S$S_max_t
  Smax_vec=c(Smax_vec,S_max)
  print(wk)
  wk=wk+1
}
# change point estimator
ind_max=which(res_S$S_new==S_max)
tau_est=tau[ind_max[1]]
print(paste('change point=',tau_est,sep=""))
# plot CUSUM
if(flagCreatePNG){png("cusumAdj3.png", width = 7, height =5, units = 'in', res = 600)}
par(mfrow=c(1,1),mar=c(3,3,0.5,0)+0.1,mgp=c(2,1,0),cex.lab=1,cex.main=1,cex.axis=1.2,cex=1)
plot(seq(wk_begin,wk),Smax_vec,type='b',xlab='Week',ylab='CUSUM',xaxt='none')
axis(1,at=seq(wk_begin,wk),labels=format(uniqueWeeksMonday[wk_begin:wk],"%m/%d"))
abline(h=h_s,lty=2)
if(flagCreatePNG){dev.off()}
# plot identified cluster - zero pop ZCs removed 
clust_center=as.character(fl_tracts$GEOID[indNonZeroTract[res_S$region_max]])
print(paste("cluster center is at tract",clust_center,sep=""))
rad_cluster=r[res_S$radius_max]
# second and third largest clusters
reg_max2=res_S$RadiusRegionTop3[2,2]; reg_max3=res_S$RadiusRegionTop3[2,3]
rad_max2=r[res_S$RadiusRegionTop3[1,2]]; rad_max3=r[res_S$RadiusRegionTop3[1,3]]
# find tracts included in cluster- zero pop ZCs not removed
ind_InCluster=which(ODmatrix[indNonZeroTract[res_S$region_max],]<=rad_cluster)
# find tracts included in the second and third largest clusters
ind_InCluster2=which(ODmatrix[indNonZeroTract[reg_max2],]<=rad_max2)
ind_InCluster3=which(ODmatrix[indNonZeroTract[reg_max3],]<=rad_max3)
colorTracts=rep("NA",ns)
colorTracts[indNonZeroTract]='navajowhite2'
colorTracts[ind_InCluster2]="red"  # comment to not show second and third largest clusters
colorTracts[ind_InCluster3]="blue" # comment to not show second and third largest clusters
colorTracts[ind_InCluster]="seagreen4"
if(flagCreatePNG){png("ClusterAdj0.png", width = 7, height =5, units = 'in', res = 600)}
plot(panhandleTracts,border='grey',main= "")# paste('Clusters identified in wk ',wk))
plot(fl_tracts,border=NA,add=TRUE,col=alpha(colorTracts,0.5),
     axes = FALSE)
text(coordinates(fl_tracts)[indNonZeroTract[res_S$region_max],1]+0.1,coordinates(fl_tracts)[indNonZeroTract[res_S$region_max],2]+0.1,"1")
text(coordinates(fl_tracts)[indNonZeroTract[reg_max2],1]+0.1,coordinates(fl_tracts)[indNonZeroTract[reg_max2],2],"2")
text(coordinates(fl_tracts)[indNonZeroTract[reg_max3],1]+0.1,coordinates(fl_tracts)[indNonZeroTract[reg_max3],2]-.1,"3")
if(flagCreatePNG){dev.off()}

