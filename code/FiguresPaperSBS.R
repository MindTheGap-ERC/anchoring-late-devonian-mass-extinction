library(astrochron)
library(tidyverse)
library(ggpubr)
setwd('D:/Users/Anne-Christine/Documents/ANNECHRI2/POST/UTRECHT/Spectral ANALYSIS/SteinBruchSchmidt/Final_End2018') #Set Working Directory
setwd('I:/ANNECHRI2/POST/UTRECHT/Spectral ANALYSIS/SteinBruchSchmidt/Final_End2018') #Set Working Directory
shell.exec(getwd())

# Figure 2 - Vertical plot of data -----------------------------------------------------------------
data=read.csv("SbS_XRF_forfactor3.csv",header = T, sep=";")
Height=data[,c(1)]

MS<-data[,c(1,23)]

Ca<-data[,c(1,9)]

Ti<-data[, c(1, 10)]
Til=logT(Ti, opt=2)
colnames(Til) <- c("Altitude","LogTi") # put new column names

C<-data[,c(1,21)] 
C<-C[complete.cases(C),] # we remove the lines with no data 

plotMS<-ggplot(MS, aes(Altitude, MS)) + 
  geom_line(color='darkgreen') 
plotMS<-plotMS + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotLogTi<-ggplot(Ti, aes(Altitude, LogTi)) + 
  geom_line(color='red3') + theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.title.y = element_text("LogTi"))
plotC<-ggplot(C, aes(Altitude, d13C)) + 
  geom_line(color='slateblue4')   + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggarrange(plotMS, plotLogTi, plotC, ncol=1, nrow=3, align='v')
library(plotly)
plot_ly(Ca, x = ~Altitude, y = ~Ca, type = 'scatter', mode = 'lines')

# eTimeOpt Results--------------------------

data=read.csv("SbS_XRF_forfactor3.csv",header = T, sep=";")
Height=data[,c(1)]
mydata=data[,c(2:23)]
mydata[is.na(mydata)] <- 0
MS<-data[,c(1,23)]
MS=na.omit(MS)
MSi=linterp(MS, 0.02)

Ti<-data[, c(1, 10)]
Til=logT(Ti, opt=2)
Tili=linterp(Til) # we create a LogTi record interpolated at 0.02
Tii=linterp(Ti, 0.02) # we create a Ti record interpolated at 0.02

C<-data[,c(1,21)] 
C<-C[complete.cases(C),] # we remove the lines with no data 
Ci=linterp(C, 0.02)

# we set up the target for the estimation of amplitude modulation 
targetE=c(130.7, 123.8, 98.9, 94.9)
bergerPeriods(372)
targetP=c(19.98, 16.87)

# Amplitude modulation of precession by short eccentricity 
etimeOptLogTi1=eTimeOpt(Tili,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.55, numsed=100,linLog=1,
                        limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                        detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)
etimeOptTi1=eTimeOpt(Tii,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
                        limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                        detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)
etimeOptMS1=eTimeOpt(MSi,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
                     limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                     detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)
etimeOptC1=eTimeOpt(Ci,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
                    limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                    detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)

# Amplitude modulation of short eccentricity by long eccentricity 
etimeOptLogTi2=eTimeOpt(Tili,win=0.02*100,step=0.02*10,sedmin=0.01,sedmax=3,  numsed=100,linLog=1,
                      limit=T,fit=2,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                      detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)
etimeOptTi2=eTimeOpt(Tii,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6,  numsed=100,linLog=1,
                      limit=T,fit=2,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                      detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)
etimeOptMS2=eTimeOpt(MSi,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
                     limit=T,fit=2,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                     detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)
etimeOptC2=eTimeOpt(Ci,win=0.02*100,step=0.02*10,sedmin=0.1,sedmax=0.6, numsed=100,linLog=1,
                    limit=T,fit=2,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                    detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)

# extract the optimal fits for the envelope*power optimization
sedratesLogTi1=eTimeOptTrack(etimeOptLogTi1[3])
sedratesTi1=eTimeOptTrack(etimeOptTi1[3])
sedratesMS1=eTimeOptTrack(etimeOptMS1[3])
sedratesC1=eTimeOptTrack(etimeOptC1[3])

sedratesLogTi2=eTimeOptTrack(etimeOptLogTi2[3])
sedratesTi2=eTimeOptTrack(etimeOptTi2[3])
sedratesMS2=eTimeOptTrack(etimeOptMS2[3])
sedratesC2=eTimeOptTrack(etimeOptC2[3])



# We tune the different time series 
timeLogTi1=sedrate2time(sedratesLogTi1)
sedrateTunedLogTi1=tune(Tili, timeLogTi1, extrapolate=T)

timeTi1=sedrate2time(sedratesTi1)
sedrateTunedTi1=tune(Tii, timeTi1, extrapolate=T)

timeMS1=sedrate2time(sedratesMS1)
sedrateTunedMS1=tune(MS, timeMS1, extrapolate=T)

timeC1=sedrate2time(sedratesC1)
sedrateTunedC1=tune(Ci, timeC1, extrapolate=T)

timeLogTi2=sedrate2time(sedratesLogTi2)
sedrateTunedLogTi2=tune(Tili, timeLogTi2, extrapolate=T)

timeTi2=sedrate2time(sedratesTi2)
sedrateTunedTi2=tune(Tii, timeTi2, extrapolate=T)

timeMS2=sedrate2time(sedratesMS2)
sedrateTunedMS2=tune(MS, timeMS2, extrapolate=T)

timeC2=sedrate2time(sedratesC2)
sedrateTunedC2=tune(C, timeC2, extrapolate=T)

# Tuning to the of ash bed
ash=idPts(SBS_MS_Tune)
SBS_MS_eTimeOpt_anchored=anchorTime(sedrateTunedMS1,time=ash$x,age=372360, timeDir=2)



# ORTA Results--------------------------

# First we run the monte Carlo Simulation with the Tie Points we d 
# Start from "rough" tie-point ages. These are the ages listed in BLACK on Figures 5 and 6
Pointers_original=c(2750,2530,2430,2280,2130,1930,1700,1315,860,700,500,400,300,150,100,-20,-180,-450)


#The time difference between the tie-point ages in the line above
offset=c(200,100,150,150,200,230,385,455,160,200,100,100,150,50,120,160,270)

# Initiate some variables
astropower_sum=c()
Pointers_sum_best=c()
astropower_sum_best=1
offset_matrix=matrix(0,10000,17)
offset_thisloop=c()
start_i=1

# we disturb the time difference between two consecutibe tie points with 2 percent 
sd_disturbance=0.03

# we create  the "rough" tie-point ages for different sections stratigraphic levels
meter=c(0.5, 1.5, 3.1, 3.9, 4.5)
age=c(700, 500, 300, 100, -20)
Pointers_SBS<-cbind(meter, age)

# we read the file corresponding to the different records in the distance domain 
SBS_d13C_carb=read.csv("C_SBS_withoutOutliers.csv",header=T, sep=";")
SBS_MS=linterp(SBS_d13C_carb)

SBS_MS=read.csv("SbS_MS.csv",header=T, sep=";")
SBS_MS=linterp(SBS_MS)

# Monte Carlo loop

# what is the offset ?
#the loop is done 1000 - it takes the pointers as in the files and we use these pointers to caclculate the spetrum
# then we disturb a bit and check the spectrum and then if better we keep going into this direction 
# if not, we trash it. So from the starting point, it goes into one dir.
# with these numbers eevery 400 we force it to go back to 0 and to try the other dir 

# CHECK POINTERS USED
for (i in 1:5000) {
  if (i %in% c(401,801,1201,1601,2001,2401,2801,3201,3601,4001,4401,4801,5201,5601,6001,6401,6801,7201,7601,8001,8401,8801,9201,9601)==TRUE) {
    start_i=i
    offset=c(220,110,200,100,270,250,300,450,150,200,200,100,130,60,80,180,270)
  }
  Pointers=c()
  Pointers[1]=Pointers_original[1]
  for (j in 1:17) {
    Pointers[j+1]=Pointers[j]-offset[j]*rnorm(1,1,sd_disturbance) # Disturb the time-differences between consecutive tie-points
  }
  
  
  # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  colnames(Pointers_SBS)=c('Height',"Relative Age (kyr)")
  # with pointers 2 and 1
  # Pointers_SBS[,2]=Pointers[c(10:16)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  #Pointers_SBS[,2]=Pointers[c(10:16)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  Pointers_SBS[,2]=Pointers[c(10, 11, 13, 15, 16)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  
  # Calibrate the proxy records using the "tune" function in astrochron
  # SBS
  SBS_MS_T1=tune(SBS_MS,Pointers_SBS,extrapolate=T, genplot=F) # Convert depth to time according to the "disturbed" tie-point ages
  
  # SBS_MS
  SBS_MS_T1_i=linterp(SBS_MS_T1, 4, genplot=F)
  SBS_MS_T1_i=detrend(SBS_MS_T1_i, genplot=F)
  SBS_MS_T1_MTM=mtm(SBS_MS_T1_i,tbw=2,padfac = 5, genplot=F,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07) # Calculate MTM power spectrum for Section C, according to the "disturbed" tie-point ages
  
  idx1=which(SBS_MS_T1_MTM$Frequency>0.008 & SBS_MS_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(SBS_MS_T1_MTM$Frequency>0.025 & SBS_MS_T1_MTM$Frequency<0.035) # Obliquity band
  idx3=which(SBS_MS_T1_MTM$Frequency>0.045 & SBS_MS_T1_MTM$Frequency<0.065) # Precession band
  
  SBS_MS_100peak_idx=which(SBS_MS_T1_MTM$Power==max(SBS_MS_T1_MTM$Power[idx1]))
  SBS_MS_100_misfit=abs(SBS_MS_T1_MTM$Frequency[SBS_MS_100peak_idx]-0.0105)/0.0105
  SBS_MS_100_height=(100-SBS_MS_T1_MTM$AR1_CL[SBS_MS_100peak_idx])/100
  SBS_MS_100peak_Fidx=which(SBS_MS_T1_MTM$Harmonic_CL==max(SBS_MS_T1_MTM$Harmonic_CL[idx1]))
  SBS_MS_100_Fmisfit=abs(SBS_MS_T1_MTM$Frequency[SBS_MS_100peak_Fidx]-0.0105)/0.0105
  SBS_MS_100_Fheight=(100-SBS_MS_T1_MTM$Harmonic_CL[SBS_MS_100peak_Fidx])/100
  
  SBS_MS_obpeak_idx=which(SBS_MS_T1_MTM$Power==max(SBS_MS_T1_MTM$Power[idx2]))
  SBS_MS_ob_misfit=abs(SBS_MS_T1_MTM$Frequency[SBS_MS_obpeak_idx]-0.031)/0.031
  SBS_MS_ob_height=(100-SBS_MS_T1_MTM$AR1_CL[SBS_MS_obpeak_idx])/100
  SBS_MS_obpeak_Fidx=which(SBS_MS_T1_MTM$Harmonic_CL==max(SBS_MS_T1_MTM$Harmonic_CL[idx2]))
  SBS_MS_ob_Fmisfit=abs(SBS_MS_T1_MTM$Frequency[SBS_MS_obpeak_Fidx]-0.031)/0.031
  SBS_MS_ob_Fheight=(100-SBS_MS_T1_MTM$Harmonic_CL[SBS_MS_obpeak_Fidx])/100
  
  SBS_MS_P1peak_idx=which(SBS_MS_T1_MTM$Power==max(SBS_MS_T1_MTM$Power[idx3]))
  SBS_MS_P1_misfit=abs(SBS_MS_T1_MTM$Frequency[SBS_MS_P1peak_idx]-0.055)/0.055
  SBS_MS_P1_height=(100-SBS_MS_T1_MTM$AR1_CL[SBS_MS_P1peak_idx])/100
  SBS_MS_P1peak_Fidx=which(SBS_MS_T1_MTM$Harmonic_CL==max(SBS_MS_T1_MTM$Harmonic_CL[idx3]))
  SBS_MS_P1_Fmisfit=abs(SBS_MS_T1_MTM$Frequency[SBS_MS_P1peak_Fidx]-0.055)/0.055
  SBS_MS_P1_Fheight=(100-SBS_MS_T1_MTM$Harmonic_CL[SBS_MS_P1peak_Fidx])/100
  
  SBS_MS_misfit=mean(c(SBS_MS_100_misfit,SBS_MS_ob_misfit,SBS_MS_P1_misfit,SBS_MS_100_height,SBS_MS_ob_height,SBS_MS_P1_height))
  
  
  # End average of "tied-in" sections (H32, CG1, Sinsin and Fuhe + SBS)
  astropower_sum[i]=mean(c(SBS_MS_misfit))
  for (k in 1:17) {offset_thisloop[k]=Pointers[k]-Pointers[k+1]}
  offset_matrix[i,c(1:17)]=offset_thisloop
  sd_disturbance=astropower_sum[i]^1.5
  
  if (astropower_sum[i]<astropower_sum_best){
    astropower_sum_best=astropower_sum[i]
    Pointers_sum_best=Pointers
    save.image(file = "Sum_best.Rdata")  }
  
  if (astropower_sum[i]<=min(astropower_sum[c(start_i:i)])){
    for (k in 1:17) {offset[k]=Pointers[k]-Pointers[k+1]}
  }
  
}



# Calibrate the proxy records using the New tie points
# if we have to go back again to the pointers
Pointers_optimal=read.csv("Pointers_Best_Sum_5000_Pointers_SBS_27_2.csv",header=T, sep=";")
Pointers_optimal=c(2750, 2517.561387, 2393.783308, 2204.039405, 2104.823618, 1807.241569, 1574.014651, 1232.296504, 885.6331306, 719.8978336, 496.4931781, 304.5426969, 222.3660329, 92.35163494, 34.46326485, -75.05072288, -263.5613211, 463.7788785)

meter=c(0.5, 1.5, 3.1, 3.9, 4.5)
age=c(700, 500, 300, 100, -20)
Pointers_SBS<-cbind(meter, age)
colnames(Pointers_SBS)=c('Height',"Relative Age (kyr)")
Pointers_SBS[,2]=Pointers_optimal[c(10, 11, 13, 15, 16)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages



data=read.csv("SbS_XRF_forfactor3.csv",header = T, sep=";")
Height=data[,c(1)]

MS<-data[,c(1,23)]
MS=na.omit(MS)

C<-data[,c(1,21)] 
Ti<-data[, c(1, 10)]
Til=logT(Ti, opt=2)
colnames(Til) <- c("Altitude","LogTi") # put new column names

SBS_MS_Tune=tune(MS, Pointers_SBS,extrapolate=T) # Convert depth to time according to the "disturbed" tie-point ages
colnames(SBS_MS_Tune) <- c("kyr","MS") # put new column names
SBS_C_Tune=tune(C, Pointers_SBS,extrapolate=T) 
SBS_Ti_Tune=tune(Til,Pointers_SBS,extrapolate=T) # Convert depth to time according to the "disturbed" tie-point ages

# position of ash bed
ash=idPts(SBS_MS_Tune)
# we create the new tuned time scale 
SBS_MS_anchored=anchorTime(SBS_MS_Tune,time=ash$x,age=372360, timeDir=2)
SBS_C_anchored=anchorTime(SBS_C_Tune,time=ash$x,age=372360, timeDir=2)
SBS_Ti_anchored=anchorTime(SBS_Ti_Tune,time=ash$x,age=372360, timeDir=2)



# Figure 3 Wavelet and sedimentation rate evolution-----------------------------------------------------------------
data=read.csv("SbS_XRF_forfactor3.csv",header = T, sep=";")
Height=data[,c(1)]
MS<-data[,c(1,23)]
MS<-MS[complete.cases(MS),] # we remove the lines with no data 
MSi=linterp(MS, 0.02)
MSL=noLow(MSi, smooth=0.2)
wtLS3<-wt(MSL, dt=1, dj=1/12)
dev.off()
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot.biwavelet(wtLS3, plot.sig = FALSE, type="power.corr.norm", main="Bias-corrected", plot.cb=TRUE, zlim=c(1,10), lwd.coi=1, col.coi = 'black')

# eha  
eha_SBS <- eha(MSi, detrend=T, win=0.60, step=0.015, output=4,pl=1,pad=5000,genplot=3,fmax=20, 
               xlab="Frequency (cycles/ka)",ylab="Age (ka)")

# Implementation of Sedimentation rate  from ORTA 
# we use the sedimentation rate extracted from orta:
SBS_MS_T1=read.csv("SBS_MS_T1_Pointers27_1.csv",header = T, sep=",")
SBS_MS=read.csv("SbS_MS.csv",header = T, sep=";")
SBS_MS_T1$X1 <- SBS_MS_T1$X1+334.5

library(dplyr)

SBS_MST=cbind(SBS_MS, SBS_MS_T1$X1)
colnames(SBS_MST) <- c('Alt_m', 'MS', 'Alt_kyr')

# Instantaneous sedimentation rate at each point 
SBS_MST2 <- SBS_MST %>%
  mutate(diff_m = Alt_m - lag(Alt_m, default = first(Alt_m)))
SBS_MST3 <- SBS_MST2 %>%
  mutate(diff_kyr = -(Alt_kyr - lag(Alt_kyr, default = first(Alt_kyr))))
SBS_MST4 <- SBS_MST3 %>%
  mutate(SRInst = diff_m/diff_kyr*100)
SBS_MST5=SBS_MST4[-1,] # remove first line (no previous sample, it's NA)

# Implementation of the sedimentation rate from TimeOpt on the MS signal 
etimeOptMS1=eTimeOpt(MSL,win=0.02*100,step=0.02*10,sedmin=0.15,sedmax=0.55, numsed=100,linLog=1,
                     limit=T,fit=1,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                     detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)
etimeOptMS2=eTimeOpt(MSi,win=0.02*100,step=0.02*10,sedmin=0.15,sedmax=0.55, numsed=100,linLog=1,
                     limit=T,fit=2,fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=targetE,targetP=targetP,
                     detrend=T,ydir=1,output=1,genplot=T,check=T,verbose=1)

# extract the optimal fits for the envelope*power optimization
sedratesMS1=eTimeOptTrack(etimeOptMS1[3])
sedratesMS2=eTimeOptTrack(etimeOptMS2[3])


plot(SBS_MST5$Alt_m, SBS_MST5$SRInst, type='l', ylim=c(0.3, 0.6))  
lines(sedratesMS1$Location, sedratesMS1$Sedrate, col='red')
lines(sedratesMS2$Location, sedratesMS2$Sedrate, col='darkgreen')



# Figure 6 Frasnian-Famennian Time Scale -----------------------------------------------------------------

colnames(SBS_MS_anchored) <- c("time","MS")
colnames(SBS_C_anchored) <- c("time","C")
colnames(SBS_Ti_anchored) <- c("time","Ti")

SBS_MS_anchoredi<-linterp(SBS_MS_anchored)
filterMS100 <- bandpass(SBS_MS_anchoredi, flow=1/80, fhigh=1/140)
filterMS400 <- bandpass(SBS_MS_anchoredi, flow=1/450, fhigh=1/350)
SBS_C_anchoredi<-linterp(SBS_C_anchored)
filterC100 <- bandpass(SBS_C_anchoredi, flow=1/80, fhigh=1/150)
filterC400 <- bandpass(SBS_C_anchoredi, flow=1/450, fhigh=1/350)
SBS_Ti_anchoredi<-linterp(SBS_Ti_anchored)
filterTi100 <- bandpass(SBS_Ti_anchoredi, flow=1/80, fhigh=1/150)
filterTi400 <- bandpass(SBS_Ti_anchoredi, flow=1/450, fhigh=1/350)

plot(SBS_MS_anchored, type='l', col='darkgreen')
lines(filterMS100$time, filterMS100$MS, col='darkolivegreen2', type='l')
lines(filterMS400$time, filterMS400$MS, col='seagreen2', type='l')

plot(SBS_C_anchored, type='l', col='darkgreen')
lines(filterC100$time, filterC100$C, col='darkolivegreen2', type='l')
lines(filterC400$time, filterC400$C, col='seagreen2', type='l')

plot(SBS_Ti_anchored, type='l', col='darkgreen')
plot(filterTi100$time, filterTi100$Ti, col='darkolivegreen2', type='l')
lines(filterTi400$time, filterTi400$Ti, col='seagreen2', type='l')


# Figure 7 Spectral Power -----------------------------------------------------------------
SBS_MS_anchored=read.csv("SBS_MS_anchored_Pointers_27_2.csv",header = T, sep=";")
SBS_Ti_anchored=read.csv("SBS_TiLog_anchored_Pointers_27_2.csv",header = T, sep=",")
SBS_C_anchored=read.csv("SBS_d13C_T1_anchored_Pointers_27_2.csv",header = T, sep=";")
SBS_C_anchored=read.csv("SBS_C_anchored_eTimeOpt.csv", sep=",",header = T)


SBS_C_nl=noLow(SBS_C_anchored, smooth= 0.2)
SBS_MS_nl=noLow(SBS_MS_anchored, smooth=0.2)
SBS_Ti_nl=noLow(SBS_Ti_anchored, smooth=0.2)

SBS_Clog=logT(SBS_C_nl, opt=2, c=2)
SBS_MSlog=logT(SBS_MS_nl, opt=2, c=2)
SBS_Tilog=logT(SBS_Ti_nl, opt=2, c=2)

SBS_C_nli=linterp(SBS_Clog, dt=3)
SBS_MS_nli=linterp(SBS_MSlog, dt=3)
SBS_Ti_nli=linterp(SBS_Tilog, dt=3)

pwrC=eha(SBS_C_nli, win=150, step=2, output=2,pl=1,pad=5000,genplot=3,
        xlab="Frequency (cycles/ka)",ylab="Age (ka)", fmax=0.1);
pwrMS=eha(SBS_MS_nli, win=150, step=2, output=2,pl=1,pad=5000,genplot=3,
            xlab="Frequency (cycles/ka)",ylab="Age (ka)", fmax=0.1);
pwrTi=eha(SBS_Ti_nli, win=150, step=2, output=2,pl=1,pad=5000,genplot=3,
           xlab="Frequency (cycles/ka)",ylab="Age (ka)", fmax=0.1);

integrate_obl_C=integratePower(pwrC,flow = 0.025, fhigh = 0.036, npts=366,pad=5000,ln=T)
integrate_eccS_C=integratePower(pwrC,flow = 1/150, fhigh = 1/80, npts=346,pad=5000,ln=T)

integrate_obl_MS=integratePower(pwrMS,flow = 0.025, fhigh = 0.048, npts=366,pad=5000,ln=T)
integrate_eccS_MS=integratePower(pwrMS,flow = 1/160, fhigh = 1/70,  npts=346,pad=5000,ln=T)

integrate_obl_Ti=integratePower(pwrTi,flow = 0.025, fhigh = 0.048, npts=366,pad=5000,ln=T)
integrate_eccS_Ti=integratePower(pwrTi,flow = 1/160, fhigh = 1/70,  npts=346,pad=5000,ln=T)


