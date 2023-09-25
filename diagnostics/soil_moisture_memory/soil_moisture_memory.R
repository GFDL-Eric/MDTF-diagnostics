################################
# .libPaths("/home/water2/ab5/R/x86_64-redhat-linux-gnu-library/3.2")
library(colorRamps)
library(maps)
library(fields)
library(ncdf4)

POD_HOME  <- Sys.getenv("POD_HOME")
WK_DIR <- Sys.getenv("WK_DIR")
OBS_DATA <- Sys.getenv("OBS_DATA")
DATADIR <- Sys.getenv("DATADIR")
CASENAME <- Sys.getenv("CASENAME")
yr1 <- Sys.getenv("FIRSTYR")
yr2 <- Sys.getenv("LASTYR")
source(paste0(POD_HOME,"/lib_LoCo.R"))

MRSOS_FILE <- Sys.getenv("MRSOS_FILE")

#soil moisture content [kg.m-2]
name.mrsos  <- Sys.getenv("mrsos_var")

##########################################
print("Taking soil moisture from model")
print(MRSOS_FILE)
data_mrsos <- nc_open(MRSOS_FILE)
lons.i  <- ncvar_get(data_mrsos, "lon")
lats.i  <- ncvar_get(data_mrsos, "lat")
#[kg.m-2]
sm.i  <- ncvar_get(data_mrsos, name.mrsos)
nc_close(data_mrsos) 

print("PP: checking the longitude system.")
#check if lons ranges from [0,360]
if (max(lons.i) > 250)
{
    lons.i  <-  lons.i - 180
    sm.i    <-  ShiftLon(sm.i)
}

##########################################
#threshold to have significant correlation
thr.sm  <-  1/exp(1)
#threshold of minimum length of time series
thr.len <-  20

##########################################
print("Auto-correlation of soil moisture")
lag.sm  <-  array(NA,dim=dim(sm.i)[1:2])
for (nx in 1:dim(sm.i)[1])#lon
for (ny in 1:dim(sm.i)[2])#lat
{
    ts.tmp  <-  sm.i[nx,ny,]
    #check the length
    if (is.na(mean(ts.tmp)) | length(ts.tmp) < thr.len)
        next

    acf.tmp <-  acf(ts.tmp,lag.max=length(ts.tmp)-1,plot=F)$acf
    lag.sm[nx,ny] <-  min(seq(1,length(acf.tmp))[acf.tmp < thr.sm])
}

##########################################
print("Plotting soil moisture memory")
png(paste0(WK_DIR,"/model/SMM.png",sep=""), width=700)

lag.sm[!is.na(lag.sm) & (lag.sm > 11)] <-  11
lag.sm[!is.na(lag.sm) & (lag.sm <= 0)] <-  NA

col.br  <-  seq(0,11,1)

image.plot(lons.i,lats.i,lag.sm,zlim=range(col.br),
           breaks=col.br, col=tim.colors(11),
           xlab="", ylab="",
           main=paste0("Soil Moisture Memory (months)"))
plot(coastsCoarse,add=T)

dev.off()

############################################################
print("Normal End of Soil_Moisture_Memory.R")
############################################################
