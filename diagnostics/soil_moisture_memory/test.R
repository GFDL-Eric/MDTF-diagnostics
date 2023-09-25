################################
source("lib_LoCo.R")

name.user  <- Sys.getenv("USER")
#air specific humidity profile [kg.kg-1]
name.mrsos  <- 'mrsos'

path.i  <- paste0('/work/',name.user,'/zun_mdtf/inputdata/model/GFDL.CM4.c96L32.amip/mon/')

########################################
#print("Taking specific humidity profile from model")
#filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.mrsos.mon.nc')
#data  <- nc_open(filei)
#lons.i <- ncvar_get(data, "lon")
#lats.i <- ncvar_get(data, "lat")
#sm.i  <- ncvar_get(data, name.mrsos)
#nc_close(data) 

#if (max(lons.i) > 250)
#{
#    lons.i  <-  lons.i - 180
#    sm.i    <-  ShiftLon(sm.i)
#}

##########################################
#threshold to have significant correlation
thr.sm  <-  1/exp(1)
#threshold of minimum length of time series
thr.len <-  20

##########################################
#print("Auto-correlation of soil moisture")
#lag.sm  <-  array(NA,dim=dim(sm.i)[1:2])
#for (nx in 1:dim(sm.i)[1])#lon
#for (ny in 1:dim(sm.i)[2])#lat
#{
#    ts.tmp  <-  sm.i[nx,ny,]
#    #check the length
#    if (is.na(mean(ts.tmp)) | length(ts.tmp) < thr.len)
#        next
#
#    acf.tmp <-  acf(ts.tmp,lag.max=length(ts.tmp)-1,plot=F)$acf
#    lag.sm[nx,ny] <-  min(seq(1,length(acf.tmp))[acf.tmp < thr.sm])
#}

##########################################
print("Plotting soil moisture memory")
png(paste0("SMM.png",sep=""), width=700)

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
print("Normal End of Convect Triggering Potential.R")
############################################################
