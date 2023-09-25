################################
source("lib_LoCo.R")

name.user  <- Sys.getenv("USER")
#air specific humidity profile [kg.kg-1]
Lv  <-  2.66e6

path.i  <- paste0('/work/',name.user,'/zun_mdtf/inputdata/model/GFDL.ESM4.piControl/daily/')

#######################################
#filei <-  paste0(path.i,'GFDL.ESM4.piControl.pr.daily.nc')
#data  <- nc_open(filei)
#lons.i <- ncvar_get(data, "lon")
#lats.i <- ncvar_get(data, "lat")
#pr.i  <- ncvar_get(data, 'pr')
#nc_close(data) 
#
#filei <-  paste0(path.i,'GFDL.ESM4.piControl.mrsos.daily.nc')
#data  <- nc_open(filei)
#sm.i  <- ncvar_get(data, 'mrsos')
#nc_close(data) 
#
#filei <-  paste0(path.i,'GFDL.ESM4.piControl.tasmax.daily.nc')
#data  <- nc_open(filei)
#tsmax.i  <- ncvar_get(data, 'tasmax')
#nc_close(data) 
#
#filei <-  paste0(path.i,'GFDL.ESM4.piControl.hfls.daily.nc')
#data  <- nc_open(filei)
#lhf.i  <- ncvar_get(data, 'hfls')
#nc_close(data) 
#
#filei <-  paste0(path.i,'GFDL.ESM4.piControl.rlds.daily.nc')
#data  <- nc_open(filei)
#rlds.i  <- ncvar_get(data, 'rlds')
#nc_close(data) 
#
#filei <-  paste0(path.i,'GFDL.ESM4.piControl.rlus.daily.nc')
#data  <- nc_open(filei)
#rlus.i  <- ncvar_get(data, 'rlus')
#nc_close(data) 
#
#filei <-  paste0(path.i,'GFDL.ESM4.piControl.rsds.daily.nc')
#data  <- nc_open(filei)
#rsds.i  <- ncvar_get(data, 'rsds')
#nc_close(data) 
#
#filei <-  paste0(path.i,'GFDL.ESM4.piControl.rsus.daily.nc')
#data  <- nc_open(filei)
#rsus.i  <- ncvar_get(data, 'rsus')
#nc_close(data) 
#
#
#if (max(lons.i) > 250)
#{
#    lons.i  <-  lons.i - 180
#    pr.i    <-  ShiftLon(pr.i)
#    sm.i    <-  ShiftLon(sm.i)
#    tsmax.i <-  ShiftLon(tsmax.i)
#    lhf.i   <-  ShiftLon(lhf.i)
#    rlds.i  <-  ShiftLon(rlds.i)
#    rlus.i  <-  ShiftLon(rlus.i)
#    rsds.i  <-  ShiftLon(rsds.i)
#    rsus.i  <-  ShiftLon(rsus.i)
#}
#
#
#rnet.i  <-  rsds.i+rlds.i-rsus.i-rlus.i
#rnet.t  <-  apply(rnet.i,c(1,2),mean,na.rm=T)
#pr.t    <-  apply(pr.i,c(1,2),mean,na.rm=T)
#
#lsm.i <-  apply(sm.i,c(1,2),sd)
#lsm.i[lsm.i > 0] <-  1
#lsm.i[lsm.i == 0] <-  0
#
#ai.i  <-  .8*rnet.t/(Lv*pr.t)
#
#dff.o <-  GetBinMat(sm.i,sm.i,ai.i,lsm.i)

png(paste0("DFF.png",sep=""), width=500, height=450)

par(mgp=c(3,.5,0),
    oma=c(.2,.2,.2,.2),mar=c(2,2,.5,.5),cex=1.5)

n.bin <-  100/dim(dff.o)[1]

image.plot(seq(n.bin,100,n.bin),seq(n.bin,100,n.bin),
           dff.o, col=col.dff(50),xlab='',ylab='',
           bigplot=c(.1,.85,.1,.95),
           smallplot=c(.87,.9,.1,.95))

mtext(1,text='soil moisture percentile',line=2.3,cex=1.5)
mtext(2,text='Climatological AI percentile',line=2.3,cex=1.5)

dev.off()

############################################################
print("Normal End of DFF.R")
############################################################
