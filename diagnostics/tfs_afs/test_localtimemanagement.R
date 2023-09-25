################################
source("lib_LoCo.R")

name.user  <- Sys.getenv("USER")
#air specific humidity profile [kg.kg-1]
name.hus  <- 'hus'
#air temperature profile [K]
name.ta   <- 'ta'
#surface pressure [Pa]
name.ps   <- 'ps'
#sensible heat flux [W.m-2]
name.shf  <- 'hfss'
#latent heat flux [W.m-2]
name.lhf  <- 'hfls'
#precipitation [kg.m-2.s-1]
name.pre  <- 'pr'

print("##############################################")

path.i  <- paste0('/work/',name.user,'/zun_mdtf/inputdata/model/GFDL.CM4.c96L32.amip/6hr/')

#########################################
#print("Taking specific humidity profile from model")
#filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.hus.6hr.nc')
#data  <- nc_open(filei)
#times.utc <-  nc.get.time.series(data)
#lons.i <- ncvar_get(data, "lon")
#lats.i <- ncvar_get(data, "lat")
#lev.i  <- ncvar_get(data, "lev")
#ap.i  <- ncvar_get(data, "ap")
#b.i  <- ncvar_get(data, "b")
#hus.i <- ncvar_get(data, name.hus)
#ps.i  <- ncvar_get(data, name.ps)
#nc_close(data) 
#
#print("Taking air temperature profile from model")
#filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.ta.6hr.nc')
#data <- nc_open(filei)
### Here we read all years
#ta.i <- ncvar_get(data, name.ta)
#nc_close(data) 
#
#print("Taking sensible heat flux from model")
#filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.hfss.6hr.nc')
#data <- nc_open(filei)
### Here we read all years
#shf.i <- ncvar_get(data, name.shf)
#nc_close(data) 
#
#print("Taking latent heat flux from model")
#filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.hfls.6hr.nc')
#data <- nc_open(filei)
### Here we read all years
#lhf.i <- ncvar_get(data, name.lhf)
#nc_close(data) 
#
#print("Taking precipitation from model")
#filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.pr.6hr.nc')
#data <- nc_open(filei)
### Here we read all years
#pre.i <- ncvar_get(data, name.pre)
#nc_close(data) 
#
#print("Taking land mask from model")
#filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.sftlf.6hr.nc')
#data  <- nc_open(filei)
#msk.i <- ncvar_get(data,"sftlf")
#nc_close(data) 
#
####################
##data preprocessing
####################
##only consider regions between [-lat.lim,lat.lim]
#ef.i  <-  lhf.i/(lhf.i+shf.i)
#
#lat.lim <-  60
#print("PP: checking study regions")
#lats.t  <-  lats.i
#lats.t[lats.t < -lat.lim]    <-  1e3
##index of the south point
#lat.s.ind <-  which.min(abs(lats.t + lat.lim))
#lats.t  <-  lats.i
#lats.t[lats.t > lat.lim]    <-  -1e3
##index of the north point
#lat.n.ind <-  which.min(abs(lats.t - lat.lim))
#
#if (lat.n.ind <= lat.s.ind)
#        stop('Error: study area does not cover [-60S--60N]')
##corp regions for [-60,60]
#lats.i  <-  lats.i[lat.s.ind:lat.n.ind]
#ta.i    <-  ta.i[,lat.s.ind:lat.n.ind,,]
#hus.i   <-  hus.i[,lat.s.ind:lat.n.ind,,]
#ps.i    <-  ps.i[,lat.s.ind:lat.n.ind,]
#ef.i    <-  ef.i[,lat.s.ind:lat.n.ind,]
#pre.i   <-  pre.i[,lat.s.ind:lat.n.ind,]
#msk.i   <-  msk.i[,lat.s.ind:lat.n.ind]
#
#print("PP: checking the longitude system.")
##check if lons ranges from [0,360]
#if (max(lons.i) > 250)
#{
#    lons.i  <-  lons.i - 180
#    ta.i    <-  ShiftLon(ta.i)
#    hus.i   <-  ShiftLon(hus.i)
#    ps.i    <-  ShiftLon(ps.i)
#    pre.i   <-  ShiftLon(pre.i)
#    ef.i    <-  ShiftLon(ef.i)
#    msk.i    <-  ShiftLon(msk.i)
#}
#
#
#print("PP: collecting the time information.")
#dates.utc <-  unique(format(times.utc,'%Y-%m-%d'))
#n.date  <-  length(dates.utc)
#n.lon   <-  length(lons.i)
#n.lat   <-  length(lats.i)


###########################################
#print("Start calculate CTP--Hi_low")
###########################################
#ctp.o <-  array(NA,dim=c(n.lon,n.lat,n.date))
#hi.o  <-  array(NA,dim=c(n.lon,n.lat,n.date))
#ef.o  <-  array(NA,dim=c(n.lon,n.lat,n.date))
#pre.o <-  array(NA,dim=c(n.lon,n.lat,n.date))
for (nx in 1:n.lon)
{
    #Get local time and dates
    times.lst <-  UTC2LST(times.utc,lons.i[nx])
    timed.lst <-  format(times.lst,'%Y-%m-%d')
    dates.lst <-  unique(format(times.lst,'%Y-%m-%d'))
    print(paste("nx=",nx))
    for (ny in 1:n.lat)
    {
        #mask grid cells with land less than half
        if (msk.i[nx,ny] < 50)
            next

        #calculate the sunrise and sunset local time
        for (nt in 1:n.date)
        {
            srs.lst <-  GetSunRiseSet(as.Date(dates.lst[nt]),
                                      lons.i[nx],lats.i[ny]) 

            times.today <-  times.lst[timed.lst ==
                            format(srs.lst[1],'%Y-%m-%d')]
            ind.today   <-  seq(1,length(times.lst))[timed.lst == dates.lst[nt]]

            #current date doesn't contain data in required period
            #(morning and afternoon)
            hf.hr <-  1800
            ind.morning <-  (times.today > srs.lst[1] - hf.hr) &
                            (times.today < mean(srs.lst)-hf.hr)
            ind.midday  <-  (times.today > mean(srs.lst)-hf.hr) &
                            (times.today < mean(srs.lst)+hf.hr)
            ind.afternoon <-  (times.today > mean(srs.lst)+hf.hr) &
                              (times.today < srs.lst[2]+hf.hr)
            if (!(any(ind.morning) &
                  any(ind.afternoon)))
                next

            if (any(ind.midday))
            {
                idx.m <-  min(ind.today) +
                          which.max(ind.morning) - 1
                idx.a <-  min(ind.today) +
                          which.max(ind.midday) - 1
                idx.e <-  min(ind.today) +
                          which.max(ind.afternoon) - 1
                ctp.o[nx,ny,nt]  <-  CalCTP(ps.i[nx,ny,idx.m]*b.i+ap.i,
                                            ta.i[nx,ny,,idx.m],
                                            hus.i[nx,ny,,idx.m],
                                            ps.i[nx,ny,idx.m])
                hi.o[nx,ny,nt]   <-  CalHilow(ps.i[nx,ny,idx.m]*b.i+ap.i,
                                              ta.i[nx,ny,,idx.m],
                                              hus.i[nx,ny,,idx.m],
                                              ps.i[nx,ny,idx.m])
                ef.o[nx,ny,nt]  <-  mean(ef.i[nx,ny,idx.m:idx.a])
                pre.o[nx,ny,nt] <-  mean(pre.i[nx,ny,idx.a:idx.e])
            } else
            {
                idx.m <-  min(ind.today) +
                          seq(1,length(times.today))[ind.morning] - 1
                idx.e <-  min(ind.today) +
                          seq(1,length(times.today))[ind.afternoon] - 1
                idx.mm  <-  min(idx.m)
                ctp.o[nx,ny,nt]  <-  CalCTP(ps.i[nx,ny,idx.mm]*b.i+ap.i,
                                            ta.i[nx,ny,,idx.mm],
                                            hus.i[nx,ny,,idx.mm],
                                            ps.i[nx,ny,idx.mm])
                hi.o[nx,ny,nt]   <-  CalHilow(ps.i[nx,ny,idx.mm]*b.i+ap.i,
                                              ta.i[nx,ny,,idx.mm],
                                              hus.i[nx,ny,,idx.mm],
                                              ps.i[nx,ny,idx.mm])
                ef.o[nx,ny,nt]  <-  mean(ef.i[nx,ny,idx.m])
                pre.o[nx,ny,nt] <-  mean(pre.i[nx,ny,idx.e])
            }
        }
    }
}

###########################################
#print("Plotting daily mean of CTP and Hi_low")
###########################################
#ctp.t <-  apply(ctp.o,c(1,2),mean,na.rm=T)
#hi.t  <-  apply(hi.o,c(1,2),mean,na.rm=T)
#
#ctp.br  <-  seq(-400,400,50)
#ctp.tmp <-  ctp.t
#ctp.tmp[!is.na(ctp.tmp) & (ctp.tmp < -400)] <-  -400
#ctp.tmp[!is.na(ctp.tmp) & (ctp.tmp >  400)] <-  400
#png(paste0("day_mean_ctp.png"), width=700)
#image.plot(lons.i,lats.i,ctp.tmp,zlim=range(ctp.br),breaks=ctp.br,
#           col=col.df(length(ctp.br)-1),xlab='',ylab='')
#plot(coastsCoarse,add=T)
#dev.off()
#
###################
#hi.br  <-  seq(0,50,5)
#hi.tmp <-  hi.t
#hi.tmp[!is.na(hi.tmp) & (hi.tmp < 0)] <-  0
#hi.tmp[!is.na(hi.tmp) & (hi.tmp > 50)] <-  50
#png(paste0("day_mean_hi_low.png"), width=700)
#image.plot(lons.i,lats.i,hi.tmp,zlim=range(hi.br),breaks=hi.br,
#           col=col.val(length(hi.br)-1),xlab='',ylab='')
#plot(coastsCoarse,add=T)
#dev.off()
#
#
###########################################
#print("Plotting the diagnosis based on the CTP-Hi_low")
###########################################
#diag.t  <-  DiagCTPHilow(ctp.t,hi.t)
#png(paste0("diag_hi_low.png"), width=700)
#image.plot(lons.i,lats.i,diag.t,
#           col=col.ind(12)[c(10,1,12,4)],
#           axis.args=list(at=seq(1,4),
#                          labels=c('Atmosphere','Wet',
#                                   'Transition','Dry'),
#                          las=3),xlab='',ylab='')
#plot(coastsCoarse,add=T)
#dev.off()
#
###########################################
#print("Write data out into NetCDF")
###########################################
#fileo <-  paste0('MDTF_CTP_Hilow.nc')
#x.dim <-  ncdim_def("lon","longitude",lons.i)
#y.dim <-  ncdim_def("lat","latitude",lats.i)
#t.dim <-  ncdim_def("time",paste0("days since ",
#                    dates.utc[1]),seq(.5,n.date-.5,1),unlim=T)
#
#ctp.v <-  ncvar_def("ctp","J.kg-1",list(x.dim,y.dim,t.dim),
#                    NA,"Convective Triggering Potential")
#hi.v  <-  ncvar_def("hi_low","C",list(x.dim,y.dim,t.dim),
#                    NA,"Low level humidity Index")
#
#nco <-  nc_create(fileo,list(ctp.v,hi.v))
#
#ncvar_put(nco,ctp.v,ctp.o)
#ncvar_put(nco,hi.v,hi.o)
#
#nc_close(nco)
#
#############################################################
#print("Normal End of Convect Triggering Potential.R")
#############################################################
