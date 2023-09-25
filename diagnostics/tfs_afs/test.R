################################
source("lib_LoCo.R")

name.user  <- Sys.getenv("USER")
#mm
pre.thr <-  1

#sensible heat flux [W.m-2]
name.shf  <- 'hfss'
#latent heat flux [W.m-2]
name.lhf  <- 'hfls'
#precipitation [kg.m-2.s-1]
name.pre  <- 'pr'

print("##############################################")

path.i  <- paste0('/work/',name.user,'/zun_mdtf/inputdata/model/GFDL.CM4.c96L32.amip/6hr/')

##############################################
print("Taking sensible heat flux from model")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.hfss.6hr.nc')
data <- nc_open(filei)
times.utc <-  nc.get.time.series(data)
cal.i <-  ncatt_get(data,"time","calendar")$value
lons.i <- ncvar_get(data, "lon")
lats.i <- ncvar_get(data, "lat")
## Here we read all years
shf.i <- ncvar_get(data, name.shf)
nc_close(data) 

print("Taking latent heat flux from model")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.hfls.6hr.nc')
data <- nc_open(filei)
## Here we read all years
lhf.i <- ncvar_get(data, name.lhf)
nc_close(data) 

print("Taking precipitation from model")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.pr.6hr.nc')
data <- nc_open(filei)
## Here we read all years
pre.i <- ncvar_get(data, name.pre)
nc_close(data) 

print("Taking land mask from model")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.sftlf.6hr.nc')
data  <- nc_open(filei)
msk.i <- ncvar_get(data,"sftlf")
nc_close(data) 

###################
#data preprocessing
###################
#only consider regions between [-lat.lim,lat.lim]

lat.lim <-  60
print("PP: checking study regions")
lats.t  <-  lats.i
lats.t[lats.t < -lat.lim]    <-  1e3
#index of the south point
lat.s.ind <-  which.min(abs(lats.t + lat.lim))
lats.t  <-  lats.i
lats.t[lats.t > lat.lim]    <-  -1e3
#index of the north point
lat.n.ind <-  which.min(abs(lats.t - lat.lim))

if (lat.n.ind <= lat.s.ind)
        stop('Error: study area does not cover [-60S--60N]')
#corp regions for [-60,60]
lats.i  <-  lats.i[lat.s.ind:lat.n.ind]
shf.i    <-  shf.i[,lat.s.ind:lat.n.ind,]
lhf.i    <-  lhf.i[,lat.s.ind:lat.n.ind,]
pre.i   <-  pre.i[,lat.s.ind:lat.n.ind,]
msk.i   <-  msk.i[,lat.s.ind:lat.n.ind]
ef.i  <-  lhf.i/(lhf.i+shf.i)

print("PP: checking the longitude system.")
#check if lons ranges from [0,360]
if (max(lons.i) > 250)
{
    lons.i  <-  lons.i - 180
    pre.i   <-  ShiftLon(pre.i)
    ef.i    <-  ShiftLon(ef.i)
    shf.i   <-  ShiftLon(shf.i)
    lhf.i   <-  ShiftLon(lhf.i)
    msk.i    <-  ShiftLon(msk.i)
}


print("PP: collecting the time information.")
dates.utc <-  unique(format(times.utc,'%Y-%m-%d'))
n.date  <-  length(dates.utc)
n.lon   <-  length(lons.i)
n.lat   <-  length(lats.i)


##########################################
print("Start calculate CTP--Hi_low")
##########################################
ctp.o <-  array(NA,dim=c(n.lon,n.lat,n.date))
hi.o  <-  array(NA,dim=c(n.lon,n.lat,n.date))
ef.o  <-  array(NA,dim=c(n.lon,n.lat,n.date))
pre.am.o <-  array(NA,dim=c(n.lon,n.lat,n.date))
pre.pm.o <-  array(NA,dim=c(n.lon,n.lat,n.date))
for (nx in 1:n.lon)
{
    #Get local time and dates
    times.lst <-  UTC2LST(times.utc,lons.i[nx])
    timed.lst <-  format(times.lst,'%Y-%m-%d')
    dates.lst <-  unique(format(times.lst,'%Y-%m-%d'))
    print(paste("nx=",nx))

    for (nt in 1:n.date)
    {
        times.today <-  times.lst[timed.lst ==
                                  dates.utc[nt]]
        ind.today   <-  seq(1,length(times.lst))[timed.lst == dates.utc[nt]]

        srs.a <-  as.PCICt(paste(dates.utc[nt],"06:00:00"),
                           cal=cal.i)
        srs.m <-  as.PCICt(paste(dates.utc[nt],"12:00:00"),
                           cal=cal.i)
        srs.p <-  as.PCICt(paste(dates.utc[nt],"18:00:00"),
                           cal=cal.i)

        #current date doesn't contain data in required period
        #(morning and afternoon)
        hf.hr <-  3600
        ind.am <-  (times.today >= srs.a - hf.hr) &
                   (times.today < srs.m - hf.hr)
        ind.md <-  (times.today >= srs.m - hf.hr) &
                   (times.today < srs.m + hf.hr)
        ind.pm <-  (times.today >= srs.m + hf.hr) &
                   (times.today <= srs.p + hf.hr)

        if (!(any(ind.am) &
          any(ind.pm)))
        next

        if (any(ind.md))
        {
            idx.a <-  min(ind.today) +
                      seq(1,length(times.today))[ind.am] - 1
            idx.a.min <-  min(idx.a)
            idx.m <-  min(ind.today) +
                      seq(1,length(times.today))[ind.md] - 1
            idx.p <-  min(ind.today) +
                      seq(1,length(times.today))[ind.pm] - 1

            for (ny in 1:n.lat)
            {
                if (msk.i[nx,ny] < 50)
                    next
                ef.o[nx,ny,nt]  <-  ef.i[nx,ny,idx.m]
                pre.am.o[nx,ny,nt] <-  mean(pre.i[nx,ny,c(idx.a,idx.m)])
                pre.pm.o[nx,ny,nt] <-  mean(pre.i[nx,ny,c(idx.m,idx.p)])
            }
        } else
        {
            idx.a <-  min(ind.today) +
                      seq(1,length(times.today))[ind.am] - 1
            idx.a.min  <-  min(idx.a)
            idx.a.max  <-  max(idx.a)
            idx.p <-  min(ind.today) +
                      seq(1,length(times.today))[ind.pm] - 1
            for (ny in 1:n.lat)
            {
                if (msk.i[nx,ny] < 50)
                    next
                ef.o[nx,ny,nt]  <-  ef.i[nx,ny,idx.a.max]
                pre.am.o[nx,ny,nt] <-  mean(pre.i[nx,ny,idx.a])
                pre.pm.o[nx,ny,nt] <-  mean(pre.i[nx,ny,idx.p])
            }
        }
    }
}

#convert the unit from [kg.m-2.s-1] to [mm]
pre.am.o <-  pre.am.o*3600*6
pre.pm.o <-  pre.pm.o*3600*6

save.image(file='tt.RData')


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
