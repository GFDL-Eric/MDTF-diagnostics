################################
# .libPaths("/home/water2/ab5/R/x86_64-redhat-linux-gnu-library/3.2")

POD_HOME  <- Sys.getenv("POD_HOME")
WK_DIR <- Sys.getenv("WK_DIR")
OBS_DATA <- Sys.getenv("OBS_DATA")
DATADIR  <- Sys.getenv("DATADIR")
CASENAME <- Sys.getenv("CASENAME")
yr1 <- Sys.getenv("FIRSTYR")
yr2 <- Sys.getenv("LASTYR")
source(paste0(POD_HOME,"/lib_LoCo.R"))

#air specific humidity profile [kg.kg-1]
name.hus  <- Sys.getenv("hus_var")
#air temperature profile [K]
name.ta   <- Sys.getenv("ta_var")
#surface pressure [Pa]
name.ps   <- Sys.getenv("ps_var")

print("##############################################")

HUS_FILE <- Sys.getenv("HUS_FILE")
TA_FILE  <- Sys.getenv("TA_FILE")
MSK_FILE <- Sys.getenv("SFTLF_FILE")

##########################################
print("Taking specific humidity profile from model")
data  <- nc_open(HUS_FILE)
times.utc <-  nc.get.time.series(data)
lons.i <- ncvar_get(data, "lon")
lats.i <- ncvar_get(data, "lat")
lev.i  <- ncvar_get(data, "lev")
ap.i  <- ncvar_get(data, "ap")
b.i  <- ncvar_get(data, "b")
hus.i <- ncvar_get(data, name.hus)
ps.i  <- ncvar_get(data, name.ps)
nc_close(data) 

print("Taking air temperature profile from model")
data <- nc_open(TA_FILE)
## Here we read all years
ta.i <- ncvar_get(data, name.ta)
nc_close(data) 

print("Taking land mask from model")
data  <- nc_open(MSK_FILE)
msk.i <- ncvar_get(data,"sftlf")
nc_close(data) 

#################
#data preprocessing
#################
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
ta.i    <-  ta.i[,lat.s.ind:lat.n.ind,,]
hus.i   <-  hus.i[,lat.s.ind:lat.n.ind,,]
ps.i    <-  ps.i[,lat.s.ind:lat.n.ind,]
msk.i   <-  msk.i[,lat.s.ind:lat.n.ind]

print("PP: checking the longitude system.")
#check if lons ranges from [0,360]
if (max(lons.i) > 250)
{
    lons.i  <-  lons.i - 180
    ta.i    <-  ShiftLon(ta.i)
    hus.i   <-  ShiftLon(hus.i)
    ps.i    <-  ShiftLon(ps.i)
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
for (nx in 1:n.lon)
{
    #Get local time and dates
    times.lst <-  UTC2LST(times.utc,lons.i[nx])
    dates.lst <-  unique(format(times.lst,'%Y-%m-%d'))
    for (ny in 1:n.lat)
    {
        #mask grid cells with land less than half
        if (msk.i[nx,ny] < 50)
            next

        #calculate the sunrise and sunset local time
        for (nt in 1:n.date)
        {
            srs.lst <-  GetSunRiseSet(as.Date(dates.utc[nt]),
                                      lons.i[nx],lats.i[ny]) 

            dt.t  <-  as.numeric(times.lst - srs.lst[1])

            #current date doesn't contain morning
            if (max(dt.t) <= 0)
                next

            dt.t[dt.t <= 0]  <-  NA
            idx.t <-  which.min(dt.t)

            if (times.lst[idx.t] > mean(srs.lst))
                next

            ctp.o[nx,ny,nt]  <-  CalCTP(ps.i[nx,ny,idx.t]*b.i+ap.i,
                                        ta.i[nx,ny,,idx.t],
                                        hus.i[nx,ny,,idx.t],
                                        ps.i[nx,ny,idx.t])
            hi.o[nx,ny,nt]   <-  CalHilow(ps.i[nx,ny,idx.t]*b.i+ap.i,
                                          ta.i[nx,ny,,idx.t],
                                          hus.i[nx,ny,,idx.t],
                                          ps.i[nx,ny,idx.t])
        }
    }
}

##########################################
print("Plotting daily mean of CTP and Hi_low")
##########################################
ctp.t <-  apply(ctp.o,c(1,2),mean,na.rm=T)
hi.t  <-  apply(hi.o,c(1,2),mean,na.rm=T)

ctp.br  <-  seq(-400,400,50)
ctp.tmp <-  ctp.t
ctp.tmp[!is.na(ctp.tmp) & (ctp.tmp < -400)] <-  -400
ctp.tmp[!is.na(ctp.tmp) & (ctp.tmp >  400)] <-  400
png(paste0(WK_DIR,"/model/day_mean_ctp.png"), width=700)
image.plot(lons.i,lats.i,ctp.tmp,zlim=range(ctp.br),breaks=ctp.br,
           col=col.df(length(ctp.br)-1),xlab='',ylab='')
plot(coastsCoarse,add=T)
dev.off()

##################
hi.br  <-  seq(0,50,5)
hi.tmp <-  hi.t
hi.tmp[!is.na(hi.tmp) & (hi.tmp < 0)] <-  0
hi.tmp[!is.na(hi.tmp) & (hi.tmp > 50)] <-  50
png(paste0(WK_DIR,"/model/day_mean_hi_low.png"), width=700)
image.plot(lons.i,lats.i,hi.tmp,zlim=range(hi.br),breaks=hi.br,
           col=col.val(length(hi.br)-1),xlab='',ylab='')
plot(coastsCoarse,add=T)
dev.off()


##########################################
print("Plotting the diagnosis based on the CTP-Hi_low")
##########################################
diag.t  <-  DiagCTPHilow(ctp.t,hi.t)

png(paste0(WK_DIR,"/model/diag_hi_low.png"), width=700)
image.plot(lons.i,lats.i,diag.t,
           col=col.ind(12)[c(10,1,12,4)],
           axis.args=list(at=seq(1,4),
                          labels=c('Atmosphere','Wet',
                                   'Transition','Dry'),
                          las=3),xlab='',ylab='')
plot(coastsCoarse,add=T)
dev.off()

##########################################
print("Write data out into NetCDF")
##########################################
fileo <-  paste0(WK_DIR,'/model/MDTF_CTP_Hilow_',CASENAME,'.nc')
x.dim <-  ncdim_def("lon","longitude",lons.i)
y.dim <-  ncdim_def("lat","latitude",lats.i)
t.dim <-  ncdim_def("time",paste0("days since ",
                    dates.utc[1]),seq(.5,n.date-.5,1),unlim=T)

ctp.v <-  ncvar_def("ctp","J.kg-1",list(x.dim,y.dim,t.dim),
                    NA,"Convective Triggering Potential")
hi.v  <-  ncvar_def("hi_low","C",list(x.dim,y.dim,t.dim),
                    NA,"Low level humidity Index")

nco <-  nc_create(fileo,list(ctp.v,hi.v))

ncvar_put(nco,ctp.v,ctp.o)
ncvar_put(nco,hi.v,hi.o)

nc_close(nco)

############################################################
print("Normal End of Convect Triggering Potential.R")
############################################################
