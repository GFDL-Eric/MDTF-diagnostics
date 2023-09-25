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

#sensible heat flux [W.m-2]
name.shf  <- Sys.getenv("hfss_var")
#latent heat flux [W.m-2]
name.lhf  <- Sys.getenv("hfls_var")
#precipitation [kg.m-2.s-1]
name.pre  <- Sys.getenv("pr_var")

print("##############################################")

SHF_FILE <- Sys.getenv("HFSS_FILE")
LHF_FILE <- Sys.getenv("HFLS_FILE")
PRE_FILE <- Sys.getenv("PR_FILE")
MSK_FILE <- Sys.getenv("SFTLF_FILE")

##########################################
print("Taking land mask from model")
data  <- nc_open(MSK_FILE)
msk.i <- ncvar_get(data,"sftlf")
nc_close(data) 

print("Taking sensible heat flux from model")
data <- nc_open(SHF_FILE)
times.utc <-  nc.get.time.series(data)
cal.i  <-  ncatt_get(data,"time","calendar")$value
lons.i <- ncvar_get(data, "lon")
lats.i <- ncvar_get(data, "lat")
shf.i <- ncvar_get(data, name.shf)
nc_close(data) 

print("Taking latent heat flux from model")
data <- nc_open(LHF_FILE)
lhf.i <- ncvar_get(data, name.lhf)
nc_close(data) 

print("Taking precipitation from model")
data <- nc_open(PRE_FILE)
pre.i <- ncvar_get(data, name.pre)
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
print("Start calculate TFS--AFS")
##########################################
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

#convert the unit from [kg.m-2.s-1] to [mm per 6 hour]
pre.am.o <-  pre.am.o*3600*6
pre.pm.o <-  pre.pm.o*3600*6

tfs.o <-  array(NA,dim=c(n.lon,n.lat))
afs.o <-  array(NA,dim=c(n.lon,n.lat))
for (nx in 1:n.lon)
{
    print(paste('nx =',nx))
for (ny in 1:n.lat)
{
    if (msk.i[nx,ny] < 50)
        next

    v.o  <-  CalTFSAFS(ef.o[nx,ny,],pre.am.o[nx,ny,],
                       pre.pm.o[nx,ny,])
    tfs.o[nx,ny]  <-  v.o$tfs
    afs.o[nx,ny]  <-  v.o$afs
}
}

##########################################
print("Plotting TFS and AFS")
##########################################
thr.t <-  .99
t.thr <-  quantile(abs(tfs.o),probs=thr.t,na.rm=T)
a.thr <-  quantile(abs(afs.o),probs=thr.t,na.rm=T)

png(paste0(WK_DIR,"/model/tfs.png"), width=700,height=450)
par(oma=rep(.5,4),mar=c(2,2,0,0),cex=1.2)
plt.t <-  tfs.o
plt.t[plt.t > t.thr]  <-  t.thr
plt.t[plt.t < -t.thr]  <-  -t.thr
image.plot(lons.i,lats.i,plt.t,zlim=c(-t.thr,t.thr),
           col=col.df(50),xlab='',ylab='',
           bigplot=c(.05,.87,.1,.9),
           smallplot=c(.9,.93,.1,.85))
plot(coastsCoarse,add=T)
dev.off()

png(paste0(WK_DIR,"/model/afs.png"), width=700,height=450)
par(oma=rep(.5,4),mar=c(2,2,0,0),cex=1.2)
plt.t <-  afs.o
plt.t[plt.t > t.thr]  <-  t.thr
plt.t[plt.t < -t.thr]  <-  -t.thr
image.plot(lons.i,lats.i,plt.t,zlim=c(-t.thr,t.thr),
           col=col.df(50),xlab='',ylab='',
           bigplot=c(.05,.87,.1,.9),
           smallplot=c(.9,.93,.1,.85))
plot(coastsCoarse,add=T)
dev.off()

##########################################
print("Write data out into NetCDF")
##########################################
fileo <-  paste0(WK_DIR,'/model/TFS_AFS_',CASENAME,'.nc')
x.dim <-  ncdim_def("lon","longitude",lons.i)
y.dim <-  ncdim_def("lat","latitude",lats.i)

tfs.v <-  ncvar_def("tfs","1",list(x.dim,y.dim),
                    NA,"Triggering feedback strength")
afs.v <-  ncvar_def("afs","1",list(x.dim,y.dim),
                    NA,"Amplification feedback strength")

nco <-  nc_create(fileo,list(tfs.v,afs.v))

ncvar_put(nco,tfs.v,tfs.o)
ncvar_put(nco,afs.v,afs.o)

nc_close(nco)

#############################################################
print("Normal End of Convect Triggering Potential.R")
#############################################################
