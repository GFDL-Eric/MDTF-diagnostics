################################
source("lib_LoCo.R")

name.user  <- Sys.getenv("USER")
#air specific humidity profile [kg.kg-1]
name.q2m <- 'huss'
#air temperature profile [K]
name.t2m <- 'tas'
#surface pressure [Pa]
name.ps   <- 'ps'
name.shf  <- 'hfss'
name.sm   <- 'mrsos'

print("##############################################")

path.i  <- paste0('/work/',name.user,'/zun_mdtf/inputdata/model/GFDL.CM4.c96L32.amip/mon/')

#######################################
print("Taking 2-m specific humidity")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.huss.mon.nc')
data  <- nc_open(filei)
times.utc <-  nc.get.time.series(data)
lons.i <- ncvar_get(data, "lon")
lats.i <- ncvar_get(data, "lat")
q2m.i <- ncvar_get(data, name.q2m)
nc_close(data) 

print("Taking air temperature")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.tas.mon.nc')
data  <- nc_open(filei)
t2m.i <- ncvar_get(data, name.t2m)
nc_close(data) 

print("Taking soil moisture")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.mrsos.mon.nc')
data  <- nc_open(filei)
sm.i  <- ncvar_get(data, name.sm)
nc_close(data) 

print("Taking sensible heat flux")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.hfss.mon.nc')
data  <- nc_open(filei)
shf.i <- ncvar_get(data, name.shf)
nc_close(data) 

print("Taking surface pressure")
filei <-  paste0(path.i,'GFDL.CM4.c96L32.amip.ps.mon.nc')
data  <- nc_open(filei)
ps.i  <- ncvar_get(data, name.ps)
nc_close(data) 

n.lon <-  length(lons.i)
n.lat <-  length(lats.i)


#calculate PLCL
d2m.i <-  DewT(q2m.i,ps.i)
plcl.i  <-  ps.i - ps.i*((t2m.i-d2m.i)/223.15+1)^(-3.5)

msk.i <-  apply(sm.i,c(1,2),sd)
msk.i[msk.i > 0]  <-  1

sm.i[msk.i < .5]  <-  NA
plcl.i[msk.i < .5]  <-  NA
shf.i[msk.i < .5]  <-  NA


###################
#data preprocessing
###################

print("PP: checking the longitude system.")
#check if lons ranges from [0,360]
if (max(lons.i) > 250)
{
    lons.i  <-  lons.i - 180
    sm.i    <-  ShiftLon(sm.i)
    shf.i   <-  ShiftLon(shf.i)
    plcl.i  <-  ShiftLon(plcl.i)
    msk.i  <-  ShiftLon(msk.i)
}


print("PP: collecting the time information.")
dates.utc <-  unique(format(times.utc,'%Y-%m-%d'))
n.date  <-  length(dates.utc)

times.mon <-  as.numeric(format(times.utc,'%m'))

djf.ind <-  seq(1,length(times.mon))[times.mon < 3 | times.mon == 12]
mam.ind <-  seq(1,length(times.mon))[times.mon > 2 & times.mon < 6]
jja.ind <-  seq(1,length(times.mon))[times.mon > 5 & times.mon < 9]
son.ind <-  seq(1,length(times.mon))[times.mon > 8 & times.mon < 12]

#[lon,lat,4season]
land.leg  <-  array(NA,dim=c(n.lon,n.lat,4))
land.p    <-  array(NA,dim=c(n.lon,n.lat,4))
atmos.leg <-  array(NA,dim=c(n.lon,n.lat,4))
atmos.p   <-  array(NA,dim=c(n.lon,n.lat,4))
tot.leg   <-  array(NA,dim=c(n.lon,n.lat,4))
tot.p     <-  array(NA,dim=c(n.lon,n.lat,4))

sea.ind <-  c('djf','mam','jja','son')
for (nx in 1:n.lon)
{
    print(paste0('nx = ',nx))
for (ny in 1:n.lat)
{
        #mask grid cells with land less than half
        if (msk.i[nx,ny] < .5 | is.na(mean(sm.i[nx,ny,])))
            next

        for (ns in 1:4)
        {
            eval(parse(text=paste0('t.ind <- ',sea.ind[ns],'.ind')))
            sm.ts <-  CircRemove(sm.i[nx,ny,t.ind],times.utc[t.ind])
            shf.ts <-  CircRemove(shf.i[nx,ny,t.ind],times.utc[t.ind])
            plcl.ts <-  CircRemove(plcl.i[nx,ny,t.ind],
                                   times.utc[t.ind])
            cor.l <-  cor.test(sm.ts,shf.ts)
            cor.a <-  cor.test(plcl.ts,shf.ts)
            land.leg[nx,ny,ns]  <-  cor.l$estimate*sd(shf.ts)
            land.p[nx,ny,ns]    <-  cor.l$p.value
            atmos.leg[nx,ny,ns]  <-  cor.a$estimate*sd(plcl.ts)
            atmos.p[nx,ny,ns]   <-  cor.a$p.value
            tot.leg[nx,ny,ns]  <-  cor.l$estimate*
                                   atmos.leg[nx,ny,ns]
            tot.p[nx,ny,ns] <-  max(cor.l$p.value,cor.a$p.value)
        }
}
}

##########################################
print("Plotting the two-legged metrix")
##########################################

if (n.lat%%2 == 0)
{
    north.lat  <-  n.lat/2+1
    south.lat  <-  n.lat/2
} else
{
    north.lat  <-  (n.lat+1)/2
    south.lat  <-  (n.lat-1)/2
}

thr.t <-  .99
l.thr <-  quantile(abs(land.leg),probs=thr.t,na.rm=T)
a.thr <-  quantile(abs(atmos.leg),probs=thr.t,na.rm=T)
t.thr <-  quantile(abs(tot.leg),probs=thr.t,na.rm=T)

png(paste0("land_leg.png"), width=700,height=450)
par(oma=rep(.5,4),mar=c(2,2,0,0),cex=1.2)
leg.t <-  land.leg[,,3]
leg.t[,1:south.lat] <-  land.leg[,1:south.lat,1]
leg.t[leg.t > l.thr] <-  l.thr
leg.t[leg.t < -l.thr] <-  -l.thr
image.plot(lons.i,lats.i,leg.t,zlim=c(-l.thr,l.thr),
           col=col.df(50),xlab='',ylab='',
           bigplot=c(.05,.87,.1,.9),
           smallplot=c(.9,.93,.1,.85),
           legend.args=list(text='[W.m-2]',cex=1.3,
                            line=.5))
plot(coastsCoarse,add=T)
dev.off()

png(paste0("atmos_leg.png"), width=700,height=450)
par(oma=rep(.5,4),mar=c(2,2,0,0),cex=1.2)
leg.t <-  atmos.leg[,,3]
leg.t[,1:south.lat] <-  atmos.leg[,1:south.lat,1]
leg.t[leg.t > a.thr] <-  a.thr
leg.t[leg.t < -a.thr] <-  -a.thr
image.plot(lons.i,lats.i,leg.t,zlim=c(-a.thr,a.thr),
           col=col.df(50),xlab='',ylab='',
           bigplot=c(.05,.87,.1,.9),
           smallplot=c(.9,.93,.1,.85),
           legend.args=list(text='[Pa]',cex=1.3,
                            line=.5))
plot(coastsCoarse,add=T)
dev.off()

png(paste0("total_leg.png"), width=700,height=450)
par(oma=rep(.5,4),mar=c(2,2,0,0),cex=1.2)
leg.t <-  tot.leg[,,3]
leg.t[,1:south.lat] <-  tot.leg[,1:south.lat,1]
leg.t[leg.t > t.thr] <-  t.thr
leg.t[leg.t < -t.thr] <-  -t.thr
image.plot(lons.i,lats.i,leg.t,zlim=c(-t.thr,t.thr),
           col=col.df(50),xlab='',ylab='',
           bigplot=c(.05,.87,.1,.9),
           smallplot=c(.9,.93,.1,.85),
           legend.args=list(text='[Pa]',cex=1.3,
                            line=.5))
plot(coastsCoarse,add=T)
dev.off()


#############################################################
print("Write data out into NetCDF")
#############################################################
fileo  <-  "tt.nc"
x.dim <-  ncdim_def("lon","degrees_east",lons.i)
y.dim <-  ncdim_def("lat","degrees_north",lats.i)
z.dim <-  ncdim_def("season","1",seq(1,4),
                    longname="DJF_MAM_JJA_SON")

ll.v  <-  ncvar_def("land_leg","W.m-2",list(x.dim,y.dim,z.dim),
                    NA,"coupling strength between the land and the flux")
ll.p  <-  ncvar_def("land_p","-",list(x.dim,y.dim,z.dim),
                    NA,"coupling significance between the land and the flux")
aa.v  <-  ncvar_def("atmos_leg","Pa",list(x.dim,y.dim,z.dim),
                    NA,"coupling strength between the flux and the atmosphere")
aa.p  <-  ncvar_def("atmos_p","-",list(x.dim,y.dim,z.dim),
                    NA,"coupling significance between the flux and the atmosphere")
tt.v  <-  ncvar_def("tot_leg","Pa",list(x.dim,y.dim,z.dim),
                    NA,"coupling strength between the land and the atmosphere")
tt.p  <-  ncvar_def("tot_p","-",list(x.dim,y.dim,z.dim),
                    NA,"coupling significance between the land and the atmosphere")

nco <-  nc_create(fileo,list(ll.v,ll.p,aa.v,aa.p,tt.v,tt.p))
ncvar_put(nco,ll.v,land.leg)
ncvar_put(nco,ll.p,land.p)
ncvar_put(nco,aa.v,atmos.leg)
ncvar_put(nco,aa.p,atmos.p)
ncvar_put(nco,tt.v,tot.leg)
ncvar_put(nco,tt.p,tot.p)

nc_close(nco)
#############################################################
#print("Normal End of Land Atmosphere Coupling.R")
#############################################################
