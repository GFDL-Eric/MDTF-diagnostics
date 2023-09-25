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

if (!exists('dates.utc'))
    load('tt.RData')

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



##############################################
print("Plotting daily mean of CTP and Hi_low")
##############################################
thr.t <-  .99
t.thr <-  quantile(abs(tfs.o),probs=thr.t,na.rm=T)
a.thr <-  quantile(abs(afs.o),probs=thr.t,na.rm=T)

png(paste0("tfs.png"), width=700,height=450)
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

png(paste0("afs.png"), width=700,height=450)
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
fileo <-  paste0('TFS_AFS.nc')
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
#print("Normal End of Convect Triggering Potential.R")
#############################################################
