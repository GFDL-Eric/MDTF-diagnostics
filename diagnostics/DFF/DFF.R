################################
# .libPaths("/home/water2/ab5/R/x86_64-redhat-linux-gnu-library/3.2")
POD_HOME  <- Sys.getenv("POD_HOME")
WK_DIR <- Sys.getenv("WK_DIR")
OBS_DATA <- Sys.getenv("OBS_DATA")
DATADIR <- Sys.getenv("DATADIR")
CASENAME <- Sys.getenv("CASENAME")
yr1 <- Sys.getenv("FIRSTYR")
yr2 <- Sys.getenv("LASTYR")
source(paste0(POD_HOME,"/lib_LoCo.R"))

Lv  <-  2.66e6

PR_FILE <- Sys.getenv("PR_FILE")
MRSOS_FILE <- Sys.getenv("MRSOS_FILE")
TASMAX_FILE <- Sys.getenv("TASMAX_FILE")
HFLS_FILE <- Sys.getenv("HFLS_FILE")
RLDS_FILE <- Sys.getenv("RLDS_FILE")
RLUS_FILE <- Sys.getenv("RLUS_FILE")
RSDS_FILE <- Sys.getenv("RSDS_FILE")
RSUS_FILE <- Sys.getenv("RSUS_FILE")

#precipitation [kg.m-2.s-1]
name.pr   <- Sys.getenv("pr_var")
#soil moisture content [kg.m-2]
name.mrsos  <- Sys.getenv("mrsos_var")
#maximum air temperature [K]
name.tasmax <- Sys.getenv("tasmax_var")
#latent heat flux [W.m-2]
name.hfls <- Sys.getenv("hfls_var")
#radiation terms [W.m-2]
name.rlds <- Sys.getenv("rlds_var")
name.rlus <- Sys.getenv("rlus_var")
name.rsds <- Sys.getenv("rsds_var")
name.rsus <- Sys.getenv("rsus_var")

##########################################
print("Read data")
data <- nc_open(PR_FILE)
lons.i <- ncvar_get(data, "lon")
lats.i <- ncvar_get(data, "lat")
pr.i  <-  ncvar_get(data,name.pr)
nc_close(data) 

data <- nc_open(MRSOS_FILE)
sm.i  <-  ncvar_get(data,name.mrsos)
nc_close(data) 

data <- nc_open(TASMAX_FILE)
tsmax.i  <-  ncvar_get(data,name.tasmax)
nc_close(data) 

data <- nc_open(HFLS_FILE)
lhf.i <-  ncvar_get(data,name.hfls)
nc_close(data) 

data <- nc_open(RLDS_FILE)
rlds.i <-  ncvar_get(data,name.rlds)
nc_close(data) 

data <- nc_open(RLUS_FILE)
rlus.i <-  ncvar_get(data,name.rlus)
nc_close(data) 

data <- nc_open(RSDS_FILE)
rsds.i <-  ncvar_get(data,name.rsds)
nc_close(data) 

data <- nc_open(RSUS_FILE)
rsus.i <-  ncvar_get(data,name.rsus)
nc_close(data) 

print("PP: checking the longitude system.")
#check if lons ranges from [0,360]
if (max(lons.i) > 250)
{
    lons.i  <-  lons.i - 180
    pr.i    <-  ShiftLon(pr.i)
    sm.i    <-  ShiftLon(sm.i)
    tsmax.i <-  ShiftLon(tsmax.i)
    lhf.i   <-  ShiftLon(lhf.i)
    rlds.i  <-  ShiftLon(rlds.i)
    rlus.i  <-  ShiftLon(rlus.i)
    rsds.i  <-  ShiftLon(rsds.i)
    rsus.i  <-  ShiftLon(rsus.i)
}

print("calculate the arid index.")
rnet.i  <-  rsds.i+rlds.i-rsus.i-rlus.i
rnet.t  <-  apply(rnet.i,c(1,2),mean,na.rm=T)
pr.t    <-  apply(pr.i,c(1,2),mean,na.rm=T)

lsm.i <-  apply(sm.i,c(1,2),sd)
lsm.i[lsm.i > 0] <-  1
lsm.i[lsm.i == 0] <-  0

ai.i  <-  .8*rnet.t/(Lv*pr.t)

############################################################
print("Calculate DFF")
############################################################
dff.o <-  GetBinMat(sm.i,sm.i,ai.i,lsm.i)

############################################################
print("Plot DFF")
############################################################
png(paste0(WK_DIR,"/model/DFF.png",sep=""), width=500, height=450)

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

#write the data out as a table
fileo <-  paste0(WK_DIR,'/model/DFF_',CASENAME,'.txt')
write.table(dff.o,file=fileo,row.names=F,col.names=F,
            sep='\t')

############################################################
print("Normal End of DFF.R")
############################################################
