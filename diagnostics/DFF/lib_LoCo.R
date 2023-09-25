library(RColorBrewer)
library(maps)
library(fields)
library(PCICt)
library(ncdf4.helpers)
library(suncalc)
library(ncdf4)
library(humidity)
library(foreach)
library(rworldmap)
library(doParallel)
library(abind)
library(plyr)

col.df  <-  colorRampPalette(brewer.pal(11,'RdBu'))
col.val <-  colorRampPalette(brewer.pal(9,'YlGnBu'))
col.ind <-  colorRampPalette(brewer.pal(12,'Set3'))
col.dff <-  colorRampPalette(brewer.pal(11,'RdYlGn'))

#####################
#functions
#####################
UTC2LST <-  function(ts.i,lon.i)
{
    #ts.i is the input time based on UTC
    #lon.i is the input longitude
    ts.o    <-  ts.i + lon.i*3600/15
    return(ts.o)
}

LST2UTC <-  function(ts.i,lon.i)
{
    #ts.i is the input time based on UTC
    #lon.i is the input longitude
    ts.o    <-  ts.i - lon.i*3600/15
    return(ts.o)
}

GetSunRiseSet  <-  function(date.i,lon.i,lat.i)
{
    v.tmp   <-  getSunlightTimes(date.i,lat=lat.i,lon=lon.i,
                                 keep=c("sunrise","sunset"))
    #sunrise
    srs  <-  UTC2LST(as.PCICt(v.tmp$sunrise,cal="365_day"),lon.i)
    srs  <-  c(srs,UTC2LST(as.PCICt(v.tmp$sunset,cal="365_day"),lon.i))
    return(srs)
}

CalCTP  <-  function(pa.i,ta.i,qa.i,ps.i,nseg=20){
  #calculate the CTP at the given time step
  #pa.i, pressure profile [Pa]
  #ta.i, temperature profile [K]
  #qa.i, specific humidity profile [kg/kg]
  #nseg, number of segments for integration

  #parameters
  grav  <-  9.81
  Rd  <-  287.04
  cp  <-  1005.7
  Rcp <-  Rd/cp
  Lv  <-  2.5e6
  Lv_cp <-  Lv/cp
  Rv  <-  461.5
  ep  <-  0.622

  #check if profile is upside down
  if (pa.i[1] < pa.i[length(pa.i)])
  {
    pa.i  <-  rev(pa.i)
    ta.i  <-  rev(ta.i)
    qa.i  <-  rev(qa.i)
  }

  #get the increment of pressure
  d.p <-  (3e4 - 1e4)/nseg

  #initialize the variables
  v.o <-  0
  t.par.old <-  LocLin(ps.i-1e4,pa.i,ta.i)
  
  for (i in 1:nseg)
  {
    #the pressure of the bottom and top
    p.bot <-  ps.i -1e4 - d.p*(i-1)
    #p.mid <-  ps.i -1e4 - d.p*(2*i-1)/2
    p.top <-  ps.i -1e4 - d.p*i

    #find the corresponding t and q
    #[K]
    t.bot <-  LocLin(p.bot,pa.i,ta.i)
    t.top <-  LocLin(p.top,pa.i,ta.i)
    #[Pa]
    q.bot <-  LocLin(p.bot,pa.i,qa.i)
    q.top <-  LocLin(p.top,pa.i,qa.i)

    #----------------------------------------------------
    #--- Get moist adiabatic lapse rate [K/m] and 
    #--- the depth of the layer from lower to upper level
    #----------------------------------------------------
    p.mid <-  PrsIntPol(p.bot,p.top,p.bot,p.top)
    t.mid <-  PrsIntPol(t.bot,t.top,p.bot,p.top)
    q.mid <-  PrsIntPol(q.bot,q.top,p.bot,p.top)
    dz    <-  (p.bot - p.top)/(grav * p.mid/(Rd*t.mid*(1+(q.mid/ep))/(1+q.mid)))
    #calculate saturated specific humidity [kg.kg-1]
    q.sat  <-  SSH(p.mid,t.mid)

    moist.lapse <-  (grav/cp) *
                    (1 + (Lv * q.sat)/(Rd*t.mid)) /
                    (1 + (Lv**2 * q.sat)/(cp*Rv*t.mid**2))

    #----------------------------------------------------
    #--- Get the parcel temperature
    #----------------------------------------------------
    t.par <-  t.par.old - moist.lapse*dz

    #----------------------------------------------------
    #--- Get mid-point temps from environment and parcel
    #----------------------------------------------------
    t.par.mid <-  0.5 * (t.par + t.par.old)
    t.mid     <-  0.5 * (t.bot + t.top)

    #----------------------------------------------------
    #--- Integrate from old increment to increment level
    #----------------------------------------------------
    v.o <-  v.o + (Rd*(t.par.mid - t.mid))*log(p.bot/p.top)

    #----------------------------------------------------
    #--- Update the last increment values 
    #----------------------------------------------------
    t.par.old <-  t.par
    t.bot     <-  t.top
    q.bot     <-  q.top
    p.bot     <-  p.top
  }
  return(v.o)
}

RH2SH <-  function(q.i,t.i,p.i){
  #covert relative humidity to specific humidity
  #require 'humidity' package
  #q.i: humidity [%]
  #t.i: temperature [K]
  #p.i: pressure [hPa]
  es.t  <-  SVP(t.i)
  e.t <-  WVP2(q.i,es.t)
  q.o <-  SH(e.t, p = p.i)

  return(q.o)
}

LocLin  <-  function(r.i,r.a,v.a){
  #get linear interpolate according to reference value and arrays
  #r.i: input reference
  #r.a: reference arrays
  #v.a: value arrays

  if ((r.i > max(r.a)) | (r.i < min(r.a)))
    stop("LocLin: reference input is out of the reference array")

  #find the index besides the r.i
  ind.a <-  which.min(r.a[r.a >= r.i])
  ind.b <-  which.max(r.a[r.a < r.i])

  if (r.a[ind.a] == r.i)
  {
    v.o <-  v.a[ind.a]
  } else
    v.o <-  v.a[ind.a] + (v.a[ind.b] - v.a[ind.a])*
            (r.i - r.a[ind.a])/(r.a[ind.b] - r.a[ind.a])

  return(v.o)
}

PrsIntPol <-  function(v.b,v.t,p.b,p.t)
{
  #pressure based interpolation
  #unit should be [Pa]
  #v.b: value of bottom
  #v.t: value of top
  #p.b: pressure of bottom
  #p.t: pressure of top
  v.o <-  (v.b*log(p.b) + v.t*log(p.t))/
          log(p.b*p.t)

  return(v.o)
}

SSH <-  function(p.i,t.i){
  #calculate saturated specific humidity
  #needs 'humidity' package
  #p.i: input pressure [Pa]
  #t.i: input temperature [K]
  
  es.t  <-  SVP(t.i)
  e.t <-  WVP2(100,es.t)
  v.o <-  SH(e.t,p.i)
  return(v.o)
}


CalPres <-  function(sq.i,rq.i,t.i)
{
  #calculate atmospheric pressure [Pa]
  #sq.i: specific humidity [kg.kg-1]
  #rq.i: relative humidity [%]
  #t.i: temperature [K]
  v.o <-  rq.i*SVP.ClaCla(t.i)*(0.622+0.378*sq.i)/sq.i
  return(v.o)
}

DewT  <-  function(q.i,p.i){
  #dew point temperature [K]
  #q.i: specific humidity [kg.kg-1]
  #p.i: pressure [Pa]

  #some parameters
  a.t <-  610.8
  b.t <-  237.3
  c.t <-  17.2693882

  #--------------------------------------------
  #--- Vapor pressure and convert to Pa
  #--------------------------------------------
  e.t <-  q.i*p.i/(0.622+0.378*q.i)

  #--------------------------------------------
  #--- Vapor pressure and convert to Pa
  #--------------------------------------------
  v.o <-  (log(e.t/a.t)*b.t) / (c.t-log(e.t/a.t)) + 273.15

  return(v.o)
}

CalHilow  <-  function(pa.i,ta.i,qa.i,ps.i){
  #calculate the Hi_low at the given time step
  #pa.i, pressure profile [Pa]
  #ta.i, temperatuure profile [K]
  #qa.i, specific humidity profile [kg/kg]


  #check if profile is upside down
  if (pa.i[1] < pa.i[length(pa.i)])
  {
    pa.i  <-  rev(pa.i)
    ta.i  <-  rev(ta.i)
    qa.i  <-  rev(qa.i)
  }

  #calculate t50, t150, q50, q150
  t50   <-  LocLin(ps.i-5e3,pa.i,ta.i)
  t150  <-  LocLin(ps.i-15e3,pa.i,ta.i)
  q50   <-  LocLin(ps.i-5e3,pa.i,qa.i)
  q150  <-  LocLin(ps.i-15e3,pa.i,qa.i)

  #calculate the dew point temperature (K)
  t.dew50   <-  DewT(q50,ps.i-5e3)
  t.dew150  <-  DewT(q150,ps.i-15e3)

  v.o <-  (t50 - t.dew50) + (t150 - t.dew150)
  return(v.o)
}

acomb   <-  function(...) abind(...,along=4)

ShiftLon    <-  function(v.i){
    d.v <-  dim(v.i)
    demi.lon    <-  d.v[1]/2

    if (length(d.v) == 2)
    {
        v.o <-  v.i
        v.o[1:demi.lon,]   <-  v.i[(demi.lon+1):d.v[1],]
        v.o[(demi.lon+1):d.v[1],]   <-  v.i[1:demi.lon,]
    } else if (length(d.v) == 3)
    {
        v.o <-  v.i
        v.o[1:demi.lon,,]   <-  v.i[(demi.lon+1):d.v[1],,]
        v.o[(demi.lon+1):d.v[1],,]   <-  v.i[1:demi.lon,,]
    } else
    {
        v.o <-  v.i
        v.o[1:demi.lon,,,]   <-  v.i[(demi.lon+1):d.v[1],,,]
        v.o[(demi.lon+1):d.v[1],,,]   <-  v.i[1:demi.lon,,,]
    }
    return(v.o)
}


DiagCTPHilow  <-  function(ctp.i,hi.i){
    n.lon <-  dim(ctp.i)[1]
    n.lat <-  dim(ctp.i)[2]

    v.o <-  ctp.i
    for (nx in 1:n.lon)
    for (ny in 1:n.lat)
    {
      if (is.na(ctp.i[nx,ny]) | is.na(hi.i[nx,ny]))
          next

      if (ctp.i[nx,ny] < 0)
      {
        v.o[nx,ny] <-  1
      } else if ((hi.i[nx,ny] < 5) | (hi.i[nx,ny] > 15))
      {
        v.o[nx,ny] <-  1
      } else if (hi.i[nx,ny] <= 10)
      {
        v.o[nx,ny] <-  2
      } else if (ctp.i[nx,ny] < 180)
      {
        v.o[nx,ny] <-  3
      } else
        v.o[nx,ny] <-  4
    }

    return(v.o)
}


DFF <-  function(sm.pc,v.pc,ai.pc,sm.4x,v.4x,ai.4x,lsm,n.bin=50){
#sm [lon,lat,time]
#vi [lon,lat,time]
#ai [lon,lat]
    bin.pc <-  GetBinMat(sm.pc,v.pc,ai.pc,lsm,n.bin=n.bin)
    bin.4x <-  GetBinMat(sm.4x,v.4x,ai.pc,lsm,n.bin=n.bin)

    return(bin.4x - bin.pc)
}


GetBinMat  <-  function(sm.i,v.i,ai.i,lsm.i,n.bin=50){
    n.lon <-  dim(sm.i)[1]
    n.lat <-  dim(sm.i)[2]
    n.t   <-  dim(sm.i)[3]

    v.o <-  array(NA,dim=c(n.bin,n.bin))

    sm.t <-  array(NA,dim=c(n.lon,n.lat,n.bin))
    vv.t <-  array(NA,dim=c(n.lon,n.lat,n.bin))

    for (i in 1:n.lon)
    for (j in 1:n.lat)
    {
        if (lsm.i[i,j] < .5)
            next

        v.t <-  CalBin(sm.i[i,j,],v.i[i,j,],n.bin=n.bin)
        sm.t[i,j,] <-  v.t$sm
        vv.t[i,j,] <-  v.t$v
    }

    sm.tt <-  matrix(sm.t,prod(dim(sm.t)[1:2]),dim(sm.t)[3])
    vv.tt <-  matrix(vv.t,prod(dim(vv.t)[1:2]),dim(vv.t)[3])
    ai.tt <-  matrix(ai.i,prod(dim(ai.i)[1:2]))

    sm.tt <-  sm.tt[lsm.i > .5,]
    vv.tt <-  vv.tt[lsm.i > .5,]
    ai.tt <-  ai.tt[lsm.i > .5]

    ind.ai  <-  as.numeric(CutQuantile(ai.tt,n.bin,
                                       labels=T))

    for (i in 1:n.bin)
        v.o[,i] <-  apply(vv.tt[ind.ai == i,],2,mean,na.rm=T)


    return(v.o)
}

CalBin  <-  function(sm.i,v.i,n.bin=50){
    sm.o  <-  as.vector(tapply(sm.i,CutQuantile(sm.i,n.bin),mean)) 
    v.o   <-  as.vector(tapply(v.i,CutQuantile(sm.i,n.bin),mean)) 

    v.oo  <-  list(sm=sm.o,v=v.o)
    return(v.oo)
}

CutQuantile <-  function(v.i,n.bin=50,labels=F){
    v.max <-  max(abs(v.i),na.rm=T)
    v.t <-  v.i + rnorm(length(v.i),0,sd=v.max*1e-5) 
    v.o <-  cut(v.t, breaks=c(quantile(v.t,
                              probs = seq(0,1,by=1/n.bin),
                              na.rm = T)),
                include.lowest=TRUE,labels=seq(1,n.bin))
    return(v.o)
}
