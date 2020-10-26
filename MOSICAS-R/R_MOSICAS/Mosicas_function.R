########################################################################################################
########################################################################################################
###############################  MOSICAS FUNCTION ######################################################
########################################################################################################
########################################################################################################

########################################################################################################
#   	Subroutines	Definition
########################################################################################################
#   	emerg		emergence
#   	laiglob	Green leaf area index
#   	intercep	Interception efficiency
#   	convert	Conversion into total dry mass
#   	partdmt	Partition of total dry mass into above ground dry mass
#   	partstem	Partition of above ground dry mass into stem dry mass
#   	partsugar	Partition of stem dry mass into sucrose and structure
#   	humstem	stem sucrose and water contents
#   	humcan		Above ground biomass water content
#   	elong		Stalks: growth of Tvd height
#   	flor		Flowering
#     Astrono
#     counter
########################################################################################################

########################################################################################################
###################### emerg #########################################################################
########################################################################################################
# talcas alive stalks number calculation 
# Tempmn et tempmx: minimum and maximum temperature
# emergtb: base temperature for emergence  (plant parameter)
# emergtt: TT (termal time) to reach emergence (plant parameter)
# emergsdd: TT of emergence
# xdd: degree days of the day
# nbstem: Alive stem per m2 # only 1 or 0 in this version (see Moscas-Fortran for tillering)
########################################################################################################
emerg <- function(tempmn,tempmx,emergsdd,nbstem, emergtb, 
                   emergtt){
  
  xdd=max(c(0.,(tempmn+tempmx)/2-emergtb))
  emergsdd=emergsdd+xdd
  
  if (emergsdd>=emergtt) { nbstem=1 }
  else{nbstem=0.}
  
  OUT <- list(tempmn,tempmx,emergsdd,nbstem, emergtb, 
              emergtt)
  names(OUT) <- c("tempmn","tempmx","emergsdd","nbstem", "emergtb", 
                  "emergtt")
  return(OUT)
}
# END talcas
########################################################################################################


########################################################################################################
###################### laiglob #########################################################################
########################################################################################################
# Green leaf area index calculation
# lai & glai: daily growth in lai ; slai = lai mortality
# tempmn,tempmx
# laitb: lai base temperature
# laigrowth: rate of lai growth
# laiwksen: sensitivity of lai to water stress
# swdflai: soil water stress factor for lai
########################################################################################################
laiglob <- function(tempmn,tempmx,laitb,laigrowth,laiwksen,laimax,swdflai,sddlai,lai,nbstem){
  
  glai = 0 ; slai = 0 ;xddlai=0 
  if (nbstem>0) {
    if (lai==0) { lai=0.05 }
    
    xddlai = (tempmx+tempmn)/2 - laitb ;     xddlai=max(0,xddlai) ; sddlai=sddlai+xddlai

    # Daily LAI growth taking into account LAI competition (gompertz derivative)
    glai=(laigrowth)*lai*(1-lai/laimax)*xddlai #previous version (gompertz derivative)
    # water stress effect on leaf growth
    glai=glai*swdflai
    # leaf mortality due to water stress
    slai=max(c(0,laiwksen*(1-swdflai)*lai))
        
    lai=max(c(0,lai+glai-slai))
  }
  
OUT <-list(tempmn,tempmx,laitb,laigrowth,laiwksen,laimax,swdflai,sddlai,lai,nbstem)
names(OUT) <- c("tempmn","tempmx","laitb","laigrowth","laiwksen","laimax","swdflai","sddlai","lai","nbstem")
  
return(OUT)
}

#end function laiglob0
########################################################################################################

########################################################################################################
###################### intercep ########################################################################
########################################################################################################
# calculation of interception efficiency (ei)
# calculation of intercepted photosynthetical active radiation (pari)
# ke: parameter, extenction coefficient
# solrad, par: incident global radiation and photosynthetical radiation
# lai: green leaf area index
########################################################################################################
intercep <- function(isobei=0,obei=NULL,ke,lai,solrad,ei,par,pari) {
  
  ei = 1 - exp(-(ke * lai))
  
  # if ei is forced, not implemented yet
  if(isobei==1){ei=obei; lai=-log(1-ei)/ke}
  
  par=0.5*solrad; pari=0.5*solrad*ei
  
  OUT <- list(isobei,obei,ke,lai,solrad,ei,par,pari)
  names(OUT) <- c("isobei","obei","ke","lai","solrad","ei","par","pari")
  return (OUT)
}
# End function
########################################################################################################

########################################################################################################
###################### convert ########################################################################
########################################################################################################
# Conversion of daily intercepted photosynthetic radiation (pari) into
# daily total dry mass accumulation (gmst)
# astrgx : Global extraterrestrial  radiation (MJ/m2/d)
# pari: daily intercepted photosynthetic radiation (MJ/m2/d)
# tempmn,tempmx : daily min and max temperatures
# ruemax: conversion coefficient of intercepted photosynthetic radiation into total dry mass (gr/MJ)
# ruetopt: Optimum temperature for conversion
# ruetk: Effect of temperature on conversion
# solrad: Global Incident radiation (MJ/m2/d)
# swdfrue: water stress index for dry mass accumulation (0-1)
# gdmt: daily total dry mass accumulation (g/m2/d)
# dmt: total dry mass (g/m2)
# p01: coefficient for maintenance effect on conversion
# p04: coefficient for diffused radiation effect on conversion
########################################################################################################
convert <- function(ruemax,ruetk,ruetopt,tempmn,                
                     tempmx,pari,solrad,astrgx,swdf1,gdmt,dmt,p01,p04){
  
  gdmt=0; kage=1; ktemp=1
  
  # Average temperature during daylight
  tdiu=0.25*tempmn+0.75*tempmx
  # Temperature effect on rue
  ktemp = 1-ruetk*(abs(tdiu-ruetopt))^2#0.5   attention dans mosicas serait 0.5 ?!? CHANGER
  ktemp= max(c(0,min(c(1,ktemp))))
  # Diffuse radiation effect
  rap=min(c(1.0, solrad/astrgx))
  ruex=ruemax
  pari=pari*(1-p04*(rap-0.75)) 
  # Maintenance effect (kage)
  kage=1-(dmt*p01/10000)^3
  kage=max(c(0,min(1,kage)))
  # Daily increase in total dry mass
  gdmt=ruex*kage*pari*swdfrue*ktemp
  dmt=dmt+gdmt
  
  OUT <- list(ruemax,ruetk,ruetopt,tempmn,                
              tempmx,pari,solrad,astrgx,swdfrue,gdmt,dmt,p01,p04)
  names(OUT) <- c("ruemax","ruetk","ruetopt","tempmn",                
                     "tempmx","pari","solrad","astrgx","swdfrue","gdmt","dmt","p01","p04")
  return(OUT)
}
# END function
########################################################################################################

########################################################################################################
###################### partdmt ########################################################################
########################################################################################################
# Partition of total dry mass into above ground biomass
# ratoon: crop cycle (0: plant crop, 1: first ratoon, ..)
# gdmt: daily total dry mass accumulation (g/m2/d)
# dmt: total dry mass (g/m2)
# sdd: Thermal time since planting or previous harvesting (base temperature=0)
# gdmroot: daily root dry mass accumulation (g/m2/d) (g/m2/d)
# dmroot: root dry mass (g/m2)
# gdmbaer: daily above ground dry mass accumulation (g/m2/d) (g/m2/d)
# dmbaer: above ground dry mass (g/m2)
# prootini : initial allocation to roots
# prootend : final allocation to roots
# prootdec : attenuation coefficient for allocation to roots
########################################################################################################
partdmt <- function(ratoon,gdmt,sdd, dmt,                       
                     gdmroot,dmroot,gdmbaer,dmbaer,
                    prootini,prootdec,prootend){

  gdmrac=0; gdmbaer=0
  
  if (ratoon>=1) {kroot=prootend}
  
  if (ratoon <1) {
    kroot=max(c(prootend,(prootini-prootdec*sdd)))
  }
  gdmroot=kroot*gdmt
  dmroot=dmroot+gdmroot
    
  if ((dmt-dmroot)>500) { # allocation to above ground mass
    gdmbaer=(gdmt-gdmroot)}
  if ((dmt-dmroot)<=500) {
    gdmbaer=(0.7+(dmt-dmroot)/500*0.3)*(gdmt-gdmroot)}

  dmbaer=dmbaer+gdmbaer
  
  OUT <- list(ratoon,gdmt,sdd, dmt, gdmroot,dmroot,gdmbaer,dmbaer,prootini,prootdec,prootend)
  names(OUT) <- c("ratoon","gdmt","sdd", "dmt", "gdmroot",
                  "dmroot","gdmbaer","dmbaer","prootini","prootdec","prootend")
  return(OUT)
}
# END partmst
########################################################################################################

########################################################################################################
###################### partstem #######################################################################
########################################################################################################
# Partition of above ground dry mass into millable stalk dry mass
# gdmbaer: daily above ground dry mass accumulation (g/m2/d) (g/m2/d)
# dmbaer: above ground dry mass (g/m2)
# jat: age when stem appearance
# pstemini:Beginning of stem dry mass appearance (Plant parameter (g/m2)
# pstemdec: Extinction coefficient of  daily fraction of aboveground dry mass
#             allocated to stem  (Plant parameter)
# pstemend:	Final  daily fraction of aboveground dry mass allocated to
#              stem (Plant parameter)
# gdmstm: daily increase in stem dry mass (g/m2/d)
# dmstm: stem dry mass (g/m2)
########################################################################################################
partstem <- function (pstemini,pstemend,pstemdec,dmbaer,gdmbaer,  
                       jat,gdmstm,dmstm){
  
  if (dmbaer>pstemini) { 
    jat=jat+1
    gdmstm=pstemend*(1-exp(-pstemdec/1000*(dmbaer-pstemini)))*gdmbaer}
  else {
    jat=0;gdmstm=0
  }
  dmstm=dmstm+gdmstm
  
  OUT <- list(pstemini,pstemend,pstemdec,dmbaer,gdmbaer,  
              jat,gdmstm,dmstm)
  names(OUT) <- c("pstemini","pstemend","pstemdec","dmbaer","gdmbaer",  
                     "jat","gdmstm","dmstm")
  return(OUT)
}
# END parttige
########################################################################################################

########################################################################################################
###################### partsugar ######################################################################
########################################################################################################
# Partition of millable stalk dry mass into sucrose and structure
# gdmstm: daily increase in stem dry mass (g/m2/d)
# dmstm: stem dry mass (g/m2)
# tempmn,tempmx : daily min and max temperatures
# pstrudec	Extinction coefficient of  daily fraction of stem dry mass
#             allocated to structures   Plant parameter
# psugend	Final daily fraction of stem dry mass allocated to sucrose
#             Plant parameter  g/g
# pstrutb	Temperature treshold from which fraction of stem dry mass
#             allocated to structures is decreasing.  Plant parameter
# pstrutgrowth	Temperature effect on daily fraction of stem dry mass
#             allocated to structures  /°C.	  Plant parameter
# gdmsug: daily increase in succrose mass (g/m2/d)
# dmsug: succrose mass (g/m2)
# swdfrue: water stress index for dry mass accumulation (0-1)
# swdflai: water stress index for growth (0-1)
########################################################################################################
partsugar <- function(tempmn, tempmx,gdmstm,dmstm,swdfrue,swdflai,   
                       pstrutb,pstrutgrowth,psugend,pstrudec,gdmsug,dmsug){
  
  gdmstst=gdmstm
  tstru=1
  tmo=(tempmn+tempmx)/2
  
  if (tmo>pstrutb) {
    tstru=1+pstrutgrowth*(tmo-pstrutb)
    }  
  
  gdmstst=(1-psugend*(1-exp(-pstrudec/1000*dmstm)))*gdmstm 
  gdmstst=gdmstst*tstru*(1-(swdfrue-swdflai))

  if (gdmstst>gdmstm) {gdmstst=gdmstm}
  gdmsug=gdmstm-gdmstst
  dmsug=dmsug+gdmsug
  
  OUT <- list(tempmn, tempmx,gdmstm,dmstm,swdfrue,swdflai,   
              pstrutb,pstrutgrowth,psugend,pstrudec,gdmsug,dmsug)
  names(OUT) <- c("tempmn", "tempmx","gdmstm","dmstm","swdfrue","swdflai",   
                  "pstrutb","pstrutgrowth","psugend","pstrudec","gdmsug","dmsug")
  return(OUT)
}
# END SUBROUTINE
########################################################################################################

########################################################################################################
###################### humstem #######################################################################
########################################################################################################
# Calculation of stem fresh mass, sucrose and water contents
# dmstm: stem dry mass (g/m2)
# hum_stm millable stalk water content (%)
# dmsug: succrose mass (g/m2)
# yldcan, Yield or millable stalk fresh mass (T/Ha)
# sug_stm, millable stalk fresh mass sucrose content (%)
# jat: age when millable stalk appearance
########################################################################################################
humstem <- function(dmstm,hum_stm,dmsug,fmstm,sug_stm,jat,humstemini,humstemdec){
  
  if (dmstm>0) {
    hum_stm = humstemini  - humstemdec * jat 
    #hum_stm = 88 - 0.070374 * jat  #rmse=8.14 Hum%tige=85.77-0.603*agetu (R570, Rep, CalReu) 121107o  (jfm121109)
    hum_stm=min(c(hum_stm,90)); hum_stm=max(c(hum_stm,68))
    fmstm=dmstm/(1-hum_stm/100)
    
    if (fmstm>0.1) {sug_stm=dmsug/fmstm}
  }
  
  OUT <- list(dmstm,hum_stm,dmsug,fmstm,sug_stm,jat,humstemini,humstemdec)
  names(OUT) <- c("dmstm","hum_stm","dmsug","fmstm","sug_stm","jat","humstemini","humstemdec")
  return(OUT)
}  
# END hustem
########################################################################################################

########################################################################################################
###################### humcan #########################################################################
########################################################################################################
# Above ground biomass water content
# humaerini: initial ground biomass water content % Plant parameter
# humaerend: begining in decrease of ground biomass water content  Plant parameter
# humaertb: Base temperature for ground biomass water content  Plant parameter
# humaerdec: Rate of decrease in ground biomass water content  Plant parameter
# hum_aer: above ground biomass water content
# tempmn,tempmx : daily min and max temperatures
# humaersdd: Thermal time for ground biomass water content
# dmbaer: above ground dry mass (g/m2)
########################################################################################################

humcan <- function(humaerini, humaerend, humaertb, humaerdec,       
                   hum_aer,tempmx,tempmn,humaersdd,dmbaer){
  
  if (dmbaer<=0) {humaersdd=0; hum_aer=humaerini}
  else {
    dd=max(c(0,(tempmx+tempmn)/2-humaertb))
    humaersdd=humaersdd + dd
    
    if (humaersdd<humaerend) {hum_aer=humaerini}
    else {hum_aer=hum_aer - dd * humaerdec}
    hum_aer=max(50.,hum_aer)
  }
  
  OUT <- list(humaerini, humaerend, humaertb, humaerdec,       
              hum_aer,tempmx,tempmn,humaersdd,dmbaer)
  names(OUT) <- c("humaerini", "humaerend", "humaertb", "humaerdec",       
                  "hum_aer","tempmx","tempmn","humaersdd","dmbaer")
  return(OUT)
}
# END humcan
########################################################################################################

########################################################################################################
###################### Astrono ##########################################################################
########################################################################################################

Astrono <- function(jcal,lati,longi,alti,astlatrad,astdr,
                    astdecl,astomeg,astrgo,astrgx,astdurjour,astdurjourp){
  
  astlatrad= lati * 3.14159265358979 / 180                
  astdr=1+0.033*cos(0.0172*jcal)                            
  astdecl=0.409*sin(0.0172*jcal-1.39)                     
  astomeg=acos(-tan(astlatrad)*tan(astdecl))                
  astrgo=cos(astdecl)*cos(astlatrad)*sin(astomeg)         
  astrgo= 37.6*astdr*(astomeg*sin(astdecl)*sin(astlatrad)+astrgo)
  astdurjour=7.64*astomeg                                 
  astrgx=(0.75+0.00002*alti)*astrgo                     
  
  OUT <- list(jcal,lati,longi,alti,astlatrad,astdr,
              astdecl,astomeg,astrgo,astrgx,astdurjour,astdurjourp)
  names(OUT) <- c("jcal","lati","longi","alti","astlatrad","astdr",
                  "astdecl","astomeg","astrgo","astrgx","astdurjour","astdurjourp")
  return(OUT)
}
########################################################################################################



counter <- function (i) {
  if (i%%10 == 0) 
    cat(i)
  else cat(".")
  if (i%%50 == 0) 
    cat("\n")
  flush.console()
}


