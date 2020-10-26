##########################################
#### Mosicas water balance
##########################################



if (Input_read$WATMODEL %in% c('MOSICASWATBAL')){
  
  ### using the Ceres2wf submodel
  ### (after Gabrielle et al, SSSAJ 59: 1403-1412, 1995)
  ### avec modIFications pour la canne (routine de Transpiration)
  ### version avec une seule couche pour les propriétés hydro-dynamiques

  swc2=swc
  ### initialization
  stocki <- sum(swc*dlayer*10) #initiale water stock at the beginning of the day
  
  
  ### Modification of soil characteristics for a specific soil, based on Battie-Laclau et ., XXXX
#  if(!is.na(Soilgen_input$soiltype)){
#  if (Soilgen_input$soiltype=="latosol") {
#    id <- which((swc-swcfc) > 0.2*(swcsat-swcfc))
#    swc[id] <- swcfc[id]+0.2*(swcsat[id]-swcfc[id])  
#  }}
  
  drain=0.
  deepdrain=0
  runoff=0.
  rain <- Weather_Irrig_data$rr[iday]
  irrdose <- Weather_Irrig_data$irrig[iday]
  
  # runoff, currently assumed 0
  # runoff = rain + irrdose - (swcsat_mm - stocki) ; runoff <- max(c(0,runoff))
  raininf = rain + irrdose - runoff
  
  ###################################
  ####### water infiltration ########
  # assuming all incoming water can infiltrate
  # and excess of SAT is by-pass flow
  
  # water (mm) infiltrating each layer i from layer i-1 above, initialization
  fluxin <- rep(0,nlayer) #drl1
  # water quantity (mm) that can infiltrate in each layer
  limflux <- swcsat*dlayer*10-swc*dlayer*10 ; limflux[which(limflux<0)] <- 0
  # cumulative
  cumlimflux <- sapply(1:nlayer,function(i) sum(limflux[1:i]))
  # defining depth of filled layers
  id <- which(cumlimflux > raininf)
  if(length(id)==0){
    fluxin = limflux
    deepdrain = deepdrain + raininf - sum(limflux)
  }
  if(length(id)>0){
    if(id[1]==1){fluxin[1] = raininf} # only first layer filled
    if(id[1] >1){
      fluxin[id[1]] <- raininf - sum(cumlimflux[id[1]-1]) 
      fluxin[1:(id[1]-1)] <- limflux[1:(id[1]-1)]
    }
    deepdrain=deepdrain+0
  }
  
  swc <- swc + fluxin/(dlayer*10)
  
  
  ####################################
  ####### gravimetric drainage #######
  # computing drain using Gardner's approximation in Gabrielle et al., 1995. Drain occurs under gravimetry forces only
  # and each layer loses water above field capacity independently (within a day) from what may drip from the upper layer
  # K = hydraulic conductivity (cm/d)
  # Ks saturated hydraulic conductivity (cm/d)
  # A dimensionless parameter
  K = Ks * exp(A*(swc-swcfc))
  # No downwards flux when swc > field capacity
  K[which(swc < swcfc)] <- 0
  # downwards flux, mm/d
  fluxout <- K*10 # downwards flux in mm
  # Flux from layer i can't be higher that swc-swcfc of layer i 
  limflux <- (swc-swcfc)*dlayer*10 ; limflux[which(limflux<0)] <- 0
  id <- which(fluxout>limflux)
  fluxout[id] <- limflux[id]
  # Flux from layer i can't be higher than swcsat - swc of layer i+1
  # there is no constraints on the deepest layer
  limflux <- c(((swcsat-swc)*dlayer*10)[2:nlayer],10^6) ; limflux[which(limflux<0)] <- 0
  id <- which(fluxout>limflux)
  fluxout[id] <- limflux[id]
  # new soil water content
  swc <- swc - fluxout/(dlayer*10)
  swc[2:nlayer] <- swc[2:nlayer] + fluxout[1:(nlayer-1)]/(dlayer[2:nlayer]*10)
  # deep drainage
  deepdrain = deepdrain + fluxout[nlayer]
  

  ##########################################
  ####### potential soil evaporation ############
  etp <- Weather_Irrig_data$etp[iday]
  # eos = potential soil evaporation
  if (lai <= 1){eos = etp*(1-0.43*lai)} #ceres original
  if (lai > 1){ eos = etp/1.1 * exp(-0.4*lai)} #ceres original, in STICS it's etp * exp(-(ke-0.2)*LAI)
  
  
  ############ actual soil evaporation #####################
  # Based on Ritchie 1972
  # three phase: 1st phase es = potential ; transition phase the soil start to dry (es reduced) ; 2nd phase the soil is dry (es reduced)
  # where define the transition limit between the 1st and 3rd phase
  #U=5.0 # mm
  winf=raininf
  es=0
  
  # phase 1 wet soil and transition
  if((sumes1<U)|((sumes1>=U)&(winf>=sumes2))){

    if(sumes1<U) {  sumes1 = max(c(0,sumes1-winf)) } # the soil is wet (sumes1<U)
    if((sumes1>=U)&(winf>=sumes2)){ sumes1 = max(c(0,U - (winf - sumes2)))} # the soil was dry but became wet with rain
    
    # the soil is drying
    sumes1=sumes1 + eos
    
    if(sumes1 <= U){ 
      es = eos
      sumes2=0} # the soil is still wet
    if(sumes1 > U) {  # the soil reached the dry state this specific day 
      es = eos-0.4*(sumes1-U)
      sumes2=0.6*(sumes1-U)
      t=(sumes2/alphaes)^2
      }
  }
  
  # phase 2, dry soil
  if((sumes1>=U)&(winf<sumes2)){ # the soil is dry and remain dry
    t=t+1
    es= alphaes*sqrt(t)-sumes2
    
    if(winf>0){ # particular case when it rains a little
      esx=0.8*winf
      if(esx<es){esx=es+winf}
      es=min(esx,eos)}
    if(winf==0){es=min(es,eos)}
    
    sumes2=sumes2+es-winf
    t=(sumes2/alphaes)^2
  }
  

    # new soil water content taking into account the fact that the swc of the first layer can go below swcwp
  swc[1] = swc[1] - es/(dlayer[1]*10)
  
  if(swc[1] >= swcwp[1]*swef){es = es}
  if(swc[1] < swcwp[1]*swef){ # very dry soil
    es = es - (swcwp[1]*swef-swc[1])*dlayer[1]*10
    swc[1] = swcwp[1]*swef
  }
  
  
  
  ######## Unsaturated flow below drained upper limit  # ALL in Gabrielle et al., 1995
  ########   k-theta relation ### CAPILLARY RISE !!!!!!!!!!!!!!!!!!!!
  ######## and van Keulen and Wolf (eds, 1986), for the
  ######## h-theta relation       (dbar in cm/d)
  ########       ....in order to apply a coherent Darcy-Richards'law...
  Flow <- rep(0,nlayer) # upwards capillary rising to layer i
  swctemp <- swc

  ini=1 ; if(dlayer[1]<5){ini=2} # avoid issue with very small layer
  thet1 <- swc[ini:(nlayer-1)]-swcsat[ini:(nlayer-1)] 
  thet2 <- swc[(ini+1):nlayer]-swcsat[(ini+1):nlayer]    
  K = Ks * exp(A*(thet1+thet2)/2) ; 
  K[ which(K>Ks)] <- Ks
  dlayerAvg <- (dlayer[ini:(nlayer-1)]+dlayer[(ini+1):nlayer])/2
  psi1 <- exp(sqrt(abs(log(swcsat[ini:(nlayer-1)]/swc[ini:(nlayer-1)]))/gam)) # cm
  psi2 <- exp(sqrt(abs(log(swcsat[(ini+1):nlayer]/swc[(ini+1):nlayer]))/gam)) # cm
  
  Flow[ini:(nlayer-1)] <- K*(psi1-psi2-dlayerAvg)/dlayerAvg # cm
  Flow[1] <- min(c(Flow[1],(swcfc[1]-swc[1])*dlayer[1]))
  Flow[which(Flow<0)] <- 0
  
  swctemp = swctemp+Flow/dlayer
  swctemp[2:nlayer] = swctemp[2:nlayer] - Flow[1:(nlayer-1)]/dlayer[2:nlayer]
  
  # deep capillary rise
  # to activate if there is a water table
  #Limit <-c(250,500,1000,2500,5000,10^10)
  #FlowWaterTable <- c(0,0.03,0.05,0.06,0.08,0.09) # cm
  #psi1 <- exp(sqrt(log(swcsat[nlayer]/abs(swc[nlayer]))/gam)) # cm
  #Flow[nlayer] <- FlowWaterTable[which(psi1<Limit)[1]]
  #swctemp[nlayer] <- swctemp[nlayer] + Flow[nlayer]/dlayer[nlayer]
  
  # endding upwards flux
  swc <- swctemp
  
  #CES=CES+ES ?
  
  ###################################################
  ###### Transpiration ##############################
  
  # no transpiration after maturity KEEP ?
  #if(istage>4){ep=0}
  
  ep = etp *PARAM$kcmax * sqrt(min(c(1,lai/3)))
  #ep = etp *(1+(PARAM$kcmax-1)/(1+exp(-1.5*(lai-3)))) # stics model alternative but less convincing
  
#  if(!is.na(Soilgen_input$soiltype)){
#  if (Soilgen_input$soiltype=="latosol") {
#    ep=etp*(PARAM$kcmax*1.03)*(1-exp(-0.5*lai)) # à creuser
#  }}
  
  if(lai>0){ ep = min(c(PARAM$kcmax*etp-es),ep)  }

  # maximum transpiration in mm/d
  tmp = ep
  # maximum evapotranspiration in mm/d
  etm=ep+es
  
  # root growth and depth

  # Daily increase in root front
  if(tmo >roottb){ rootfront = rootfront+ (tmo-roottb)*rootgrowth} # root front in cm
  if(rootfront > rootfrontmax){ rootfront = rootfrontmax}
  
  # daily growth in root density (rlv, cm/cm3) = drlv
  drlv = rep(1.6*10^-2 * rldmax * (max(c(0,tmo-roottb))*rootgrowth), nlayer)
  
#  if(!is.na(Soilgen_input$soiltype)){
#  if (Soilgen_input$soiltype=="latosol") {
#    drlv = rep(0.8*10^-2 * rldmax * ((tmo-roottb)*rootgrowth), nlayer)
#  }}
  
  # no growth below root front
  dep <- sapply(1:nlayer, function(i) sum(dlayer[1:i])-dlayer[i])
  drlv[which(dep>rootfront)] <- 0
  
  # soil water stress
  esw <- swcsat-swcwp
  swdf <- rep(1,nlayer)
  idswdf <- which((swc-swcwp) < 0.25*esw)
  swdf[idswdf] <-  4*(swc[idswdf]-swcwp[idswdf])/esw[idswdf]
  swdf <- replace(swdf,swdf<0,0)
  swdf <- replace(swdf,swdf>1,1)
  
  # new root density rlv
  rlv <- rlv + drlv*swdf
  rlv <- replace(rlv,rlv>rldmax, rldmax)
  
  #### root water uptake 
  # rwu = root water uptake in cm water / cm soil / cm roots      
  rwu <- 2.67*10^-3*exp(62.*(swc-swcwp))/(6.68-log(rlv))
  # rwumax = maximum root water uptake
  rwu <- replace(rwu, rwu > rwumax,rwumax)
  rwu <- replace(rwu, swc < swcwp,0)
  
  # quantitative water uptake per layer
  rwu_cm = rwu*dlayer*rlv
  rwu_cm[which(rwu_cm > dlayer*(swc-swcwp))] <- (dlayer*(swc-swcwp))[which(rwu_cm > dlayer*(swc-swcwp))]
  rwu_cm <- replace(rwu_cm,rwu_cm<0,0)
  trwu = sum(rwu_cm)*10 # in mm/d
  trwu = max(c(trwu,0.001)) # to avoid bug when high soil evaporation in the 1st layer
  
  # water use and soil water deficit factor
  # reduction of water use when ep < trwu
  wuf = min(c(ep/trwu,1))
  # actualisation of rwu
  rwu_cm = rwu_cm*wuf
  # new soil water content
  swc <- swc - rwu_cm/dlayer
  # total soil water content, mm
  #tsw <- sum(swc*dlayer*10)
  #pesw <- tsw - sum(swcwp*dlayer*10) 
  
  # soil water deficit factor
  # swdf2 for LAI
  # swdf1 for RUE
  swdflai = 1 ; swdfrue = 1
  if(ep>0){
    if(trwu/ep < 1.5){swdflai = 0.67*trwu/ep}
    swdfrue = min(c(1,trwu/ep))
  }
  # actual plant transpiration, mm
  ep = trwu*wuf
  etr = es + ep # actual evapotranspiration
  trp = ep # actual plant transpiration
  
  # Check of water balance
  # finale stock
  stockf= sum(swc*dlayer*10)
  
  error=stockf-stocki+(etr+deepdrain-raininf)
  if(abs(error)>0.1){
    print(paste(iday,"Warning: error in water balance",sep=' '))
    deepdrain = deepdrain-error
    error=stockf-stocki+(etr+deepdrain-raininf)
  }
  
  # cumulative values
  depstockfr <- replace(rlv,rlv>0,1)
  id <- which(depstockfr==0)[1]-1
  depstockfr[id] <- (rootfront-dep[id])/dlayer[id]
  
  swcwpr <- sum(swcwp*dlayer*10*depstockfr)
  swcsatr <- sum(swcsat*dlayer*10*depstockfr)
  swcfcr <- sum(swcfc*dlayer*10*depstockfr)
  stockfr<- sum(swc*dlayer*10*depstockfr) - swcwpr
  
  #available stock, mm
  astock = stockf - Calcsoilwatercapacities$swcwp_mm
  astockfr = stockfr
  
  # global soil water deficit
  swdef = min(c(1,stockfr/(swcfcr-swcwpr))) 
  swdef = max(c(0.001,swdef)) # avoid bug later

  
}


if (Input_read$WATMODEL %in% c('MOSICASFAO')){}

if ((Input_read$WATMODEL != 'MOSICASFAO') & (Input_read$WATMODEL != 'MOSICASWATBAL')){
  print('Warning no water balance model selected')
}

##############################################################################################
if (Input_read$IWATBAL==0){swdfrue=1;swdflai=1;swdef=1}


##############################################################################################
# actualisation of water stress factor based on Brisson 1992 and  Slabber 1982  8.5  (Ozier 12)
if(lai>0){
  srue=max(c(0,(0.94-0.26*psilb/tmp)))
  swdfrue=1 ; swdf2 = lai
  if(swdef <= 1.5*srue){swdflai = swdef/(1.5*srue)}
  if(swdef <= srue){swdfrue = swdef/(srue)}
  
  swdflai = swdflai^PARAM$sthydlai
  swdfrue = swdfrue^PARAM$sthydrue
}


