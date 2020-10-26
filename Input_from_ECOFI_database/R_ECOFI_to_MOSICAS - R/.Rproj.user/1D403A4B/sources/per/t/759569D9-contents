######## trait input file ###########
INPUTTRAITMOS <- function(TREATCODE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE=''){
  TREATCODE <- sort(as.vector(TREATCODE))
  MATPC <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_plotcycle.txt',header = T,sep='\t')
  MATVC <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_varcycle.txt',header = T,sep='\t')
  MATP <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_plot.txt',header = T,sep='\t')
  MATT <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_trial.txt',header = T,sep='\t')
  MATWT <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_wstrial.txt',header = T,sep='\t')
  MATGS <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_genesoil.txt',header = T,sep='\t')
  head(MATVC)
  
  irricode <- sapply(1:length(TREATCODE), function(i) MATPC$irricode[MATPC$plotcyclecode == TREATCODE[i]])
  irricode <- as.vector(irricode)
  inisoilvol <- sapply(1:length(TREATCODE), function(i) MATPC$inisoilvol[MATPC$plotcyclecode == TREATCODE[i]])
  cycle <- sapply(1:length(TREATCODE), function(i) MATPC$cycle[MATPC$plotcyclecode == TREATCODE[i]])
  varcode <- sapply(1:length(TREATCODE), function(i) MATVC$varcode[MATVC$plotcyclecode==TREATCODE[i]])
  
  plotcode <- as.vector(sapply(1:length(TREATCODE), function(i) MATPC$plotcode[MATPC$plotcyclecode == TREATCODE[i]]))
  
  
  nbbgi <- sapply(1:length(TREATCODE), function(i) MATP$nbbgi[MATP$plotcode==plotcode[i]])
  rowspace <- sapply(1:length(TREATCODE), function(i) MATP$rowspace[MATP$plotcode==plotcode[i]])
  
  trialcode <- as.vector(sapply(1:length(TREATCODE), function(i) MATP$trialcode[MATP$plotcode==plotcode[i]]))
  
  trialalt <- sapply(1:length(TREATCODE), function(i) MATT$trialalt[MATT$trialcode==trialcode[i]])
  trialat <- sapply(1:length(TREATCODE), function(i) MATT$trialat[MATT$trialcode==trialcode[i]])
  trialong <- sapply(1:length(TREATCODE), function(i) MATT$trialong[MATT$trialcode==trialcode[i]])
  soilcode <- as.vector(sapply(1:length(TREATCODE), function(i) MATT$soilcode[MATT$trialcode==trialcode[i]]))
  rootdepth <- sapply(1:length(TREATCODE), function(i) MATGS$rootdepth[MATGS$soilcode==soilcode[i]][1])
    
  wstacode <- unlist(sapply(1:length(TREATCODE), function(i) MATWT$wscode[MATWT$trialcode==trialcode[i]][1]))
  
  begdate <- sapply(1:length(TREATCODE), function(i) MATVC$cyclestartingdate[MATVC$plotcyclecode==TREATCODE[i]])
  enddate <- sapply(1:length(TREATCODE), function(i) MATVC$cyclendingdate[MATVC$plotcyclecode==TREATCODE[i]])
  begdate <- format(as.Date(begdate,format='%d/%m/%Y'), "%d-%m-%Y")
  enddate <- format(as.Date(enddate,format='%d/%m/%Y'), "%d-%m-%Y")
  
  treatcode <- TREATCODE
  
  MATOUT <- data.frame(nbbgi,treatcode,cycle,varcode,rowspace,begdate,enddate,trialat,trialong,
                       trialalt,soilcode,inisoilvol,rootdepth,irricode,wstacode)  
  write.table(MATOUT,file = paste(SETOUT,'/TREATMENT_Mosicas',IDCODE,'.csv',sep=''),col.names=T,
              row.names=F,sep=';',quote=F,dec='.')
}

######## weather station input file ###########
INPUTWSMOS <- function(TREATCODE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE=''){
  TREATCODE <- sort(as.vector(TREATCODE))
  MATWD <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_wdataday.txt',header = T,sep ="\t")
  MATWDec <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_wdatadec.txt',header = T,sep ="\t")
  MATWS <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_ws.txt',header = T,sep ="\t")
  MATTRAIT <- read.table(paste(SETOUT,'/TREATMENT_Mosicas',IDCODE,'.csv',sep=''),header=T,stringsAsFactors = F,sep=';',dec='.')
  Wscode <- levels(as.factor(MATTRAIT$wstacode))
  
  MATWS <- MATWS[which(MATWS$wscode %in% Wscode),]
  MATWSOUT <- data.frame(wstacode=MATWS$wscode,
                         wsalt=MATWS$wsalt,wslat=MATWS$wslat,wslong=MATWS$wslong,
                         TAV=rep(NA,dim(MATWS)[1]),AMP=rep(NA,dim(MATWS)[1]),
                         countrycode=MATWS$countrycode,
                         wstype=MATWS$wstype,
                         wsname=MATWS$wsname)
  write.table(MATWSOUT,file = paste(SETOUT,'/WEATHER_STATION_Mosicas.csv',sep=''),col.names=T,
              row.names=F,sep=';',quote=F,dec='.')
  
  IDLACK <-c()
  for (i in 1:length(Wscode)){
    MATTEMP <- MATWD[MATWD$wscode==Wscode[i],]
    if(dim(MATTEMP)[1]==0){
      MATTEMPB <- MATWDec[MATWDec$wscode==Wscode[i],]
        if(dim(MATTEMPB)[1]==0){
          IDLACK <- c(IDLACK,i)
          MATTEMPB <- data.frame(decade=1:36,tmind=rep(0,36),tmaxd=rep(0,36),
                                 tmeand=rep(0,36),rhmind=rep(NA,36),rhmaxd=rep(NA,36),
                                 windmeand=rep(0,36),windmaxd=rep(0,36),radmeand=rep(0,36),
                                 etomeand=rep(0,36),rainfallmeand=rep(0,36),wscode=rep(Wscode[i],36))
          }
      
      YEARMIN <- min(as.numeric(format(as.Date(MATTRAIT[MATTRAIT$wstacode==Wscode[i],]$begdate,format='%d-%m-%Y'),'%Y')))
      YEARMAX <- max(as.numeric(format(as.Date(MATTRAIT[MATTRAIT$wstacode==Wscode[i],]$enddate,format='%d-%m-%Y'),'%Y')))
      #1988 à 2010
      MATTEMPT <- MATTEMPB
      for (j in 1:(YEARMAX-YEARMIN+1)){MATTEMPT <- rbind(MATTEMPT,MATTEMPB)}

      intday <- ((1:dim(MATTEMPT)[1])-1)*10+1
      
      tntemp <- approx(intday,MATTEMPT$tmind,xout=1:max(intday))$y
      txtemp <- approx(intday,MATTEMPT$tmaxd,xout=1:max(intday))$y
      rgtemp <- approx(intday,MATTEMPT$radmeand,xout=1:max(intday))$y
      tntemp <- approx(intday,MATTEMPT$tmind,xout=1:max(intday))$y
      etptemp <- approx(intday,MATTEMPT$etomeand,xout=1:max(intday))$y
      rrtemp <- approx(intday,MATTEMPT$rainfallmeand,xout=1:max(intday))$y
      MATTEMP <- data.frame(weatherdate=as.Date(paste('01/01/',YEARMIN,sep=''),format='%d/%m/%Y')+(1:max(intday))-1,
                            eto=etptemp,tmax=txtemp,tmin=tntemp,grad=rgtemp,
                            rhmax=rep(NA,max(intday)),rhmin=rep(NA,max(intday)),
                            windtot=rep(NA,max(intday)),tmoy=(tntemp+txtemp)/2,
                            rainfall=rrtemp)
    }
    
    date <- MATTEMP$weatherdate ; date <- format(as.Date(date,format='%d/%m/%Y'), "%d-%m-%Y")
    rr <- MATTEMP$rainfall ; 
    etp <- MATTEMP$eto ; 
    tx <- MATTEMP$tmax ; 
    tn <- MATTEMP$tmin ; 
    rg <- MATTEMP$grad ; 
    ux <- MATTEMP$rhmax ;
    un <- MATTEMP$rhmin ;
    vt <- MATTEMP$windtot ;
    tm <- MATTEMP$tmoy ; 
    MATOUT <- data.frame(date,rr,etp,tx,tn,rg,ux,un,vt,tm)
    MATOUT <- MATOUT[order(as.Date(MATOUT$date,format= "%d-%m-%Y")),]
    
    # is there some missing days ?
    #dim(MATTEMP)
    #NDAY = as.numeric(as.Date(MATOUT$date[length(date)],format= "%d-%m-%Y")-as.Date(MATOUT$date[1],format= "%d-%m-%Y"))+1
    
    # cleaning data, remove NA to replace by the average the same day in all years
    MATOUTb <- MATOUT
    TEST <- c('rr','etp','tx','tn','rg')
    for (itest in 1:length(TEST)){
      ID <- which(is.na(MATOUT[,TEST[itest]]))
      if(length(ID)>0){
        inidate <- sapply(1:length(date), function(j) paste0('01-01-',format(as.Date(date[j],format='%d-%m-%Y'),'%Y')))
        juldate <- as.numeric(as.Date(date,format='%d-%m-%Y')-as.Date(inidate,format='%d-%m-%Y'))+1
        MATOUTb[ID,TEST[itest]] <- sapply(1:length(ID), function(j) mean(MATOUT[which(juldate==juldate[ID[j]]),TEST[itest]],na.rm=T))
        #plot(MATOUT[ID,TEST[itest]])
      }
    }
    
    
    write.table(MATOUTb,file = paste(SETOUT,'/',Wscode[i],'.txt',sep=''),col.names=T,
                row.names=F,sep='\t',quote=F)
    
  }
  print(Wscode[IDLACK])
}

######## Irrigation input file ###########
INPUTIRRIGMOS <- function(TREATCODE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE=''){

  MATTRAIT <- read.table(paste(SETOUT,'/TREATMENT_Mosicas',IDCODE,'.csv',sep=''),header=T,stringsAsFactors = F,sep=';',dec='.')
    irricode <- unique(as.vector(MATTRAIT$irricode)) ; irricode <- irricode[-which(irricode=='RAINFED')]
  MATIRRIGATION <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_irrigation.txt',header = T,sep ="\t")
    MATTEMP <- data.frame(irricode,efficiency=MATIRRIGATION$efficiency[which(MATIRRIGATION$irricode %in% irricode)])
  
  MATIRRDATA <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_irridata.txt',header = T,sep ="\t")
  MATIRRDATA <- MATIRRDATA[which(MATIRRDATA$irricode %in% irricode),]
  
  dateirri <- as.Date(MATIRRDATA$irridate,format='%d/%m/%Y')
  dateirri <- format(dateirri, "%d-%m-%Y")
  
  MATOUT <- data.frame(irricode=MATIRRDATA$irricode,
                       irritype = rep(1,dim(MATIRRDATA)[1]),
                       irridose = MATIRRDATA$irridose,
                       dateirri=dateirri)
  MATOUT <- merge(MATOUT,MATTEMP,by='irricode',all.x=T)
  MATOUT <- MATOUT[,c(1,2,5,3,4)]
  
  write.table(MATOUT,file = paste(SETOUT,'/irrig.txt',sep=''),col.names=T,
              row.names=F,sep='\t',quote=F)
  
}

######## soil input file ###########
INPUTSOILMOS <- function(TREATCODE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE=''){

  MATTRAIT <- read.table(paste(SETOUT,'/TREATMENT_Mosicas',IDCODE,'.csv',sep=''),header=T,stringsAsFactors = F,dec='.',sep=';')
    soilcode <- unique(as.vector(MATTRAIT$soilcode)) 
  
  MATGS <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_genesoil.txt',header = T,sep ="\t")
  MATTEMP <- MATGS[which(MATGS$soilcode %in% soilcode),]
  MATOUT <- data.frame(soilcode=MATTEMP$soilcode,
                       soiltype=MATTEMP$typesoil,
                       nblayer=MATTEMP$nblayer,
                       p0=MATTEMP$p0,
                       ru=MATTEMP$ru,
                       soildepth=MATTEMP$depthsoil)
  write.table(MATOUT,file = paste(SETOUT,'/SOILGEN_Mosicas.csv',sep=''),col.names=T,
              row.names=F,sep=';',quote=F,dec='.')

  MATLS <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_layersoil.txt',header = T,sep ="\t")
  MATLS <- MATLS[which(MATLS$soilcode %in% soilcode),]
  
  head(MATLS)
  MATOUT <- data.frame(soilcode=NA,layernum=NA,thickness=NA,swcwp=NA,
                        swcfc=NA,swcsat=NA,bdens=NA)
  for (i in 1:length(soilcode)){
    MATTEMP <- MATLS[which(MATLS$soilcode==soilcode[i]),]
    MATTEMP2 <- MATTEMP[MATTEMP$layerlib=='epaiss',]
    MATTEMP2[order(MATTEMP2$layernum),]$layerval
    MATTEMP3 <- data.frame(soilcode=MATTEMP2$soilcode,
                           layernum=MATTEMP2$layernum,
                           thickness=MATTEMP2[order(MATTEMP2$layernum),]$layerval)
    
    MATTEMP2 <- MATTEMP[MATTEMP$layerlib=='hpf',]
    MATTEMP3$swcwp <- MATTEMP2[order(MATTEMP2$layernum),]$layerval
    MATTEMP2 <- MATTEMP[MATTEMP$layerlib=='hcc',]
    MATTEMP3$swcfc <- MATTEMP2[order(MATTEMP2$layernum),]$layerval
    MATTEMP2 <- MATTEMP[MATTEMP$layerlib=='hsat',]
    MATTEMP3$swcsat <- MATTEMP2[order(MATTEMP2$layernum),]$layerval
    MATTEMP2 <- MATTEMP[MATTEMP$layerlib=='bdens',]
    MATTEMP3$bdens <- MATTEMP2[order(MATTEMP2$layernum),]$layerval
    
    MATOUT <- rbind(MATOUT,MATTEMP3)
  }
  
  MATOUT <- MATOUT[-1,]
  write.table(MATOUT,file = paste(SETOUT,'/SOILLAYER_Mosicas.csv',sep=''),col.names=T,
              row.names=F,sep=';',quote=F,dec='.')
  
  MATTRAIT <- read.table(paste(SETOUT,'/TREATMENT_Mosicas',IDCODE,'.csv',sep=''),header=T,stringsAsFactors = F,sep=';',dec='.')
  soilcode <- unique(as.vector(MATTRAIT$soilcode)) 
  
  MATIS <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_inisoil.txt',header = T,sep ="\t")
  MATTEMP <- MATIS[which(MATIS$plotcyclecode %in% MATTRAIT$treatcode),]
  MATOUT <- data.frame(treatcode=MATTEMP$plotcyclecode,layernum=MATTEMP$layer,swc=MATTEMP$hvsoil)
  
  write.table(MATOUT,file = paste(SETOUT,'/SOILINI_Mosicas.csv',sep=''),col.names=T,
             row.names=F,sep=';',quote=F,dec='.')
  
  
}

####### plant and soil observations
INPUTOBSMOS <- function(TREATCODE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE=''){
  
  TREATCODE <- sort(as.vector(TREATCODE))
  MATOBS <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_obsplant.txt',header = T,sep ="\t")
  MATOBS <- MATOBS[which(MATOBS$varcyclecode %in% TREATCODE),]
  head(MATOBS)
  
  MATstagnb <- MATOBS[which(MATOBS$obsplantlib %in% 'nbtv'),]
  MATstadnb <- MATOBS[which(MATOBS$obsplantlib %in% 'nbtm'),]
  MATstafnb <- MATOBS[which(MATOBS$obsplantlib %in% 'nbtfl'),]
#  MATstafnb <- data.frame('nbtfl','CVA1-0-0-AAA-R579','2003-03-24',NA) ; names(MATstafnb) <- names(MATstagnb)
  MATstahtvd <- MATOBS[which(MATOBS$obsplantlib %in% 'htm'),]
  MATblatnb <- MATOBS[which(MATOBS$obsplantlib %in% 'nblts'),]
  MATblagnb <- MATOBS[which(MATOBS$obsplantlib %in% 'nblvs'),]
  MATblastvd <- MATOBS[which(MATOBS$obsplantlib %in% 'ltvdsurfs'),]
  MATglai <- MATOBS[which(MATOBS$obsplantlib %in% 'lai'),]
  MATei <- MATOBS[which(MATOBS$obsplantlib %in% 'ei'),]
  MATstamsu <- MATOBS[which(MATOBS$obsplantlib %in% 'mssucrtu'),]
  MATstamfm <- MATOBS[which(MATOBS$obsplantlib %in% 'mftu'),]
  MATstamdm <- MATOBS[which(MATOBS$obsplantlib %in% 'mstu'),]
  MATstamwc <- MATOBS[which(MATOBS$obsplantlib %in% 'humtig'),]
  MATstamsuc <- MATOBS[which(MATOBS$obsplantlib %in% 'brixjt'),]
  MATstamfbc <- MATOBS[which(MATOBS$obsplantlib %in% 'fibrt'),]
  MATstampur <- MATOBS[which(MATOBS$obsplantlib %in% 'brixjt'),] ######### problème pas bon
  MATblagdm <- MATOBS[which(MATOBS$obsplantlib %in% 'mslv'),]
  MATbladm <- MATOBS[which(MATOBS$obsplantlib %in% 'mslt'),]
  #MATbladm <- data.frame('mslt','CVA1-0-0-AAA-R579','2003-03-24',NA) ; names(MATbladm) <- names(MATstagnb)
  MATagrfm <- MATOBS[which(MATOBS$obsplantlib %in% 'mfa'),]
  MATagrdm <- MATOBS[which(MATOBS$obsplantlib %in% 'msa'),]
  MATagrwc <- MATOBS[which(MATOBS$obsplantlib %in% 'humaer'),]
  MATsendm <- MATOBS[which(MATOBS$obsplantlib %in% 'msasen'),]
  
  NAMESLIST <- c('stagnb','stadnb','stafnb','stahtvd','blatnb','blagnb','blastvd','glai','ei','stamsu','stamfm','stamdm',
                 'stamwc','stamsuc','stamfbc','stambjc','stampur','blagdm','bladm','agrfm','agrdm','agrwc','sendm')
  MATLIST <- list(MATstagnb,MATstadnb,MATstafnb,MATstahtvd,MATblatnb,MATblagnb,MATblastvd,MATglai,MATei,
                  MATstamsu,MATstamfm,MATstamdm,MATstamwc,MATstamsuc,MATstamfbc,MATstampur,MATstampur,MATblagdm,
                  MATbladm,MATagrfm,MATagrdm,MATagrwc,MATsendm)
  
  MATOBSPLANT <- merge(MATLIST[[1]][,c('obsplantdate','obsplantval','varcyclecode')],
                       MATLIST[[2]][,c('obsplantdate','obsplantval','varcyclecode')],by=c('varcyclecode','obsplantdate'),all=T)
  names(MATOBSPLANT)[3:4] <- c(NAMESLIST[1],NAMESLIST[2])
  for (i in 3:23){
    MATOBSPLANT <- merge(MATOBSPLANT,MATLIST[[i]][,c('obsplantdate','obsplantval','varcyclecode')],
                         by=c('varcyclecode','obsplantdate'),all=T)
    names(MATOBSPLANT)[i+2] <- c(NAMESLIST[i])
  }
#  head(MATOBSPLANT)
  
  #MATOBSPLANT <- replace(MATOBSPLANT,is.na(MATOBSPLANT), '-999')
  calage <- rep(1,dim(MATOBSPLANT)[1])
  MATOBSOUT <- data.frame(MATOBSPLANT[,1:2],calage,MATOBSPLANT[,3:25])
  OBSDATE <- MATOBSOUT$obsplantdate
  OBSDATE <- format(as.Date(OBSDATE,format='%d/%m/%Y'), "%d-%m-%Y")
  MATOBSOUT$obsplantdate <- OBSDATE
  
  names(MATOBSOUT)[1:2] <- c('treatcode','obsdate')
  
  write.table(MATOBSOUT,file = paste(SETOUT,'/OBS_plant',IDCODE,'.csv',sep=''),col.names=T,
              row.names=F,sep=';',quote=F,dec='.')
  
  # soil observation
  MATOBS <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_obsoil.txt',header = T,sep ="\t")
  MATOBS <- MATOBS[which(MATOBS$obsoilib=='hvsol'),]
  
  MATVC <- read.table(file = 'DATA/SAUV_ECOFI/ecofi_varcycle.txt',header = T,sep='\t')
  head(MATVC)
  MATTEMP <- MATVC[which(MATVC$varcyclecode %in% TREATCODE),1:2]
  
  MATOBS <- merge(MATTEMP,MATOBS,by='plotcyclecode')
  head(MATOBS)

  MATOUT <- data.frame(treatcode=NA,obsdate=as.Date('2000-01-01',format='%d/%m/%Y'),swc1=NA,swc2=NA,swc3=NA,swc4=NA,swc5=NA,swc6=NA)  
  
  for (i in 1:length(TREATCODE)){
  
  ID <- which(MATOBS$varcyclecode %in% TREATCODE[i])
  if(sum(ID) !=0){
    MATTEMP <- MATOBS[ID,-1]
    MATTEMP$obsoildate <- format(as.Date(MATTEMP$obsoildate,format='%d/%m/%Y'),format='%d/%m/%Y')
    
    MATTEMP <- MATTEMP[order(MATTEMP[,c('layer')]),]
    head(MATTEMP)
    
    MATLIST <- list()
    for (ilay in 1:max(MATTEMP$layer,na.rm=T)){
    MATTEMP2 <- MATTEMP[which(MATTEMP$layer==ilay),]
    MATTEMP2 <- MATTEMP2[order(MATTEMP2$obsoildate),]
    MATLIST[[ilay]] <- MATTEMP2}
    

    MATOUTTREAT <- MATLIST[[1]][,c('varcyclecode','obsoildate','obsoilval')] 
    if(max(MATTEMP$layer,na.rm=T)>1){
      for (ilay in 2:max(MATTEMP$layer,na.rm=T)){
        MATOUTTREAT <- merge(MATOUTTREAT,
                             MATLIST[[ilay]][,c('varcyclecode','obsoildate','obsoilval')],
                             by=c('varcyclecode','obsoildate'))
      }
    }
    IDLACK <- 8-dim(MATOUTTREAT)[2]
    if(IDLACK >0){MATOUTTREAT[,9-IDLACK] <-NA}
    MATOUTTREAT <- MATOUTTREAT[,1:8]
    names(MATOUTTREAT)[1:8] <- c('treatcode','obsdate','swc1','swc2','swc3','swc4','swc5','swc6')
    head(MATOUTTREAT)
    
    MATOUT <- rbind(MATOUT,MATOUTTREAT)
    
  }
  }
  MATOUT <- MATOUT[-1,]
  #MATOUT <- replace(MATOUT,is.na(MATOUT),-999)

  write.table(MATOUT,file = paste(SETOUT,'/OBS_soil',IDCODE,'.csv',sep=''),col.names=T,
              row.names=F,sep=';',quote=F,dec='.')
  
  
}

