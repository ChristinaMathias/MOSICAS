INPUT_READ <- function(FILE='SIM_Mosicas.csv',SEP=';',SetInput=setInput,setOutput,
                           OUTNAME,WATMODEL,IFREQ,IOBS,IWATBAL,INBAL,NBRUN,RUN,
                            IRRIGMODE,IRRIGFREQ,IRRIGMM){
  # read run input
  MATTEMP <- read.table(file = paste(SetInput,FILE,sep='/'),header = T,sep = SEP,stringsAsFactors=F,dec='.')
  OUTNAME <- MATTEMP$Input[which(MATTEMP$Input.name=='OUTNAME')]
  WATMODEL <- MATTEMP$Input[which(MATTEMP$Input.name=='WATMODEL')]
  IFREQ <- as.numeric(MATTEMP$Input[which(MATTEMP$Input.name=='IFREQ')])
  IOBS <- as.numeric(MATTEMP$Input[which(MATTEMP$Input.name=='IOBS')])
  IWATBAL <- as.numeric(MATTEMP$Input[which(MATTEMP$Input.name=='IWATBAL')])
  INBAL <- as.numeric(MATTEMP$Input[which(MATTEMP$Input.name=='INBAL')])
  IRRIGMODE <- as.numeric(MATTEMP$Input[which(MATTEMP$Input.name=='IRRIGMODE')])
  IRRIGFREQ <- as.numeric(MATTEMP$Input[which(MATTEMP$Input.name=='IRRIGFREQ')])
  IRRIGMM <- as.numeric(MATTEMP$Input[which(MATTEMP$Input.name=='IRRIGMM')])
  NBRUN <- as.numeric(MATTEMP$Input[which(MATTEMP$Input.name=='NBRUN')])
  RUN <- MATTEMP$Input[which(MATTEMP$Input.name=='RUN'):(which(MATTEMP$Input.name=='RUN')+NBRUN-1)]
  
  return(list(OUTNAME=OUTNAME,WATMODEL=WATMODEL,IFREQ=IFREQ,IOBS=IOBS,IWATBAL=IWATBAL,INBAL=INBAL,NBRUN=NBRUN,RUN=RUN,
              IRRIGMODE=IRRIGMODE,IRRIGFREQ=IRRIGFREQ,IRRIGMM=IRRIGMM))
}


TREATMENT_READ <- function(FILE='TREATMENT_Mosicas.csv',SEP=';',SetInput=setInput,
                           Treatcode=Input_read$RUN){
  # read run input
  MATTEMP <- read.table(file = paste(SetInput,FILE,sep='/'),header = T,sep = SEP,stringsAsFactors=F,dec='.')
  MATTEMP <- MATTEMP[which(MATTEMP$treatcode %in% Treatcode),]
  # create list output
  LISTTEMP <- list() ;   length(LISTTEMP) <-length(names(MATTEMP))
  LISTTEMP <- lapply(1:length(LISTTEMP), function(i) MATTEMP[,i])
  names(LISTTEMP) <-names(MATTEMP)
  return(LISTTEMP)
}


TREATMENT_ID <-function(TREATMENT_read=Treatment_read,ITREAT=Itreat){
  LISTTEMP <- lapply(1:length(TREATMENT_read),function(i) TREATMENT_read[[i]][ITREAT])
  names(LISTTEMP) <- names(TREATMENT_read)
  return(LISTTEMP)
}
  
WEATHERSTATION_INPUT <- function(SetInput=setInput,FILE='WEATHER_STATION_Mosicas.csv',SEP=";",
                                 WSATCODE =Treatment_ID$wstacode){
  MATTEMP <- read.table(file = paste(SetInput,FILE,sep='/'),header = T,sep = SEP,stringsAsFactors=F,dec='.')
  MATTEMP <- MATTEMP[which(MATTEMP$wstacode == WSATCODE),]
  # create list output
  LISTTEMP <- list() ;   length(LISTTEMP) <-length(names(MATTEMP))
  LISTTEMP <- lapply(1:length(LISTTEMP), function(i) MATTEMP[,i])
  names(LISTTEMP) <-names(MATTEMP)
  print('Read WEATHER_STATION_Mosicas')
  return(LISTTEMP)
}  

WEATHER_DATA <- function(WSATCODE=Treatment_ID$wstacode,begdate=Treatment_ID$begdate,enddate=Treatment_ID$enddate,
                         SetInput=setInput){
  MATTEMP <- read.table(file=paste(SetInput,'/Weather and irrigation data/',WSATCODE,'.txt',sep=''),header=T)
  DATINI <- as.Date(begdate,format='%d/%m/%Y')
  DATEND <- as.Date(enddate,format='%d/%m/%Y')
  MATTEMP$date = as.Date(MATTEMP$date,format='%d-%m-%Y')
  if(min(MATTEMP$date) > DATINI){print('Warning beginning date < than weather data')}
  if(max(MATTEMP$date) < DATEND){print('Warning endding date > than weather data')}
  IDINI <- which(MATTEMP$date==DATINI)
  IDEND <- which(MATTEMP$date==DATEND)
  MATTEMP <- MATTEMP[IDINI:IDEND,]
  MATTEMP <- replace(MATTEMP,MATTEMP==-999,NA)
  head(MATTEMP)
  IDNA <- sapply(2:5, function(i) length(which(is.na(MATTEMP[,i]))))
  if(sum(IDNA)>0){print('Warning missing daily weather data in rr, etp, tx, tn or rg')}
    MATTEMP$juldatmet <- as.numeric(MATTEMP$date-DATINI)+1
  
  return(MATTEMP)

}

IRRIG_DATA <- function(IRRIMODE=Input_read$IRRIGMODE,IRRICODE=Treatment_ID$irricode,
                       IRRIFREQ=Input_read$IRRIGFREQ,IRRIMM=Input_read$IRRIGMM,
                       begdate=Treatment_ID$begdate,enddate=Treatment_ID$enddate,SetInput=setInput,
                       Weatherdata=Weather_data){
  
  # no irrigation
  if((IRRIMODE==0)|(IRRICODE=='RAINFED')){MATTEMP <- NULL;print('No irrigation')}
  
  # irrigation read in file
  if((IRRIMODE==1)&(IRRICODE!='RAINFED')){
    MATTEMP <- read.table(file=paste(SetInput,'/Weather and irrigation data/irrig.txt',sep=''),header=T)
    MATTEMP <- MATTEMP[which(MATTEMP$irricode==IRRICODE),]
  
    if (dim(MATTEMP)[1]==0){print('Warning irrigation data not found, assumed automatic') ; IRRIMODE =2}
    if (dim(MATTEMP)[1]>=1){
      DATINI <- as.Date(begdate,format='%d/%m/%Y')
      juldayirri <- as.numeric(as.Date(MATTEMP$dateirri,format='%d-%m-%Y')-DATINI)+1
      MATTEMP$juldayirri = juldayirri
      MATTEMP$dateirri = as.Date(MATTEMP$dateirri,format='%d-%m-%Y')
    }
  }
  # automatic irrigation
  if((IRRIMODE==2)&(IRRICODE!='RAINFED')){
    print('Automatic irrigation')
    DATINI <- as.Date(begdate,format='%d/%m/%Y')
    DATEND <- as.Date(enddate,format='%d/%m/%Y')
    # duration of the crop cycle
    LENDAY <- as.numeric(DATEND-DATINI)
    # number of irrigation applicaton
    NIRR <- round(LENDAY/IRRIFREQ,0)
    # date of application in julian day
    juldayirri <- (0:(NIRR-1))*IRRIFREQ+1
    # day of application
    dateirri <- DATINI + juldayirri-1
    # out
    MATTEMP <- data.frame(irricode=rep('Auto',length(dateirri)),
                          irritype=rep(1,length(dateirri)),
                          efficiency=rep(1,length(dateirri)),
                          irridose=rep(IRRIMM,length(dateirri)),
                          dateirri=dateirri,juldayirri=juldayirri)
  }

  if((IRRIMODE==0)|(IRRICODE=='RAINFED')){
    Weatherdata$irrig=0 ; 
    MATTEMP2=Weatherdata
    names(MATTEMP2)[which(names(MATTEMP2)=="juldatmet")] <- 'julday'}
  if((IRRIMODE>0)&(IRRICODE!='RAINFED')){
    names(Weatherdata)[which(names(Weatherdata)=="juldatmet")] <- 'julday'
    MATTEMP$irrig = MATTEMP$efficiency*MATTEMP$irridose
    names(MATTEMP)[which(names(MATTEMP)=="juldayirri")] <- 'julday'
    MATTEMP2 <- merge(Weatherdata,MATTEMP[,c('julday','irrig')],by='julday',all.x = T)
    MATTEMP2$irrig[which(is.na(MATTEMP2$irrig))] <- 0
  }
    
  return(MATTEMP2)
}

SOILGEN_INPUT <- function(SetInput=setInput,FILE='SOILGEN_Mosicas.csv',SEP=";",
                                 SOILCODE =Treatment_ID$soilcode){
  MATTEMP <- read.table(file = paste(SetInput,FILE,sep='/'),header = T,sep = SEP,stringsAsFactors=F,dec='.')
  MATTEMP <- MATTEMP[which(MATTEMP$soilcode == SOILCODE),]
  # create list output
  LISTTEMP <- list() ;   length(LISTTEMP) <-length(names(MATTEMP))
  LISTTEMP <- lapply(1:length(LISTTEMP), function(i) MATTEMP[,i])
  names(LISTTEMP) <-names(MATTEMP)
  print('Read SOILGEN_Mosicas')
  return(LISTTEMP)
}  

SOILLAYER_INPUT <- function(SetInput=setInput,FILE='SOILLAYER_Mosicas.csv',SEP=";",
                          SOILCODE =Treatment_ID$soilcode,Soilgeninput=Soilgen_input){
  if(Soilgen_input$nblayer>1){
    MATTEMP <- read.table(file = paste(SetInput,FILE,sep='/'),header = T,sep = SEP,stringsAsFactors=F,dec='.')
    MATTEMP <- MATTEMP[which(MATTEMP$soilcode == SOILCODE),]
    # create list output
    LISTTEMP <- list() ;   length(LISTTEMP) <-length(names(MATTEMP))
    LISTTEMP <- lapply(1:length(LISTTEMP), function(i) MATTEMP[,i])
    names(LISTTEMP) <-names(MATTEMP)
  }
  if(Soilgeninput$nblayer==1){
    Soilgeninput$nblayer=6
    LISTTEMP <- list()
    LISTTEMP$soilcode = rep('Hypothetical',6)
    LISTTEMP$layernum = 1:6
    LISTTEMP$thickness = c(5,rep(round((Soilgeninput$soildepth-5)/5,0),5))
    LISTTEMP$swcwp = rep(0.3,6)
    LISTTEMP$swcfc = rep(0.3+Soilgeninput$ru/1000,6)
    LISTTEMP$swcsat = LISTTEMP$swcfc+0.05
    LISTTEMP$bdens=rep(1,6)
  }
  print('Read SOILLAYER_Mosicas')
  return(LISTTEMP)
}  

SOILINI_INPUT <- function(SetInput=setInput,FILE='SOILINI_Mosicas.csv',SEP=";",
                            TREATCODE =Treatment_ID$treatcode,Soillayerinput=Soillayer_input,
                          TreatmentID=Treatment_ID){
  
  MATTEMP <- read.table(file = paste(SetInput,FILE,sep='/'),header = T,sep = SEP,stringsAsFactors=F,dec='.')
  MATTEMP <- MATTEMP[which(MATTEMP$treatcode == TREATCODE),]
  
  #in case there more observed layer than define layers
  if(dim(MATTEMP)[1]> 0){MATTEMP <- MATTEMP[1:max(Soillayerinput$layernum),]}
  
  ### if no value we calculate using the inisoilvol
  if(dim(MATTEMP)[1]==0){
    MATTEMP <- data.frame(as.vector(rep(TREATCODE,max(Soillayerinput$layernum))),
                          1:max(Soillayerinput$layernum),
                          (Soillayerinput$swcwp +(Soillayerinput$swcfc-Soillayerinput$swcwp)*TreatmentID$inisoilvol))
    names(MATTEMP) <- c('treatcode','layernum','swc')
  }
  # create list output
  LISTTEMP <- list() ;   length(LISTTEMP) <-length(names(MATTEMP))
  LISTTEMP <- lapply(1:length(LISTTEMP), function(i) MATTEMP[,i])
  names(LISTTEMP) <-names(MATTEMP)
  print('Read SOILINI_Mosicas')
  return(LISTTEMP)
}  

CALCSOILWATERCAPACITIES <- function(NLAYER=Soilgen_input$nblayer,SOILLAYER_input=Soillayer_input,
                                    SOILINI_input=Soilini_input){

    if (Soilgen_input$nblayer==1){NLAYER=6} # default value in Soillayer_input
      
    swcsat_mm <- sum(SOILLAYER_input$swcsat*SOILLAYER_input$thickness*10)
    swcfc_mm <- sum(SOILLAYER_input$swcfc*SOILLAYER_input$thickness*10)
    swcwp_mm <- sum(SOILLAYER_input$swcwp*SOILLAYER_input$thickness*10)
    swcini_mm <- sum(SOILINI_input$swc*SOILLAYER_input$thickness*10)

#    tsw=0.;tswcwp=0.; esw <-c()
    
    # attention ?? vérifier si sert
#    for (i in 1:NLAYER){
#      esw <- c(esw,SOILLAYER_input$swcfc[i]-SOILLAYER_input$swcwp[i])
#      if(esw[i]<=0){print('Error in soil layer water capacity')}
#    tswcwp <- tswcwp + SOILLAYER_input$swcwp[i]*SOILLAYER_input$thickness[i]
#    tsw <- tsw + (SOILINI_input$swc[i]-SOILLAYER_input$swcwp[i])/esw[i]
#    }
    return(list(swcsat_mm=swcsat_mm,swcfc_mm=swcfc_mm,swcwp_mm=swcwp_mm,
                #tsw=tsw,tswcwp=tswcwp,swcini_mm=swcini_mm,esw=esw))
                swcini_mm=swcini_mm))

    # semblerait que tswwp ne son pas important à supprimer !
  
}


OPEN_OUTPUT <- function(SetOutput=setOutput,OUTNAME=Input_read$OUTNAME,INIDATE,ENDDATE,
                        TREATCODE=Input_read$RUN, IFREQ,IOBS,ISENS=0,rsens=0){
  
  MATOUT <- list()
  # nombre of day for each simulation
  duration <- sapply(1:length(INIDATE), function(i) 
    as.numeric(as.Date(ENDDATE[i],format='%d/%m/%Y')-as.Date(INIDATE[i],format='%d/%m/%Y'))+1)
  # names of output
  NAMESOUT <- c('nosimul','treatcode','lati','longi','alti','dates','age','agestem','glai','ei',
                'stemnb','stemfm','stemdm','stemwc','stemsu','stemsuc',
                'rue','agrfm','agrdm','agrwc','bladm','rootdm',
                'hv1','hv2','hv3','hv4','hv5','hv6','stock','stockfr','runoff','deepdrain',
                'rootfront','droot1','droot2','droot3','droot4','droot5','droot6','precip',
                'tmin','tmax','tmean','solrad','etp','etm','etr','tmp','trp','eos','es','swdef','swdfrue','swdflai','swan',
                'pari','stmean','spar','spari','spr','sirr','setp','setm','setr','stmp','strp',
                'swdefm','swdfruem','swdflaim','sdrain','srunoff','sevap','durjour','srgx') 
  if(ISENS !=1){  
  # output only at harvest
  if(IFREQ==0){
    for (i in 1:length(INIDATE)){
      MATOUT[[i]] <- data.frame(matrix(rep(NA,length(NAMESOUT)*1),nrow=1))
      names(MATOUT[[i]]) <- NAMESOUT
      }
  }
  # output every day
  if(IFREQ==1){
    for (i in 1:length(INIDATE)){
      MATOUT[[i]] <- data.frame(matrix(rep(NA,length(NAMESOUT)*duration[i]),nrow=duration[i]))
      names(MATOUT[[i]]) <- NAMESOUT
    }
  }
    names(MATOUT) <- TREATCODE
  }
  
  if(ISENS==1){
    for (i in 1:(rsens)){
      MATOUT[[i]] <- data.frame(matrix(rep(NA,length(NAMESOUT)*1),nrow=1))
      names(MATOUT[[i]]) <- NAMESOUT
    }
    names(MATOUT) <- sapply(1:rsens,function(i) paste('Sim',i,sep='_'))
  }
  

  return(MATOUT)
}
#

Plant_param_input <- function(SetInput=setInput,SEP=';',VARCODE=Treatment_ID$varcode,CYCLE=Treatment_ID$cycle){
  MATTEMP <- read.table(file = paste(SetInput,'Plant_parameter.csv',sep='/'),header = T,sep = SEP,stringsAsFactors=F,dec='.')
  MATTEMP <- MATTEMP[which(MATTEMP$varcode==VARCODE),]
  LISTTEMP <- list()
  for (i in 1:dim(MATTEMP)[1]){
    if(CYCLE==0){LISTTEMP[[i]] <- MATTEMP$plantcrop[[i]]}
    if(CYCLE>0){LISTTEMP[[i]] <- MATTEMP$ratoon[[i]]}
  }
  names(LISTTEMP) <- MATTEMP$paramname
  return(LISTTEMP)
}

OBS_READ <- function(FILE='OBS_plant.csv',SEP=';',SetInput=setInput,
                           Treatcode=Input_read$RUN){
  # read run input
  MATTEMP <- read.table(file = paste(SetInput,FILE,sep='/'),header = T,sep = SEP,stringsAsFactors=F,dec='.')
  MATTEMP <- MATTEMP[which(MATTEMP$treatcode %in% Treatcode),]
  names(MATTEMP)[4:dim(MATTEMP)[2]] <-  sapply(4:(dim(MATTEMP)[2]),function(i) paste('o',names(MATTEMP)[[i]],sep=''))

  return(MATTEMP)
}
