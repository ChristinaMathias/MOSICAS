#########################
#### sensitivity and calibration function for Mosicas



#SensMosicas.fun <- function(Newpar) {
  
  
  ############### Reading input and opening output
  source(paste(setIni,'R_MOSICAS/Mosicas_input_output.R',sep='/'))
  ############### growth function
  source(paste(setIni,'R_MOSICAS/Mosicas_function.R',sep='/'))
  
  #### Read simulation input
  Input_read <- INPUT_READ(FILE='SIM_Mosicas.csv',SEP=';',SetInput=setInput)
  # whatever happen, only the first simulation units (RUN) is used
  Input_read$NBRUN=1 ; Input_read$RUN = Input_read$RUN[1]
  Input_read$IFREQ = 0 # only values at harvest used
  
  #### Read simulations run
  Treatment_read <- TREATMENT_READ(FILE='TREATMENT_Mosicas.csv',SEP=';',SetInput=setInput,Treatcode=Input_read$RUN)
  
  #### Open output (list of data.frame)
  MATOUT <- OPEN_OUTPUT(SetOutput=setOutput,OUTNAME=Input_read$OUTNAME,INIDATE=Treatment_read$begdate,
                        ENDDATE=Treatment_read$enddate, TREATCODE=Input_read$RUN,
                        IFREQ=0,IOBS=Input_read$IOBS,ISENS = 1,rsens = dim(Newpar)[1])
  
  ##################################################
  #### Unchanged input
  ##################################################
  Itreat=1 ; isobei=0
  
  ################ Choice of treatment
  Treatment_ID <- TREATMENT_ID(TREATMENT_read=Treatment_read,ITREAT=Itreat)
  
  ############### read weather station information
  Weatherstation_input <- WEATHERSTATION_INPUT (WSATCODE =Treatment_ID$wstacode)
  
  ############### read weather data
  Weather_data <- WEATHER_DATA(WSATCODE=Treatment_ID$wstacode,begdate=Treatment_ID$begdate,enddate=Treatment_ID$enddate)
  
  ############### irrigation data
  Weather_Irrig_data <- IRRIG_DATA(IRRIMODE=Input_read$IRRIGMODE,IRRICODE=Treatment_ID$irricode,
                                   IRRIFREQ=Input_read$IRRIGFREQ,IRRIMM=Input_read$IRRIGMM,
                                   begdate=Treatment_ID$begdate,enddate=Treatment_ID$enddate,SetInput=setInput)
  NDAYSIMUL = max(Weather_Irrig_data$julday)
  DATESIMUL <- Weather_Irrig_data$date
  
  ############### Read soil generic input
  Soilgen_input <- SOILGEN_INPUT (SOILCODE =Treatment_ID$soilcode)
  
  ############### Read soil layer input
  Soillayer_input <- SOILLAYER_INPUT (SOILCODE =Treatment_ID$soilcode,Soilgeninput=Soilgen_input)
  
  ############### Read soil initialization input
  Soilini_input <- SOILINI_INPUT (TREATCODE =Treatment_ID$treatcode,Soillayerinput=Soillayer_input)
  
  ############### Total soil water capacities
  Calcsoilwatercapacities <- CALCSOILWATERCAPACITIES(NLAYER=Soilgen_input$nblayer,
                                                     SOILLAYER_input=Soillayer_input,SOILINI_input=Soilini_input)
  
  ############### plant parameter
  PARAM <- Plant_param_input(VARCODE=Treatment_ID$varcode,CYCLE=Treatment_ID$cycle)
  
  
  ##################################################
  #### Beginning of sensitivity loop
  ##################################################
  
  for (Isim in 1:dim(Newpar)[1]) {
    
    for(i in 1:dim(Newpar)[2]){
      ID <- which(names(PARAM)==names(as.data.frame(Newpar))[i])
      PARAM[[ID]] <- Newpar[Isim,i]
    }
    
    ############### Initialisation variable
    source(paste(setIni,'R_MOSICAS/Mosicas_Initialisation_variable.R',sep='/'))
    
    
    ######################################
    ######## begin daily loop ############
    
    for (iday in 1:NDAYSIMUL){
      
      ########## daily growth model
      source(paste(setIni,'R_MOSICAS/Mosicas_daily_growth.R',sep='/'))

      ########## Writing daily output    
      ISENS = 1 # activation of output
      source(paste(setIni,'R_MOSICAS/Mosicas_daily_output.R',sep='/'))

      
      
    }
    print(paste('Sim',Isim,'on',dim(Newpar)[1],sep=' '))
  }
  
  MATOUTTEMP <- do.call("rbind", MATOUT)
  #return(MATOUTTEMP)
#}


