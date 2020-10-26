############################################################
######### Mosicas run program ##############################
############################################################
remove(list = ls())
#library(profvis)
#T1<-Sys.time()


############################################################
#### Input and output folder have to be in the same place than "RUN_Mosica.R"
#require(rstudioapi)
setIni <- setwd(c(dirname(rstudioapi::getActiveDocumentContext()$path)))
setInput <- paste(c(dirname(rstudioapi::getActiveDocumentContext()$path)),"/input",sep='') 
setOutput <- paste(c(dirname(rstudioapi::getActiveDocumentContext()$path)),"/output",sep='') 


############### Reading input and opening output
source(paste(setIni,'R_MOSICAS/Mosicas_input_output.R',sep='/'))
############### growth function
source(paste(setIni,'R_MOSICAS/Mosicas_function.R',sep='/'))


#### output folder creation, if not existing
dir.create('output')

#### Read simulation input
Input_read <- INPUT_READ(FILE='SIM_Mosicas.csv',SEP=';',SetInput=setInput)

#### Read simulations run
Treatment_read <- TREATMENT_READ(FILE='TREATMENT_Mosicas.csv',SEP=';',SetInput=setInput,Treatcode=Input_read$RUN)

#### Open output (list of data.frame)
MATOUT <- OPEN_OUTPUT(SetOutput=setOutput,OUTNAME=Input_read$OUTNAME,INIDATE=Treatment_read$begdate,
                      ENDDATE=Treatment_read$enddate, TREATCODE=Input_read$RUN,
                      IFREQ=Input_read$IFREQ,IOBS=Input_read$IOBS)


##################################################
#### Beginning of treatment loop
##################################################

for (Itreat in 365:length(Treatment_read$treatcode)){
  
  print('Beginning treatment loop')
  
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
  #NDAYSIMUL = max(Weather_Irrig_data$julday)
  NDAYSIMUL = as.numeric(as.Date(Treatment_ID$enddate,format='%d/%m/%Y')-as.Date(Treatment_ID$begdate,format='%d/%m/%Y'))+1
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
  
  ############### Forcing interception, to do
  isobei=0
  
  ############### plant parameter
  PARAM <- Plant_param_input(VARCODE=Treatment_ID$varcode,CYCLE=Treatment_ID$cycle)
  
  ############### Initialisation variable
  source('R_MOSICAS/Mosicas_Initialisation_variable.R')
  
  
  ######################################
  ######## begin daily loop ############
  
  for (iday in 716:NDAYSIMUL){


    ########## daily growth model
    source(paste(setIni,'R_MOSICAS/Mosicas_daily_growth.R',sep='/'))
    
    ########## Writing daily output    
    source(paste(setIni,'R_MOSICAS/Mosicas_daily_output.R',sep='/'))
    
    counter(iday)
    

  }
  print(Itreat)
}

save.image(paste('output/RUN_Mosicas_',Input_read$OUTNAME,sep=''))


#T2<-Sys.time()

#Tdiff= difftime(T2, T1)


plot(c(MATOUT[[1]]$glai))
plot(c(MATOUT[[1]]$swdef))

plot(MATOUT[[1]]$agrdm,col=1,type='l',lwd=2)
points(MATOUT[[1]]$stemdm,col=2,type='l',lwd=2)

#plot(MATOUT[[1]]$eos)
#points(MATOUT[[1]]$es,col=2)
#plot(MATOUT[[1]]$agrdm,col=1,type='l',lwd=2)
#points(MATOUT[[1]]$stamdm,col=2,type='l',lwd=2)

#plot(MATOUT[[1]]$swdef)
#plot(MATOUT[[1]]$swdf1)
#plot(MATOUT[[1]]$swdf2)


plot(MATOUT[[1]]$rootdm/MATOUT[[1]]$agrdm,xlim=c(1,300),type='l')

