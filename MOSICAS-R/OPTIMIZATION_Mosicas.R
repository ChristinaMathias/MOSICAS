###################################################
########## mosicas parameter optimization #########
###################################################
rm(list=ls())

############################################################
#### Input and output folder have to be in the same place than "SENSITIVITY_Mosicas.R"
setini <- setwd(c(dirname(rstudioapi::getActiveDocumentContext()$path)))
setinput <- paste(c(dirname(rstudioapi::getActiveDocumentContext()$path)),"/input",sep='') 
setoutput <- paste(c(dirname(rstudioapi::getActiveDocumentContext()$path)),"/output",sep='') 

# list of parameter to be used in the sensitivity analysis
MATPAR <- read.csv2(paste(setinput,'/Plant_parameter.csv',sep=''),header=T,stringsAsFactors = F,dec='.')
# all parameters
ParnameALL <- MATPAR$paramname[which(MATPAR$varcode=='R570')]
# specific parameters
Parname <- c("ke")
Nbpar <- length(Parname)

# initial value of parameter and bounds values
Valpar <- MATPAR$plantcrop[which(MATPAR$varcode=='R570')] # plant crop or ratoon crop
Valparini <- Valpar[ParnameALL %in% Parname]

Bsup <- Valparini*1.5
Binf <- Valparini*0.5
bounds <- matrix(c(Binf,Bsup),nrow=Nbpar,byrow=F)

# Output variable used
Outputname='agrdm'


# Optimization using Rgenoud

library(rgenoud)

OPTIM <- genoud(MosicasOptimFunction,  # objective function, 1st argument has to be a vector of parameter value
              namepar=Parname, # second argument of the function, with the parameter names
              outputname=Outputname, # third argument of the function with the output variable used
              setIni = setini,
              nvars = Nbpar,
              Domains = bounds,
#              default.domains= X,   # 1  inivalue +/- default.domain
              starting.values= Valparini,
              boundary.enforcement = 2, # no possibilities to go outside the boundaries
              print.level = 2,
              max = FALSE,
              pop.size = 3,#50
              max.generations = 10, # 200
              wait.generations = 1, #◙ 3
              hard.generation.limit = TRUE,
              MemoryMatrix = TRUE,
              BFGS = FALSE#,
              #unif.seed = 737,
              #int.seed = 747
)


MosicasOptimFunction <- function(optimParameterValue,namepar,outputname,setIni=setini){
  
  
  ############### Reading input and opening output
  source(paste(setIni,'R_MOSICAS/Mosicas_input_output.R',sep='/'))
  ############### growth function
  source(paste(setIni,'R_MOSICAS/Mosicas_function.R',sep='/'))
  setInput <- paste(setIni,"/input",sep='') 
  setOutput <- paste(setIni,"/output",sep='') 
  
  #### Read simulation input
  Input_read <- INPUT_READ(FILE='SIM_Mosicas.csv',SEP=';',SetInput=setInput)
  Input_read$IFREQ = 1 
  
  #### Read simulations run
  Treatment_read <- TREATMENT_READ(FILE='TREATMENT_Mosicas.csv',SEP=';',SetInput=setInput,Treatcode=Input_read$RUN)
  
  #### Open output (list of data.frame)
  MATOUT <- OPEN_OUTPUT(SetOutput=setOutput,OUTNAME=Input_read$OUTNAME,INIDATE=Treatment_read$begdate,
                        ENDDATE=Treatment_read$enddate, TREATCODE=Input_read$RUN,
                        IFREQ=Input_read$IFREQ,IOBS=Input_read$IOBS)
  
  
  ##################################################
  #### Beginning of treatment loop
  ##################################################
  
  for (Itreat in 1:length(Treatment_read$treatcode)){
    
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
    
    ############### Forcing
    isobei=0
    
    ############### plant parameter
    PARAM <- Plant_param_input(VARCODE=Treatment_ID$varcode,CYCLE=Treatment_ID$cycle)
    PARAM[namepar] <- optimParameterValue
    

    ############### Initialisation variable
    source(paste(setIni,'R_MOSICAS/Mosicas_Initialisation_variable.R',sep='/'))
    
    
    ######################################
    ######## begin daily loop ############
    
    for (iday in 1:NDAYSIMUL){
      
      ########## daily growth model
      source(paste(setIni,'R_MOSICAS/Mosicas_daily_growth.R',sep='/'))
      
      ########## Writing daily output    
      source(paste(setIni,'R_MOSICAS/Mosicas_daily_output.R',sep='/'))
      

      
      
    }
  }
  
  # data.frame with all simulations
  MATOUTTEMP <- do.call("rbind", MATOUT)
  
  # data.frame with observed values
  MATOBS <- OBS_READ(FILE='OBS_plant.csv',SEP=';',SetInput=setInput,
                       Treatcode=unique(MATOUTTEMP$treatcode))
  
  # merging according to observation date and trial
  MATOUTTEMP$mergedates = as.Date(MATOUTTEMP$dates, format='%Y-%m-%d')
  MATOBS$mergedates = as.Date(MATOBS$obsdate, format='%d/%m/%Y')
  MATMERGE <- merge(MATOUTTEMP,MATOBS,by=c('treatcode','mergedates'),all.x = F,all.y=T)
  
  # MSE values
  SIMVAL <- MATMERGE[,outputname]
  OBSVAL <- MATMERGE[,paste('o',outputname,sep='')]
  MSE <- sum((SIMVAL-OBSVAL)^2,na.rm=T)/length(OBSVAL[!is.na(OBSVAL)])
  
  return(MSE)
}

