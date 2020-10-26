####### Run parameter optimization

### identification of available units of simulations (USM) in trait.txt
MATTRAIT <- read.table('input/trait.txt',header=T)
print(MATTRAIT$treatcode)

# selection of USMs that will be used for parameter optimization
# for all USMs, observed variables are given in obsplant.txt 
# if you want only ratoon or plot crop, check the crop cycle
Treats <- as.vector(MATTRAIT$treatcode[which(MATTRAIT$cycle>0)])
# specific USMs
#Treats <- c('BCSU-0-0-IRR-R570','BCSU-0-0-PLU-R570')

###############################################################
############ parameter identification and values ##############
## based on an example
## observed variable used for calibration
var='glai'
## parameter to calibrate
param=c('laicroi','laiwksen','sthydcroi')
## initial value of parameters
central=c(0.00415,0.0187,1)
## range of iteration +/- XX%
iter=c(95,95,95)
## step of the recuit analysis
pas=c(1,1,1)

# Writing calibration inputs
if(length(param)==1){CALIB <- list(NA,NA,NA,NA)}
if(length(param)==2){CALIB <- list(NA,NA,NA,NA,NA)}
if(length(param)==3){CALIB <- list(NA,NA,NA,NA,NA,NA)}
if(length(param)==4){CALIB <- list(NA,NA,NA,NA,NA,NA,NA)}
if(length(param)==5){CALIB <- list(NA,NA,NA,NA,NA,NA,NA,NA)}
CALIB[[1]] <- 1
CALIB[[2]] <- c(var,100)
CALIB[[3]] <- length(param)
for (i in 1:length(param)){CALIB[[i+3]] <- c(param[i],central[i],iter[i],pas[i])}

writeLines(unlist(lapply(CALIB, paste, collapse=" ")), 'parcal.txt')
writeLines(unlist(lapply(CALIB, paste, collapse=" ")), 'parcalini.txt')


###########################################################
##### Modification of sim.txt according to simulation option ######
##### carefull to adjust the number of calibrated parameters

# open sim.txt
LINES <- readLines('sim.txt')
Ntreatment=length(Treats)

if(Ntreatment==1){LINES <- LINES[1:19]}
if(Ntreatment>1){LINES <- c(LINES[1:19],rep(NA,Ntreatment-1))}

# change input and output folders
LINES[2] <- paste(getwd(),'/input',sep='')
LINES[4] <- paste(getwd(),'/output',sep='')

# localisation of USMs file: "obs" = read trait.txt, "str" = read trstr.txt
LINES[6] <- "obs" #"str"

# On the same line :
#  - NAMERUN = simulation name
NAMERUN <- 'Test'
#  - Choice of water balance model: MOSICASCERESb (recommended), MOSICASCERES, MOSICASCERESFAO56 (when p0 is used)
#  - output frequency: 0 = end of the simulation, 1 every days 
#  - output at the same date that measurements: 0 = no, 1 = yes
#  - water stress activation: 0 = no, 1 = yes
#  - nigroten stress activation: 0 = no, 1 = yes
LINES[8] <- paste(NAMERUN,' 1 ','MOSICASCERESb',' 0 0 1 1',sep='')

# type of simulation:
#  - first number: 0 simple run, 1 sensitivity analysis, 2 = parameter optimization
#  - second number: number of parameter used for sensitivity or optimization
LINES[10] <- "2 3"

# type irrigation: 0 = no irrigation, 1 = read irrigation (if any) in irrig.txt
LINES[16] <- 1

# List of simulations
LINES[17] <- '==TRAITEMENTS'
LINES[18] <- Ntreatment
LINES[19:(19+Ntreatment-1)] <- Treats
LINES[1:20]

# change sim.txt in the main folder
writeLines(LINES, 'sim.txt')

##################################################
##### run MOSICAS model #####
shell('mosicas.exe')

##################################################
##### read optimization results ##################

# home made function
RETNUM <- function(x){  as.numeric(strsplit(x," ")[[1]])[!is.na(as.numeric(strsplit(x," ")[[1]]))]}

# read results in recuit.txt and create a table
CALIB_VAL <- readLines('recuit.txt')
OptimVal <- c();
for (i in 1:length(param)){
  OptimVal <- c(OptimVal,RETNUM(CALIB_VAL[(i)+(length(param)+2)*4-1]))
}
MAT <- data.frame(param,OptimVal)

print(MAT)

######

