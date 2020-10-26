####### Run the mosicas.exe model

### identification of available units of simulations (USM) in trait.txt
MATTRAIT <- read.table('input/trait.txt',header=T)
print(MATTRAIT$treatcode)

# selection of USMs that will be simulated
# all of them
Treats <- as.vector(MATTRAIT$treatcode)
# specific ones
#Treats <- c('BCSU-0-0-IRR-R570','BCSU-0-0-PLU-R570')

# Modification of sim.txt according to simulation option
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
LINES[10] <- "0 1"
  
# type irrigation: 0 = no irrigation, 1 = read irrigation (if any) in irrig.txt
LINES[16] <- 1

# List of simulations
LINES[17] <- '==TRAITEMENTS'
LINES[18] <- Ntreatment
LINES[19:(19+Ntreatment-1)] <- Treats
LINES[1:20]
  
# change sim.txt in the main folder
writeLines(LINES, 'sim.txt')
  
# run MOSICAS model
shell('mosicas.exe')

  
  