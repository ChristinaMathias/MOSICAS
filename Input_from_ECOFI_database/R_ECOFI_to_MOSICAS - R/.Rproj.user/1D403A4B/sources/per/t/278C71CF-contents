## Create MOSICAS input file from all trials with a comon variety in ECOFI

###############
## List of all trials with choosen cultivar, ex R570
#VAR <- c("R570")
MATVARCYCLE <- read.table(file='DATA/SAUV_ECOFY/ecofi_varcycle.txt',header=T)
#MATVARCYCLE <- MATVARCYCLE[which(MATVARCYCLE$varcode %in% VAR),]
head(MATVARCYCLE)
LISTVARCYCLE <- MATVARCYCLE$varcyclecode

###########################################
#### Création input trait
source('SCRIPT/R_function_mosicas_input.R')

# trait input file
INPUTTRAITMOS(TREATCODE=LISTVARCYCLE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE='')
# weather station and data input files
INPUTWSMOS(TREATCODE=LISTVARCYCLE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE='')
# irrigation input files
INPUTIRRIGMOS(TREATCODE=LISTVARCYCLE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE='')
# Soil inputs files
INPUTSOILMOS(TREATCODE=LISTVARCYCLE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE='')
# plant and soil observation inputs files
INPUTOBSMOS(TREATCODE=LISTVARCYCLE,SETOUT='DATA/MOSICAS_INPUT/',IDCODE='')




