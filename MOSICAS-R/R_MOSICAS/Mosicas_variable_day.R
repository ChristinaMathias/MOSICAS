########################################################
######## Daily variable calculation ####################
########################################################

# Minimum temperature
tempmn <- Weather_Irrig_data$tn[iday]
# Maximum temperature
tempmx <- Weather_Irrig_data$tx[iday]
# Average temperature
tmo=(tempmn+tempmx)*0.5
# sum of degree day (base 0)
stmo = stmo + tmo

# Global radiation
solrad <- Weather_Irrig_data$rg[iday]
# Photosynthetic active radiation
par=0.5*solrad;  
# sum of incident active radiation
spar=spar+par

# Julian day for astronomic calculation
jcal = as.numeric(format(DATESIMUL[iday],"%j"))















