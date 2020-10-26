########################################################
######## Daily output ####################
########################################################
#### 1) sum of variables / indicators ###
#### 2) re-initialising value for the next day
#### 3) writing daily output or only at harvest


######## Sum of variables ########
spari <- INTERCEP$pari+ spari
spr <- Weather_Irrig_data$rr[iday] + spr
setp <- Weather_Irrig_data$etp[iday] + setp
sirr <- Weather_Irrig_data$irrig[iday] + sirr
sswdef <- sswdef + swdef
sswdfrue <- swdfrue+ sswdfrue
sswdflai <- swdflai+ sswdflai
srgx <- ASTRONO$astrgx +srgx
setm <- etm+setm
setr <- etr+setr
stmp <- tmp+stmp
strp <- trp+strp
srunoff <- runoff+srunoff
sevap <- es+sevap
sdrain <- deepdrain+sdrain

######## Re-initialising value for the next day ########
emergsdd <- EMERG$emergsdd
nbstem <- EMERG$nbstem
sddlai <- LAIGLOB$sddlai
lai <- LAIGLOB$lai
ei <- INTERCEP$ei
pari <- INTERCEP$pari
dmt <- CONVERT$dmt 
dmroot <- PARTDMT$dmroot 
dmbaer <- PARTDMT$dmbaer 
dmstm <- PARTSTEM$dmstm 
dmsug <- PARTSUGAR$dmsug 
fmstm <- HUMSTEM$fmstm
sug_stm <- HUMSTEM$sug_stm
hum_stm <- HUMSTEM$hum_stm
hum_aer <- HUMCAN$hum_aer
humaersdd <- HUMCAN$humaersdd
jat <- PARTSTEM$jat


######## Writing output, case IFREQ = 1 #######
if (Input_read$IFREQ==1){ILIST=Itreat; LINE=iday ; WRITE=T}

######## Writing output, case IFREQ = 0 #######
ENDDAY <- as.numeric(as.Date(Treatment_ID$enddate,format='%d/%m/%Y') -
                       as.Date(Treatment_ID$begdate,format='%d/%m/%Y'))+1

if((Input_read$IFREQ==0)&(iday==ENDDAY)&(ISENS != 1)){ILIST=Itreat;LINE=1; WRITE=T}
  
######## Writing output, case sensitivity analysis #######
if((Input_read$IFREQ==0)&(iday==ENDDAY)&(ISENS==1)){ILIST=Isim ; LINE=1; WRITE=T}

if(WRITE==T){
MATOUT[[ILIST]]$nosimul[LINE] <- Itreat
MATOUT[[ILIST]]$treatcode[LINE] <- Treatment_ID$treatcode
MATOUT[[ILIST]]$lati[LINE] <- Treatment_ID$trialat
MATOUT[[ILIST]]$longi[LINE] <- Treatment_ID$trialong
MATOUT[[ILIST]]$alti[LINE] <- Treatment_ID$trialalt
MATOUT[[ILIST]]$dates[LINE] <- as.character(as.Date(Treatment_ID$begdate,format='%d-%m-%Y')+iday-1)
MATOUT[[ILIST]]$age[LINE] <- iday-1
MATOUT[[ILIST]]$agestem[LINE] <- PARTSTEM$jat
MATOUT[[ILIST]]$glai[LINE] <- LAIGLOB$lai
MATOUT[[ILIST]]$ei[LINE] <- INTERCEP$ei
MATOUT[[ILIST]]$stemnb[LINE] <- EMERG$nbstem
MATOUT[[ILIST]]$stemfm[LINE] <-  HUMSTEM$fmstm / 1000 * 10 # conversion from g/m2 to t/ha
MATOUT[[ILIST]]$stemdm[LINE] <- PARTSTEM$dmstm / 1000 * 10 # conversion from g/m2 to t/ha
MATOUT[[ILIST]]$stemwc[LINE] <- HUMSTEM$hum_stm
MATOUT[[ILIST]]$stemsu[LINE] <- PARTSUGAR$dmsug / 1000 * 10 # conversion from g/m2 to t/ha
MATOUT[[ILIST]]$stemsuc[LINE] <- HUMSTEM$sug_stm
MATOUT[[ILIST]]$rue[LINE] <- CONVERT$gdmt/CONVERT$pari
MATOUT[[ILIST]]$agrfm[LINE] <- PARTDMT$dmbaer / (1-HUMCAN$hum_aer/100) / 1000 * 10 # conversion from g/m2 to t/ha 
MATOUT[[ILIST]]$agrdm[LINE] <- PARTDMT$dmbaer / 1000 * 10 # conversion from g/m2 to t/ha
MATOUT[[ILIST]]$agrwc[LINE] <- HUMCAN$hum_aer
MATOUT[[ILIST]]$bladm[LINE] <- (PARTDMT$dmbaer - PARTSTEM$dmstm)  / 1000 * 10 # conversion from g/m2 to t/ha
MATOUT[[ILIST]]$rootdm[LINE] <- PARTDMT$dmroot / 1000 * 10 # conversion from g/m2 to t/ha
MATOUT[[ILIST]]$hv1[LINE] <- swc[1]  
MATOUT[[ILIST]]$hv2[LINE] <- swc[2]  
MATOUT[[ILIST]]$hv3[LINE] <- swc[3]  
MATOUT[[ILIST]]$hv4[LINE] <- swc[4]  
MATOUT[[ILIST]]$hv5[LINE] <- swc[51]  
MATOUT[[ILIST]]$stock[LINE] <- astock  
MATOUT[[ILIST]]$stockfr[LINE] <- astockfr  
MATOUT[[ILIST]]$rootfront[LINE] <- rootfront  
MATOUT[[ILIST]]$droot1[LINE] <- rlv[1]  
MATOUT[[ILIST]]$droot2[LINE] <- rlv[2]  
MATOUT[[ILIST]]$droot3[LINE] <- rlv[3]  
MATOUT[[ILIST]]$droot4[LINE] <- rlv[4]  
MATOUT[[ILIST]]$droot5[LINE] <- rlv[5]  
MATOUT[[ILIST]]$droot6[LINE] <- rlv[6]  
MATOUT[[ILIST]]$precip[LINE] <-  Weather_Irrig_data$rr[iday]
MATOUT[[ILIST]]$tmin[LINE] <-  Weather_Irrig_data$tn[iday]
MATOUT[[ILIST]]$tmax[LINE] <-  Weather_Irrig_data$tx[iday]
MATOUT[[ILIST]]$tmean[LINE] <-  tmo
MATOUT[[ILIST]]$solrad[LINE] <- Weather_Irrig_data$rg[iday]
MATOUT[[ILIST]]$etp[LINE] <- Weather_Irrig_data$etp[iday]
MATOUT[[ILIST]]$etm[LINE] <- etm
MATOUT[[ILIST]]$etr[LINE] <- etr
MATOUT[[ILIST]]$tmp[LINE] <- tmp
MATOUT[[ILIST]]$trp[LINE] <- trp
MATOUT[[ILIST]]$es[LINE] <- es
MATOUT[[ILIST]]$eos[LINE] <- eos
MATOUT[[ILIST]]$runoff[LINE] <- runoff
MATOUT[[ILIST]]$deepdrain[LINE] <- deepdrain
MATOUT[[ILIST]]$swdef[LINE] <- swdef
MATOUT[[ILIST]]$swdfrue[LINE] <- swdfrue
MATOUT[[ILIST]]$swdflai[LINE] <- swdflai
MATOUT[[ILIST]]$pari[LINE] <- INTERCEP$pari
MATOUT[[ILIST]]$stmean[LINE] <- stmo
MATOUT[[ILIST]]$spar[LINE] <- spar
MATOUT[[ILIST]]$durjour[LINE] <- ASTRONO$astdurjour
MATOUT[[ILIST]]$srgx[LINE] <- srgx
MATOUT[[ILIST]]$spari[LINE] <- spari
MATOUT[[ILIST]]$spr[LINE] <- spr
MATOUT[[ILIST]]$sirr[LINE] <- sirr
MATOUT[[ILIST]]$setp[LINE] <- setp
MATOUT[[ILIST]]$setm[LINE] <- setm
MATOUT[[ILIST]]$setr[LINE] <- setr
MATOUT[[ILIST]]$stmp[LINE] <- stmp
MATOUT[[ILIST]]$strp[LINE] <- strp
MATOUT[[ILIST]]$swdefm[LINE] <- sswdef/iday
MATOUT[[ILIST]]$swdfruem[LINE] <- sswdfrue/iday
MATOUT[[ILIST]]$swdflaim[LINE] <- sswdflai/iday
MATOUT[[ILIST]]$sdrain[LINE] <- sdrain
MATOUT[[ILIST]]$srunoff[LINE] <- srunoff
MATOUT[[ILIST]]$sevap[LINE] <- sevap

}
