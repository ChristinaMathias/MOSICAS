################################################
#####  daily MOSIWAS growth model for one simulatio ##
################################################
## Daily variable calculation 
source('R_MOSICAS/Mosicas_variable_day.R')

## astronomic calculation
ASTRONO <- Astrono(jcal,lati=Treatment_ID$trialat,longi=Treatment_ID$trialong,alti=Treatment_ID$trialalt,
                   astlatrad,astdr,astdecl,astomeg,astrgo,astrgx,astdurjour,astdurjourp)

## Soil water balance
source('R_MOSICAS/Mosicas_waterbalance.R')


################################################
##### Growth ######################
################################################

### Stalks appearance and senescence
EMERG <- emerg(tempmn,tempmx,emergsdd,nbstem, emergtb=PARAM$emergtb, 
                 emergtt=PARAM$emergtt)

### Calculation of Green leaf area index and interception
LAIGLOB <- laiglob(tempmn,tempmx,laitb=PARAM$laitb,laigrowth=PARAM$laigrowth,laiwksen=PARAM$laiwksen,
                   laimax=PARAM$laimax,
                   swdflai,sddlai,lai,nbstem=EMERG$nbstem)

INTERCEP <-  intercep(isobei,obei=NULL,ke=PARAM$ke,lai=LAIGLOB$lai,solrad,ei,par,pari)

### Convertion into above ground dry mass (convertion)
CONVERT <- convert(ruemax=PARAM$ruemax,ruetk=PARAM$ruetk,ruetopt=PARAM$ruetopt,
                   tempmn,tempmx,pari=INTERCEP$pari,
                   solrad,astrgx=ASTRONO$astrgx,swdfrue,gdmt,dmt,p01=PARAM$p01,p04=PARAM$p04)

### partition into aboveground dry mass
PARTDMT <-  partdmt(Treatment_ID$cycle,gdmt=CONVERT$gdmt,sdd=stmo, dmt=CONVERT$dmt,
                    gdmroot,dmroot,gdmbaer,dmbaer,
                    prootini=PARAM$prootini,prootdec=PARAM$prootdec,prootend=PARAM$prootend)

### Partition of above ground into stem dry mass
PARTSTEM <- partstem (pstemini=PARAM$pstemini,pstemend=PARAM$pstemend,pstemdec=PARAM$pstemdec,
                      dmbaer=PARTDMT$dmbaer,gdmbaer=PARTDMT$gdmbaer,
                      jat,gdmstm,dmstm)

### Partition of millable stalks into sugar
PARTSUGAR <-  partsugar(tempmn, tempmx,gdmstm=PARTSTEM$gdmstm,dmstm=PARTSTEM$dmstm,swdfrue,swdflai,
                        pstrutb=PARAM$pstrutb,pstrutgrowth=PARAM$pstrutgrowt,
                        psugend=PARAM$psugend,pstrudec=PARAM$pstrudec,
                        gdmsug,dmsug)

### Yield and stem water content calculation
HUMSTEM <- humstem(dmstm=PARTSTEM$dmstm,hum_stm,dmsug=PARTSUGAR$dmsug,fmstm,sug_stm,jat,
                         humstemini=PARAM$humstemini,humstemdec=PARAM$humstemdec)

### Above ground water content calculation
HUMCAN <- humcan(humaerini=PARAM$humaerini, humaerend=PARAM$humaerend, humaertb=PARAM$humaertb, 
                 humaerdec=PARAM$humaerdec,
                 hum_aer,tempmx,tempmn,humaersdd,dmbaer=PARTDMT$dmbaer)

