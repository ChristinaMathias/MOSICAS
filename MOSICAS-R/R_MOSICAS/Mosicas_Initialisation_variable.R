#####################################
#### Initialisation variables #######

# output
ISENS=0 # 1= sensitivity analysis
WRITE=F # output flag

## calculated output variables #####
spari = 0 ; spr = 0 ; sirr = 0
setp =0;sswdef =0;sswdfrue =0;sswdflai =0
srgx =0;setm =0;setr =0;stmp =0
strp =0;srunoff =0;sevap =0;sdrain =0
stmo =0;spar=0

## Plant growth

## thermal time
sdd=0
emergsdd=0;sddlai=0.;
humaersdd=0

### blades

### stalks

nbstem=0.;
astdurjourp=12.
hum_stm=0.
sug_stm=0

### lai
lai=0.;ei=0.
### stalk height
htvd=0

### biomasses
dmt=0.
dmstru=0;#mst=0;
dmbaer=0.;
dmstm=0;dmsug=0
dmstst=0;dmtbl=0;dmgbl=0;dmroot=0;
fmstm=0

## soil ####
swc <- Soilini_input$swc
swcwp <- Soillayer_input$swcwp
swcfc <- Soillayer_input$swcfc # dul
swcsat <- Soillayer_input$swcsat # sat
swcwp <- Soillayer_input$swcwp
swcsat_mm <- Calcsoilwatercapacities$swcsat_mm 
dlayer <- Soillayer_input$thickness

## Soil evaporation and soil parameters- RItchie 1972 ####
sumes1=100.*(1-sum((swc-swcwp)/(swcfc-swcwp)*dlayer)/sum(dlayer))#Calcsoilwatercapacities$tsw) 
sumes2=0.
swef=0.9-0.00038*(Soillayer_input$thickness[1]-30)**2
t=0
U=5.0
alphaes=3.5
gam=0.019;
Ks=4.5;A=82.5


# plant transpiration and water stress
psilb=PARAM$psilb
swdfrue=1;swdflai=1;swdef=1

### Vertisol parameters
Wver=0.3; K2ver=0.6  

###  divers,général

nage=0 #initialisation days counter
jat=0; forc1eff=0

############# Initialisation of soil water balance parameters 
#if(!is.na(Soilgen_input$soiltype)){
#  if (Soilgen_input$soiltype=="latosol") {rootgrowth=0.03; roottb=0.5}}


####### root profile initialisation  #####
nlayer <- max(Soillayer_input$layernum)
dlayer = Soillayer_input$thickness
rldmax = PARAM$rldmax
rwumax = PARAM$rwumax
roottb = PARAM$roottb # 0° by default

### limitation of profrac (root depth) down to soil depth
# total soil depth
dep <- sum(dlayer[1:nlayer])
# root front depth, cm
rootfront=0 
# root front growth velocity, cm/°d
rootgrowth = PARAM$rootgrowth #0.08 valeur dans Ceres
# maximum root front, cm
rootfrontmax = min(c(dep,Treatment_ID$rootdepth))
## initialisation of root length density distribution, cm/cm3
rlv <- rep(0,nlayer)

## initialisation rootfront and rlv when ratoons
if (Treatment_ID$cycle>0) {
  dep <- 0 ; for (i in 2:nlayer){dep <- c(dep,dlayer[i-1]+dep[i-1])}
  rootfront=rootfrontmax
  
  for (i in 1:nlayer){
    if(rootfront>dep[i]){
      if (rootfront >= (dep[i]+dlayer[i])) {
        v1=384649.* (1/(dep[i]+30))^(3.4163)
        v2=384649.* (1/(dep[i]+dlayer[i]+30))^(3.4163)
        rlv[i]=(v1+v2)/2}
      if (rootfront < (dep[i]+dlayer[i])){
        v1=384649.* (1/(dep[i]+30))^(3.4163)
        v2=384649.* (1/(rootfront+30))^(3.4163)
        rlv[i]=(v1+v2)/2}
      }}
}




