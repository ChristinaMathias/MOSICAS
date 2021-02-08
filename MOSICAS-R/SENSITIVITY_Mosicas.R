###################################################
########## mosicas sensitivity analysis ###########
###################################################
rm(list=ls())
library(sensitivity)

############################################################
#### Input and output folder have to be in the same place than "SENSITIVITY_Mosicas.R"
setIni <- setwd(c(dirname(rstudioapi::getActiveDocumentContext()$path)))
setInput <- paste(c(dirname(rstudioapi::getActiveDocumentContext()$path)),"/input",sep='') 
setOutput <- paste(c(dirname(rstudioapi::getActiveDocumentContext()$path)),"/output",sep='') 

# list of parameter to be used in the sensitivity analysis
MATPAR <- read.csv2(paste(setInput,'/Plant_parameter.csv',sep=''),header=T,stringsAsFactors = F,dec='.')
# all parameters
ParnameALL <- MATPAR$paramname[which(MATPAR$varcode=='R570')]
# specific parameters
Parname <- c("humaerdeb","humaerdec","humaerfin","humaertb", "humstemdeb","humstemdec", 
             "kcmax","ke","laicroi","laimax","laitb","laiwksen",
             "p01","p04","pracdeb","pracdec","pracfin", 
             "psilb", "pstrudec", "pstrufin", "pstrutb","pstrutcroi","pstemdeb","pstemdec",
             "pstemfin", "ruemax","ruetk", "ruetopt","sthydbio", "sthydcroi",  
             "taldebtt")
Parname <- ParnameALL
Nbpar <- length(Parname)

# initial value of parameter
Valpar <- MATPAR$plantcrop[which(MATPAR$varcode=='R570')] # plant crop or ratoon crop
Valpar <- Valpar[ParnameALL %in% Parname]

# limits, here +/-5%
Valpar[which(Valpar==0)] <- 0.001 
Bsup <- Valpar*1.05
Binf <- Valpar*0.95
bounds <- apply(cbind(Binf, Bsup), 1, function(x){list(min = x[1], max = x[2])})

# name of the output
ID_SENS = 'SENSI_PLU_R570_SYP'

###################################################
######### sampling design #########################
### choice between: Moriss, Fast or Sobol methods
Method = 'Fast'  # alternative Moriss / Fast / Variance based analysis
rsensM=100 #default 100
rsensF=350 #default 450
rsensV=2000 #default 1000

if(Method=='Moriss'){
  sa<-morris(model=NULL,
            factors = Parname, 
            r = rsensM, #default 100 
            binf=Binf, 
            bsup=Bsup, 
            design = list(type = "oat", levels = 10, grid.jump = 5))

  # plot sampling
  plot(sa$X)
}

if(Method=='Fast'){
  sa <- fast99(model = NULL,
                factors = Parname,
                n = rsensF,
                q = rep("qunif", Nbpar),
                q.arg=bounds)
  plot(sa$X[,1:2])
}

if(Method=='Variance'){
  library(lhs)
  # we will test only 2nd order 
  Nblimit = factorial(Nbpar)/(factorial(2)*factorial(Nbpar-2))
  if(rsensV < 2*Nblimit){print('Warning, not enough sammpling combinations for interaction analysis')}
  sa <- list()
  #  sa$X <- improvedLHS(rsensV,Nbpar)
  sa$X <- randomLHS(rsensV,Nbpar)
  for(i in 1:Nbpar){ sa$X[,i] <- sa$X[,i]*(Bsup[i]-Binf[i]) + Binf[i]}
  sa$X <- as.data.frame(sa$X)
  names(sa$X) <- Parname
  plot(sa$X[,1:2])
}

###################################################
########## run simulations

Newpar = sa$X
# sensitivity function
source(paste(setIni,'R_MOSICAS/Mosicas_function_sensitivity.R',sep='/'))
# output
y <- MATOUTTEMP

# save
save.image(paste(setOutput,'/SENS_Mosicas_',Method,'_',ID_SENS,sep=''))



###################################################
########## sensitivity 

##### Variance based methods ######################
 ### exemple

Method='Variance'
load(paste(setOutput,'/SENS_Mosicas_',Method,'_',ID_SENS,sep=''))

Varout=c('spari','agrfm','stemfm','stemsu')
PARSI <- colnames(SAFAST$X)
SIP <- list()
SIT <- list()

for(j in 1:length(Varout)){

  SAFAST <- sa
  SAFAST$y <- y[,Varout[j]]
  tell(SAFAST)
  SIP[[j]] <- (SAFAST$D1/SAFAST$V)
  SIT[[j]] <- 1-SAFAST$Dt/SAFAST$V
  MATIS[[j]] <- rbind(SIT[[j]],SIT[[j]]-SIP[[j]])

}

#x11()
PARSI[ORDER]

ORDER <- rev(c(31,32,9,10,10,12,8,14,26,27,28,13,
           15,16,17,23,24,25,19,20,21,22,1:6,
           7,18,29,30, 33,34,35,36))

par(mfrow=c(1,4),oma=c(4,6,1,1),mar=c(0,1,0,0))
barplot(MATIS[[1]][,ORDER],names.arg = PARSI[ORDER],las=2,horiz=T,mgp=c(2,0.5,0),col=gray(c(0.8,0.1)),
        border=NA,xlim=c(0,0.5),width=1,space=0.5,ylim=c(1,53),tck=-0.02)
axis(2,at=(1:length(PARSI))*1.5-0.5,labels=rep('',length(PARSI)),tck=-0.02)
box();text(0.25,53,labels=expression(Sigma~'iPAR'),cex=1.2)
abline(h=1.5*c(3,7,13,17,23,27,29,33)+1+0.5+0.25,lty=2)

barplot(MATIS[[2]][,ORDER],names.arg = rep('',length(PARSI)),las=2,horiz=T,mgp=c(2,0.5,0),col=gray(c(0.8,0.1)),
        border=NA,xlim=c(0,0.7),width=1,space=0.5,ylim=c(1,53),tck=-0.02)
axis(2,at=(1:length(PARSI))*1.5-0.5,labels=rep('',length(PARSI)),tck=-0.02)
box();text(0.35,53,labels=expression(ABV[FM]),cex=1.2)
abline(h=1.5*c(3,7,13,17,23,27,29,33)+1+0.5+0.25,lty=2)

barplot(MATIS[[3]][,ORDER],names.arg = rep('',length(PARSI)),las=2,horiz=T,mgp=c(2,0.5,0),col=gray(c(0.8,0.1)),
        border=NA,xlim=c(0,0.4),width=1,space=0.5,ylim=c(1,53),tck=-0.02)
axis(2,at=(1:length(PARSI))*1.5-0.5,labels=rep('',length(PARSI)),tck=-0.02)
box();text(0.2,53,labels=expression(Sta[FM]),cex=1.2)
abline(h=1.5*c(3,7,13,17,23,27,29,33)+1+0.5+0.25,lty=2)

barplot(MATIS[[4]][,ORDER],names.arg = rep('',length(PARSI)),las=2,horiz=T,mgp=c(2,0.5,0),col=gray(c(0.8,0.1)),
        border=NA,xlim=c(0,0.4),width=1,space=0.5,ylim=c(1,53),tck=-0.02)
axis(2,at=(1:length(PARSI))*1.5-0.5,labels=rep('',length(PARSI)),tck=-0.02)
box();text(0.2,53,labels=expression(Sug[DM]),cex=1.2)
abline(h=1.5*c(3,7,13,17,23,27,29,33)+1+0.5+0.25,lty=2)
mtext('Fast99 Sensitivity Index',side = 1,outer=T,cex=1.2,line=3)



# fast indice calculation
if(Method=='Fast'){
Ifast1 <- sa$D1/sa$V
Ifastt <- 1-sa$Dt/sa$V
}


